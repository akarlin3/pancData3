#!/usr/bin/env python3
"""Async batch graph analysis using Google Gemini vision.

Scans the most recent ``saved_files_*`` timestamped folder for all PNG/JPG
images produced by the MATLAB pipeline, sends each to the Gemini Flash
vision model, and extracts structured metadata:

- Axis labels, units, and ranges
- Observed trends (direction, series name, description)
- Inflection points with approximate (x, y) coordinates

Results are validated via strict Pydantic schemas and written to a flat CSV
(``graph_analysis_results.csv``) inside the output folder.  The CSV is then
consumed by downstream analysis scripts (``generate_report.py``,
``cross_reference_dwi.py``, ``statistical_relevance.py``, etc.).

When the ``GEMINI_API_KEY`` is not set or the ``google-genai`` package is
not installed, the script automatically falls back to a **local
filename-based heuristic** that infers graph type, axes, and comparison
type from the structured filenames produced by the MATLAB pipeline.  This
fallback produces lower-fidelity results (no visual analysis) but keeps the
downstream CSV pipeline functional.

Usage:
    python batch_graph_analysis.py                       # auto-detect folder
    python batch_graph_analysis.py /path/to/saved_files  # explicit folder
    python batch_graph_analysis.py --local /path/to/saved_files  # force local
"""

from __future__ import annotations

import asyncio
import base64
import csv
import json
import os
import re
import sys
from pathlib import Path
from typing import Optional

from tqdm import tqdm  # type: ignore

# Ensure analysis/ root is on sys.path so 'shared' is importable when run as subprocess.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import get_config, resolve_folder, setup_utf8_stdout  # type: ignore

# Ensure emoji and special characters print correctly on Windows consoles.
setup_utf8_stdout()

# Lazy-import google.genai so that the module can be imported (and tested)
# without the google-genai package installed.  The actual imports happen
# inside ``analyze_image`` and ``main``, which are only called when running
# the vision pipeline.
genai = None  # type: ignore[assignment]
types = None  # type: ignore[assignment]

from pydantic import BaseModel, Field, ValidationError  # type: ignore


def _ensure_genai():
    """Import google.genai lazily; raise ImportError if unavailable."""
    global genai, types
    if genai is None:
        from google import genai as _genai  # type: ignore
        from google.genai import types as _types  # type: ignore
        genai = _genai
        types = _types


# ── Pydantic schema ──────────────────────────────────────────────────────────
# These models define the expected JSON structure returned by the vision model.
# Pydantic validates the raw JSON and provides clear error messages when the
# model response deviates from the schema.

class Axis(BaseModel):
    """Description of a single axis (x, y, or colour-bar)."""
    label: str = Field(..., description="Axis label text as shown on the graph")
    units: Optional[str] = Field(None, description="Units (e.g., mm\u00b2/s, Gy, days)")
    range_min: Optional[float] = Field(None, description="Minimum value on this axis")
    range_max: Optional[float] = Field(None, description="Maximum value on this axis")
    scale_type: Optional[str] = Field(None, description="Scale type: linear, log, categorical")
    tick_count: Optional[int] = Field(None, description="Number of major tick marks visible on this axis")
    tick_labels: list[str] = Field(default_factory=list, description="Tick labels for categorical axes (e.g., ['BL', 'Fx1', 'Fx10', 'W2'])")


class Trend(BaseModel):
    """A single trend observed in the graph (one per data series)."""
    series: Optional[str] = Field(None, description="Name of the data series, if labeled")
    direction: str = Field(..., description="increasing, decreasing, flat, non-monotonic, U-shaped, etc.")
    description: str = Field(..., description="Brief description of the trend")
    magnitude: Optional[str] = Field(None, description="Approximate magnitude of change (e.g., '~20% decline', '0.001 to 0.002')")
    statistical_significance: Optional[str] = Field(None, description="P-value or significance level noted for this trend, if any")
    confidence_band: Optional[str] = Field(None, description="Confidence/uncertainty band description if visible (e.g., '95% CI shaded region')")
    start_value: Optional[float] = Field(None, description="Approximate starting y-value of this data series")
    end_value: Optional[float] = Field(None, description="Approximate ending y-value of this data series")


class InflectionPoint(BaseModel):
    """An inflection point or notable change in curvature / direction.

    Primarily relevant for line and trajectory plots; for other graph types
    the vision model is instructed to return an empty list.
    """
    approximate_x: Optional[float] = Field(None, description="Approximate x-coordinate")
    approximate_y: Optional[float] = Field(None, description="Approximate y-coordinate")
    description: str = Field(..., description="What happens at this inflection point")
    magnitude: Optional[float] = Field(None, description="Approximate magnitude of the change at this point")


class StatisticalTest(BaseModel):
    """A statistical test result visible on or described in the graph."""
    test_name: str = Field(..., description="Name of the test (e.g., Wilcoxon, t-test, ANOVA, log-rank)")
    statistic_value: Optional[float] = Field(None, description="Test statistic value if shown")
    p_value: Optional[float] = Field(None, description="P-value if shown")
    comparison_groups: Optional[str] = Field(None, description="Groups being compared (e.g., 'LF vs LC')")


class Outlier(BaseModel):
    """A notable outlier data point visible in the graph."""
    approximate_x: Optional[float] = Field(None, description="Approximate x-coordinate of the outlier")
    approximate_y: Optional[float] = Field(None, description="Approximate y-coordinate of the outlier")
    series: Optional[str] = Field(None, description="Data series the outlier belongs to, if identifiable")
    description: str = Field(..., description="Brief description (e.g., 'extreme high ADC value')")


class ReferenceLine(BaseModel):
    """A reference line, threshold, or boundary shown on the graph."""
    orientation: str = Field(..., description="horizontal, vertical, or diagonal")
    value: Optional[float] = Field(None, description="Value where the line is placed (y for horizontal, x for vertical)")
    label: Optional[str] = Field(None, description="Label text for the reference line (e.g., 'threshold = 0.001')")
    style: Optional[str] = Field(None, description="Line style: solid, dashed, dotted")


class GraphAnalysis(BaseModel):
    """Full structured analysis of a single graph image.

    This is the top-level Pydantic model that the vision model's JSON
    response is parsed into.  Each instance is later flattened into a CSV
    row by :func:`flatten`.
    """
    file_path: str = Field(..., description="Path to the image file")
    graph_title: Optional[str] = Field(None, description="Title shown on the graph")
    graph_type: str = Field(..., description="Type of graph: line, bar, scatter, heatmap, box, histogram, parameter_map, etc.")
    x_axis: Optional[Axis] = Field(None, description="X-axis details")
    y_axis: Optional[Axis] = Field(None, description="Y-axis details")
    color_axis: Optional[Axis] = Field(None, description="Colorbar / color-axis if present")
    trends: list[Trend] = Field(default_factory=list, description="List of observed trends")
    inflection_points: list[InflectionPoint] = Field(default_factory=list, description="Notable inflection points")
    statistical_tests: list[StatisticalTest] = Field(default_factory=list, description="Statistical tests visible or referenced in the graph")
    outliers: list[Outlier] = Field(default_factory=list, description="Notable outlier data points visible in the graph")
    reference_lines: list[ReferenceLine] = Field(default_factory=list, description="Reference lines, thresholds, or boundaries shown on the graph")
    issues: list[str] = Field(default_factory=list, description="List of detected quality issues with the graph")
    summary: str = Field(..., description="One-paragraph plain-English summary of the graph")
    sample_size: Optional[str] = Field(None, description="Sample size annotation if visible (e.g., 'n=45', 'N=120')")
    data_series_count: Optional[int] = Field(None, description="Number of distinct data series visible in the graph")
    error_bars: Optional[str] = Field(None, description="Type of error bars if present: SD, SEM, CI, IQR, or null")
    annotations: list[str] = Field(default_factory=list, description="Text annotations visible on the graph (e.g., significance brackets, data labels)")
    clinical_relevance: Optional[str] = Field(None, description="Brief clinical interpretation of the findings")
    data_density: Optional[str] = Field(None, description="Data density: sparse, moderate, dense")
    spatial_pattern: Optional[str] = Field(None, description="For parameter maps: uniform, gradient, focal, diffuse, heterogeneous")
    legend_items: list[str] = Field(default_factory=list, description="Legend entry labels visible on the graph")
    subpanel_count: Optional[int] = Field(None, description="Number of subpanels/subplots in a multi-panel figure (null if single panel)")
    comparison_type: Optional[str] = Field(None, description="Type of comparison: paired, unpaired, longitudinal, cross-sectional, dose-response")
    figure_quality: Optional[str] = Field(None, description="Overall rendering quality: high, medium, low")


# ── Helpers ──────────────────────────────────────────────────────────────────


def collect_images(folder: Path) -> list[Path]:
    """Recursively collect all PNG/JPG images under *folder*.

    Parameters
    ----------
    folder : Path
        Root directory to search (typically a ``saved_files_*`` folder).

    Returns
    -------
    list[Path]
        Sorted list of image file paths found in the folder tree.
    """
    exts = {".png", ".jpg", ".jpeg"}
    images = []
    for f in sorted(folder.rglob("*")):
        if f.suffix.lower() in exts and f.is_file():
            images.append(f)
    return images


def image_to_base64(path: Path) -> str:
    """Read an image file and return its base64-encoded string.

    The base64 string is used to embed the image in API request payloads.

    Parameters
    ----------
    path : Path
        Path to the image file.

    Returns
    -------
    str
        Base64-encoded file contents (ASCII string).
    """
    return base64.standard_b64encode(path.read_bytes()).decode("utf-8")


def media_type_for(path: Path) -> str:
    """Return the MIME type string for an image file.

    Parameters
    ----------
    path : Path
        Path to the image file.

    Returns
    -------
    str
        MIME type (e.g. ``"image/png"`` or ``"image/jpeg"``).
        Defaults to ``"image/png"`` for unrecognised extensions.
    """
    ext = path.suffix.lower()
    return {
        ".png": "image/png",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
    }.get(ext, "image/png")


# ── Local filename-based fallback analyzer ───────────────────────────────────
# When the Gemini API is unavailable, this function infers graph metadata
# from the structured filenames produced by the MATLAB pipeline.  The
# results are lower-fidelity (no visual analysis) but keep the downstream
# CSV pipeline functional for report generation and cross-DWI comparison.

# Filename keywords → graph_type mapping.  Order matters: first match wins.
_GRAPH_TYPE_KEYWORDS: list[tuple[list[str], str]] = [
    (["heatmap", "dice_heatmap", "hausdorff_heatmap"], "heatmap"),
    (["parameter_map", "param_map", "overlay"], "parameter_map"),
    (["histogram", "hist"], "histogram"),
    (["boxplot", "box_plot", "box"], "box"),
    (["scatter", "correlation"], "scatter"),
    (["violin"], "violin"),
    (["bar", "volume_comparison"], "bar"),
    (["longitudinal", "trajectory", "timeseries", "time_series"], "line"),
    (["kaplan", "survival", "km_curve"], "line"),
    (["roc", "auc"], "line"),
    (["dose_vs", "dvh"], "line"),
]

# Filename keywords → axis label / units hints.
_AXIS_HINTS: dict[str, dict[str, str | None]] = {
    "adc": {"label": "ADC", "units": "mm\u00b2/s"},
    "ivim": {"label": "IVIM Parameter", "units": None},
    "d_mean": {"label": "D (mean)", "units": "mm\u00b2/s"},
    "f_mean": {"label": "f (mean)", "units": None},
    "dstar": {"label": "D* (mean)", "units": "mm\u00b2/s"},
    "dose": {"label": "Dose", "units": "Gy"},
    "time": {"label": "Time", "units": "days"},
    "longitudinal": {"label": "Timepoint", "units": None},
    "survival": {"label": "Time", "units": "days"},
    "kaplan": {"label": "Time", "units": "days"},
    "dice": {"label": "Method", "units": None},
    "hausdorff": {"label": "Method", "units": None},
    "volume": {"label": "Method", "units": None},
}

# Filename keywords → comparison_type hints.
_COMPARISON_HINTS: dict[str, str] = {
    "longitudinal": "longitudinal",
    "trajectory": "longitudinal",
    "byoutcome": "unpaired",
    "lf_vs_lc": "unpaired",
    "dose_vs": "dose-response",
    "dvh": "dose-response",
    "repeat": "paired",
    "baseline": "cross-sectional",
}


def _infer_graph_type(name_lower: str) -> str:
    """Infer graph_type from a lowercased filename stem."""
    for keywords, gtype in _GRAPH_TYPE_KEYWORDS:
        for kw in keywords:
            if kw in name_lower:
                return gtype
    return "unknown"


def _infer_axes(name_lower: str) -> tuple[Axis | None, Axis | None]:
    """Infer x/y axis labels from filename keywords."""
    x_axis = None
    y_axis = None
    for keyword, hints in _AXIS_HINTS.items():
        if keyword in name_lower:
            # For time-related x-axes, set x_axis
            if keyword in ("time", "longitudinal", "survival", "kaplan"):
                x_axis = Axis(
                    label=hints["label"],
                    units=hints.get("units"),
                )
            else:
                # For metric keywords, set as y_axis
                y_axis = Axis(
                    label=hints["label"],
                    units=hints.get("units"),
                )
    return x_axis, y_axis


def _infer_comparison_type(name_lower: str) -> str | None:
    """Infer comparison_type from filename keywords."""
    for keyword, ctype in _COMPARISON_HINTS.items():
        if keyword in name_lower:
            return ctype
    return None


def analyze_image_local(image_path: Path) -> GraphAnalysis:
    """Produce a best-effort ``GraphAnalysis`` from filename heuristics.

    This function does **not** inspect the image pixels.  It uses the
    structured naming conventions of the MATLAB pipeline
    (e.g. ``Longitudinal_Mean_Metrics_Standard.png``,
    ``core_method_dice_heatmap.png``) to infer graph type, axis labels,
    and comparison type.

    Parameters
    ----------
    image_path : Path
        Path to the graph image.

    Returns
    -------
    GraphAnalysis
        A structured analysis with ``graph_type``, axis hints, and a
        summary noting that the result was produced by the local fallback.
    """
    rel_path = str(image_path)
    stem = image_path.stem  # filename without extension
    name_lower = stem.lower()

    graph_type = _infer_graph_type(name_lower)
    x_axis, y_axis = _infer_axes(name_lower)
    comparison_type = _infer_comparison_type(name_lower)

    # Build a human-readable summary from the filename.
    pretty_name = stem.replace("_", " ")
    summary = (
        f"Local fallback analysis of '{pretty_name}' "
        f"(inferred type: {graph_type}). "
        "No visual analysis was performed; metadata was inferred from the "
        "filename. Re-run with GEMINI_API_KEY set for full vision analysis."
    )

    # Detect DWI type from path for clinical context.
    dwi_types = {"Standard", "dnCNN", "IVIMnet"}
    clinical_note = None
    for part in image_path.parts:
        if part in dwi_types:
            clinical_note = f"Graph from {part} DWI processing pipeline."
            break

    return GraphAnalysis(
        file_path=rel_path,
        graph_title=pretty_name,
        graph_type=graph_type,
        x_axis=x_axis,
        y_axis=y_axis,
        summary=summary,
        comparison_type=comparison_type,
        clinical_relevance=clinical_note,
        figure_quality="unknown",
    )


SYSTEM_PROMPT = """\
You are an expert scientific figure analyst specializing in medical physics and \
MRI diffusion-weighted imaging. You will receive a single graph image from a \
pancreatic DWI analysis pipeline.

Your task: extract comprehensive structured information about the graph.

Respond with ONLY a valid JSON object matching this exact schema (no markdown, \
no commentary, no code fences):

{
  "graph_title": "<string or null>",
  "graph_type": "<line|bar|scatter|heatmap|box|histogram|parameter_map|violin|other>",
  "x_axis": {"label": "...", "units": "...", "range_min": ..., "range_max": ..., "scale_type": "linear|log|categorical", "tick_count": ..., "tick_labels": []} or null,
  "y_axis": {"label": "...", "units": "...", "range_min": ..., "range_max": ..., "scale_type": "linear|log|categorical", "tick_count": ..., "tick_labels": []} or null,
  "color_axis": {"label": "...", "units": "...", "range_min": ..., "range_max": ..., "scale_type": "...", "tick_count": ..., "tick_labels": []} or null,
  "trends": [
    {"series": "...", "direction": "...", "description": "...", "magnitude": "...", "statistical_significance": "...", "confidence_band": "...", "start_value": ..., "end_value": ...}
  ],
  "inflection_points": [
    {"approximate_x": ..., "approximate_y": ..., "description": "...", "magnitude": ...}
  ],
  "statistical_tests": [
    {"test_name": "...", "statistic_value": ..., "p_value": ..., "comparison_groups": "..."}
  ],
  "outliers": [
    {"approximate_x": ..., "approximate_y": ..., "series": "...", "description": "..."}
  ],
  "reference_lines": [
    {"orientation": "horizontal|vertical|diagonal", "value": ..., "label": "...", "style": "solid|dashed|dotted"}
  ],
  "issues": ["..."],
  "summary": "...",
  "sample_size": "<string like 'n=45' or null>",
  "data_series_count": <integer or null>,
  "error_bars": "<SD|SEM|CI|IQR or null>",
  "annotations": ["text annotation 1", "..."],
  "clinical_relevance": "<brief clinical interpretation or null>",
  "data_density": "<sparse|moderate|dense or null>",
  "spatial_pattern": "<uniform|gradient|focal|diffuse|heterogeneous or null>",
  "legend_items": ["item 1", "..."],
  "subpanel_count": <integer or null>,
  "comparison_type": "<paired|unpaired|longitudinal|cross-sectional|dose-response or null>",
  "figure_quality": "<high|medium|low or null>"
}

Rules:
- Use null for any field you cannot determine. Use empty lists [] when no items apply.
- For axis objects:
  - scale_type: "linear" (default), "log" for logarithmic, "categorical" for discrete categories.
  - tick_count: number of major tick marks visible on the axis.
  - tick_labels: for categorical axes, list the category labels in order (e.g., ["BL", "Fx1", "Fx10", "W2"]).
- For trends:
  - magnitude: approximate description of the change (e.g., "~20% decline", "0.001 to 0.002 mm²/s").
  - statistical_significance: p-value or significance level noted for this specific trend.
  - confidence_band: describe any shaded confidence/uncertainty regions (e.g., "95% CI shaded").
  - start_value/end_value: approximate starting and ending y-values if readable from the graph.
- For inflection_points: magnitude is the approximate size of the change/jump at that point. \
Inflection points apply mainly to line/trajectory plots; for other types, return an empty list.
- For statistical_tests: extract any statistical tests visible on the graph or mentioned in annotations \
(e.g., Wilcoxon rank-sum, t-test, ANOVA, log-rank, chi-squared). Include the test statistic, p-value, \
and comparison groups when visible.
- For outliers: report notable outlier data points that stand far from the main distribution or trend. \
Include their approximate coordinates and which data series they belong to.
- For reference_lines: report any horizontal, vertical, or diagonal reference lines, thresholds, or \
boundaries visible on the graph (e.g., dashed threshold lines, identity lines in scatter plots). \
Include the line's orientation, value, label if shown, and line style.
- For sample_size: report any "n=..." or "N=..." annotations visible on the graph.
- For data_series_count: count the number of distinct data series, groups, or categories shown.
- For error_bars: identify the type if labeled or inferable (SD, SEM, CI for confidence intervals, IQR).
- For annotations: list any text annotations overlaid on the data area (significance brackets like \
"*", "**", "***", "ns"; data point labels; threshold markers; correlation or regression equations).
- For clinical_relevance: provide a 1-sentence interpretation in the context of pancreatic cancer \
radiotherapy response monitoring using DWI metrics.
- For data_density: estimate whether the data is sparse (<10 points), moderate (10-50), or dense (>50).
- For parameter maps (spatial images), set graph_type to "parameter_map" and set spatial_pattern to one of: \
"uniform" (homogeneous), "gradient" (directional change), "focal" (localised hotspot), \
"diffuse" (widely spread), or "heterogeneous" (mixed pattern).
- For legend_items: list all legend entry labels in order.
- For subpanel_count: count the number of subpanels/subplots if this is a multi-panel figure (null for single panel).
- For comparison_type: classify the study design visible in the graph.
- For figure_quality: assess the overall rendering quality (resolution, label readability, colour clarity).
- For heatmaps (Dice, Hausdorff), extract axis labels as the methods being compared.
- Keep the summary concise (2-4 sentences) and include key quantitative findings.
- For "issues", report any quality problems you observe in the graph. Common issues include:
  missing or unreadable axis labels, truncated or clipped data, overlapping text or legends,
  axis scale problems (e.g., misleading origin), empty or nearly empty plots, outliers that
  distort the scale, low sample size warnings visible on the graph, rendering artefacts,
  and missing titles or legends. Return an empty list if no issues are found.
"""


# ── Gemini model configuration ──────────────────────────────────────────────
# Loaded from the centralised analysis config (analysis_config.json /
# shared._DEFAULTS).  Values can be overridden without editing source code.
_vision_cfg = get_config()["vision"]

# Model to use for vision analysis.
GEMINI_MODEL = _vision_cfg["gemini_model"]

# ── Async API call ───────────────────────────────────────────────────────────

# Concurrency limit for API requests.  Kept low (2) to stay well within
# Gemini's rate limits.
SEM_LIMIT = _vision_cfg["max_concurrent_requests"]

# Number of retries on rate-limit (429) errors.  Uses exponential
# backoff: 15s, 30s, 60s, 120s.
MAX_RETRIES = _vision_cfg["max_retries"]

# Base seconds for exponential backoff on rate-limit retries.
BACKOFF_BASE = _vision_cfg["backoff_base_seconds"]

# Maximum output tokens for the vision model response.
MAX_OUTPUT_TOKENS = _vision_cfg["max_output_tokens"]

# Per-request timeout in seconds.  If the Gemini API does not respond
# within this window the request is cancelled and treated as a failure.
REQUEST_TIMEOUT = _vision_cfg["request_timeout_seconds"]


class _RateLimitCoordinator:
    """Coordinates rate-limit backoff across all workers.

    When any worker hits a 429/rate-limit error it calls :meth:`signal`,
    which sets a shared cooldown timestamp.  Before each API call workers
    call :meth:`wait_if_cooling_down` to honour the cooldown.  This
    prevents the thundering-herd problem where multiple workers slam the
    API simultaneously after individual backoffs expire.
    """

    def __init__(self):
        self._resume_at: float = 0.0  # monotonic timestamp
        self._lock = asyncio.Lock()
        self._consecutive_hits: int = 0

    async def signal(self):
        """Record a rate-limit hit and extend the global cooldown."""
        async with self._lock:
            self._consecutive_hits += 1
            # Exponential cooldown: base * 2^(hits-1), capped at 4× base.
            wait = min(BACKOFF_BASE * (2 ** (self._consecutive_hits - 1)),
                       BACKOFF_BASE * 4)
            now = asyncio.get_event_loop().time()
            new_resume = now + wait
            if new_resume > self._resume_at:
                self._resume_at = new_resume
                print(f"  [RATE-LIMIT] Global cooldown: {wait:.0f}s "
                      f"(consecutive hits: {self._consecutive_hits})", flush=True)

    async def wait_if_cooling_down(self):
        """Sleep until the global cooldown expires (no-op if not active)."""
        now = asyncio.get_event_loop().time()
        if now < self._resume_at:
            delay = self._resume_at - now
            await asyncio.sleep(delay)

    async def record_success(self):
        """Reset the consecutive-hit counter after a successful request."""
        async with self._lock:
            self._consecutive_hits = 0


def _repair_truncated_json(raw: str) -> dict | None:
    """Attempt to repair JSON truncated by token limits.

    Progressively closes open strings, arrays, and objects from the end
    of the raw text.  Returns the parsed dict on success, or ``None``
    if the response is beyond salvage.
    """
    text = raw.rstrip()
    # Try closing open structures (up to a handful of attempts).
    closers = ['"', "]", "}", "]", "}"]
    for i in range(1, len(closers) + 1):
        candidate = text + "".join(closers[:i])  # type: ignore
        try:
            data = json.loads(candidate)
            if isinstance(data, dict):
                return data
        except json.JSONDecodeError:
            continue
    # Brute-force: strip back to last complete key-value, close object.
    # Find the last complete value (ends with , or a closing bracket).
    for trim in range(1, min(300, len(text))):
        stub = text[:-trim].rstrip().rstrip(",")  # type: ignore
        for suffix in ["}", "]}", "\"}", "\"]}", "\"]}",
                        "]}", "\"]}", "\"]}}"]:
            try:
                data = json.loads(stub + suffix)
                if isinstance(data, dict):
                    return data
            except json.JSONDecodeError:
                continue
    return None


def _is_rate_limit_error(exc: Exception) -> bool:
    """Check whether *exc* is a rate-limit / quota error."""
    err_str = str(exc).lower()
    return ("429" in err_str or "rate" in err_str
            or "resource" in err_str or "quota" in err_str)


async def analyze_image(
    client,
    image_path: Path,
    rate_limiter: _RateLimitCoordinator,
    progress_bar: tqdm | None = None,
) -> GraphAnalysis:
    """Send one image to Gemini and parse the structured JSON response.

    The function handles rate-limit retries (coordinated via a shared
    :class:`_RateLimitCoordinator`), markdown fence stripping, JSON
    parsing, and Pydantic validation.  On any non-recoverable error a
    fallback ``GraphAnalysis`` with ``graph_type="unknown"`` or ``"error"``
    is returned so that the batch can continue.

    Parameters
    ----------
    client : genai.Client
        Authenticated Google Gen AI client.
    image_path : Path
        Path to the graph image to analyse.
    rate_limiter : _RateLimitCoordinator
        Shared coordinator that pauses all workers on 429 errors.

    Returns
    -------
    GraphAnalysis
        Validated structured analysis, or a fallback object on failure.
    """
    image_bytes = await asyncio.to_thread(image_path.read_bytes)  # type: ignore
    mime = media_type_for(image_path)
    rel_path = str(image_path)

    if progress_bar is not None:
        progress_bar.set_postfix_str(f"analyzing {image_path.name}", refresh=True)

    # Retry loop with exponential backoff for rate-limit errors.
    response = None
    for attempt in range(MAX_RETRIES + 1):
        # Honour any global cooldown before sending a request.
        await rate_limiter.wait_if_cooling_down()

        try:
            response = await asyncio.wait_for(
                client.aio.models.generate_content(
                    model=GEMINI_MODEL,
                    contents=[
                        types.Part.from_bytes(  # type: ignore
                            data=image_bytes,
                            mime_type=mime,
                        ),
                        f"Analyze this graph image. File: {image_path.name}\n"
                        "Return ONLY the JSON object described in your instructions.",
                    ],
                    config=types.GenerateContentConfig(  # type: ignore
                        system_instruction=SYSTEM_PROMPT,
                        max_output_tokens=MAX_OUTPUT_TOKENS,
                    ),
                ),
                timeout=REQUEST_TIMEOUT,
            )
            await rate_limiter.record_success()
            break  # Success: exit the retry loop.
        except asyncio.TimeoutError:
            if attempt < MAX_RETRIES:
                wait = 2 ** attempt * BACKOFF_BASE
                print(f"  [TIMEOUT] {image_path.name}: retry {attempt+1}/{MAX_RETRIES} in {wait}s", flush=True)
                await asyncio.sleep(wait)
            else:
                raise  # Give up after exhausting all retries.
        except Exception as e:
            if _is_rate_limit_error(e) and attempt < MAX_RETRIES:
                await rate_limiter.signal()
                # Worker-local backoff on top of global cooldown.
                wait = 2 ** attempt * BACKOFF_BASE
                print(f"  [RATE-LIMIT] {image_path.name}: retry {attempt+1}/{MAX_RETRIES} in {wait}s", flush=True)
                await asyncio.sleep(wait)
            elif _is_rate_limit_error(e):
                raise  # Give up after exhausting all retries.
            else:
                raise  # Non-rate-limit errors are re-raised immediately.

    # Guard against the retry loop exhausting without a successful response.
    if response is None:
        raise RuntimeError(f"No response received for {image_path.name} after {MAX_RETRIES + 1} attempts")

    raw_text = response.text.strip()  # type: ignore

    # Strip markdown code fences if the model wraps its JSON response
    # (e.g. "```json\n{...}\n```").
    if raw_text.startswith("```"):
        raw_text = re.sub(r"^```(?:json)?\s*", "", raw_text)
        raw_text = re.sub(r"\s*```$", "", raw_text)

    # ── JSON parsing (with truncation repair) ──
    try:
        data = json.loads(raw_text)
    except json.JSONDecodeError:
        # Attempt to repair truncated JSON from token-limited responses.
        data = _repair_truncated_json(raw_text)
        if data is None:
            if progress_bar is not None:
                progress_bar.update(1)  # type: ignore
            else:
                print(f"  \u274c  JSON parse error for {image_path.name}", flush=True)
            return GraphAnalysis(
                file_path=rel_path,  # type: ignore
                graph_type="unknown",  # type: ignore
                summary=f"JSON parse error (unrepairable). Raw response: {raw_text[:200]}",  # type: ignore
            )

    # Inject the file path (not returned by the model).
    data["file_path"] = rel_path  # type: ignore

    # ── Pydantic validation ──
    try:
        analysis = GraphAnalysis(**data)
    except ValidationError as e:
        if progress_bar is not None:
            progress_bar.update(1)
        else:
            print(f"  \u26a0\ufe0f  Validation warning for {image_path.name}: {e}", flush=True)
        # Fallback: preserve whatever fields parsed successfully.
        return GraphAnalysis(
            file_path=rel_path,  # type: ignore
            graph_type=data.get("graph_type", "unknown"),  # type: ignore
            summary=data.get("summary", f"Validation error: {e}"),  # type: ignore
        )

    if progress_bar is not None:
        progress_bar.set_postfix_str(image_path.name, refresh=True)
        progress_bar.update(1)
    return analysis


# ── Flatten to CSV rows ──────────────────────────────────────────────────────
# The nested Pydantic models need to be flattened into a single-level dict
# for CSV serialisation.  Axis fields are prefixed (``x_axis_``, ``y_axis_``,
# ``color_axis_``) and complex sub-structures (trends, inflection points)
# are serialised as JSON strings.

CSV_COLUMNS = [
    "file_path",
    "graph_title",
    "graph_type",
    # Axis fields (5 per axis: label, units, range_min, range_max, scale_type)
    "x_axis_label",
    "x_axis_units",
    "x_axis_range_min",
    "x_axis_range_max",
    "x_axis_scale_type",
    "x_axis_tick_count",
    "x_axis_tick_labels_json",
    "y_axis_label",
    "y_axis_units",
    "y_axis_range_min",
    "y_axis_range_max",
    "y_axis_scale_type",
    "y_axis_tick_count",
    "y_axis_tick_labels_json",
    "color_axis_label",
    "color_axis_units",
    "color_axis_range_min",
    "color_axis_range_max",
    "color_axis_scale_type",
    "color_axis_tick_count",
    "color_axis_tick_labels_json",
    # Complex nested fields
    "num_trends",
    "trends_json",
    "num_inflection_points",
    "inflection_points_json",
    "num_statistical_tests",
    "statistical_tests_json",
    "num_outliers",
    "outliers_json",
    "num_reference_lines",
    "reference_lines_json",
    "num_issues",
    "issues_json",
    # Scalar metadata fields
    "summary",
    "sample_size",
    "data_series_count",
    "error_bars",
    "num_annotations",
    "annotations_json",
    "clinical_relevance",
    "data_density",
    "spatial_pattern",
    "num_legend_items",
    "legend_items_json",
    "subpanel_count",
    "comparison_type",
    "figure_quality",
]


def flatten(a: GraphAnalysis) -> dict:
    """Flatten a :class:`GraphAnalysis` into a flat dict for CSV writing.

    Nested ``Axis`` objects are expanded into prefixed scalar columns.
    ``Trend`` and ``InflectionPoint`` lists are serialised as JSON strings
    so the CSV remains rectangular.

    Parameters
    ----------
    a : GraphAnalysis
        Validated analysis result for a single graph image.

    Returns
    -------
    dict
        Flat dictionary whose keys match :data:`CSV_COLUMNS`.
    """
    def axis_fields(ax: Optional[Axis], prefix: str) -> dict:
        """Expand an Axis model into prefixed scalar fields."""
        if ax is None:
            return {f"{prefix}_label": "", f"{prefix}_units": "",
                    f"{prefix}_range_min": "", f"{prefix}_range_max": "",
                    f"{prefix}_scale_type": "",
                    f"{prefix}_tick_count": "",
                    f"{prefix}_tick_labels_json": "[]"}
        return {
            f"{prefix}_label": ax.label or "",
            f"{prefix}_units": ax.units or "",
            f"{prefix}_range_min": ax.range_min if ax.range_min is not None else "",
            f"{prefix}_range_max": ax.range_max if ax.range_max is not None else "",
            f"{prefix}_scale_type": ax.scale_type or "",
            f"{prefix}_tick_count": ax.tick_count if ax.tick_count is not None else "",
            f"{prefix}_tick_labels_json": json.dumps(ax.tick_labels),
        }

    row = {
        "file_path": a.file_path,
        "graph_title": a.graph_title or "",
        "graph_type": a.graph_type,
    }
    row.update(axis_fields(a.x_axis, "x_axis"))
    row.update(axis_fields(a.y_axis, "y_axis"))
    row.update(axis_fields(a.color_axis, "color_axis"))
    row["num_trends"] = len(a.trends)  # type: ignore
    # Serialise trend/inflection lists as JSON so downstream scripts can
    # parse them back with json.loads().
    row["trends_json"] = json.dumps([t.model_dump() for t in a.trends])
    row["num_inflection_points"] = len(a.inflection_points)  # type: ignore
    row["inflection_points_json"] = json.dumps([ip.model_dump() for ip in a.inflection_points])
    row["num_statistical_tests"] = len(a.statistical_tests)  # type: ignore
    row["statistical_tests_json"] = json.dumps([st.model_dump() for st in a.statistical_tests])
    row["num_outliers"] = len(a.outliers)  # type: ignore
    row["outliers_json"] = json.dumps([o.model_dump() for o in a.outliers])
    row["num_reference_lines"] = len(a.reference_lines)  # type: ignore
    row["reference_lines_json"] = json.dumps([rl.model_dump() for rl in a.reference_lines])
    row["num_issues"] = len(a.issues)  # type: ignore
    row["issues_json"] = json.dumps(a.issues)
    row["summary"] = a.summary
    row["sample_size"] = a.sample_size or ""
    row["data_series_count"] = a.data_series_count if a.data_series_count is not None else ""  # type: ignore
    row["error_bars"] = a.error_bars or ""
    row["num_annotations"] = len(a.annotations)  # type: ignore
    row["annotations_json"] = json.dumps(a.annotations)
    row["clinical_relevance"] = a.clinical_relevance or ""
    row["data_density"] = a.data_density or ""
    row["spatial_pattern"] = a.spatial_pattern or ""
    row["num_legend_items"] = len(a.legend_items)  # type: ignore
    row["legend_items_json"] = json.dumps(a.legend_items)
    row["subpanel_count"] = a.subpanel_count if a.subpanel_count is not None else ""  # type: ignore
    row["comparison_type"] = a.comparison_type or ""
    row["figure_quality"] = a.figure_quality or ""
    return row


# ── Main ─────────────────────────────────────────────────────────────────────

def _should_use_local(argv: list[str]) -> bool:
    """Check whether ``--local`` was passed on the command line."""
    return "--local" in argv


def _write_results_csv(
    folder: Path,
    final_results: list[GraphAnalysis],
    errors: int,
    images: list[Path],
    mode_label: str,
) -> int:
    """Write the results CSV and print a summary.

    Parameters
    ----------
    folder : Path
        Output folder.
    final_results : list[GraphAnalysis]
        Analysed results (one per image).
    errors : int
        Number of images that failed analysis.
    images : list[Path]
        Original image list (for the error count denominator).
    mode_label : str
        Label for the analysis mode (e.g. "Gemini API" or "local fallback").

    Returns
    -------
    int
        Exit code: 1 if any errors, 0 otherwise.
    """
    if errors:
        print(f"\n  {errors}/{len(images)} images failed (see 'error' rows in CSV)")

    out_csv = folder / "graph_analysis_results.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for analysis in final_results:
            writer.writerow(flatten(analysis))

    print()
    print(f"\U0001f4c1 Results saved to: {out_csv}")
    print(f"   Total graphs analyzed: {len(final_results)} ({mode_label})")

    type_counts: dict[str, int] = {}
    for r in final_results:
        type_counts[r.graph_type] = type_counts.get(r.graph_type, 0) + 1
    print(f"   Graph types: {json.dumps(type_counts, indent=2)}")

    return 1 if errors else 0


def _run_local_fallback(folder: Path, images: list[Path]) -> int:
    """Run the local filename-based fallback analyser for all images.

    Parameters
    ----------
    folder : Path
        Output folder.
    images : list[Path]
        Image files to analyse.

    Returns
    -------
    int
        Exit code (always 0 for local fallback).
    """
    print(f"\U0001f680 Found {len(images)} graph images in {folder.name}")
    print("   Mode: local fallback (filename heuristics)")
    print()

    pbar = tqdm(
        total=len(images),
        desc="Analyzing graphs (local)",
        unit="img",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}] {postfix}",
    )

    final_results: list[GraphAnalysis] = []
    for img in images:
        pbar.set_postfix_str(img.name, refresh=True)
        final_results.append(analyze_image_local(img))
        pbar.update(1)
    pbar.close()

    return _write_results_csv(folder, final_results, 0, images, "local fallback")


async def main():
    """Discover images, analyse them, and write CSV.

    This is the async entry point.  It uses a **worker pool** pattern
    instead of ``asyncio.gather`` + semaphore: a fixed number of worker
    coroutines pull images from an ``asyncio.Queue`` one at a time.

    This avoids the freeze caused by the old architecture, where all tasks
    were created upfront and backoff sleeps held semaphore slots, starving
    queued tasks and causing cascading rate-limit failures.

    **Fallback behaviour:**

    - If ``--local`` is passed, the local filename-based fallback is used
      unconditionally.
    - If ``GEMINI_API_KEY`` is not set or the ``google-genai`` package is
      not installed, the script falls back to local analysis automatically
      (with a warning) instead of exiting with an error.
    - If ``vision.fallback_to_local`` is ``false`` in the analysis config,
      missing API credentials will still cause a hard exit (preserving the
      original behaviour for CI environments that require Gemini).

    Steps:
    1. Resolves the output folder.
    2. Collects all PNG/JPG images recursively.
    3. Checks for ``--local`` flag, API key, and ``google-genai`` package.
    4. If Gemini is available: spawns async workers for API analysis.
    5. If Gemini is unavailable: runs local filename-based fallback.
    6. Writes ``graph_analysis_results.csv`` to the output folder.
    """
    force_local = _should_use_local(sys.argv)

    folder = resolve_folder(sys.argv)
    images = collect_images(folder)

    if not images:
        sys.exit(f"ERROR: No image files found in {folder}")

    # ── Determine analysis mode ──
    fallback_enabled = _vision_cfg.get("fallback_to_local", True)
    use_local = force_local

    if not force_local:
        # Check whether the google-genai package is installed.
        try:
            _ensure_genai()
        except ImportError:
            if fallback_enabled:
                print(
                    "\u26a0\ufe0f  google-genai package not installed. "
                    "Falling back to local filename-based analysis.\n"
                    "   Install with: pip install google-genai\n"
                )
                use_local = True
            else:
                sys.exit(
                    "ERROR: google-genai package not installed.\n"
                    "  Install with: pip install google-genai\n"
                    "  Or set vision.fallback_to_local=true in analysis_config.json "
                    "to allow local fallback."
                )

        if not use_local:
            api_key = os.environ.get("GEMINI_API_KEY")
            if not api_key:
                if fallback_enabled:
                    print(
                        "\u26a0\ufe0f  GEMINI_API_KEY not set. "
                        "Falling back to local filename-based analysis.\n"
                        "   For full vision analysis, set the key via:\n"
                        "   1. Running 'python run_analysis.py' (interactive prompt)\n"
                        "   2. Creating 'analysis/.env' with GEMINI_API_KEY=your-key\n"
                        "   3. export GEMINI_API_KEY='your-key'\n"
                    )
                    use_local = True
                else:
                    sys.exit(
                        "ERROR: GEMINI_API_KEY environment variable not set.\n"
                        "  You can set it by:\n"
                        "  1. Running 'python run_analysis.py' to be prompted interactively.\n"
                        "  2. Creating an 'analysis/.env' file with GEMINI_API_KEY=your-api-key\n"
                        "  3. Exporting it in your shell: export GEMINI_API_KEY='your-api-key'\n"
                        "  Or set vision.fallback_to_local=true in analysis_config.json "
                        "to allow local fallback."
                    )

    # ── Local fallback path ──
    if use_local:
        return _run_local_fallback(folder, images)

    # ── Gemini API path ──
    api_key = os.environ.get("GEMINI_API_KEY")
    print(f"\U0001f680 Found {len(images)} graph images in {folder.name}")
    print(f"   Model: {GEMINI_MODEL} (Google Gemini)")
    print(f"   Concurrency: {SEM_LIMIT}")
    print()

    client = genai.Client(api_key=api_key)  # type: ignore
    rate_limiter = _RateLimitCoordinator()

    # ── Worker pool ──
    # Instead of creating N tasks (one per image) behind a semaphore, we
    # create exactly SEM_LIMIT workers that pull from a shared queue.
    # This ensures only SEM_LIMIT images are ever in-flight, backoff
    # sleeps don't block other images from progressing, and rate-limit
    # coordination is handled globally.
    queue: asyncio.Queue[tuple[int, Path]] = asyncio.Queue()
    for idx, img in enumerate(images):
        queue.put_nowait((idx, img))

    # Pre-allocate result slots (preserves original image order).
    results: list[GraphAnalysis | None] = [None] * len(images)

    pbar = tqdm(
        total=len(images),
        desc="Analyzing graphs",
        unit="img",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}] {postfix}",
    )

    async def _worker(worker_id: int):
        """Pull images from the queue and analyse them one at a time."""
        while True:
            try:
                idx, img = queue.get_nowait()
            except asyncio.QueueEmpty:
                return  # No more work — worker exits.
            try:
                result = await analyze_image(client, img, rate_limiter, pbar)
                results[idx] = result
            except Exception as e:
                print(f"  [ERROR] {img.name}: {e}", flush=True)
                results[idx] = GraphAnalysis(
                    file_path=str(img),  # type: ignore
                    graph_type="error",  # type: ignore
                    summary=f"API error: {e}",  # type: ignore
                )
                if pbar is not None:
                    pbar.update(1)
            queue.task_done()

    # Spawn workers and wait for all images to be processed.
    workers = [asyncio.create_task(_worker(i)) for i in range(SEM_LIMIT)]
    await asyncio.gather(*workers)
    pbar.close()

    # ── Collect results ──
    final_results: list[GraphAnalysis] = []
    errors = 0
    for img, res in zip(images, results):
        if res is None:
            errors += 1
            final_results.append(GraphAnalysis(
                file_path=str(img),  # type: ignore
                graph_type="error",  # type: ignore
                summary="Worker did not produce a result",  # type: ignore
            ))
        elif res.graph_type == "error":  # type: ignore
            errors += 1
            final_results.append(res)
        else:
            final_results.append(res)

    return _write_results_csv(folder, final_results, errors, images, "Gemini API")


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    if exit_code:
        sys.exit(exit_code)  # type: ignore
