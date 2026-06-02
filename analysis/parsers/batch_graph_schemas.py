#!/usr/bin/env python3
"""Pydantic schemas and image helpers for batch graph analysis.

This module defines the strict Pydantic models that the vision model's JSON
response is validated against, the system prompt sent to the vision model,
basic image-collection / encoding helpers, and the CSV flattening logic that
turns a :class:`GraphAnalysis` into a flat dict for serialisation.

It is imported by :mod:`parsers.batch_graph_analysis` (which re-exports every
public name for backwards compatibility) and by the sibling
``batch_graph_local`` and ``batch_graph_vision`` modules.
"""

from __future__ import annotations

import base64
import json
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field  # type: ignore


# ── Pydantic schema ──────────────────────────────────────────────────────────
# These models define the expected JSON structure returned by the vision model.
# Pydantic validates the raw JSON and provides clear error messages when the
# model response deviates from the schema.

class Axis(BaseModel):
    """Description of a single axis (x, y, or colour-bar)."""
    label: str = Field(..., description="Axis label text as shown on the graph")
    units: Optional[str] = Field(None, description="Units (e.g., mm²/s, Gy, days)")
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
