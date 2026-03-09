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

**Requires** the ``GEMINI_API_KEY`` environment variable to be set.

Usage:
    python batch_graph_analysis.py                       # auto-detect folder
    python batch_graph_analysis.py /path/to/saved_files  # explicit folder
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

from shared import resolve_folder, setup_utf8_stdout

# Ensure emoji and special characters print correctly on Windows consoles.
setup_utf8_stdout()

from google import genai
from google.genai import types
from pydantic import BaseModel, Field, ValidationError


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


class Trend(BaseModel):
    """A single trend observed in the graph (one per data series)."""
    series: Optional[str] = Field(None, description="Name of the data series, if labeled")
    direction: str = Field(..., description="increasing, decreasing, flat, non-monotonic, U-shaped, etc.")
    description: str = Field(..., description="Brief description of the trend")


class InflectionPoint(BaseModel):
    """An inflection point or notable change in curvature / direction.

    Primarily relevant for line and trajectory plots; for other graph types
    the vision model is instructed to return an empty list.
    """
    approximate_x: Optional[float] = Field(None, description="Approximate x-coordinate")
    approximate_y: Optional[float] = Field(None, description="Approximate y-coordinate")
    description: str = Field(..., description="What happens at this inflection point")


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
    summary: str = Field(..., description="One-paragraph plain-English summary of the graph")


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

Your task: extract structured information about the graph.

Respond with ONLY a valid JSON object matching this exact schema (no markdown, \
no commentary, no code fences):

{
  "graph_title": "<string or null>",
  "graph_type": "<line|bar|scatter|heatmap|box|histogram|parameter_map|violin|other>",
  "x_axis": {"label": "...", "units": "...", "range_min": ..., "range_max": ...} or null,
  "y_axis": {"label": "...", "units": "...", "range_min": ..., "range_max": ...} or null,
  "color_axis": {"label": "...", "units": "...", "range_min": ..., "range_max": ...} or null,
  "trends": [
    {"series": "...", "direction": "...", "description": "..."}
  ],
  "inflection_points": [
    {"approximate_x": ..., "approximate_y": ..., "description": "..."}
  ],
  "summary": "..."
}

Rules:
- Use null for any field you cannot determine.
- For parameter maps (spatial images), set graph_type to "parameter_map" and describe the spatial pattern in summary.
- For heatmaps (Dice, Hausdorff), extract axis labels as the methods being compared.
- Inflection points apply mainly to line/trajectory plots; for other types, return an empty list.
- Keep the summary concise (2-4 sentences).
"""


# ── Gemini model configuration ──────────────────────────────────────────────

# Model to use for vision analysis.
GEMINI_MODEL = "gemini-2.5-flash"

# ── Async API call ───────────────────────────────────────────────────────────

# Concurrency limit for API requests.  Kept low (2) to stay well within
# Gemini's rate limits.
SEM_LIMIT = 2

# Number of retries on rate-limit (429) errors.  Uses exponential
# backoff: 15s, 30s, 60s, 120s.
MAX_RETRIES = 4


async def analyze_image(
    client: genai.Client,
    image_path: Path,
    semaphore: asyncio.Semaphore,
) -> GraphAnalysis:
    """Send one image to Gemini and parse the structured JSON response.

    The function handles rate-limit retries, markdown fence stripping, JSON
    parsing, and Pydantic validation.  On any non-recoverable error a
    fallback ``GraphAnalysis`` with ``graph_type="unknown"`` or ``"error"``
    is returned so that the batch can continue.

    Parameters
    ----------
    client : genai.Client
        Authenticated Google Gen AI client.
    image_path : Path
        Path to the graph image to analyse.
    semaphore : asyncio.Semaphore
        Controls concurrent API request count.

    Returns
    -------
    GraphAnalysis
        Validated structured analysis, or a fallback object on failure.
    """
    async with semaphore:
        image_bytes = image_path.read_bytes()
        mime = media_type_for(image_path)
        rel_path = str(image_path)

        print(f"  \u2699\ufe0f  Analyzing: {image_path.name} ...", flush=True)

        # Retry loop with exponential backoff for rate-limit errors.
        for attempt in range(MAX_RETRIES + 1):
            try:
                response = await client.aio.models.generate_content(
                    model=GEMINI_MODEL,
                    contents=[
                        types.Part.from_bytes(
                            data=image_bytes,
                            mime_type=mime,
                        ),
                        f"Analyze this graph image. File: {image_path.name}\n"
                        "Return ONLY the JSON object described in your instructions.",
                    ],
                    config=types.GenerateContentConfig(
                        system_instruction=SYSTEM_PROMPT,
                        max_output_tokens=2048,
                    ),
                )
                break  # Success: exit the retry loop.
            except Exception as e:
                # Check for rate-limit errors (HTTP 429 or resource exhausted).
                err_str = str(e).lower()
                is_rate_limit = "429" in err_str or "rate" in err_str or "resource" in err_str or "quota" in err_str
                if is_rate_limit and attempt < MAX_RETRIES:
                    # Exponential backoff: 15s, 30s, 60s, 120s
                    wait = 2 ** attempt * 15
                    print(f"  [RATE-LIMIT] {image_path.name}: retry {attempt+1}/{MAX_RETRIES} in {wait}s", flush=True)
                    await asyncio.sleep(wait)
                elif is_rate_limit:
                    raise  # Give up after exhausting all retries.
                else:
                    raise  # Non-rate-limit errors are re-raised immediately.

        raw_text = response.text.strip()

        # Strip markdown code fences if the model wraps its JSON response
        # (e.g. "```json\n{...}\n```").
        if raw_text.startswith("```"):
            raw_text = re.sub(r"^```(?:json)?\s*", "", raw_text)
            raw_text = re.sub(r"\s*```$", "", raw_text)

        # ── JSON parsing ──
        try:
            data = json.loads(raw_text)
        except json.JSONDecodeError as e:
            print(f"  \u274c  JSON parse error for {image_path.name}: {e}", flush=True)
            return GraphAnalysis(
                file_path=rel_path,
                graph_type="unknown",
                summary=f"JSON parse error: {e}. Raw response: {raw_text[:200]}",
            )

        # Inject the file path (not returned by the model).
        data["file_path"] = rel_path

        # ── Pydantic validation ──
        try:
            analysis = GraphAnalysis(**data)
        except ValidationError as e:
            print(f"  \u26a0\ufe0f  Validation warning for {image_path.name}: {e}", flush=True)
            # Fallback: preserve whatever fields parsed successfully.
            return GraphAnalysis(
                file_path=rel_path,
                graph_type=data.get("graph_type", "unknown"),
                summary=data.get("summary", f"Validation error: {e}"),
            )

        print(f"  \u2705  Done: {image_path.name} ({analysis.graph_type})", flush=True)
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
    "x_axis_label",
    "x_axis_units",
    "x_axis_range_min",
    "x_axis_range_max",
    "y_axis_label",
    "y_axis_units",
    "y_axis_range_min",
    "y_axis_range_max",
    "color_axis_label",
    "color_axis_units",
    "color_axis_range_min",
    "color_axis_range_max",
    "num_trends",
    "trends_json",
    "num_inflection_points",
    "inflection_points_json",
    "summary",
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
        """Expand an Axis model into four prefixed scalar fields."""
        if ax is None:
            return {f"{prefix}_label": "", f"{prefix}_units": "",
                    f"{prefix}_range_min": "", f"{prefix}_range_max": ""}
        return {
            f"{prefix}_label": ax.label or "",
            f"{prefix}_units": ax.units or "",
            f"{prefix}_range_min": ax.range_min if ax.range_min is not None else "",
            f"{prefix}_range_max": ax.range_max if ax.range_max is not None else "",
        }

    row = {
        "file_path": a.file_path,
        "graph_title": a.graph_title or "",
        "graph_type": a.graph_type,
    }
    row.update(axis_fields(a.x_axis, "x_axis"))
    row.update(axis_fields(a.y_axis, "y_axis"))
    row.update(axis_fields(a.color_axis, "color_axis"))
    row["num_trends"] = len(a.trends)
    # Serialise trend/inflection lists as JSON so downstream scripts can
    # parse them back with json.loads().
    row["trends_json"] = json.dumps([t.model_dump() for t in a.trends])
    row["num_inflection_points"] = len(a.inflection_points)
    row["inflection_points_json"] = json.dumps([ip.model_dump() for ip in a.inflection_points])
    row["summary"] = a.summary
    return row


# ── Main ─────────────────────────────────────────────────────────────────────

async def main():
    """Discover images, send them to the Gemini vision API, and write CSV.

    This is the async entry point.  It:
    1. Resolves the output folder.
    2. Collects all PNG/JPG images recursively.
    3. Validates the ``GEMINI_API_KEY`` environment variable.
    4. Launches concurrent analysis tasks (bounded by :data:`SEM_LIMIT`).
    5. Collects results, substituting error placeholders for failures.
    6. Writes ``graph_analysis_results.csv`` to the output folder.
    """
    folder = resolve_folder(sys.argv)
    images = collect_images(folder)

    if not images:
        sys.exit(f"ERROR: No image files found in {folder}")

    print(f"\U0001f680 Found {len(images)} graph images in {folder.name}")
    print(f"   Model: {GEMINI_MODEL} (Google Gemini)")
    print(f"   Concurrency: {SEM_LIMIT}")
    print()

    # ── API key validation ──
    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        sys.exit(
            "ERROR: GEMINI_API_KEY environment variable not set.\n"
            "  Set it with: export GEMINI_API_KEY='your-api-key'"
        )

    client = genai.Client(api_key=api_key)
    semaphore = asyncio.Semaphore(SEM_LIMIT)

    # Launch all analysis tasks concurrently.
    # ``return_exceptions=True`` ensures one failed image does not cancel the
    # entire batch -- failed tasks return their exception object instead.
    tasks = [analyze_image(client, img, semaphore) for img in images]
    raw_results = await asyncio.gather(*tasks, return_exceptions=True)

    # ── Collect results, replacing exceptions with error placeholders ──
    results: list[GraphAnalysis] = []
    errors = 0
    for img, res in zip(images, raw_results):
        if isinstance(res, Exception):
            errors += 1
            print(f"  [ERROR] {img.name}: {res}", flush=True)
            results.append(GraphAnalysis(
                file_path=str(img),
                graph_type="error",
                summary=f"API error: {res}",
            ))
        else:
            results.append(res)

    if errors:
        print(f"\n  {errors}/{len(images)} images failed (see 'error' rows in CSV)")

    # ── Write output CSV ──
    out_csv = folder / "graph_analysis_results.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for analysis in results:
            writer.writerow(flatten(analysis))

    print()
    print(f"\U0001f4c1 Results saved to: {out_csv}")
    print(f"   Total graphs analyzed: {len(results)}")

    # Print a quick breakdown of graph types found.
    type_counts: dict[str, int] = {}
    for r in results:
        type_counts[r.graph_type] = type_counts.get(r.graph_type, 0) + 1
    print(f"   Graph types: {json.dumps(type_counts, indent=2)}")


if __name__ == "__main__":
    asyncio.run(main())
