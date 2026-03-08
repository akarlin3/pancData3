#!/usr/bin/env python3
"""
batch_graph_analysis.py — Async batch graph analysis using Claude 4.6 Opus vision.

Scans the most recent saved_files_* timestamped folder for all PNG images,
sends each to the Claude 4.6 Opus vision model to extract:
  - axis labels, units, and ranges
  - observed trends
  - inflection points (with approximate coordinates)

Results are validated via a strict Pydantic schema and saved to CSV.

Usage:
    python batch_graph_analysis.py
"""

from __future__ import annotations

import asyncio
import base64
import csv
import glob
import io
import json
import os
import re
import sys
from pathlib import Path
from typing import Optional

# Force UTF-8 stdout on Windows to avoid cp1252 emoji crashes
if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import anthropic
from pydantic import BaseModel, Field, ValidationError


# ── Pydantic schema ──────────────────────────────────────────────────────────

class Axis(BaseModel):
    """Description of a single axis."""
    label: str = Field(..., description="Axis label text as shown on the graph")
    units: Optional[str] = Field(None, description="Units (e.g., mm²/s, Gy, days)")
    range_min: Optional[float] = Field(None, description="Minimum value on this axis")
    range_max: Optional[float] = Field(None, description="Maximum value on this axis")


class Trend(BaseModel):
    """A trend observed in the graph."""
    series: Optional[str] = Field(None, description="Name of the data series, if labeled")
    direction: str = Field(..., description="increasing, decreasing, flat, non-monotonic, U-shaped, etc.")
    description: str = Field(..., description="Brief description of the trend")


class InflectionPoint(BaseModel):
    """An inflection point or notable change in curvature / direction."""
    approximate_x: Optional[float] = Field(None, description="Approximate x-coordinate")
    approximate_y: Optional[float] = Field(None, description="Approximate y-coordinate")
    description: str = Field(..., description="What happens at this inflection point")


class GraphAnalysis(BaseModel):
    """Full structured analysis of a single graph image."""
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

def find_latest_saved_folder(base_dir: str) -> Path:
    """Return the most recently modified saved_files_* directory."""
    pattern = os.path.join(base_dir, "saved_files_*")
    dirs = sorted(glob.glob(pattern), reverse=True)
    if not dirs:
        sys.exit("ERROR: No saved_files_* folders found in " + base_dir)
    # Pick the lexicographically last (most recent timestamp)
    return Path(dirs[0])


def collect_images(folder: Path) -> list[Path]:
    """Recursively collect all PNG/JPG images."""
    exts = {".png", ".jpg", ".jpeg"}
    images = []
    for f in sorted(folder.rglob("*")):
        if f.suffix.lower() in exts and f.is_file():
            images.append(f)
    return images


def image_to_base64(path: Path) -> str:
    """Read an image file and return its base64 encoding."""
    return base64.standard_b64encode(path.read_bytes()).decode("utf-8")


def media_type_for(path: Path) -> str:
    """Return MIME type for the image."""
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


# ── Async API call ───────────────────────────────────────────────────────────

SEM_LIMIT = 2  # max concurrent API requests (low to respect rate limits)
MAX_RETRIES = 4  # retry on rate-limit errors with exponential backoff


async def analyze_image(
    client: anthropic.AsyncAnthropic,
    image_path: Path,
    semaphore: asyncio.Semaphore,
) -> GraphAnalysis:
    """Send one image to Claude Sonnet 4 and parse the structured response."""
    async with semaphore:
        b64 = image_to_base64(image_path)
        mime = media_type_for(image_path)
        rel_path = str(image_path)

        print(f"  ⚙️  Analyzing: {image_path.name} ...", flush=True)

        # Retry loop for rate-limit errors
        for attempt in range(MAX_RETRIES + 1):
            try:
                response = await client.messages.create(
                    model="claude-sonnet-4-20250514",
                    max_tokens=2048,
                    system=SYSTEM_PROMPT,
                    messages=[
                        {
                            "role": "user",
                            "content": [
                                {
                                    "type": "image",
                                    "source": {
                                        "type": "base64",
                                        "media_type": mime,
                                        "data": b64,
                                    },
                                },
                                {
                                    "type": "text",
                                    "text": (
                                        f"Analyze this graph image. File: {image_path.name}\n"
                                        "Return ONLY the JSON object described in your instructions."
                                    ),
                                },
                            ],
                        }
                    ],
                )
                break  # success
            except anthropic.RateLimitError:
                if attempt < MAX_RETRIES:
                    wait = 2 ** attempt * 15  # 15, 30, 60, 120 seconds
                    print(f"  [RATE-LIMIT] {image_path.name}: retry {attempt+1}/{MAX_RETRIES} in {wait}s", flush=True)
                    await asyncio.sleep(wait)
                else:
                    raise  # give up after all retries

        raw_text = response.content[0].text.strip()

        # Strip markdown fences if the model wraps its response
        if raw_text.startswith("```"):
            raw_text = re.sub(r"^```(?:json)?\s*", "", raw_text)
            raw_text = re.sub(r"\s*```$", "", raw_text)

        try:
            data = json.loads(raw_text)
        except json.JSONDecodeError as e:
            print(f"  ❌  JSON parse error for {image_path.name}: {e}", flush=True)
            return GraphAnalysis(
                file_path=rel_path,
                graph_type="unknown",
                summary=f"JSON parse error: {e}. Raw response: {raw_text[:200]}",
            )

        data["file_path"] = rel_path

        try:
            analysis = GraphAnalysis(**data)
        except ValidationError as e:
            print(f"  ⚠️  Validation warning for {image_path.name}: {e}", flush=True)
            # Fallback: store what we can
            return GraphAnalysis(
                file_path=rel_path,
                graph_type=data.get("graph_type", "unknown"),
                summary=data.get("summary", f"Validation error: {e}"),
            )

        print(f"  ✅  Done: {image_path.name} ({analysis.graph_type})", flush=True)
        return analysis


# ── Flatten to CSV rows ──────────────────────────────────────────────────────

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
    """Flatten a GraphAnalysis into a dict suitable for CSV writing."""
    def axis_fields(ax: Optional[Axis], prefix: str) -> dict:
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
    row["trends_json"] = json.dumps([t.model_dump() for t in a.trends])
    row["num_inflection_points"] = len(a.inflection_points)
    row["inflection_points_json"] = json.dumps([ip.model_dump() for ip in a.inflection_points])
    row["summary"] = a.summary
    return row


# ── Main ─────────────────────────────────────────────────────────────────────

async def main():
    base_dir = str(Path(__file__).resolve().parent.parent)
    folder = find_latest_saved_folder(base_dir)
    images = collect_images(folder)

    if not images:
        sys.exit(f"ERROR: No image files found in {folder}")

    print(f"🚀 Found {len(images)} graph images in {folder.name}")
    print(f"   Model: claude-sonnet-4-20250514 (Claude Sonnet 4)")
    print(f"   Concurrency: {SEM_LIMIT}")
    print()

    # Ensure ANTHROPIC_API_KEY is set
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        sys.exit(
            "ERROR: ANTHROPIC_API_KEY environment variable not set.\n"
            "  Set it with: export ANTHROPIC_API_KEY='sk-ant-...'"
        )

    client = anthropic.AsyncAnthropic(api_key=api_key)
    semaphore = asyncio.Semaphore(SEM_LIMIT)

    # Launch all tasks — return_exceptions so one failure doesn't kill the batch
    tasks = [analyze_image(client, img, semaphore) for img in images]
    raw_results = await asyncio.gather(*tasks, return_exceptions=True)

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

    # Write CSV
    out_csv = folder / "graph_analysis_results.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for analysis in results:
            writer.writerow(flatten(analysis))

    print()
    print(f"📁 Results saved to: {out_csv}")
    print(f"   Total graphs analyzed: {len(results)}")

    # Quick summary
    types = {}
    for r in results:
        types[r.graph_type] = types.get(r.graph_type, 0) + 1
    print(f"   Graph types: {json.dumps(types, indent=2)}")


if __name__ == "__main__":
    asyncio.run(main())
