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

from tqdm import tqdm

from shared import get_config, resolve_folder, setup_utf8_stdout

# Ensure emoji and special characters print correctly on Windows consoles.
setup_utf8_stdout()

# Lazy-import google.genai so that the module can be imported (and tested)
# without the google-genai package installed.  The actual imports happen
# inside ``analyze_image`` and ``main``, which are only called when running
# the vision pipeline.
genai = None  # type: ignore[assignment]
types = None  # type: ignore[assignment]

from pydantic import BaseModel, Field, ValidationError


def _ensure_genai():
    """Import google.genai lazily; raise ImportError if unavailable."""
    global genai, types
    if genai is None:
        from google import genai as _genai
        from google.genai import types as _types
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
    issues: list[str] = Field(default_factory=list, description="List of detected quality issues with the graph")
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
  "issues": ["description of issue 1", "description of issue 2"],
  "summary": "..."
}

Rules:
- Use null for any field you cannot determine.
- For parameter maps (spatial images), set graph_type to "parameter_map" and describe the spatial pattern in summary.
- For heatmaps (Dice, Hausdorff), extract axis labels as the methods being compared.
- Inflection points apply mainly to line/trajectory plots; for other types, return an empty list.
- Keep the summary concise (2-4 sentences).
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
        candidate = text + "".join(closers[:i])
        try:
            data = json.loads(candidate)
            if isinstance(data, dict):
                return data
        except json.JSONDecodeError:
            continue
    # Brute-force: strip back to last complete key-value, close object.
    # Find the last complete value (ends with , or a closing bracket).
    for trim in range(1, min(300, len(text))):
        stub = text[:-trim].rstrip().rstrip(",")
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
    image_bytes = await asyncio.to_thread(image_path.read_bytes)
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
                        types.Part.from_bytes(
                            data=image_bytes,
                            mime_type=mime,
                        ),
                        f"Analyze this graph image. File: {image_path.name}\n"
                        "Return ONLY the JSON object described in your instructions.",
                    ],
                    config=types.GenerateContentConfig(
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

    raw_text = response.text.strip()

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
                progress_bar.update(1)
            else:
                print(f"  \u274c  JSON parse error for {image_path.name}", flush=True)
            return GraphAnalysis(
                file_path=rel_path,
                graph_type="unknown",
                summary=f"JSON parse error (unrepairable). Raw response: {raw_text[:200]}",
            )

    # Inject the file path (not returned by the model).
    data["file_path"] = rel_path

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
            file_path=rel_path,
            graph_type=data.get("graph_type", "unknown"),
            summary=data.get("summary", f"Validation error: {e}"),
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
    "num_issues",
    "issues_json",
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
    row["num_issues"] = len(a.issues)
    row["issues_json"] = json.dumps(a.issues)
    row["summary"] = a.summary
    return row


# ── Main ─────────────────────────────────────────────────────────────────────

async def main():
    """Discover images, send them to the Gemini vision API, and write CSV.

    This is the async entry point.  It uses a **worker pool** pattern
    instead of ``asyncio.gather`` + semaphore: a fixed number of worker
    coroutines pull images from an ``asyncio.Queue`` one at a time.

    This avoids the freeze caused by the old architecture, where all tasks
    were created upfront and backoff sleeps held semaphore slots, starving
    queued tasks and causing cascading rate-limit failures.

    Steps:
    1. Resolves the output folder.
    2. Collects all PNG/JPG images recursively.
    3. Validates the ``GEMINI_API_KEY`` environment variable.
    4. Spawns ``SEM_LIMIT`` worker coroutines pulling from a queue.
    5. Collects results, substituting error placeholders for failures.
    6. Writes ``graph_analysis_results.csv`` to the output folder.
    """
    _ensure_genai()
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
                    file_path=str(img),
                    graph_type="error",
                    summary=f"API error: {e}",
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
                file_path=str(img),
                graph_type="error",
                summary="Worker did not produce a result",
            ))
        elif res.graph_type == "error":
            errors += 1
            final_results.append(res)
        else:
            final_results.append(res)

    if errors:
        print(f"\n  {errors}/{len(images)} images failed (see 'error' rows in CSV)")

    # ── Write output CSV ──
    out_csv = folder / "graph_analysis_results.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for analysis in final_results:
            writer.writerow(flatten(analysis))

    print()
    print(f"\U0001f4c1 Results saved to: {out_csv}")
    print(f"   Total graphs analyzed: {len(final_results)}")

    # Print a quick breakdown of graph types found.
    type_counts: dict[str, int] = {}
    for r in final_results:
        type_counts[r.graph_type] = type_counts.get(r.graph_type, 0) + 1
    print(f"   Graph types: {json.dumps(type_counts, indent=2)}")

    # Return non-zero exit code if any images failed (e.g. timeouts).
    return 1 if errors else 0


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    if exit_code:
        sys.exit(exit_code)
