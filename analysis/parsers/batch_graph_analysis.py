#!/usr/bin/env python3
"""Async batch graph analysis using Google Gemini and/or Anthropic Claude vision.

Scans the most recent ``saved_files_*`` timestamped folder for all PNG/JPG
images produced by the MATLAB pipeline, sends each to the configured vision
model, and extracts structured metadata:

- Axis labels, units, and ranges
- Observed trends (direction, series name, description)
- Inflection points with approximate (x, y) coordinates

Results are validated via strict Pydantic schemas and written to a flat CSV
(``graph_analysis_results.csv``) inside the output folder.  The CSV is then
consumed by downstream analysis scripts (``generate_report.py``,
``cross_reference_dwi.py``, ``statistical_relevance.py``, etc.).

**Provider selection** is controlled by the ``vision.provider`` config key
(or the ``--provider`` CLI flag):

- ``"gemini"`` — use Google Gemini (default)
- ``"claude"`` — use Anthropic Claude
- ``"both"`` — run **both** APIs on every image and write an additional
  ``graph_analysis_comparison.csv`` noting per-image differences

When the required API key is not set or the SDK package is not installed,
the script automatically falls back to a **local filename-based heuristic**
that infers graph type, axes, and comparison type from the structured
filenames produced by the MATLAB pipeline.  This fallback produces
lower-fidelity results (no visual analysis) but keeps the downstream CSV
pipeline functional.

This module is the orchestrator: it discovers images, picks a provider, runs
the worker pool, and writes the result / comparison CSVs.  The Pydantic
schemas and CSV helpers live in :mod:`parsers.batch_graph_schemas`, the local
filename fallback in :mod:`parsers.batch_graph_local`, and the vision-API
coroutines in :mod:`parsers.batch_graph_vision`.  Every public name from those
modules is re-imported here so that existing
``from parsers.batch_graph_analysis import ...`` call sites keep working.

Usage:
    python batch_graph_analysis.py                       # auto-detect folder
    python batch_graph_analysis.py /path/to/saved_files  # explicit folder
    python batch_graph_analysis.py --local /path/to/saved_files  # force local
    python batch_graph_analysis.py --provider claude      # use Claude API
    python batch_graph_analysis.py --provider both        # run both, compare
"""

from __future__ import annotations

import asyncio
import csv
import json
import sys
from pathlib import Path

from tqdm import tqdm  # type: ignore

# Ensure analysis/ root is on sys.path so 'shared' is importable when run as subprocess.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import get_api_key, get_config, resolve_folder, setup_utf8_stdout  # type: ignore

# Ensure emoji and special characters print correctly on Windows consoles.
setup_utf8_stdout()

# ── Re-exports for backwards compatibility ───────────────────────────────────
# Pydantic schemas, image helpers, the system prompt, and CSV flattening live
# in batch_graph_schemas; the local fallback lives in batch_graph_local; the
# vision-API machinery (constants, rate limiter, analyze_* coroutines) lives in
# batch_graph_vision.  These are re-imported here so that
# ``from parsers.batch_graph_analysis import ...`` continues to expose them.
from parsers.batch_graph_schemas import (  # noqa: F401
    Axis,
    Trend,
    InflectionPoint,
    StatisticalTest,
    Outlier,
    ReferenceLine,
    GraphAnalysis,
    collect_images,
    image_to_base64,
    media_type_for,
    SYSTEM_PROMPT,
    CSV_COLUMNS,
    flatten,
)
from parsers.batch_graph_local import (  # noqa: F401
    _GRAPH_TYPE_KEYWORDS,
    _AXIS_HINTS,
    _COMPARISON_HINTS,
    _infer_graph_type,
    _infer_axes,
    _infer_comparison_type,
    analyze_image_local,
)
from parsers.batch_graph_vision import (  # noqa: F401
    _ensure_genai,
    _ensure_anthropic,
    GEMINI_MODEL,
    CLAUDE_MODEL,
    VISION_PROVIDER,
    SEM_LIMIT,
    MAX_RETRIES,
    BACKOFF_BASE,
    MAX_OUTPUT_TOKENS,
    REQUEST_TIMEOUT,
    PER_IMAGE_TIMEOUT,
    _RateLimitCoordinator,
    _repair_truncated_json,
    _is_rate_limit_error,
    analyze_image,
    _parse_vision_response,
    _is_claude_rate_limit_error,
    analyze_image_claude,
)
import parsers.batch_graph_vision as _vision
from parsers.batch_graph_compare import (  # noqa: F401
    COMPARISON_COLUMNS,
    _compare_results,
    _write_comparison_csv,
)

# Shared vision config dict (used for fallback flags below).
_vision_cfg = get_config()["vision"]


# ── Main ─────────────────────────────────────────────────────────────────────

def _should_use_local(argv: list[str]) -> bool:
    """Check whether ``--local`` was passed on the command line."""
    return "--local" in argv


def _get_provider(argv: list[str]) -> str:
    """Extract the ``--provider`` value from *argv*, falling back to config."""
    for i, arg in enumerate(argv):
        if arg == "--provider" and i + 1 < len(argv):
            return argv[i + 1].lower()
    return VISION_PROVIDER


def _write_results_csv(
    folder: Path,
    final_results: list[GraphAnalysis],
    errors: int,
    images: list[Path],
    mode_label: str,
    csv_name: str = "graph_analysis_results.csv",
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
    csv_name : str
        Output CSV filename (default ``graph_analysis_results.csv``).

    Returns
    -------
    int
        Exit code: 1 if any errors, 0 otherwise.
    """
    if errors:
        print(f"\n  {errors}/{len(images)} images failed (see 'error' rows in CSV)")

    out_csv = folder / csv_name
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


# ── Provider availability checks ────────────────────────────────────────────


def _check_gemini_available(fallback_enabled: bool) -> bool:
    """Return True if Gemini API is usable, False otherwise.

    Prints warnings or exits depending on *fallback_enabled*.
    """
    try:
        _ensure_genai()
    except ImportError:
        if fallback_enabled:
            print(
                "⚠️  google-genai package not installed. "
                "Falling back to local filename-based analysis.\n"
                "   Install with: pip install google-genai\n"
            )
        else:
            sys.exit(
                "ERROR: google-genai package not installed.\n"
                "  Install with: pip install google-genai\n"
                "  Or set vision.fallback_to_local=true in analysis_config.json "
                "to allow local fallback."
            )
        return False

    if not get_api_key("gemini"):
        if fallback_enabled:
            print(
                "⚠️  Gemini API key not set. "
                "Falling back to local filename-based analysis.\n"
                "   Set vision.gemini_api_key in analysis_config.json\n"
                "   or export GEMINI_API_KEY in your shell.\n"
            )
        else:
            sys.exit(
                "ERROR: Gemini API key not set.\n"
                "  Set vision.gemini_api_key in analysis_config.json\n"
                "  or export GEMINI_API_KEY in your shell.\n"
                "  Or set vision.fallback_to_local=true in analysis_config.json "
                "to allow local fallback."
            )
        return False

    return True


def _check_claude_available(fallback_enabled: bool) -> bool:
    """Return True if Claude API is usable, False otherwise."""
    try:
        _ensure_anthropic()
    except ImportError:
        if fallback_enabled:
            print(
                "⚠️  anthropic package not installed. "
                "Falling back to local filename-based analysis.\n"
                "   Install with: pip install anthropic\n"
            )
        else:
            sys.exit(
                "ERROR: anthropic package not installed.\n"
                "  Install with: pip install anthropic\n"
                "  Or set vision.fallback_to_local=true in analysis_config.json "
                "to allow local fallback."
            )
        return False

    if not get_api_key("claude"):
        if fallback_enabled:
            print(
                "⚠️  Anthropic API key not set. "
                "Falling back to local filename-based analysis.\n"
                "   Set vision.anthropic_api_key in analysis_config.json\n"
                "   or export ANTHROPIC_API_KEY in your shell.\n"
            )
        else:
            sys.exit(
                "ERROR: Anthropic API key not set.\n"
                "  Set vision.anthropic_api_key in analysis_config.json\n"
                "  or export ANTHROPIC_API_KEY in your shell.\n"
                "  Or set vision.fallback_to_local=true in analysis_config.json "
                "to allow local fallback."
            )
        return False

    return True


# ── Worker pool helper ──────────────────────────────────────────────────────


async def _run_worker_pool(
    images: list[Path],
    analyze_fn,
    client,
    rate_limiter: _RateLimitCoordinator,
    desc: str = "Analyzing graphs",
) -> tuple[list[GraphAnalysis], int]:
    """Run a worker pool to analyse images with the given function.

    Parameters
    ----------
    images : list[Path]
        Images to analyse.
    analyze_fn : callable
        Async function with signature ``(client, path, limiter, pbar) -> GraphAnalysis``.
    client
        API client to pass to *analyze_fn*.
    rate_limiter : _RateLimitCoordinator
        Shared rate-limit coordinator.
    desc : str
        Progress bar description.

    Returns
    -------
    tuple[list[GraphAnalysis], int]
        (final_results, error_count).
    """
    queue: asyncio.Queue[tuple[int, Path]] = asyncio.Queue()
    for idx, img in enumerate(images):
        queue.put_nowait((idx, img))

    results: list[GraphAnalysis | None] = [None] * len(images)

    pbar = tqdm(
        total=len(images),
        desc=desc,
        unit="img",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}] {postfix}",
    )

    async def _worker(worker_id: int):
        while True:
            try:
                idx, img = queue.get_nowait()
            except asyncio.QueueEmpty:
                return
            try:
                result = await asyncio.wait_for(
                    analyze_fn(client, img, rate_limiter, pbar),
                    timeout=PER_IMAGE_TIMEOUT,
                )
                results[idx] = result
            except asyncio.TimeoutError:
                print(
                    f"  [TIMEOUT] {img.name}: exceeded per-image limit "
                    f"({PER_IMAGE_TIMEOUT}s) — skipping",
                    flush=True,
                )
                results[idx] = GraphAnalysis(
                    file_path=str(img),
                    graph_type="error",
                    summary=f"Per-image timeout ({PER_IMAGE_TIMEOUT}s exceeded)",
                )
                if pbar is not None:
                    pbar.update(1)
            except (KeyboardInterrupt, asyncio.CancelledError):
                raise
            except Exception as e:
                print(f"  [ERROR] {img.name}: {e}", flush=True)
                results[idx] = GraphAnalysis(
                    file_path=str(img),
                    graph_type="error",
                    summary=f"API error: {e}",
                )
                if pbar is not None:
                    pbar.update(1)
            finally:
                queue.task_done()

    workers = [asyncio.create_task(_worker(i)) for i in range(SEM_LIMIT)]
    await asyncio.gather(*workers)
    pbar.close()

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

    return final_results, errors


async def main():
    """Discover images, analyse them via configured provider(s), and write CSV.

    This is the async entry point.  It uses a **worker pool** pattern:
    a fixed number of worker coroutines pull images from an
    ``asyncio.Queue`` one at a time.

    **Provider modes:**

    - ``"gemini"`` — Gemini API only (default, original behaviour).
    - ``"claude"`` — Claude API only.
    - ``"both"`` — run both APIs sequentially, write the primary CSV from
      Gemini, a second CSV from Claude, and a comparison CSV noting
      per-image differences between the two.

    **Fallback behaviour:**

    - If ``--local`` is passed, the local filename-based fallback is used
      unconditionally.
    - If the required API key or package is missing, the script falls back
      to local analysis automatically (unless ``fallback_to_local=false``).
    """
    force_local = _should_use_local(sys.argv)
    provider = _get_provider(sys.argv)

    folder = resolve_folder(sys.argv)
    images = collect_images(folder)

    if not images:
        sys.exit(f"ERROR: No image files found in {folder}")

    # ── Determine analysis mode ──
    fallback_enabled = _vision_cfg.get("fallback_to_local", True)
    use_local = force_local

    if not force_local and provider in ("gemini", "both"):
        if not _check_gemini_available(fallback_enabled):
            if provider == "both":
                print("  Gemini unavailable; falling back to Claude-only mode.")
                provider = "claude"
            else:
                use_local = True

    if not force_local and provider in ("claude", "both"):
        if not _check_claude_available(fallback_enabled):
            if provider == "both":
                print("  Claude unavailable; falling back to Gemini-only mode.")
                provider = "gemini"
            else:
                use_local = True

    # ── Local fallback path ──
    if use_local:
        return _run_local_fallback(folder, images)

    print(f"\U0001f680 Found {len(images)} graph images in {folder.name}")
    print(f"   Provider: {provider}")
    print(f"   Concurrency: {SEM_LIMIT}")

    # ── Single-provider: Gemini ──
    if provider == "gemini":
        print(f"   Model: {GEMINI_MODEL} (Google Gemini)")
        print()
        client = _vision.genai.Client(api_key=get_api_key("gemini"))  # type: ignore
        rate_limiter = _RateLimitCoordinator()
        final_results, errors = await _run_worker_pool(
            images, analyze_image, client, rate_limiter, "Analyzing (Gemini)")
        return _write_results_csv(folder, final_results, errors, images, "Gemini API")

    # ── Single-provider: Claude ──
    if provider == "claude":
        print(f"   Model: {CLAUDE_MODEL} (Anthropic Claude)")
        print()
        client = _vision.anthropic_mod.AsyncAnthropic(api_key=get_api_key("claude"))  # type: ignore
        rate_limiter = _RateLimitCoordinator()
        final_results, errors = await _run_worker_pool(
            images, analyze_image_claude, client, rate_limiter, "Analyzing (Claude)")
        return _write_results_csv(folder, final_results, errors, images, "Claude API")

    # ── Dual-provider: both ──
    assert provider == "both"
    print(f"   Gemini model: {GEMINI_MODEL}")
    print(f"   Claude model: {CLAUDE_MODEL}")
    print()

    # Run Gemini first
    print("── Phase 1/2: Gemini API ──")
    gemini_client = _vision.genai.Client(api_key=get_api_key("gemini"))  # type: ignore
    gemini_limiter = _RateLimitCoordinator()
    gemini_results, gemini_errors = await _run_worker_pool(
        images, analyze_image, gemini_client, gemini_limiter, "Analyzing (Gemini)")

    # Run Claude second
    print("\n── Phase 2/2: Claude API ──")
    claude_client = _vision.anthropic_mod.AsyncAnthropic(api_key=get_api_key("claude"))  # type: ignore
    claude_limiter = _RateLimitCoordinator()
    claude_results, claude_errors = await _run_worker_pool(
        images, analyze_image_claude, claude_client, claude_limiter, "Analyzing (Claude)")

    # Write primary CSV (Gemini results used as the canonical output)
    exit_code = _write_results_csv(
        folder, gemini_results, gemini_errors, images, "Gemini API (dual-provider mode)")

    # Write Claude CSV separately
    _write_results_csv(
        folder, claude_results, claude_errors, images, "Claude API (dual-provider mode)",
        csv_name="graph_analysis_results_claude.csv")

    # Write comparison CSV
    _write_comparison_csv(folder, gemini_results, claude_results)

    return exit_code


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    if exit_code:
        sys.exit(exit_code)  # type: ignore
