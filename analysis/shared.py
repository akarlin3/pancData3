#!/usr/bin/env python3
"""Shared utilities for the pancData3 analysis suite.

This module provides common functions used by all analysis scripts:

- **Folder resolution**: locating the latest ``saved_files_*`` output directory.
- **DWI type parsing**: extracting Standard / dnCNN / IVIMnet labels from paths.
- **Statistical extraction**: regex-based extraction of p-values and correlation
  coefficients from free-text summaries produced by the vision model.
- **CSV loading**: reading the vision analysis CSV into structured dicts.
- **Grouping**: organising graph rows by normalised graph name for cross-DWI
  comparison.
- **UTF-8 stdout setup**: prevents encoding crashes on Windows consoles when
  printing emoji characters used in pipeline logging.

The stateless configuration machinery (defaults, API-key resolution, layered
config loading, JSON helpers, and the ``DWI_TYPES`` constant) lives in the
sibling module :mod:`shared_config`.  Those names are re-exported here so that
``from shared import ...`` continues to work for every existing consumer.  The
mutable config cache and its accessors (:func:`get_config` /
:func:`reset_config_cache`) remain in this module so that in-process cache
mutation has a single owner.
"""

from __future__ import annotations

import csv
import glob
import io
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

# ── Re-exported configuration API (see shared_config.py) ─────────────────────
# These names are kept importable from ``shared`` for backward compatibility.
from shared_config import (  # noqa: F401  (re-export)
    _API_KEY_MAP,
    _DEFAULTS,
    _deep_merge,
    _is_string_well_formed,
    _normalize_dwi_type,
    _permissive_json_loads,
    _strip_json_comments,
    DWI_TYPES,
    get_api_key,
    load_analysis_config,
)

# Resolved at first access via :func:`get_config`.
_config_cache: dict | None = None


def get_config() -> dict:
    """Return the cached analysis configuration (loading it on first call).

    This is the primary entry point for scripts that need config values.
    The result is cached so that repeated calls do not re-read disk.
    """
    global _config_cache
    if _config_cache is None:
        _config_cache = load_analysis_config()
    return _config_cache


def reset_config_cache() -> None:
    """Clear the cached config (useful in tests)."""
    global _config_cache
    _config_cache = None


def setup_utf8_stdout():
    """Force UTF-8 stdout/stderr on Windows to avoid cp1252 emoji crashes.

    On Python 3.7+ the ``reconfigure`` method is available; for older
    interpreters we fall back to wrapping the underlying byte buffer.
    This is a no-op on non-Windows platforms.
    """
    if sys.platform != "win32":
        return
    # Skip when running under test frameworks — modifying stdout/stderr
    # interferes with capture mechanisms and causes teardown errors.
    if "pytest" in sys.modules or "unittest" in sys.modules:
        return
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
            sys.stderr.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass  # Captured/redirected streams may reject reconfigure
    elif (getattr(sys.stdout, "encoding", "") or "").lower() != "utf-8":
            # Wrap the raw binary buffer with a new TextIOWrapper.
            # Guard: stdout/stderr may lack .buffer when captured (e.g.,
            # pytest StringIO or piped contexts).
            if hasattr(sys.stdout, "buffer"):
                sys.stdout = io.TextIOWrapper(
                    sys.stdout.buffer, encoding="utf-8", errors="replace",
                )
            if hasattr(sys.stderr, "buffer"):
                sys.stderr = io.TextIOWrapper(
                    sys.stderr.buffer, encoding="utf-8", errors="replace",
                )


def find_latest_saved_folder(base_dir: str | None = None) -> Path:
    """Return the most recently modified ``saved_files_*`` directory.

    The timestamp suffix in the folder name (``YYYYMMDD_HHMMSS``) ensures
    that lexicographic sorting picks the newest run when sorted in reverse.

    Parameters
    ----------
    base_dir : str, optional
        Directory to search in. Defaults to the repository root
        (one level above the ``analysis/`` folder).

    Returns
    -------
    Path
        Absolute path to the latest output folder.

    Raises
    ------
    SystemExit
        If no matching directories are found.
    """
    if base_dir is None:
        # Navigate from analysis/ up to the repository root.
        base_dir = str(Path(__file__).resolve().parent.parent)
    pattern = os.path.join(base_dir, "saved_files_*")
    # Reverse-sorted so the newest (lexicographically largest) comes first.
    dirs = sorted(glob.glob(pattern), reverse=True)
    if not dirs:
        sys.exit("ERROR: No saved_files_* folders found in " + base_dir)
    return Path(dirs[0])


def resolve_folder(argv: list[str] | None = None) -> Path:
    """Resolve the output folder from CLI args or auto-detect.

    If a folder path is provided as the second element of *argv*
    (i.e. the first positional CLI argument), it is validated and returned.
    Otherwise, :func:`find_latest_saved_folder` is used for auto-detection.

    Parameters
    ----------
    argv : list[str], optional
        Typically ``sys.argv``.  ``argv[1]`` is treated as the folder path.

    Returns
    -------
    Path
        Validated path to the output folder.

    Usage in scripts::

        folder = resolve_folder(sys.argv)
    """
    if argv and len(argv) > 1:
        p = Path(argv[1])
        if p.is_dir():
            return p
        sys.exit(f"ERROR: Specified folder does not exist: {p}")
    return find_latest_saved_folder()


def parse_dwi_info(file_path: str) -> tuple[str, str]:
    """Extract DWI type and normalised base graph name from a file path.

    The function looks for a ``saved_files_*`` component in the path and
    checks whether the next component is one of the known DWI types.
    DWI-type suffixes (e.g. ``_Standard``) and the ``.png`` extension are
    stripped from the filename to produce a normalised ``base_name`` that
    can be used for cross-DWI comparisons.

    Parameters
    ----------
    file_path : str
        Absolute or relative path to a graph image file.

    Returns
    -------
    dwi_type : str
        One of ``"Standard"``, ``"dnCNN"``, ``"IVIMnet"``, or ``"Root"``
        (when the file is not inside a DWI-type subfolder).
    base_name : str
        Filename without DWI suffix or ``.png`` extension.
    """
    # Normalise path separators for cross-platform compatibility.
    fp = file_path.replace("\\", "/")
    parts = fp.split("/")
    dwi_type = "Root"
    base_name = parts[-1]
    # Walk path components to find the saved_files directory.
    for i, p in enumerate(parts):
        if "saved_files" in p and i + 1 < len(parts):
            # The directory immediately after saved_files is the DWI type.
            if parts[i + 1] in DWI_TYPES:  # type: ignore
                dwi_type = parts[i + 1]  # type: ignore
            break
    # Remove DWI-type suffixes so "Longitudinal_Mean_Metrics_Standard"
    # becomes "Longitudinal_Mean_Metrics" for cross-type matching.
    for t in ["_Standard", "_dnCNN", "_IVIMnet"]:  # type: ignore
        base_name = base_name.replace(t, "")
    return dwi_type, base_name.replace(".png", "")


def safe_text(row: dict, *keys: str) -> str:
    """Safely concatenate multiple CSV row fields into one string.

    Returns a single space-joined string.  Missing or ``None``-valued
    fields are silently replaced with empty strings so that downstream
    regex extraction never crashes on incomplete vision CSV rows.
    """
    parts = []
    for k in keys:
        v = row.get(k, "") or ""
        parts.append(v)
    return " ".join(parts)


def extract_pvalues(text: str) -> list[tuple[float, str]]:
    """Extract all p-values from a text string.

    Searches for patterns like ``p = 0.032``, ``p < 0.001``, or
    ``p-value = 1.2e-4`` (case-insensitive).

    Parameters
    ----------
    text : str
        Free-text string to search (typically a graph summary or
        trend description produced by the vision model).

    Returns
    -------
    list[tuple[float, str]]
        Each tuple contains ``(p_value, context_snippet)`` where the
        context snippet is up to 80 characters before and 40 characters
        after the match to provide surrounding context.
    """
    results: list[tuple[float, str]] = []
    patterns = [
        # Matches: p-value = 0.04, p-value<0.001 (checked first, more specific)
        r"p-value\s*[=<>]\s*([\d.]+(?:e[+-]?\d+)?)",
        # Matches: p = 0.03, p<0.001, p > 1.2e-4
        # (?<![a-zA-Z]) prevents matching "up = 2.5", "group = ...", etc.
        r"(?<![a-zA-Z])p\s*[=<>]\s*([\d.]+(?:e[+-]?\d+)?)",
    ]
    seen_spans: set[tuple[int, int]] = set()
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            # Deduplicate overlapping matches (p-value pattern is a superset).
            span = (m.start(), m.end())
            if any(s[0] <= span[0] < s[1] or s[0] < span[1] <= s[1] for s in seen_spans):
                continue
            seen_spans.add(span)
            try:
                val = float(m.group(1))
                # Skip values that are clearly not p-values (> 1.0 or negative).
                if val > 1.0:
                    continue
                # Extract surrounding context for downstream display.
                start = max(0, m.start() - 80)
                end = min(len(text), m.end() + 40)
                context = text[start:end].strip()  # type: ignore
                results.append((val, context))
            except ValueError:
                pass
    return results


def extract_correlations(text: str) -> list[tuple[float, str]]:
    """Extract correlation coefficients (r, rs, r-squared) from text.

    Searches for patterns like ``r = 0.78``, ``rs = -0.45``, or
    ``r² = 0.61`` (where ``\\xb2`` is the superscript-2 character).

    Parameters
    ----------
    text : str
        Free-text string to search.

    Returns
    -------
    list[tuple[float, str]]
        Each tuple contains ``(r_value, context_snippet)`` with up to
        60 characters of surrounding context on either side.
    """
    results: list[tuple[float, str]] = []
    patterns = [
        # Spearman rank correlation: "rs = 0.65" (checked first, more specific)
        r"(?<![a-zA-Z])rs\s*=\s*(-?[\d.]+)",
        # Coefficient of determination: "r² = 0.61"  (\xb2 = superscript 2)
        r"(?<![a-zA-Z])r\xb2\s*=\s*([\d.]+)",
        # Pearson r: "r = 0.78" or "r = -0.45"
        # (?<![a-zA-Z]) prevents matching "parameter = ...", "error = ...", etc.
        r"(?<![a-zA-Z])r\s*=\s*(-?[\d.]+)",
    ]
    seen_spans: set[tuple[int, int]] = set()
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            span = (m.start(), m.end())
            if any(s[0] <= span[0] < s[1] or s[0] < span[1] <= s[1] for s in seen_spans):
                continue
            seen_spans.add(span)
            try:
                val = float(m.group(1))
                # Skip values that are clearly not correlation coefficients.
                if abs(val) > 1.0:
                    continue
                start = max(0, m.start() - 60)
                end = min(len(text), m.end() + 60)
                context = text[start:end].strip()  # type: ignore
                results.append((val, context))
            except ValueError:
                pass
    return results


def load_graph_csv(folder: Path) -> list[dict]:
    """Load ``graph_analysis_results.csv`` from a ``saved_files`` folder.

    This CSV is produced by :mod:`batch_graph_analysis` and contains one
    row per analysed graph image with columns for axes, trends, inflection
    points, and a plain-English summary.

    Parameters
    ----------
    folder : Path
        Path to the ``saved_files_*`` output folder.

    Returns
    -------
    list[dict]
        List of row dicts (one per graph).  Returns an empty list if
        the CSV does not exist (graceful degradation).
    """
    csv_path = folder / "graph_analysis_results.csv"
    if not csv_path.exists():
        return []
    rows: list[dict] = []
    with open(csv_path, encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def group_by_graph_name(rows: list[dict]) -> dict[str, dict[str, dict]]:
    """Group CSV rows by normalised graph name, then by DWI type.

    This enables cross-DWI comparison: for the same graph
    (e.g. ``"Longitudinal_Mean_Metrics"``), all three DWI-type rows are
    gathered under a single key.

    Parameters
    ----------
    rows : list[dict]
        Rows loaded from ``graph_analysis_results.csv``.

    Returns
    -------
    dict[str, dict[str, dict]]
        Mapping of ``base_name -> {dwi_type -> row_dict}``.
    """
    groups: dict[str, dict[str, dict]] = defaultdict(dict)
    for r in rows:
        dwi_type, base_name = parse_dwi_info(r["file_path"])
        groups[base_name][dwi_type] = r  # type: ignore
    return dict(groups)
