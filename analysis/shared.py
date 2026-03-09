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

# Canonical DWI type names in the order the pipeline processes them.
# Used throughout the suite for iteration and labelling.
DWI_TYPES = ["Standard", "dnCNN", "IVIMnet"]


def setup_utf8_stdout():
    """Force UTF-8 stdout/stderr on Windows to avoid cp1252 emoji crashes.

    On Python 3.7+ the ``reconfigure`` method is available; for older
    interpreters we fall back to wrapping the underlying byte buffer.
    This is a no-op on non-Windows platforms.
    """
    if sys.platform == "win32":
        if hasattr(sys.stdout, "reconfigure"):
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
            sys.stderr.reconfigure(encoding="utf-8", errors="replace")
        elif getattr(sys.stdout, "encoding", "").lower() != "utf-8":
            # Wrap the raw binary buffer with a new TextIOWrapper.
            sys.stdout = io.TextIOWrapper(
                sys.stdout.buffer, encoding="utf-8", errors="replace",
            )
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
            if parts[i + 1] in DWI_TYPES:
                dwi_type = parts[i + 1]
            break
    # Remove DWI-type suffixes so "Longitudinal_Mean_Metrics_Standard"
    # becomes "Longitudinal_Mean_Metrics" for cross-type matching.
    for t in ["_Standard", "_dnCNN", "_IVIMnet"]:
        base_name = base_name.replace(t, "")
    return dwi_type, base_name.replace(".png", "")


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
        # Matches: p = 0.03, p<0.001, p > 1.2e-4
        r"p\s*[=<>]\s*([\d.]+(?:e[+-]?\d+)?)",
        # Matches: p-value = 0.04, p-value<0.001
        r"p-value\s*[=<>]\s*([\d.]+(?:e[+-]?\d+)?)",
    ]
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            try:
                val = float(m.group(1))
                # Extract surrounding context for downstream display.
                start = max(0, m.start() - 80)
                end = min(len(text), m.end() + 40)
                context = text[start:end].strip()
                results.append((val, context))
            except ValueError:
                pass
    return results


def extract_correlations(text: str) -> list[tuple[float, str]]:
    """Extract correlation coefficients (r, rs, r-squared) from text.

    Searches for patterns like ``r = 0.78``, ``rs = -0.45``, or
    ``r\u00b2 = 0.61`` (where ``\\xb2`` is the superscript-2 character).

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
        # Pearson r: "r = 0.78" or "r = -0.45"
        r"r\s*=\s*(-?[\d.]+)",
        # Spearman rank correlation: "rs = 0.65"
        r"rs\s*=\s*(-?[\d.]+)",
        # Coefficient of determination: "r\u00b2 = 0.61"  (\xb2 = superscript 2)
        r"r\xb2\s*=\s*([\d.]+)",
    ]
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            try:
                val = float(m.group(1))
                start = max(0, m.start() - 60)
                end = min(len(text), m.end() + 60)
                context = text[start:end].strip()
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
    with open(csv_path, encoding="utf-8") as f:
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
        groups[base_name][dwi_type] = r
    return dict(groups)
