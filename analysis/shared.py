#!/usr/bin/env python3
"""Shared utilities for the pancData3 analysis suite."""

from __future__ import annotations

import csv
import glob
import io
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

# Canonical DWI type names in pipeline order
DWI_TYPES = ["Standard", "dnCNN", "IVIMnet"]


def setup_utf8_stdout():
    """Force UTF-8 stdout/stderr on Windows to avoid cp1252 emoji crashes."""
    if sys.platform == "win32":
        sys.stdout = io.TextIOWrapper(
            sys.stdout.buffer, encoding="utf-8", errors="replace"
        )
        sys.stderr = io.TextIOWrapper(
            sys.stderr.buffer, encoding="utf-8", errors="replace"
        )


def find_latest_saved_folder(base_dir: str | None = None) -> Path:
    """Return the most recently modified saved_files_* directory.

    Parameters
    ----------
    base_dir : str, optional
        Directory to search in. Defaults to the repository root
        (one level above the analysis/ folder).
    """
    if base_dir is None:
        base_dir = str(Path(__file__).resolve().parent.parent)
    pattern = os.path.join(base_dir, "saved_files_*")
    dirs = sorted(glob.glob(pattern), reverse=True)
    if not dirs:
        sys.exit("ERROR: No saved_files_* folders found in " + base_dir)
    return Path(dirs[0])


def resolve_folder(argv: list[str] | None = None) -> Path:
    """Resolve the output folder from CLI args or auto-detect.

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

    Returns
    -------
    dwi_type : str
        One of "Standard", "dnCNN", "IVIMnet", or "Root".
    base_name : str
        Filename without DWI suffix or extension.
    """
    fp = file_path.replace("\\", "/")
    parts = fp.split("/")
    dwi_type = "Root"
    base_name = parts[-1]
    for i, p in enumerate(parts):
        if "saved_files" in p and i + 1 < len(parts):
            if parts[i + 1] in DWI_TYPES:
                dwi_type = parts[i + 1]
            break
    for t in ["_Standard", "_dnCNN", "_IVIMnet"]:
        base_name = base_name.replace(t, "")
    return dwi_type, base_name.replace(".png", "")


def extract_pvalues(text: str) -> list[tuple[float, str]]:
    """Extract all p-values from a text string.

    Returns list of (p_value, context_snippet) tuples.
    """
    results: list[tuple[float, str]] = []
    patterns = [
        r"p\s*[=<>]\s*([\d.]+(?:e[+-]?\d+)?)",
        r"p-value\s*[=<>]\s*([\d.]+(?:e[+-]?\d+)?)",
    ]
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            try:
                val = float(m.group(1))
                start = max(0, m.start() - 80)
                end = min(len(text), m.end() + 40)
                context = text[start:end].strip()
                results.append((val, context))
            except ValueError:
                pass
    return results


def extract_correlations(text: str) -> list[tuple[float, str]]:
    """Extract correlation coefficients (r, rs, r²) from text.

    Returns list of (r_value, context_snippet) tuples.
    """
    results: list[tuple[float, str]] = []
    patterns = [
        r"r\s*=\s*(-?[\d.]+)",
        r"rs\s*=\s*(-?[\d.]+)",
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
    """Load graph_analysis_results.csv from a saved_files folder.

    Returns an empty list if the CSV does not exist (graceful degradation).
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
    """Group CSV rows by normalised graph name → DWI type → row dict."""
    groups: dict[str, dict[str, dict]] = defaultdict(dict)
    for r in rows:
        dwi_type, base_name = parse_dwi_info(r["file_path"])
        groups[base_name][dwi_type] = r
    return dict(groups)
