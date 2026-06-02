"""Shared imports, constants, and helpers for the test_xref_unit* modules.

This is a NON-test helper module (it does not start with ``test_``) so pytest
does not collect it directly. The behaviour-preserving split of the original
``test_xref_unit.py`` into ``test_xref_unit.py`` and
``test_xref_unit_part2.py`` imports everything here via ``from
_xref_unit_shared import *``.
"""

from __future__ import annotations

import csv
import io
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import pytest

from conftest import GRAPH_CSV_COLUMNS, SAMPLE_GRAPH_CSV_ROWS
from shared import (
    extract_correlations,
    extract_pvalues,
    group_by_graph_name,
    parse_dwi_info,
    safe_text,
)

# Names starting with "_" are skipped by ``from _xref_unit_shared import *``
# unless they are listed in ``__all__``. Export everything the split test
# modules rely on (the underscore-prefixed helpers, the re-exported shared
# functions, conftest constants, and stdlib names used inline by tests).
__all__ = [
    # stdlib / third-party re-exports used by tests
    "csv", "io", "json", "sys", "redirect_stdout", "Path", "patch", "pytest",
    # conftest constants
    "GRAPH_CSV_COLUMNS", "SAMPLE_GRAPH_CSV_ROWS",
    # shared.py functions
    "extract_correlations", "extract_pvalues", "group_by_graph_name",
    "parse_dwi_info", "safe_text",
    # helpers
    "_build_csv_rows", "_write_csv", "_run_summary", "_run_stat_relevance",
    "_run_by_graph_type", "_run_xref_dwi",
]


def _build_csv_rows(specs):
    """Build minimal CSV rows from (dwi_type, base_name, trends_list) specs."""
    rows = []
    for dwi_type, base_name, trends in specs:
        suffix = f"_{dwi_type}" if dwi_type in ("Standard", "dnCNN", "IVIMnet") else ""
        rows.append({
            **{col: "" for col in GRAPH_CSV_COLUMNS},
            "file_path": f"saved_files_X/{dwi_type}/{base_name}{suffix}.png",
            "graph_type": "line",
            "trends_json": json.dumps(trends),
            "summary": "",
            "statistical_tests_json": "[]",
            "inflection_points_json": "[]",
            "clinical_relevance": "",
        })
    return rows


def _write_csv(folder, rows):
    """Write rows to graph_analysis_results.csv in the given folder."""
    csv_path = folder / "graph_analysis_results.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _run_summary(folder):
    """Run cross_reference_summary.main() capturing stdout."""
    import importlib
    spec = importlib.util.find_spec("cross_reference.cross_reference_summary")
    mod = importlib.util.module_from_spec(spec)
    # Patch sys.argv to point to our folder and capture output
    buf = io.StringIO()
    with patch.object(sys, "argv", ["script.py", str(folder)]), redirect_stdout(buf):
        spec.loader.exec_module(mod)
        mod.main()
    return buf.getvalue()


def _run_stat_relevance(folder):
    """Run statistical_relevance.main() capturing stdout."""
    import importlib
    spec = importlib.util.find_spec("cross_reference.statistical_relevance")
    mod = importlib.util.module_from_spec(spec)
    buf = io.StringIO()
    with patch.object(sys, "argv", ["script.py", str(folder)]), redirect_stdout(buf):
        spec.loader.exec_module(mod)
        mod.main()
    return buf.getvalue()


def _run_by_graph_type(folder):
    """Run statistical_by_graph_type.main() capturing stdout."""
    import importlib
    spec = importlib.util.find_spec("cross_reference.statistical_by_graph_type")
    mod = importlib.util.module_from_spec(spec)
    buf = io.StringIO()
    with patch.object(sys, "argv", ["script.py", str(folder)]), redirect_stdout(buf):
        spec.loader.exec_module(mod)
        mod.main()
    return buf.getvalue()


def _run_xref_dwi(folder):
    """Run cross_reference_dwi.main() capturing stdout."""
    import importlib
    spec = importlib.util.find_spec("cross_reference.cross_reference_dwi")
    mod = importlib.util.module_from_spec(spec)
    buf = io.StringIO()
    with patch.object(sys, "argv", ["script.py", str(folder)]), redirect_stdout(buf):
        spec.loader.exec_module(mod)
        mod.main()
    return buf.getvalue()
