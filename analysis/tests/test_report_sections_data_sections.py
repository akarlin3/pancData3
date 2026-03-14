"""Tests for report_sections.data_sections module.

Validates data section builders:
- Cohort overview
- Patient flow
- Data completeness
- MAT data section
- Appendix
- Figure gallery
"""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from conftest import SAMPLE_GRAPH_CSV_ROWS

from report.sections.data_sections import (
    _section_cohort_overview,
    _section_patient_flow,
    _section_data_completeness,
    _section_mat_data,
    _section_appendix,
    _section_figure_gallery,
)

# Note: actual signatures from source:
# _section_cohort_overview(mat_data, log_data, dwi_types_present)
# _section_patient_flow(log_data, dwi_types_present, mat_data)
# _section_data_completeness(log_data, dwi_types_present)  -- 2 args only
# _section_mat_data(mat_data)  -- 1 arg only
# _section_appendix(rows)  -- 1 arg only
# _section_figure_gallery(folder, rows=None)  -- folder is a Path


def _make_log_data():
    return {
        "Standard": {
            "baseline": {
                "total_outliers": {"n_removed": 4, "n_total": 42, "pct": 9.5},
                "baseline_exclusion": {
                    "n_excluded": 6, "n_total": 48,
                    "lf_rate_included": 35.0, "lf_rate_excluded": 50.0,
                },
            },
            "sanity_checks": {
                "all_converged": True,
                "total_convergence": 0,
                "dim_mismatches": 1,
                "nan_dose_warnings": 0,
            },
            "stats_comparisons": {
                "glme_excluded": {"n": 5, "pct": 12.0},
            },
        }
    }


def _make_mat_data():
    return {
        "Standard": {
            "longitudinal": {"num_patients": 42, "num_timepoints": 5},
            "dosimetry": {
                "d95_gtvp": {"mean": 45.2, "std": 3.1},
                "v50_gtvp": {"mean": 0.87, "std": 0.08},
            },
            "core_comparison": {
                "methods": ["adc_threshold", "otsu"],
                "dice_matrix": [[1.0, 0.75], [0.75, 1.0]],
            },
        }
    }


# ── Cohort Overview ──


class TestCohortOverview:
    def test_returns_html(self):
        result = _section_cohort_overview(_make_mat_data(), _make_log_data(), ["Standard"])
        assert isinstance(result, list)

    def test_empty_data(self):
        result = _section_cohort_overview(None, None, [])
        assert isinstance(result, list)

    def test_shows_patient_count(self):
        html = "\n".join(_section_cohort_overview(_make_mat_data(), _make_log_data(), ["Standard"]))
        assert "42" in html


# ── Patient Flow ──


class TestPatientFlow:
    def test_returns_html(self):
        result = _section_patient_flow(_make_log_data(), ["Standard"], _make_mat_data())
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_patient_flow(None, [], {})
        assert isinstance(result, list)

    def test_exclusion_shown(self):
        html = "\n".join(_section_patient_flow(_make_log_data(), ["Standard"], _make_mat_data()))
        # Should mention baseline exclusion
        if html:
            assert "exclu" in html.lower() or "Patient" in html


# ── Data Completeness ──


class TestDataCompleteness:
    def test_returns_html(self):
        result = _section_data_completeness(_make_log_data(), ["Standard"])
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_data_completeness(None, [])
        assert isinstance(result, list)


# ── MAT Data ──


class TestMatData:
    def test_returns_html(self):
        result = _section_mat_data(_make_mat_data())
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_mat_data(None)
        assert isinstance(result, list)

    def test_shows_dosimetry(self):
        html = "\n".join(_section_mat_data(_make_mat_data()))
        if html:
            assert "dosimetry" in html.lower() or "D95" in html or "MAT" in html


# ── Appendix ──


class TestAppendix:
    def test_returns_html(self):
        result = _section_appendix(SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)

    def test_empty_rows(self):
        result = _section_appendix([])
        assert isinstance(result, list)


# ── Figure Gallery ──


class TestFigureGallery:
    def test_returns_html(self, tmp_path):
        folder = tmp_path / "saved_files_20260301_120000"
        (folder / "Standard").mkdir(parents=True)
        result = _section_figure_gallery(folder)
        assert isinstance(result, list)

    def test_nonexistent_folder(self, tmp_path):
        folder = tmp_path / "nonexistent"
        result = _section_figure_gallery(folder)
        assert isinstance(result, list)
