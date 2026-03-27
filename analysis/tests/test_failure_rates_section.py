"""Tests for the core method failure rates report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.failure_rates import _section_failure_rates


def _make_fr_data():
    """Build synthetic failure_rates mat_data entry."""
    return {
        "method_names": [
            "adc_threshold", "d_threshold", "df_intersection",
            "otsu", "gmm", "kmeans", "region_growing", "active_contours",
            "percentile", "spectral", "fdm",
        ],
        "pipeline_names": ["Standard", "DnCNN", "IVIMNet"],
        "fallback_rate": [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.05, 0.05, 0.08],
            [0.10, 0.12, 0.15],
            [0.08, 0.10, 0.12],
            [0.30, 0.35, 0.40],
            [0.25, 0.28, 0.32],
            [0.0, 0.0, 0.0],
            [0.15, 0.18, 0.20],
            [0.50, 0.50, 0.50],
        ],
        "empty_rate": [
            [0.02, 0.02, 0.03],
            [0.03, 0.03, 0.04],
            [0.05, 0.05, 0.06],
            [0.02, 0.02, 0.03],
            [0.03, 0.04, 0.05],
            [0.02, 0.03, 0.04],
            [0.08, 0.10, 0.12],
            [0.06, 0.08, 0.10],
            [0.01, 0.01, 0.02],
            [0.04, 0.05, 0.06],
            [0.0, 0.0, 0.0],
        ],
        "insufficient_rate": [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.01, 0.01, 0.02],
            [0.0, 0.0, 0.0],
            [0.02, 0.02, 0.03],
            [0.01, 0.01, 0.02],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.03, 0.04, 0.05],
            [0.0, 0.0, 0.0],
        ],
        "all_nan_rate": [
            [0.04, 0.04, 0.04],
            [0.04, 0.04, 0.04],
            [0.04, 0.04, 0.04],
            [0.04, 0.04, 0.04],
            [0.04, 0.04, 0.04],
            [0.04, 0.04, 0.04],
            [0.04, 0.04, 0.04],
            [0.04, 0.04, 0.04],
            [0.04, 0.04, 0.04],
            [0.04, 0.04, 0.04],
            [0.04, 0.04, 0.04],
        ],
        "any_failure_rate": [
            [0.06, 0.06, 0.07],
            [0.07, 0.07, 0.08],
            [0.10, 0.10, 0.12],
            [0.11, 0.11, 0.15],
            [0.19, 0.22, 0.27],
            [0.15, 0.18, 0.22],
            [0.42, 0.49, 0.56],
            [0.35, 0.40, 0.46],
            [0.05, 0.05, 0.06],
            [0.26, 0.31, 0.35],
            [0.54, 0.54, 0.54],
        ],
        "median_core_voxels": [
            [35, 32, 30],
            [30, 28, 26],
            [25, 23, 21],
            [38, 35, 33],
            [40, 37, 34],
            [36, 33, 30],
            [20, 18, 15],
            [22, 20, 17],
            [32, 30, 28],
            [28, 25, 22],
            [None, None, None],
        ],
    }


def _make_mat_data(with_fr=True):
    """Build mat_data dict with optional failure_rates."""
    data: dict = {
        "Standard": {
            "core_method": {},
            "dosimetry": {},
            "longitudinal": {},
        }
    }
    if with_fr:
        data["Standard"]["failure_rates"] = _make_fr_data()
    return data


class TestFailureRatesSection:
    """Tests for _section_failure_rates."""

    def test_with_data_returns_html(self):
        """Section should return HTML with table when data is present."""
        mat_data = _make_mat_data(with_fr=True)
        result = _section_failure_rates(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Failure Rates" in html
        assert "<table>" in html
        assert "adc_threshold" in html
        assert "Fallback" in html

    def test_high_failure_flagged(self):
        """Cells with rate > 25% should have sig-strong class."""
        mat_data = _make_mat_data(with_fr=True)
        result = _section_failure_rates(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "sig-strong" in html

    def test_moderate_failure_flagged(self):
        """Cells with rate 10-25% should have sig-moderate class."""
        mat_data = _make_mat_data(with_fr=True)
        result = _section_failure_rates(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "sig-moderate" in html

    def test_empty_data_returns_not_available(self):
        """Missing data should produce a 'not available' message."""
        mat_data = _make_mat_data(with_fr=False)
        result = _section_failure_rates(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_none_mat_data(self):
        """None mat_data should not crash."""
        result = _section_failure_rates(None, [])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_sorted_by_total_failure(self):
        """Methods should be sorted by total failure rate descending."""
        mat_data = _make_mat_data(with_fr=True)
        result = _section_failure_rates(mat_data, ["Standard"])
        html = "\n".join(result)
        # fdm has highest total failure (54%) and should appear before adc_threshold (6%)
        fdm_pos = html.find("fdm")
        adc_pos = html.find("adc_threshold")
        # fdm should appear first in the table (after the header)
        # Find within the <tbody> section
        tbody_start = html.find("<tbody>")
        fdm_in_table = html.find("fdm", tbody_start)
        adc_in_table = html.find("adc_threshold", tbody_start)
        assert fdm_in_table < adc_in_table, \
            "fdm (higher failure) should appear before adc_threshold in sorted table"

    def test_all_11_methods_in_table(self):
        """All 11 core methods should appear in the table."""
        mat_data = _make_mat_data(with_fr=True)
        result = _section_failure_rates(mat_data, ["Standard"])
        html = "\n".join(result)
        for method in _make_fr_data()["method_names"]:
            assert method in html, f"Method {method} should appear in table"

    def test_pipeline_names_shown(self):
        """All 3 pipeline names should appear as headings."""
        mat_data = _make_mat_data(with_fr=True)
        result = _section_failure_rates(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Standard" in html
        assert "DnCNN" in html
        assert "IVIMNet" in html
