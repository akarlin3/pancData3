"""Tests for the cross-pipeline Dice report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.cross_pipeline_dice import _section_cross_pipeline_dice


def _make_cpd_data():
    """Build synthetic cross_pipeline_dice mat_data entry."""
    return {
        "methods": [
            "adc_threshold", "d_threshold", "df_intersection",
            "otsu", "gmm", "kmeans", "region_growing", "active_contours",
            "percentile", "spectral", "fdm",
        ],
        "pair_labels": ["Std vs DnCNN", "Std vs IVIMNet", "DnCNN vs IVIMNet"],
        "mean_dice": [
            [0.85, 0.72, 0.68],
            [0.80, 0.70, 0.65],
            [0.78, 0.65, 0.60],
            [0.90, 0.80, 0.75],
            [0.88, 0.78, 0.72],
            [0.87, 0.76, 0.70],
            [0.45, 0.40, 0.35],  # low Dice — should be flagged
            [0.50, 0.42, 0.38],
            [0.82, 0.74, 0.69],
            [0.75, 0.65, 0.58],
            [None, None, None],  # fDM at Fx1 — NaN
        ],
        "std_dice": [
            [0.05, 0.08, 0.10],
            [0.06, 0.09, 0.11],
            [0.07, 0.10, 0.12],
            [0.04, 0.06, 0.08],
            [0.05, 0.07, 0.09],
            [0.05, 0.08, 0.10],
            [0.15, 0.18, 0.20],
            [0.12, 0.15, 0.18],
            [0.06, 0.09, 0.11],
            [0.08, 0.12, 0.14],
            [None, None, None],
        ],
    }


def _make_mat_data(with_cpd=True):
    """Build mat_data dict with optional cross_pipeline_dice."""
    data: dict = {
        "Standard": {
            "core_method": {},
            "dosimetry": {},
            "longitudinal": {},
        }
    }
    if with_cpd:
        data["Standard"]["cross_pipeline_dice"] = _make_cpd_data()
    return data


class TestCrossPipelineDiceSection:
    """Tests for _section_cross_pipeline_dice."""

    def test_with_data_returns_html(self):
        """Section should return HTML with table when data is present."""
        mat_data = _make_mat_data(with_cpd=True)
        result = _section_cross_pipeline_dice(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Cross-Pipeline Dice" in html
        assert "<table>" in html
        assert "adc_threshold" in html
        assert "Std vs DnCNN" in html

    def test_low_dice_flagged(self):
        """Cells with mean Dice < 0.5 should have sig-strong class."""
        mat_data = _make_mat_data(with_cpd=True)
        result = _section_cross_pipeline_dice(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "sig-strong" in html

    def test_empty_data_returns_not_available(self):
        """Missing data should produce a 'not available' message."""
        mat_data = _make_mat_data(with_cpd=False)
        result = _section_cross_pipeline_dice(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_none_mat_data(self):
        """None mat_data should not crash."""
        result = _section_cross_pipeline_dice(None, [])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_empty_methods(self):
        """Empty methods list should produce 'not available' message."""
        mat_data = {"Standard": {"cross_pipeline_dice": {"methods": []}}}
        result = _section_cross_pipeline_dice(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_nan_values_show_dash(self):
        """NaN/None dice values should display as dash."""
        mat_data = _make_mat_data(with_cpd=True)
        result = _section_cross_pipeline_dice(mat_data, ["Standard"])
        html = "\n".join(result)
        # fDM row has all None values — should show dashes
        assert "&mdash;" in html

    def test_all_11_methods_in_table(self):
        """All 11 core methods should appear in the table."""
        mat_data = _make_mat_data(with_cpd=True)
        result = _section_cross_pipeline_dice(mat_data, ["Standard"])
        html = "\n".join(result)
        for method in _make_cpd_data()["methods"]:
            assert method in html, f"Method {method} should appear in table"

    def test_std_dice_shown(self):
        """Standard deviation values should appear with +/- symbol."""
        mat_data = _make_mat_data(with_cpd=True)
        result = _section_cross_pipeline_dice(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "&plusmn;" in html
