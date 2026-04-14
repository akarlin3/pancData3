"""Tests for the ADC threshold optimisation report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.threshold_optimization import _section_threshold_optimization


def _make_opt(optimal_thresh=0.0009):
    return {
        "thresholds": [0.0008, 0.0009, 0.0010, 0.0011],
        "median_dice": [0.55, 0.72, 0.68, 0.60],
        "mean_dice": [0.54, 0.70, 0.66, 0.58],
        "std_dice": [0.10, 0.08, 0.09, 0.12],
        "median_vol_frac": [0.15, 0.18, 0.22, 0.28],
        "optimal_thresh": optimal_thresh,
        "optimal_dice": 0.72,
        "optimal_vol_frac": 0.18,
        "n_patients": 12,
    }


def _make_mat_data(**kw):
    return {"Standard": {"threshold_optimization": _make_opt(**kw)}}


class TestThresholdOptimizationSection:

    def test_with_data_returns_html(self):
        result = _section_threshold_optimization(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "ADC Threshold Optimisation" in html
        assert "<table>" in html

    def test_optimal_row_highlighted(self):
        result = _section_threshold_optimization(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        # Green background + optimal marker
        assert "#d4edda" in html
        assert "Optimal" in html

    def test_comparison_box_shows_default(self):
        result = _section_threshold_optimization(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "Current default" in html
        assert "0.0010" in html

    def test_comparison_box_shows_improvement(self):
        result = _section_threshold_optimization(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "Improvement" in html

    def test_no_data_returns_empty(self):
        assert _section_threshold_optimization({}, []) == []

    def test_none_returns_empty(self):
        assert _section_threshold_optimization(None, []) == []

    def test_missing_threshold_optimization_returns_empty(self):
        mat = {"Standard": {"dosimetry": {}}}
        assert _section_threshold_optimization(mat, ["Standard"]) == []

    def test_empty_thresholds_returns_empty(self):
        mat = {"Standard": {"threshold_optimization": {"thresholds": []}}}
        assert _section_threshold_optimization(mat, ["Standard"]) == []

    def test_n_patients_shown(self):
        result = _section_threshold_optimization(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "n=12" in html

    def test_thresholds_scaled_to_e3(self):
        # 0.0008 shown as 0.80 (x10^-3)
        result = _section_threshold_optimization(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "0.80" in html or "0.8" in html
