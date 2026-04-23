"""Tests for the ADC threshold optimisation report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.threshold_optimization import (
    _section_threshold_optimization,
    _section_threshold_dice,
    _section_threshold_inflection,
    _section_threshold_significance,
)


def _make_opt(
    optimal_thresh=0.0009,
    *,
    inflection_thresh=None,
    inflection_curvature=None,
    vol_frac_curvature=None,
    significance_thresh=None,
    significance_pvalue=None,
    significance_pvalues=None,
    significance_n_lc=None,
    significance_n_lf=None,
):
    opt = {
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
    if inflection_thresh is not None:
        opt["inflection_thresh"] = inflection_thresh
        opt["inflection_curvature"] = inflection_curvature
        opt["vol_frac_curvature"] = vol_frac_curvature or [None, -0.01, -0.05, None]
    if significance_thresh is not None:
        opt["significance_thresh"] = significance_thresh
        opt["significance_pvalue"] = significance_pvalue
        opt["significance_pvalues"] = significance_pvalues or [None, 0.04, 0.20, 0.60]
        opt["significance_n_lc"] = significance_n_lc if significance_n_lc is not None else 8
        opt["significance_n_lf"] = significance_n_lf if significance_n_lf is not None else 4
        opt["significance_metric"] = "wilcoxon ranksum on vol_frac"
    return opt


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


class TestThresholdInflectionSection:
    """Tactic 2 — Volume inflection."""

    def test_renders_when_inflection_present(self):
        mat = _make_mat_data(inflection_thresh=0.0010, inflection_curvature=-0.05)
        html = "\n".join(_section_threshold_inflection(mat))
        assert "Volume Inflection (Tactic 2)" in html
        assert "Knee" in html
        assert "0.0010" in html

    def test_skipped_when_no_inflection(self):
        mat = _make_mat_data()  # no inflection_thresh provided
        assert _section_threshold_inflection(mat) == []

    def test_endpoint_curvature_dashed(self):
        mat = _make_mat_data(
            inflection_thresh=0.0010,
            inflection_curvature=-0.05,
            vol_frac_curvature=[None, -0.01, -0.05, None],
        )
        html = "\n".join(_section_threshold_inflection(mat))
        # Endpoints (None curvature) render as em-dash
        assert "&mdash;" in html

    def test_knee_row_has_distinct_highlight(self):
        # Tactic 2 uses purple (#e1d5f5), distinct from Tactic 1's green.
        mat = _make_mat_data(inflection_thresh=0.0010, inflection_curvature=-0.05)
        html = "\n".join(_section_threshold_inflection(mat))
        assert "#e1d5f5" in html
        assert "#d4edda" not in html  # not the Dice highlight


class TestThresholdSignificanceSection:
    """Tactic 3 — Outcome significance."""

    def test_renders_when_significance_present(self):
        mat = _make_mat_data(significance_thresh=0.0009, significance_pvalue=0.04)
        html = "\n".join(_section_threshold_significance(mat))
        assert "Outcome Significance (Tactic 3)" in html
        assert "Most significant" in html
        assert "0.04" in html  # p-value

    def test_skipped_when_no_significance(self):
        mat = _make_mat_data()
        assert _section_threshold_significance(mat) == []

    def test_p_below_001_shown_as_lt(self):
        mat = _make_mat_data(
            significance_thresh=0.0009,
            significance_pvalue=0.0005,
            significance_pvalues=[None, 0.0005, 0.20, 0.60],
        )
        html = "\n".join(_section_threshold_significance(mat))
        assert "<0.001" in html

    def test_n_lc_lf_shown(self):
        mat = _make_mat_data(
            significance_thresh=0.0009,
            significance_pvalue=0.04,
            significance_n_lc=10,
            significance_n_lf=5,
        )
        html = "\n".join(_section_threshold_significance(mat))
        assert "n_LC = 10" in html
        assert "n_LF = 5" in html

    def test_significance_uses_distinct_highlight(self):
        # Tactic 3 uses orange (#fde2c4), distinct from Tactic 1 and 2.
        mat = _make_mat_data(significance_thresh=0.0009, significance_pvalue=0.04)
        html = "\n".join(_section_threshold_significance(mat))
        assert "#fde2c4" in html


class TestAllThreeTacticsRender:
    """Top-level entry point renders all three sections when data is present."""

    def test_all_three_sections_appear(self):
        mat = _make_mat_data(
            inflection_thresh=0.0010,
            inflection_curvature=-0.05,
            significance_thresh=0.0009,
            significance_pvalue=0.04,
        )
        html = "\n".join(_section_threshold_optimization(mat, ["Standard"]))
        assert "Reproducibility (Tactic 1)" in html
        assert "Volume Inflection (Tactic 2)" in html
        assert "Outcome Significance (Tactic 3)" in html

    def test_only_dice_when_others_missing(self):
        mat = _make_mat_data()
        html = "\n".join(_section_threshold_optimization(mat, ["Standard"]))
        assert "Reproducibility (Tactic 1)" in html
        assert "Volume Inflection (Tactic 2)" not in html
        assert "Outcome Significance (Tactic 3)" not in html
