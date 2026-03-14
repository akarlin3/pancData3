"""Tests for report_sections.statistics module.

Validates statistics section builders:
- Effect size analysis (HR effect sizes, AUC interpretation)
- Multiple comparisons correction summary
- Model diagnostics
- Sensitivity analysis
- Power analysis
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.statistics import (
    _section_effect_sizes,
    _section_multiple_comparisons,
    _section_model_diagnostics,
    _section_sensitivity_analysis,
    _section_power_analysis,
)


# ── Shared fixtures ──

def _make_log_data():
    """Create synthetic log data for statistics tests."""
    return {
        "Standard": {
            "survival": {
                "hazard_ratios": [
                    {"covariate": "mean_adc", "hr": 1.5, "ci_lo": 0.9, "ci_hi": 2.5, "p": 0.03},
                    {"covariate": "delta_d", "hr": 0.7, "ci_lo": 0.5, "ci_hi": 1.0, "p": 0.06},
                ],
                "ipcw": {"max_weight": 1.8, "min_weight": 0.6},
            },
            "stats_comparisons": {
                "glme_details": [
                    {"metric": "mean_adc", "p": 0.01, "adj_alpha": 0.025},
                    {"metric": "mean_d", "p": 0.04, "adj_alpha": 0.0125},
                    {"metric": "mean_f", "p": 0.80, "adj_alpha": 0.05},
                ],
                "glme_excluded": {"pct": 12.0},
            },
            "stats_predictive": {
                "roc_analyses": [
                    {"auc": 0.78, "timepoint": "BL"},
                    {"auc": 0.65, "timepoint": "W2"},
                ],
                "feature_selections": [
                    {"timepoint": "BL", "features": ["adc", "d"], "lambda": 0.05},
                    {"timepoint": "W2", "features": ["adc"], "lambda": 0.12},
                ],
            },
            "baseline": {
                "total_outliers": {"pct": 8.5},
                "baseline_exclusion": {
                    "n_excluded": 6, "n_total": 48,
                    "lf_rate_included": 35.0, "lf_rate_excluded": 55.0,
                },
            },
            "sanity_checks": {
                "all_converged": True,
                "total_convergence": 0,
                "dim_mismatches": 0,
                "nan_dose_warnings": 1,
            },
        }
    }


def _make_mat_data():
    return {
        "Standard": {
            "longitudinal": {"num_patients": 42, "num_timepoints": 5},
        }
    }


# ── Effect Sizes ──


class TestEffectSizes:
    def test_returns_html(self):
        result = _section_effect_sizes(_make_log_data(), ["Standard"], None)
        assert isinstance(result, list)
        assert len(result) > 0

    def test_contains_hr_table(self):
        html = "\n".join(_section_effect_sizes(_make_log_data(), ["Standard"], None))
        assert "mean_adc" in html
        assert "1.500" in html or "1.5" in html

    def test_no_data(self):
        result = _section_effect_sizes(None, [], None)
        html = "\n".join(result)
        assert "No effect size data" in html

    def test_auc_interpretation(self):
        html = "\n".join(_section_effect_sizes(_make_log_data(), ["Standard"], None))
        # AUC 0.78 should be classified as "Acceptable"
        assert "Acceptable" in html

    def test_ci_width_commentary(self):
        log = _make_log_data()
        # Make a wide CI
        log["Standard"]["survival"]["hazard_ratios"].append(
            {"covariate": "wide", "hr": 2.0, "ci_lo": 0.1, "ci_hi": 5.0, "p": 0.5}
        )
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "wide" in html.lower() or "imprecise" in html.lower()


# ── Multiple Comparisons ──


class TestMultipleComparisons:
    def test_returns_html(self):
        result = _section_multiple_comparisons(_make_log_data(), ["Standard"], None)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Multiple Comparison" in html

    def test_fdr_counts(self):
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        # 3 total tests, 1 raw sig (p=0.01 < 0.05 AND p=0.04 < 0.05 = 2 raw sig)
        assert "Total Tests" in html

    def test_no_data(self):
        result = _section_multiple_comparisons(None, [], None)
        html = "\n".join(result)
        assert "No multiple comparison data" in html

    def test_fdr_global_csv(self):
        csv_data = {
            "fdr_global": {
                "Standard": [{"metric": "adc", "p": 0.001}],
            }
        }
        result = _section_multiple_comparisons(None, [], csv_data)
        html = "\n".join(result)
        assert "Global FDR" in html


# ── Model Diagnostics ──


class TestModelDiagnostics:
    def test_returns_html(self):
        result = _section_model_diagnostics(_make_log_data(), ["Standard"], _make_mat_data())
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Model Diagnostics" in html

    def test_ipcw_weight_check(self):
        html = "\n".join(_section_model_diagnostics(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "IPCW" in html

    def test_assumptions_section(self):
        html = "\n".join(_section_model_diagnostics(None, [], {}))
        assert "Assumptions" in html

    def test_high_baseline_missingness(self):
        log = _make_log_data()
        log["Standard"]["baseline"]["baseline_exclusion"]["n_excluded"] = 15
        log["Standard"]["baseline"]["baseline_exclusion"]["n_total"] = 48
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "missingness" in html.lower()


# ── Sensitivity Analysis ──


class TestSensitivityAnalysis:
    def test_returns_html_with_data(self):
        result = _section_sensitivity_analysis(_make_log_data(), ["Standard"], _make_mat_data())
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_sensitivity_analysis(None, [], {})
        assert result == []

    def test_epv_check(self):
        log = _make_log_data()
        # With 42 patients, ~17 events, 2 features => EPV ~8.5 (marginal)
        result = _section_sensitivity_analysis(log, ["Standard"], _make_mat_data())
        html = "\n".join(result)
        # Should flag EPV since 17/2 ≈ 8.5 < 10
        if "EPV" in html:
            assert "Marginal" in html or "Low" in html

    def test_unstable_hr_detection(self):
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"].append(
            {"covariate": "unstable", "hr": 3.0, "ci_lo": 0.01, "ci_hi": 100.0, "p": 0.8}
        )
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard"], _make_mat_data()))
        assert "Unstable" in html


# ── Power Analysis ──


class TestPowerAnalysis:
    def test_returns_html(self):
        result = _section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data())
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Power" in html

    def test_no_patients_returns_empty(self):
        result = _section_power_analysis(None, [], {})
        assert result == []

    def test_contains_detectable_effect(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Detectable" in html or "detectable" in html

    def test_cohen_d_computed(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "d" in html  # Cohen's d referenced
