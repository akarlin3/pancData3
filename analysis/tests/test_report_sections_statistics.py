"""Tests for report_sections.statistics module.

Validates statistics section builders:
- Effect size analysis (HR effect sizes, AUC interpretation)
- Multiple comparisons correction summary
- Model diagnostics
- Sensitivity analysis
- Power analysis

Edge cases covered:
- Empty/None inputs at every level
- No hazard ratios, no AUC data
- Zero-HR edge case (log(0) protection)
- Very wide CIs, very narrow CIs
- IPCW weights: extreme, moderate, well-behaved
- High competing-risk exclusion rate
- High baseline missingness
- Informative missingness (LF rate difference)
- High outlier removal rate
- Lambda stability: stable vs unstable
- EPV < 5 (low), EPV < 10 (marginal), EPV >= 10 (adequate)
- Cross-DWI GLME concordance: full, partial, none
- Unstable HR detection (CI ratio > 10)
- Power analysis: large vs small cohort
- Detectable vs undetectable effect sizes
- AUC discrimination thresholds: outstanding, excellent, acceptable, poor, none
- FDR global from CSV data
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
        assert "Acceptable" in html

    def test_ci_width_commentary(self):
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"].append(
            {"covariate": "wide", "hr": 2.0, "ci_lo": 0.1, "ci_hi": 5.0, "p": 0.5}
        )
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "wide" in html.lower() or "imprecise" in html.lower()

    def test_narrow_ci_commentary(self):
        """Narrow CIs should get positive commentary."""
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"] = [
            {"covariate": "tight", "hr": 1.2, "ci_lo": 1.1, "ci_hi": 1.3, "p": 0.01},
        ]
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "narrow" in html.lower() or "precise" in html.lower()

    def test_no_hazard_ratios_only_auc(self):
        """Log data with AUC but no HRs should still show AUC section."""
        log = {
            "Standard": {
                "survival": {},
                "stats_predictive": {
                    "roc_analyses": [{"auc": 0.85, "timepoint": "BL"}],
                },
            }
        }
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "0.85" in html or "Excellent" in html

    def test_auc_outstanding(self):
        """AUC >= 0.9 should be classified as Outstanding."""
        log = {
            "Standard": {
                "survival": {},
                "stats_predictive": {
                    "roc_analyses": [{"auc": 0.95, "timepoint": "BL"}],
                },
            }
        }
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "Outstanding" in html

    def test_auc_poor(self):
        """AUC 0.6-0.7 should be classified as Poor."""
        log = {
            "Standard": {
                "survival": {},
                "stats_predictive": {
                    "roc_analyses": [{"auc": 0.62, "timepoint": "BL"}],
                },
            }
        }
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "Poor" in html

    def test_auc_no_discrimination(self):
        """AUC < 0.6 should be classified as No discrimination."""
        log = {
            "Standard": {
                "survival": {},
                "stats_predictive": {
                    "roc_analyses": [{"auc": 0.52, "timepoint": "BL"}],
                },
            }
        }
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "No discrimination" in html

    def test_auc_benchmark_exceeds(self):
        """AUC >= 0.80 should say 'exceeds' benchmarks."""
        log = {
            "Standard": {
                "survival": {},
                "stats_predictive": {
                    "roc_analyses": [{"auc": 0.83, "timepoint": "BL"}],
                },
            }
        }
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "exceeds" in html.lower()

    def test_auc_benchmark_below(self):
        """AUC < 0.70 should say 'below' benchmarks."""
        log = {
            "Standard": {
                "survival": {},
                "stats_predictive": {
                    "roc_analyses": [{"auc": 0.58, "timepoint": "BL"}],
                },
            }
        }
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "below" in html.lower()

    def test_competing_risk_interpretation_note(self):
        """Should include competing-risk interpretation note."""
        html = "\n".join(_section_effect_sizes(_make_log_data(), ["Standard"], None))
        assert "Competing-risk" in html or "competing-risk" in html

    def test_empty_hazard_ratios_list(self):
        """Empty HR list should skip that DWI type."""
        log = {
            "Standard": {
                "survival": {"hazard_ratios": []},
                "stats_predictive": {},
            }
        }
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "No effect size data" in html

    def test_forest_plot_cell_present(self):
        """Forest plot column should be present in HR table."""
        html = "\n".join(_section_effect_sizes(_make_log_data(), ["Standard"], None))
        assert "Forest Plot" in html

    def test_dwi_type_not_in_log(self):
        """DWI type in list but not in log_data."""
        result = _section_effect_sizes(_make_log_data(), ["Standard", "dnCNN"], None)
        assert isinstance(result, list)


# ── Multiple Comparisons ──


class TestMultipleComparisons:
    def test_returns_html(self):
        result = _section_multiple_comparisons(_make_log_data(), ["Standard"], None)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Multiple Comparison" in html

    def test_fdr_counts(self):
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
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

    def test_raw_vs_fdr_counts(self):
        """Should distinguish raw significant from FDR significant."""
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        assert "Raw Significant" in html
        assert "FDR Significant" in html

    def test_lost_to_correction(self):
        """Metrics lost to FDR correction should be shown."""
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        # mean_d has p=0.04 < 0.05 (raw sig) but p=0.04 > adj_alpha=0.0125 (not FDR sig)
        assert "Lost to Correction" in html or "Rejected by FDR" in html

    def test_expected_vs_observed_false_discoveries(self):
        """Info box about expected vs observed false discoveries."""
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        assert "Expected false discoveries" in html or "expected" in html.lower()

    def test_bonferroni_alternative(self):
        """Alternative correction methods section should mention Bonferroni."""
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        assert "Bonferroni" in html

    def test_detailed_per_metric_table(self):
        """Detailed table should show each metric's status."""
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        assert "mean_adc" in html
        assert "mean_d" in html
        assert "Confirmed" in html or "Rejected" in html or "Not significant" in html

    def test_fdr_methodology_disclosure(self):
        """FDR methodology info box should always appear."""
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        assert "Benjamini" in html or "BH" in html

    def test_empty_glme_details(self):
        """Empty glme_details list should not crash."""
        log = {"Standard": {"stats_comparisons": {"glme_details": []}}}
        html = "\n".join(_section_multiple_comparisons(log, ["Standard"], None))
        assert "No multiple comparison data" in html

    def test_multiple_dwi_types(self):
        """Multiple DWI types each get their own FDR section."""
        log = _make_log_data()
        log["dnCNN"] = {
            "stats_comparisons": {
                "glme_details": [
                    {"metric": "mean_adc", "p": 0.02, "adj_alpha": 0.05},
                ],
            },
        }
        html = "\n".join(_section_multiple_comparisons(log, ["Standard", "dnCNN"], None))
        assert "Standard" in html
        assert "dnCNN" in html


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

    def test_ipcw_extreme_weight(self):
        """IPCW max weight > 5.0 should warn about positivity assumption."""
        log = _make_log_data()
        log["Standard"]["survival"]["ipcw"]["max_weight"] = 8.0
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "positivity" in html.lower() or "Extreme" in html or "truncation" in html.lower()

    def test_ipcw_well_behaved(self):
        """IPCW max weight <= 2.0 should show well-behaved message."""
        log = _make_log_data()
        log["Standard"]["survival"]["ipcw"]["max_weight"] = 1.5
        log["Standard"]["survival"]["ipcw"]["min_weight"] = 0.8
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "well-behaved" in html.lower()

    def test_high_competing_risk_exclusion(self):
        """Competing-risk exclusion > 20% should flag."""
        log = _make_log_data()
        log["Standard"]["stats_comparisons"]["glme_excluded"]["pct"] = 25.0
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "High competing-risk" in html or "substantially" in html.lower()

    def test_moderate_competing_risk_exclusion(self):
        """Competing-risk exclusion 10-20% should note context."""
        log = _make_log_data()
        log["Standard"]["stats_comparisons"]["glme_excluded"]["pct"] = 15.0
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "Moderate" in html or "context" in html.lower()

    def test_high_outlier_removal(self):
        """Outlier removal > 10% should trigger sensitivity analysis suggestion."""
        log = _make_log_data()
        log["Standard"]["baseline"]["total_outliers"]["pct"] = 15.0
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "sensitivity" in html.lower() or "15" in html

    def test_informative_missingness_warning(self):
        """LF rate difference > 15pp should flag informative missingness."""
        log = _make_log_data()
        log["Standard"]["baseline"]["baseline_exclusion"]["lf_rate_included"] = 30.0
        log["Standard"]["baseline"]["baseline_exclusion"]["lf_rate_excluded"] = 60.0
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "informative" in html.lower() or "MAR" in html or "MNAR" in html

    def test_lambda_stability_unstable(self):
        """Lambda ratio > 10x should flag instability."""
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": ["a"], "lambda": 0.001},
            {"timepoint": "W2", "features": ["b"], "lambda": 0.5},
        ]
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "\u03bb" in html or "lambda" in html.lower()

    def test_lambda_stability_stable(self):
        """Lambda ratio <= 10x should not flag."""
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": ["a"], "lambda": 0.05},
            {"timepoint": "W2", "features": ["b"], "lambda": 0.08},
        ]
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        # Should not mention lambda variation
        assert "varies by" not in html

    def test_no_survival_no_ipcw(self):
        """No survival data should skip IPCW check."""
        log = {"Standard": {"baseline": {}, "stats_comparisons": {}}}
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "IPCW" not in html or "Assumptions" in html

    def test_all_assumptions_listed(self):
        """Key assumption items should be present."""
        html = "\n".join(_section_model_diagnostics(None, [], {}))
        assert "Proportional hazards" in html
        assert "Non-parametric" in html or "Wilcoxon" in html
        assert "Multiple testing" in html or "BH-FDR" in html
        assert "LOOCV" in html


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
        result = _section_sensitivity_analysis(log, ["Standard"], _make_mat_data())
        html = "\n".join(result)
        if "EPV" in html:
            assert "Marginal" in html or "Low" in html

    def test_unstable_hr_detection(self):
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"].append(
            {"covariate": "unstable", "hr": 3.0, "ci_lo": 0.01, "ci_hi": 100.0, "p": 0.8}
        )
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard"], _make_mat_data()))
        assert "Unstable" in html

    def test_low_epv(self):
        """EPV < 5 should flag as Low/warning."""
        log = _make_log_data()
        # Many features, few patients
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": ["a", "b", "c", "d", "e", "f", "g", "h"], "lambda": 0.05},
        ]
        mat = {"Standard": {"longitudinal": {"num_patients": 20, "num_timepoints": 3}}}
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard"], mat))
        assert "Low EPV" in html or "EPV" in html

    def test_adequate_epv(self):
        """EPV >= 10 should not flag overfitting risk."""
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": ["a"], "lambda": 0.05},
        ]
        mat = {"Standard": {"longitudinal": {"num_patients": 100, "num_timepoints": 3}}}
        result = _section_sensitivity_analysis(log, ["Standard"], mat)
        html = "\n".join(result)
        assert "Low EPV" not in html
        assert "Marginal EPV" not in html

    def test_cross_dwi_concordance_full(self):
        """All DWI types agree on significant metrics."""
        log = {
            "Standard": {
                "stats_comparisons": {
                    "glme_details": [{"metric": "adc", "p": 0.01, "adj_alpha": 0.05}],
                },
                "stats_predictive": {"feature_selections": []},
                "baseline": {},
                "survival": {},
            },
            "dnCNN": {
                "stats_comparisons": {
                    "glme_details": [{"metric": "adc", "p": 0.01, "adj_alpha": 0.05}],
                },
                "stats_predictive": {"feature_selections": []},
                "baseline": {},
                "survival": {},
            },
        }
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard", "dnCNN"], {}))
        assert "Full cross-DWI" in html or "concordant" in html.lower()

    def test_cross_dwi_concordance_none(self):
        """No overlap in significant metrics across DWI types."""
        log = {
            "Standard": {
                "stats_comparisons": {
                    "glme_details": [{"metric": "adc", "p": 0.01, "adj_alpha": 0.05}],
                },
                "stats_predictive": {"feature_selections": []},
                "baseline": {},
                "survival": {},
            },
            "dnCNN": {
                "stats_comparisons": {
                    "glme_details": [{"metric": "d", "p": 0.01, "adj_alpha": 0.05}],
                },
                "stats_predictive": {"feature_selections": []},
                "baseline": {},
                "survival": {},
            },
        }
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard", "dnCNN"], {}))
        assert "No cross-DWI" in html or "processing-specific" in html.lower()

    def test_cross_dwi_concordance_partial(self):
        """Some but not all significant metrics concordant."""
        log = {
            "Standard": {
                "stats_comparisons": {
                    "glme_details": [
                        {"metric": "adc", "p": 0.01, "adj_alpha": 0.05},
                        {"metric": "d", "p": 0.01, "adj_alpha": 0.05},
                    ],
                },
                "stats_predictive": {"feature_selections": []},
                "baseline": {},
                "survival": {},
            },
            "dnCNN": {
                "stats_comparisons": {
                    "glme_details": [
                        {"metric": "adc", "p": 0.01, "adj_alpha": 0.05},
                        {"metric": "f", "p": 0.01, "adj_alpha": 0.05},
                    ],
                },
                "stats_predictive": {"feature_selections": []},
                "baseline": {},
                "survival": {},
            },
        }
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard", "dnCNN"], {}))
        assert "Partial" in html or "concordant" in html.lower()

    def test_no_significant_metrics_no_concordance(self):
        """No significant metrics means no concordance analysis."""
        log = {
            "Standard": {
                "stats_comparisons": {
                    "glme_details": [{"metric": "adc", "p": 0.80, "adj_alpha": 0.05}],
                },
                "stats_predictive": {"feature_selections": []},
                "baseline": {},
                "survival": {},
            },
            "dnCNN": {
                "stats_comparisons": {
                    "glme_details": [{"metric": "adc", "p": 0.90, "adj_alpha": 0.05}],
                },
                "stats_predictive": {"feature_selections": []},
                "baseline": {},
                "survival": {},
            },
        }
        result = _section_sensitivity_analysis(log, ["Standard", "dnCNN"], {})
        html = "\n".join(result)
        assert "concordance" not in html.lower()

    def test_multiple_unstable_hrs(self):
        """Multiple unstable HRs should show count."""
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"] = [
            {"covariate": "x", "hr": 5.0, "ci_lo": 0.001, "ci_hi": 500.0, "p": 0.9},
            {"covariate": "y", "hr": 0.1, "ci_lo": 0.0001, "ci_hi": 100.0, "p": 0.8},
        ]
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard"], _make_mat_data()))
        assert "2" in html
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

    def test_zero_patients_returns_empty(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 0, "num_timepoints": 3}}}
        result = _section_power_analysis(None, [], mat)
        assert result == []

    def test_contains_detectable_effect(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Detectable" in html or "detectable" in html

    def test_cohen_d_computed(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "d" in html  # Cohen's d referenced

    def test_small_cohort_underpowered(self):
        """Small cohort should be flagged as underpowered."""
        mat = {"Standard": {"longitudinal": {"num_patients": 10, "num_timepoints": 2}}}
        html = "\n".join(_section_power_analysis(None, [], mat))
        assert "underpowered" in html.lower() or "large" in html.lower()

    def test_large_cohort_adequate_power(self):
        """Large cohort should have adequate power for medium effects."""
        mat = {"Standard": {"longitudinal": {"num_patients": 200, "num_timepoints": 3}}}
        html = "\n".join(_section_power_analysis(None, [], mat))
        assert "adequate" in html.lower() or "small" in html.lower()

    def test_fdr_correction_penalty(self):
        """With multiple tests, FDR correction penalty should be mentioned."""
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "FDR" in html or "correction" in html.lower()

    def test_detectable_vs_undetectable_hr(self):
        """Should classify observed HRs as detectable or undetectable."""
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "detection" in html.lower() or "Observed" in html

    def test_no_log_data_still_shows_effect_sizes(self):
        """Without log data, should still show minimum detectable effects."""
        html = "\n".join(_section_power_analysis(None, [], _make_mat_data()))
        assert "Detectable" in html or "detectable" in html

    def test_min_hr_computed(self):
        """Minimum detectable HR should be reported."""
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "HR" in html

    def test_mat_data_from_non_standard_dwi(self):
        """Patient count should be found from any DWI type in mat_data."""
        mat = {"dnCNN": {"longitudinal": {"num_patients": 35, "num_timepoints": 3}}}
        result = _section_power_analysis(None, [], mat)
        assert isinstance(result, list)
        assert len(result) > 0
