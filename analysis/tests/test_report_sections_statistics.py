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

    def test_empty_log_data(self):
        result = _section_effect_sizes({}, ["Standard"], None)
        html = "\n".join(result)
        assert "No effect size data" in html

    def test_dwi_type_not_in_log(self):
        result = _section_effect_sizes(_make_log_data(), ["dnCNN"], None)
        html = "\n".join(result)
        assert "No effect size data" in html

    def test_no_hazard_ratios(self):
        log = {"Standard": {"survival": {"hazard_ratios": []}, "stats_predictive": {}}}
        result = _section_effect_sizes(log, ["Standard"], None)
        html = "\n".join(result)
        assert "No effect size data" in html

    def test_auc_interpretation_acceptable(self):
        html = "\n".join(_section_effect_sizes(_make_log_data(), ["Standard"], None))
        # AUC 0.78 should be classified as "Acceptable"
        assert "Acceptable" in html

    def test_auc_interpretation_poor(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = [
            {"auc": 0.62, "timepoint": "BL"},
        ]
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "Poor" in html

    def test_auc_interpretation_excellent(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = [
            {"auc": 0.85, "timepoint": "BL"},
        ]
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "Excellent" in html

    def test_auc_interpretation_outstanding(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = [
            {"auc": 0.95, "timepoint": "BL"},
        ]
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "Outstanding" in html

    def test_auc_interpretation_no_discrimination(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = [
            {"auc": 0.52, "timepoint": "BL"},
        ]
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "No discrimination" in html

    def test_auc_benchmark_exceeds(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = [
            {"auc": 0.85, "timepoint": "BL"},
        ]
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "exceeds" in html.lower()

    def test_auc_benchmark_within(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = [
            {"auc": 0.75, "timepoint": "BL"},
        ]
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "within" in html.lower()

    def test_auc_benchmark_below(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = [
            {"auc": 0.55, "timepoint": "BL"},
        ]
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "below" in html.lower()

    def test_ci_width_wide(self):
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"].append(
            {"covariate": "wide_ci_var", "hr": 2.0, "ci_lo": 0.1, "ci_hi": 5.0, "p": 0.5}
        )
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "wide" in html.lower() or "imprecise" in html.lower()

    def test_ci_width_narrow(self):
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"] = [
            {"covariate": "tight_var", "hr": 1.2, "ci_lo": 1.0, "ci_hi": 1.4, "p": 0.01},
        ]
        html = "\n".join(_section_effect_sizes(log, ["Standard"], None))
        assert "narrow" in html.lower() or "precise" in html.lower()

    def test_forest_plot_present(self):
        html = "\n".join(_section_effect_sizes(_make_log_data(), ["Standard"], None))
        assert "Forest Plot" in html

    def test_competing_risk_note(self):
        html = "\n".join(_section_effect_sizes(_make_log_data(), ["Standard"], None))
        assert "Competing-risk" in html or "competing-risk" in html

    def test_effect_size_classification(self):
        html = "\n".join(_section_effect_sizes(_make_log_data(), ["Standard"], None))
        # log(1.5) ≈ 0.405, which is small effect
        assert "Effect" in html

    def test_sorted_by_p_value(self):
        html = "\n".join(_section_effect_sizes(_make_log_data(), ["Standard"], None))
        # mean_adc (p=0.03) should come before delta_d (p=0.06)
        pos_adc = html.find("mean_adc")
        pos_delta = html.find("delta_d")
        assert pos_adc < pos_delta

    def test_multiple_dwi_types(self):
        log = _make_log_data()
        log["dnCNN"] = {
            "survival": {
                "hazard_ratios": [
                    {"covariate": "mean_adc", "hr": 1.3, "ci_lo": 0.8, "ci_hi": 2.0, "p": 0.04},
                ],
            },
            "stats_predictive": {"roc_analyses": []},
        }
        html = "\n".join(_section_effect_sizes(log, ["Standard", "dnCNN"], None))
        assert "Standard" in html
        assert "dnCNN" in html

    def test_hr_zero_value(self):
        """HR of 0 should not crash log calculation."""
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"] = [
            {"covariate": "zero_hr", "hr": 0.0, "ci_lo": 0.0, "ci_hi": 1.0, "p": 0.5},
        ]
        result = _section_effect_sizes(log, ["Standard"], None)
        assert isinstance(result, list)

    def test_missing_covariate_name(self):
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"] = [
            {"hr": 1.5, "ci_lo": 0.9, "ci_hi": 2.5, "p": 0.03},
        ]
        result = _section_effect_sizes(log, ["Standard"], None)
        assert isinstance(result, list)

    def test_auc_zero_excluded(self):
        """AUC of 0 should not appear in the table."""
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = [
            {"auc": 0, "timepoint": "BL"},
        ]
        result = _section_effect_sizes(log, ["Standard"], None)
        html = "\n".join(result)
        assert "Discriminative Performance" not in html


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

    def test_empty_log_data(self):
        result = _section_multiple_comparisons({}, ["Standard"], None)
        html = "\n".join(result)
        assert "No multiple comparison data" in html

    def test_fdr_methodology_disclosure(self):
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        assert "Benjamini" in html
        assert "FDR" in html

    def test_raw_vs_fdr_significant_counts(self):
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        # 2 raw sig (p=0.01 < 0.05 and p=0.04 < 0.05), 1 FDR sig (p=0.01 < adj_alpha=0.025)
        assert "Raw Significant" in html
        assert "FDR Significant" in html

    def test_lost_to_correction(self):
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        # 2 raw - 1 FDR = 1 lost
        assert "Lost to Correction" in html

    def test_expected_vs_observed_rejections(self):
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        assert "Expected false discoveries" in html
        assert "Observed pre-FDR significant" in html

    def test_bonferroni_alternative(self):
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        assert "Bonferroni" in html

    def test_detail_table_statuses(self):
        html = "\n".join(_section_multiple_comparisons(_make_log_data(), ["Standard"], None))
        # mean_adc (p=0.01 < adj_alpha=0.025) -> Confirmed
        assert "Confirmed" in html
        # mean_d (p=0.04 < 0.05 but > adj_alpha=0.0125) -> Rejected by FDR
        assert "Rejected by FDR" in html
        # mean_f (p=0.80 > 0.05) -> Not significant
        assert "Not significant" in html

    def test_fdr_global_csv(self):
        csv_data = {
            "fdr_global": {
                "Standard": [{"metric": "adc", "p": 0.001}],
            }
        }
        result = _section_multiple_comparisons(None, [], csv_data)
        html = "\n".join(result)
        assert "Global FDR" in html

    def test_fdr_global_csv_multiple_types(self):
        csv_data = {
            "fdr_global": {
                "Standard": [{"metric": "adc", "p": 0.001}],
                "dnCNN": [{"metric": "adc", "p": 0.002}, {"metric": "d", "p": 0.003}],
            }
        }
        result = _section_multiple_comparisons(None, [], csv_data)
        html = "\n".join(result)
        assert "3" in html  # total of 3 metrics

    def test_fdr_global_empty(self):
        csv_data = {"fdr_global": {}}
        result = _section_multiple_comparisons(None, [], csv_data)
        html = "\n".join(result)
        assert "No multiple comparison data" in html

    def test_no_glme_details(self):
        log = {"Standard": {"stats_comparisons": {"glme_details": []}}}
        result = _section_multiple_comparisons(log, ["Standard"], None)
        html = "\n".join(result)
        assert "No multiple comparison data" in html

    def test_all_significant_after_fdr(self):
        log = {"Standard": {"stats_comparisons": {"glme_details": [
            {"metric": "m1", "p": 0.001, "adj_alpha": 0.025},
            {"metric": "m2", "p": 0.005, "adj_alpha": 0.05},
        ]}}}
        html = "\n".join(_section_multiple_comparisons(log, ["Standard"], None))
        # No "Lost to Correction" card when all survive
        assert "Confirmed" in html

    def test_excess_positive_signal(self):
        log = {"Standard": {"stats_comparisons": {"glme_details": [
            {"metric": f"m{i}", "p": 0.001, "adj_alpha": 0.025}
            for i in range(20)
        ]}}}
        html = "\n".join(_section_multiple_comparisons(log, ["Standard"], None))
        assert "positive excess suggests genuine signal" in html

    def test_no_excess_above_chance(self):
        # 20 tests, 1 raw sig = expected ~1 by chance
        details = [{"metric": f"m{i}", "p": 0.9, "adj_alpha": 0.05} for i in range(19)]
        details.append({"metric": "m19", "p": 0.04, "adj_alpha": 0.0025})
        log = {"Standard": {"stats_comparisons": {"glme_details": details}}}
        html = "\n".join(_section_multiple_comparisons(log, ["Standard"], None))
        assert "no excess above chance" in html


# ── Model Diagnostics ──


class TestModelDiagnostics:
    def test_returns_html(self):
        result = _section_model_diagnostics(_make_log_data(), ["Standard"], _make_mat_data())
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Model Diagnostics" in html

    def test_ipcw_weight_well_behaved(self):
        html = "\n".join(_section_model_diagnostics(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "IPCW" in html
        assert "well-behaved" in html

    def test_ipcw_weight_moderate(self):
        log = _make_log_data()
        log["Standard"]["survival"]["ipcw"] = {"max_weight": 3.5, "min_weight": 0.5}
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "moderately dispersed" in html

    def test_ipcw_weight_extreme(self):
        log = _make_log_data()
        log["Standard"]["survival"]["ipcw"] = {"max_weight": 8.0, "min_weight": 0.2}
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "Extreme" in html or "extreme" in html or ">5.0" in html

    def test_assumptions_section_always_present(self):
        html = "\n".join(_section_model_diagnostics(None, [], {}))
        assert "Assumptions" in html

    def test_proportional_hazards_caveat(self):
        html = "\n".join(_section_model_diagnostics(None, [], {}))
        assert "Proportional hazards" in html

    def test_loocv_caveat(self):
        html = "\n".join(_section_model_diagnostics(None, [], {}))
        assert "LOOCV" in html

    def test_missing_data_imputation_note(self):
        html = "\n".join(_section_model_diagnostics(None, [], {}))
        assert "imputation" in html.lower()

    def test_high_baseline_missingness(self):
        log = _make_log_data()
        log["Standard"]["baseline"]["baseline_exclusion"]["n_excluded"] = 15
        log["Standard"]["baseline"]["baseline_exclusion"]["n_total"] = 48
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "missingness" in html.lower()
        assert "15/48" in html

    def test_low_baseline_missingness_no_warning(self):
        log = _make_log_data()
        log["Standard"]["baseline"]["baseline_exclusion"]["n_excluded"] = 2
        log["Standard"]["baseline"]["baseline_exclusion"]["n_total"] = 48
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        # Low exclusion rate should not trigger the "High baseline missingness" warning
        assert "High baseline missingness" not in html

    def test_lf_rate_differential(self):
        log = _make_log_data()
        # Already has 35 vs 55 = 20pp diff > 15pp threshold
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "LF rate differs" in html

    def test_lf_rate_similar_no_warning(self):
        log = _make_log_data()
        log["Standard"]["baseline"]["baseline_exclusion"]["lf_rate_included"] = 35.0
        log["Standard"]["baseline"]["baseline_exclusion"]["lf_rate_excluded"] = 40.0
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "LF rate differs" not in html

    def test_high_outlier_removal(self):
        log = _make_log_data()
        log["Standard"]["baseline"]["total_outliers"]["pct"] = 15.0
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "Outlier removal rate" in html

    def test_normal_outlier_rate_no_warning(self):
        log = _make_log_data()
        log["Standard"]["baseline"]["total_outliers"]["pct"] = 5.0
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "Outlier removal rate" not in html

    def test_lambda_instability(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": ["adc"], "lambda": 0.001},
            {"timepoint": "W2", "features": ["adc"], "lambda": 0.5},
        ]
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "λ varies" in html or "lambda" in html.lower()

    def test_lambda_stable_no_warning(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": ["adc"], "lambda": 0.05},
            {"timepoint": "W2", "features": ["adc"], "lambda": 0.08},
        ]
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "varies by" not in html

    def test_competing_risk_high_exclusion(self):
        log = _make_log_data()
        log["Standard"]["stats_comparisons"]["glme_excluded"]["pct"] = 25.0
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "High competing-risk exclusion" in html

    def test_competing_risk_moderate_exclusion(self):
        html = "\n".join(_section_model_diagnostics(_make_log_data(), ["Standard"], _make_mat_data()))
        # Default has pct=12, which is > 10 => moderate
        assert "Moderate competing-risk exclusion" in html

    def test_no_ipcw_data(self):
        log = _make_log_data()
        del log["Standard"]["survival"]["ipcw"]
        result = _section_model_diagnostics(log, ["Standard"], _make_mat_data())
        html = "\n".join(result)
        # Should still produce assumptions section
        assert "Assumptions" in html

    def test_warning_boxes_for_severe_issues(self):
        log = _make_log_data()
        log["Standard"]["baseline"]["total_outliers"]["pct"] = 15.0
        log["Standard"]["survival"]["ipcw"] = {"max_weight": 8.0, "min_weight": 0.2}
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "warn-box" in html

    def test_single_feature_selection_no_lambda_check(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": ["adc"], "lambda": 0.05},
        ]
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "varies by" not in html

    def test_elastic_net_caveat(self):
        html = "\n".join(_section_model_diagnostics(None, [], {}))
        assert "Elastic net" in html or "elastic net" in html


# ── Sensitivity Analysis ──


class TestSensitivityAnalysis:
    def test_returns_html_with_data(self):
        result = _section_sensitivity_analysis(_make_log_data(), ["Standard"], _make_mat_data())
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_sensitivity_analysis(None, [], {})
        assert result == []

    def test_epv_check_marginal(self):
        log = _make_log_data()
        # With 42 patients, n_eff = 48 - 6 = 42, events ≈ 42*0.4 = 16, 2 features => EPV = 8
        result = _section_sensitivity_analysis(log, ["Standard"], _make_mat_data())
        html = "\n".join(result)
        if "EPV" in html:
            assert "Marginal" in html or "Low" in html

    def test_epv_check_low(self):
        log = _make_log_data()
        # Many features, few patients
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": [f"f{i}" for i in range(10)], "lambda": 0.05},
        ]
        result = _section_sensitivity_analysis(log, ["Standard"], _make_mat_data())
        html = "\n".join(result)
        assert "Low EPV" in html

    def test_epv_no_features(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": [], "lambda": 0.05},
        ]
        result = _section_sensitivity_analysis(log, ["Standard"], _make_mat_data())
        html = "\n".join(result)
        # No EPV warning with 0 features
        assert "EPV" not in html or "Sensitivity" not in html

    def test_unstable_hr_detection(self):
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"].append(
            {"covariate": "unstable", "hr": 3.0, "ci_lo": 0.01, "ci_hi": 100.0, "p": 0.8}
        )
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard"], _make_mat_data()))
        assert "Unstable" in html

    def test_stable_hr_no_warning(self):
        log = _make_log_data()
        # Default HRs have CI ratio < 10
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard"], _make_mat_data()))
        assert "Unstable HR" not in html

    def test_cross_dwi_concordance_none(self):
        log = _make_log_data()
        log["dnCNN"] = {
            "stats_comparisons": {
                "glme_details": [
                    {"metric": "mean_f", "p": 0.01, "adj_alpha": 0.025},  # different metric sig
                ],
            },
            "stats_predictive": {"feature_selections": []},
            "baseline": {},
            "survival": {},
        }
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard", "dnCNN"], _make_mat_data()))
        assert "No cross-DWI GLME concordance" in html

    def test_cross_dwi_concordance_partial(self):
        log = _make_log_data()
        log["dnCNN"] = {
            "stats_comparisons": {
                "glme_details": [
                    {"metric": "mean_adc", "p": 0.01, "adj_alpha": 0.025},
                    {"metric": "extra_metric", "p": 0.005, "adj_alpha": 0.025},
                ],
            },
            "stats_predictive": {"feature_selections": []},
            "baseline": {},
            "survival": {},
        }
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard", "dnCNN"], _make_mat_data()))
        assert "Partial cross-DWI GLME concordance" in html

    def test_cross_dwi_concordance_full(self):
        log = _make_log_data()
        # Standard has mean_adc sig. Make dnCNN also only have mean_adc sig.
        log["dnCNN"] = {
            "stats_comparisons": {
                "glme_details": [
                    {"metric": "mean_adc", "p": 0.01, "adj_alpha": 0.025},
                ],
            },
            "stats_predictive": {"feature_selections": []},
            "baseline": {},
            "survival": {},
        }
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard", "dnCNN"], _make_mat_data()))
        assert "Full cross-DWI GLME concordance" in html

    def test_single_dwi_no_concordance_check(self):
        log = _make_log_data()
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard"], _make_mat_data()))
        assert "cross-DWI GLME concordance" not in html

    def test_robustness_preamble(self):
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"].append(
            {"covariate": "unstable", "hr": 3.0, "ci_lo": 0.01, "ci_hi": 100.0, "p": 0.8}
        )
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard"], _make_mat_data()))
        assert "Robustness" in html or "reliability" in html

    def test_empty_items_returns_empty(self):
        """When no issues found, section returns empty."""
        log = {"Standard": {
            "stats_comparisons": {"glme_details": []},
            "stats_predictive": {"feature_selections": []},
            "baseline": {},
            "survival": {"hazard_ratios": []},
        }}
        result = _section_sensitivity_analysis(log, ["Standard"], _make_mat_data())
        # Should be empty when no issues found
        assert result == [] or "Sensitivity" not in "\n".join(result)

    def test_multiple_unstable_hrs_truncated(self):
        log = _make_log_data()
        for i in range(5):
            log["Standard"]["survival"]["hazard_ratios"].append(
                {"covariate": f"unstable_{i}", "hr": 5.0, "ci_lo": 0.001, "ci_hi": 500.0, "p": 0.9}
            )
        html = "\n".join(_section_sensitivity_analysis(log, ["Standard"], _make_mat_data()))
        # Should show count but only list first 3 names
        assert "5 covariate(s)" in html


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

    def test_no_mat_data_returns_empty(self):
        result = _section_power_analysis(_make_log_data(), ["Standard"], {})
        assert result == []

    def test_contains_detectable_effect(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Detectable" in html or "detectable" in html

    def test_cohen_d_computed(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "d" in html

    def test_min_detectable_hr(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "HR" in html

    def test_underpowered_small_cohort(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 15, "num_timepoints": 3}}}
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], mat))
        # With n=15, n_per_group=7, min_d = 2.8/sqrt(7) ≈ 1.06 > 0.5 → underpowered
        assert "underpowered" in html.lower()

    def test_adequate_power_large_cohort(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 100, "num_timepoints": 5}}}
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], mat))
        # With n=100, n_per_group=50, min_d = 2.8/sqrt(50) ≈ 0.396 < 0.5 → adequate
        assert "adequate" in html.lower()

    def test_fdr_power_penalty(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        # Should mention FDR correction impact on power
        assert "FDR correction" in html

    def test_observed_vs_detectable_effect_sizes(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Observed vs Detectable" in html

    def test_detectable_hrs_classified(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        # mean_adc HR=1.5, log(1.5)≈0.405, should be within or below detection range
        assert "covariate(s)" in html

    def test_no_log_data_still_shows_power(self):
        """Power section should still render with mat_data but no log_data."""
        result = _section_power_analysis(None, [], _make_mat_data())
        html = "\n".join(result)
        assert "Power" in html
        assert "Detectable" in html or "detectable" in html

    def test_no_glme_tests_no_fdr_note(self):
        log = {"Standard": {"survival": {"hazard_ratios": []}, "stats_comparisons": {}}}
        result = _section_power_analysis(log, ["Standard"], _make_mat_data())
        html = "\n".join(result)
        # With 0 GLME tests, no FDR penalty note
        assert "FDR correction" not in html

    def test_post_hoc_disclaimer(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Post-hoc" in html or "post-hoc" in html or "approximate" in html.lower()

    def test_wilcoxon_power_info(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Wilcoxon" in html

    def test_cox_power_info(self):
        html = "\n".join(_section_power_analysis(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Cox" in html or "hazard ratio" in html.lower()
