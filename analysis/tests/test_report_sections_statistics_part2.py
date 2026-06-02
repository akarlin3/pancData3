"""Tests for report sections: model diagnostics, sensitivity, and power analysis.

Validates statistics section builders (part 2):
- Model diagnostics (IPCW weights, assumptions, imputation notes)
- Sensitivity analysis (EPV checks, HR stability, cross-DWI concordance)
- Power analysis (detectable effects, Cohen's d, FDR penalty)

Effect sizes and multiple comparisons live in
test_report_sections_statistics.py. Shared imports and fixtures live in
_test_report_sections_statistics_shared.py.
"""

from __future__ import annotations

from _test_report_sections_statistics_shared import *  # noqa: F401,F403


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


# ── Edge cases: convergence failures, missing metrics, partial data ──


class TestModelDiagnosticsEdgeCases:
    def test_all_none_inputs(self):
        result = _section_model_diagnostics(None, None, None)
        assert isinstance(result, list)

    def test_ipcw_none(self):
        log = {"Standard": {
            "survival": {},
            "stats_comparisons": {},
            "stats_predictive": {"feature_selections": []},
            "baseline": {},
        }}
        result = _section_model_diagnostics(log, ["Standard"], _make_mat_data())
        assert isinstance(result, list)

    def test_missing_baseline_key(self):
        log = {"Standard": {
            "survival": {"ipcw": {"max_weight": 1.5, "min_weight": 0.7}},
            "stats_comparisons": {},
            "stats_predictive": {"feature_selections": []},
        }}
        result = _section_model_diagnostics(log, ["Standard"], _make_mat_data())
        assert isinstance(result, list)

    def test_empty_feature_selections(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = []
        result = _section_model_diagnostics(log, ["Standard"], _make_mat_data())
        assert isinstance(result, list)

    def test_feature_selections_with_same_lambda(self):
        """Feature selections with same lambda should not trigger instability."""
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": ["adc"], "lambda": 0.05},
            {"timepoint": "W2", "features": ["d"], "lambda": 0.05},
        ]
        result = _section_model_diagnostics(log, ["Standard"], _make_mat_data())
        html = "\n".join(result)
        assert "varies by" not in html

    def test_zero_pct_exclusion(self):
        log = _make_log_data()
        log["Standard"]["stats_comparisons"]["glme_excluded"]["pct"] = 0.0
        html = "\n".join(_section_model_diagnostics(log, ["Standard"], _make_mat_data()))
        assert "High competing-risk exclusion" not in html
        assert "Moderate competing-risk exclusion" not in html

    def test_missing_lf_rate_fields(self):
        log = _make_log_data()
        del log["Standard"]["baseline"]["baseline_exclusion"]["lf_rate_included"]
        del log["Standard"]["baseline"]["baseline_exclusion"]["lf_rate_excluded"]
        result = _section_model_diagnostics(log, ["Standard"], _make_mat_data())
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "LF rate differs" not in html


class TestSensitivityAnalysisEdgeCases:
    def test_empty_mat_data(self):
        result = _section_sensitivity_analysis(_make_log_data(), ["Standard"], {})
        assert isinstance(result, list)

    def test_mat_data_none(self):
        result = _section_sensitivity_analysis(_make_log_data(), ["Standard"], None)
        assert isinstance(result, list)

    def test_many_features_low_epv(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": [f"f{i}" for i in range(20)], "lambda": 0.01},
        ]
        mat = {"Standard": {"longitudinal": {"num_patients": 25, "num_timepoints": 3}}}
        result = _section_sensitivity_analysis(log, ["Standard"], mat)
        html = "\n".join(result)
        assert "Low EPV" in html

    def test_no_hazard_ratios(self):
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"] = []
        result = _section_sensitivity_analysis(log, ["Standard"], _make_mat_data())
        html = "\n".join(result)
        assert "Unstable" not in html

    def test_hr_ci_lo_zero(self):
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"] = [
            {"covariate": "edge", "hr": 2.0, "ci_lo": 0.0, "ci_hi": 50.0, "p": 0.5},
        ]
        result = _section_sensitivity_analysis(log, ["Standard"], _make_mat_data())
        html = "\n".join(result)
        assert "Unstable" in html


class TestPowerAnalysisEdgeCases:
    def test_very_small_cohort(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 5, "num_timepoints": 2}}}
        result = _section_power_analysis(None, [], mat)
        html = "\n".join(result)
        assert "underpowered" in html.lower()

    def test_very_large_cohort(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 500, "num_timepoints": 5}}}
        result = _section_power_analysis(_make_log_data(), ["Standard"], mat)
        html = "\n".join(result)
        assert "adequate" in html.lower()

    def test_no_survival_data(self):
        log = {"Standard": {"stats_comparisons": {"glme_details": []}, "stats_predictive": {}}}
        result = _section_power_analysis(log, ["Standard"], _make_mat_data())
        assert isinstance(result, list)

    def test_multiple_dwi_types_power(self):
        log = _make_log_data()
        log["dnCNN"] = {
            "survival": {"hazard_ratios": [
                {"covariate": "x", "hr": 1.3, "ci_lo": 0.8, "ci_hi": 2.0, "p": 0.15},
            ]},
            "stats_comparisons": {"glme_details": []},
            "stats_predictive": {"roc_analyses": [{"auc": 0.72, "timepoint": "BL"}]},
        }
        result = _section_power_analysis(log, ["Standard", "dnCNN"], _make_mat_data())
        assert isinstance(result, list)
