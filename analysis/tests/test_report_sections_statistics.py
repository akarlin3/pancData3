"""Tests for report sections: effect_sizes and multiple_comparisons.

Validates statistics section builders (part 1):
- Effect size analysis (HR effect sizes, AUC interpretation)
- Multiple comparisons correction summary

Model diagnostics, sensitivity analysis, and power analysis live in
test_report_sections_statistics_part2.py. Shared imports and fixtures live in
_test_report_sections_statistics_shared.py.
"""

from __future__ import annotations

from _test_report_sections_statistics_shared import *  # noqa: F401,F403


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


# ── Effect Sizes edge cases ──


class TestEffectSizesEdgeCases:
    def test_missing_ci_fields(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": [
                {"covariate": "x", "hr": 1.5, "p": 0.03},
            ]},
            "stats_predictive": {"roc_analyses": []},
        }}
        result = _section_effect_sizes(log, ["Standard"], None)
        assert isinstance(result, list)

    def test_very_small_hr(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": [
                {"covariate": "protective", "hr": 0.1, "ci_lo": 0.05, "ci_hi": 0.2, "p": 0.001},
            ]},
            "stats_predictive": {"roc_analyses": []},
        }}
        result = _section_effect_sizes(log, ["Standard"], None)
        html = "\n".join(result)
        assert "0.1" in html or "protective" in html

    def test_very_large_hr(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": [
                {"covariate": "risky", "hr": 10.0, "ci_lo": 2.0, "ci_hi": 50.0, "p": 0.001},
            ]},
            "stats_predictive": {"roc_analyses": []},
        }}
        result = _section_effect_sizes(log, ["Standard"], None)
        html = "\n".join(result)
        assert "10.0" in html or "risky" in html

    def test_csv_data_with_fdr(self):
        csv_data = {"fdr_global": {"Standard": [{"metric": "adc", "p": 0.001}]}}
        result = _section_effect_sizes(_make_log_data(), ["Standard"], csv_data)
        assert isinstance(result, list)

    def test_all_auc_zero(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": [
                {"covariate": "x", "hr": 1.5, "ci_lo": 0.9, "ci_hi": 2.5, "p": 0.03},
            ]},
            "stats_predictive": {"roc_analyses": [{"auc": 0, "timepoint": "BL"}]},
        }}
        result = _section_effect_sizes(log, ["Standard"], None)
        assert isinstance(result, list)

    def test_single_hr_entry(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": [
                {"covariate": "only_var", "hr": 2.0, "ci_lo": 1.1, "ci_hi": 3.5, "p": 0.02},
            ]},
            "stats_predictive": {"roc_analyses": []},
        }}
        result = _section_effect_sizes(log, ["Standard"], None)
        html = "\n".join(result)
        assert "only_var" in html


class TestMultipleComparisonsEdgeCases:
    def test_single_glme_detail(self):
        log = {"Standard": {"stats_comparisons": {"glme_details": [
            {"metric": "m1", "p": 0.01, "adj_alpha": 0.05},
        ]}}}
        result = _section_multiple_comparisons(log, ["Standard"], None)
        html = "\n".join(result)
        assert "1" in html

    def test_all_rejected_by_fdr(self):
        log = {"Standard": {"stats_comparisons": {"glme_details": [
            {"metric": "m1", "p": 0.04, "adj_alpha": 0.01},
            {"metric": "m2", "p": 0.03, "adj_alpha": 0.01},
        ]}}}
        html = "\n".join(_section_multiple_comparisons(log, ["Standard"], None))
        assert "Rejected by FDR" in html

    def test_mixed_dwi_types(self):
        log = {
            "Standard": {"stats_comparisons": {"glme_details": [
                {"metric": "m1", "p": 0.01, "adj_alpha": 0.025},
            ]}},
            "dnCNN": {"stats_comparisons": {"glme_details": [
                {"metric": "m1", "p": 0.04, "adj_alpha": 0.05},
            ]}},
        }
        result = _section_multiple_comparisons(log, ["Standard", "dnCNN"], None)
        assert isinstance(result, list)
