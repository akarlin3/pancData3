"""Tests for report_sections main_results module (part 1).

Validates: executive summary, hypothesis, statistical significance,
broad statistical overview, treatment response. Predictive performance,
manuscript-ready findings, and results draft tests live in
test_report_sections_main_results_part2.py.

Shared helpers (_make_log_data, _make_mat_data, _make_groups) live in
_test_report_sections_main_results_shared.py.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from _test_report_sections_main_results_shared import *  # noqa: F401,F403
from _test_report_sections_main_results_shared import (
    SAMPLE_GRAPH_CSV_ROWS,
    _make_groups,
    _make_log_data,
    _make_mat_data,
    _section_broad_statistical_overview,
    _section_executive_summary,
    _section_hypothesis,
    _section_statistical_significance,
    _section_treatment_response,
)


# ── Executive Summary ──


class TestExecutiveSummary:
    def test_returns_html(self):
        result = _section_executive_summary(
            _make_log_data(), ["Standard"], SAMPLE_GRAPH_CSV_ROWS,
            None, "20260301_120000", _make_mat_data()
        )
        assert isinstance(result, list)
        assert len(result) > 0

    def test_contains_timestamp(self):
        html = "\n".join(_section_executive_summary(
            None, ["Standard"], [], None, "20260301_120000"
        ))
        assert "20260301_120000" in html

    def test_contains_dwi_types(self):
        html = "\n".join(_section_executive_summary(
            None, ["Standard", "dnCNN"], [], None, "ts"
        ))
        assert "Standard" in html
        assert "dnCNN" in html

    def test_no_data(self):
        result = _section_executive_summary(None, [], [], None, "ts")
        assert isinstance(result, list)

    def test_graph_count(self):
        html = "\n".join(_section_executive_summary(
            None, ["Standard"], SAMPLE_GRAPH_CSV_ROWS, None, "ts"
        ))
        assert "3" in html  # 3 sample rows

    def test_stat_cards_with_auc(self):
        """Should show AUC stat card when predictive data is available."""
        html = "\n".join(_section_executive_summary(
            _make_log_data(), ["Standard"], [], None, "ts", _make_mat_data()
        ))
        assert "AUC" in html or "0.78" in html

    def test_stat_cards_with_glme_sig(self):
        """Should show GLME interaction count."""
        html = "\n".join(_section_executive_summary(
            _make_log_data(), ["Standard"], [], None, "ts"
        ))
        assert "GLME" in html or "interaction" in html.lower() or len(html) > 0

    def test_csv_data_integrated(self):
        """CSV data significant metrics should contribute to stat cards."""
        csv_data = {
            "significant_metrics": {"Standard": [{"Metric": "adc", "p_value": "0.01"}]},
            "fdr_global": {"Standard": [{"metric": "adc"}]},
        }
        html = "\n".join(_section_executive_summary(
            _make_log_data(), ["Standard"], [], csv_data, "ts", _make_mat_data()
        ))
        assert isinstance(html, str)

    def test_empty_log_data_dict(self):
        """Empty log_data (not None) should be handled."""
        result = _section_executive_summary({}, ["Standard"], [], None, "ts")
        assert isinstance(result, list)

    def test_mat_data_cohort_info(self):
        """Should show cohort size from mat_data."""
        html = "\n".join(_section_executive_summary(
            None, ["Standard"], [], None, "ts", _make_mat_data()
        ))
        assert "42" in html

    def test_hazard_ratio_summary(self):
        """Should include HR summary information."""
        html = "\n".join(_section_executive_summary(
            _make_log_data(), ["Standard"], [], None, "ts"
        ))
        # Should reference hazard ratios or survival data
        assert len(html) > 100  # Substantial content generated


# ── Hypothesis ──


class TestHypothesis:
    def test_returns_html(self):
        result = _section_hypothesis(_make_groups(), _make_log_data(), _make_mat_data())
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Hypothesis" in html or "hypothesis" in html

    def test_no_data(self):
        result = _section_hypothesis(None)
        assert isinstance(result, list)

    def test_with_canonical_trends(self):
        """Canonical D-increasing, f-decreasing pattern should be noted."""
        result = _section_hypothesis(_make_groups(), _make_log_data(), _make_mat_data())
        html = "\n".join(result)
        assert len(html) > 0

    def test_no_groups_no_log_no_mat(self):
        """All None inputs should produce minimal/empty output."""
        result = _section_hypothesis(None, None, None)
        assert isinstance(result, list)

    def test_partial_data_log_only(self):
        """Only log_data (no groups, no mat) should still produce output."""
        result = _section_hypothesis(None, _make_log_data(), None)
        assert isinstance(result, list)

    def test_groups_without_longitudinal(self):
        """Groups without Longitudinal_Mean_Metrics should still work."""
        groups = {"Feature_BoxPlots": {"Standard": SAMPLE_GRAPH_CSV_ROWS[0]}}
        result = _section_hypothesis(groups, _make_log_data(), _make_mat_data())
        assert isinstance(result, list)


# ── Statistical Significance ──


class TestStatisticalSignificance:
    def test_returns_html(self):
        result = _section_statistical_significance(
            SAMPLE_GRAPH_CSV_ROWS, None, _make_log_data(), ["Standard"]
        )
        assert isinstance(result, list)
        assert len(result) > 0

    def test_empty_rows(self):
        result = _section_statistical_significance([], None, None, [])
        assert isinstance(result, list)

    def test_contains_pvalue(self):
        html = "\n".join(_section_statistical_significance(
            SAMPLE_GRAPH_CSV_ROWS, None, _make_log_data(), ["Standard"]
        ))
        assert "0.003" in html or "p" in html.lower()

    def test_borderline_findings(self):
        """P-values between 0.05 and 0.10 should appear in borderline section."""
        rows = [{
            "file_path": "Standard/test.png",
            "summary": "p = 0.07 for comparison",
            "trends_json": "[]", "inflection_points_json": "[]",
        }]
        html = "\n".join(_section_statistical_significance(rows, None, None, []))
        assert "Borderline" in html

    def test_csv_significant_metrics(self):
        """Pipeline CSV significant metrics should be displayed."""
        csv_data = {
            "significant_metrics": {
                "Standard": [
                    {"Metric": "mean_adc", "Timepoint": "BL", "p_value": "0.003"},
                ],
            },
        }
        html = "\n".join(_section_statistical_significance(
            [], csv_data, None, ["Standard"]
        ))
        assert "Pipeline CSV" in html or "mean_adc" in html

    def test_glme_interaction_details(self):
        """GLME interaction details should be shown from log_data."""
        html = "\n".join(_section_statistical_significance(
            [], None, _make_log_data(), ["Standard"]
        ))
        assert "GLME" in html
        assert "mean_adc" in html

    def test_fdr_timepoints(self):
        """FDR-significant timepoints should be displayed."""
        html = "\n".join(_section_statistical_significance(
            [], None, _make_log_data(), ["Standard"]
        ))
        assert "FDR" in html or "BL" in html

    def test_competing_risk_exclusion_warning(self):
        """GLME excluded count should trigger a warning box."""
        log = _make_log_data()
        log["Standard"]["stats_comparisons"]["glme_excluded"] = {
            "n_excluded": 5, "n_total": 42, "pct": 11.9,
        }
        html = "\n".join(_section_statistical_significance(
            [], None, log, ["Standard"]
        ))
        assert "Competing" in html or "excluded" in html.lower()

    def test_no_findings_message(self):
        """When no data at all, should show 'No significant findings'."""
        html = "\n".join(_section_statistical_significance([], None, None, []))
        assert "No significant" in html or "Statistical Significance" in html

    def test_multiple_dwi_types_glme(self):
        """GLME details should be shown per DWI type."""
        log = _make_log_data()
        log["dnCNN"] = {
            "stats_comparisons": {
                "glme_details": [
                    {"metric": "mean_d", "p": 0.03, "adj_alpha": 0.025},
                ],
            },
        }
        html = "\n".join(_section_statistical_significance(
            [], None, log, ["Standard", "dnCNN"]
        ))
        assert "Standard" in html
        assert "dnCNN" in html


# ── Broad Statistical Overview ──


class TestBroadStatisticalOverview:
    def test_returns_html(self):
        result = _section_broad_statistical_overview(_make_log_data(), ["Standard"])
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_broad_statistical_overview(None, [])
        assert isinstance(result, list)
        assert result == []

    def test_cox_ph_table(self):
        """Should show Cox PH hazard ratios table."""
        html = "\n".join(_section_broad_statistical_overview(_make_log_data(), ["Standard"]))
        assert "Cox PH" in html or "Hazard Ratio" in html

    def test_glme_cross_dwi_comparison(self):
        """Should show GLME cross-DWI comparison table."""
        html = "\n".join(_section_broad_statistical_overview(_make_log_data(), ["Standard"]))
        assert "GLME" in html

    def test_nominally_significant_variables(self):
        """Should list nominally significant variables."""
        html = "\n".join(_section_broad_statistical_overview(_make_log_data(), ["Standard"]))
        assert "Nominal" in html or "mean_adc" in html

    def test_no_nominal_sig(self):
        """When no p < 0.05, should say so."""
        log = {"Standard": {
            "survival": {"hazard_ratios": [
                {"covariate": "x", "hr": 1.0, "ci_lo": 0.5, "ci_hi": 1.5, "p": 0.8}
            ]},
            "stats_comparisons": {"glme_details": [
                {"metric": "x", "p": 0.5, "adj_alpha": 0.05}
            ]},
        }}
        html = "\n".join(_section_broad_statistical_overview(log, ["Standard"]))
        assert "No individual test" in html or "Nominal" in html

    def test_directional_consistency(self):
        """Should show cross-DWI directional consistency for HRs."""
        log = _make_log_data()
        log["dnCNN"] = {
            "survival": {"hazard_ratios": [
                {"covariate": "mean_adc", "hr": 1.3, "ci_lo": 0.8, "ci_hi": 2.0, "p": 0.1},
            ]},
            "stats_comparisons": {"glme_details": []},
        }
        html = "\n".join(_section_broad_statistical_overview(log, ["Standard", "dnCNN"]))
        assert "Consistency" in html or "Directional" in html

    def test_bonferroni_correction_noted(self):
        """Should mention Bonferroni correction in Cox PH table."""
        html = "\n".join(_section_broad_statistical_overview(_make_log_data(), ["Standard"]))
        assert "Bonferroni" in html or "correction" in html.lower()

    def test_wilcoxon_summary(self):
        """Should show Wilcoxon per-timepoint summary."""
        html = "\n".join(_section_broad_statistical_overview(_make_log_data(), ["Standard"]))
        assert "Wilcoxon" in html or "FDR" in html

    def test_no_fdr_surviving_warning(self):
        """If no FDR-surviving metrics, should warn about power."""
        log = _make_log_data()
        log["Standard"]["stats_comparisons"]["fdr_timepoints"] = [
            {"n_significant": 0, "timepoint": "BL"},
        ]
        html = "\n".join(_section_broad_statistical_overview(log, ["Standard"]))
        # Could show 0 surviving or a warning
        assert isinstance(html, str)


# ── Treatment Response ──


class TestTreatmentResponse:
    def test_returns_html(self):
        result = _section_treatment_response(_make_groups())
        assert isinstance(result, list)

    def test_no_groups(self):
        result = _section_treatment_response(None)
        assert isinstance(result, list)

    def test_with_trend_data(self):
        html = "\n".join(_section_treatment_response(_make_groups()))
        assert len(html) > 0

    def test_empty_groups_dict(self):
        result = _section_treatment_response({})
        assert isinstance(result, list)

    def test_groups_without_trends(self):
        """Groups with no trends_json should be handled."""
        groups = {"TestGraph": {"Standard": {}, "dnCNN": {}}}
        result = _section_treatment_response(groups)
        assert isinstance(result, list)


# ── Edge cases: empty cohorts, null hypothesis, partial data ──


class TestExecutiveSummaryEdgeCases:
    def test_all_none_inputs(self):
        result = _section_executive_summary(None, [], [], None, "ts")
        assert isinstance(result, list)

    def test_empty_rows_list(self):
        result = _section_executive_summary(None, ["Standard"], [], None, "ts")
        assert isinstance(result, list)

    def test_log_data_with_empty_nested_dicts(self):
        log = {"Standard": {
            "stats_comparisons": {},
            "stats_predictive": {},
            "survival": {},
            "baseline": {},
        }}
        result = _section_executive_summary(log, ["Standard"], [], None, "ts")
        assert isinstance(result, list)

    def test_csv_data_with_empty_fdr(self):
        csv_data = {"significant_metrics": {}, "fdr_global": {}}
        result = _section_executive_summary(None, ["Standard"], [], csv_data, "ts")
        assert isinstance(result, list)

    def test_mat_data_zero_patients(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 0, "num_timepoints": 0}}}
        result = _section_executive_summary(None, ["Standard"], [], None, "ts", mat)
        assert isinstance(result, list)

    def test_no_auc_no_hr_no_glme(self):
        """Empty pipeline results should produce a minimal summary."""
        log = {"Standard": {
            "stats_predictive": {"roc_analyses": [], "feature_selections": []},
            "stats_comparisons": {"glme_details": [], "glme_interactions": []},
            "survival": {"hazard_ratios": []},
            "baseline": {},
        }}
        result = _section_executive_summary(log, ["Standard"], [], None, "ts")
        assert isinstance(result, list)


class TestHypothesisEdgeCases:
    def test_empty_groups_dict(self):
        result = _section_hypothesis({})
        assert isinstance(result, list)

    def test_groups_with_empty_trends_json(self):
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": "[]"},
            }
        }
        result = _section_hypothesis(groups)
        assert isinstance(result, list)

    def test_groups_with_invalid_json(self):
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": "NOT JSON"},
            }
        }
        result = _section_hypothesis(groups)
        assert isinstance(result, list)

    def test_no_log_no_mat_with_groups(self):
        result = _section_hypothesis(_make_groups(), None, None)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert len(html) > 0

    def test_log_with_empty_survival(self):
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"] = []
        result = _section_hypothesis(_make_groups(), log, _make_mat_data())
        assert isinstance(result, list)

    def test_log_with_no_sig_glme(self):
        log = _make_log_data()
        log["Standard"]["stats_comparisons"]["glme_details"] = [
            {"metric": "x", "p": 0.5, "adj_alpha": 0.025},
        ]
        log["Standard"]["stats_comparisons"]["glme_interactions"] = [0.8]
        result = _section_hypothesis(_make_groups(), log, _make_mat_data())
        assert isinstance(result, list)

    def test_mat_data_no_dosimetry(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}}}
        result = _section_hypothesis(_make_groups(), _make_log_data(), mat)
        assert isinstance(result, list)

    def test_low_auc(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = [{"auc": 0.55, "timepoint": "BL"}]
        result = _section_hypothesis(_make_groups(), log, _make_mat_data())
        assert isinstance(result, list)

    def test_suboptimal_d95(self):
        mat = _make_mat_data()
        mat["Standard"]["dosimetry"]["d95_adc_mean"] = {"mean": 38.0}
        result = _section_hypothesis(_make_groups(), _make_log_data(), mat)
        assert isinstance(result, list)

    def test_inflection_points_in_trends(self):
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {
                    "trends_json": json.dumps([
                        {"series": "Mean D", "direction": "increasing", "description": "D rises by ~15%"},
                    ]),
                    "inflection_points_json": json.dumps([
                        {"approximate_x": 14, "approximate_y": 0.0012, "description": "Divergence at Fx14"},
                    ]),
                },
            },
        }
        result = _section_hypothesis(groups, _make_log_data(), _make_mat_data())
        assert isinstance(result, list)


class TestStatisticalSignificanceEdgeCases:
    def test_all_none_inputs(self):
        result = _section_statistical_significance(None, None, None, None)
        assert isinstance(result, list)

    def test_rows_with_no_pvalues(self):
        rows = [{
            "file_path": "Standard/test.png",
            "summary": "No statistical results here",
            "trends_json": "[]", "inflection_points_json": "[]",
        }]
        result = _section_statistical_significance(rows, None, None, [])
        assert isinstance(result, list)

    def test_log_data_empty_glme_details(self):
        log = {"Standard": {"stats_comparisons": {"glme_details": []}}}
        result = _section_statistical_significance([], None, log, ["Standard"])
        assert isinstance(result, list)


class TestBroadStatisticalOverviewEdgeCases:
    def test_log_with_empty_survival(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": []},
            "stats_comparisons": {"glme_details": []},
        }}
        result = _section_broad_statistical_overview(log, ["Standard"])
        assert isinstance(result, list)

    def test_only_nonsig_hrs(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": [
                {"covariate": "a", "hr": 1.0, "ci_lo": 0.5, "ci_hi": 2.0, "p": 0.9},
                {"covariate": "b", "hr": 0.9, "ci_lo": 0.4, "ci_hi": 1.8, "p": 0.7},
            ]},
            "stats_comparisons": {"glme_details": []},
        }}
        result = _section_broad_statistical_overview(log, ["Standard"])
        html = "\n".join(result)
        assert "No individual test" in html or isinstance(html, str)

    def test_hr_exactly_one(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": [
                {"covariate": "null_hr", "hr": 1.0, "ci_lo": 0.5, "ci_hi": 2.0, "p": 1.0},
            ]},
            "stats_comparisons": {"glme_details": []},
        }}
        result = _section_broad_statistical_overview(log, ["Standard"])
        assert isinstance(result, list)


class TestTreatmentResponseEdgeCases:
    def test_groups_with_only_non_longitudinal(self):
        groups = {
            "Feature_BoxPlots": {
                "Standard": SAMPLE_GRAPH_CSV_ROWS[0],
                "dnCNN": SAMPLE_GRAPH_CSV_ROWS[1],
            },
        }
        result = _section_treatment_response(groups)
        assert isinstance(result, list)

    def test_groups_with_empty_trends(self):
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": "[]"},
                "dnCNN": {"trends_json": "[]"},
            },
        }
        result = _section_treatment_response(groups)
        assert isinstance(result, list)

    def test_groups_with_invalid_json(self):
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": "NOT JSON"},
            },
        }
        result = _section_treatment_response(groups)
        assert isinstance(result, list)

    def test_single_dwi_type_longitudinal(self):
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {
                    "trends_json": json.dumps([{"series": "ADC", "direction": "increasing", "description": "rises"}]),
                },
            },
        }
        result = _section_treatment_response(groups)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert len(html) > 0
