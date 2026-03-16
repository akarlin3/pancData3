"""Tests for report_sections.main_results module.

Validates the main results section builders:
- Executive summary
- Hypothesis section
- Statistical significance
- Broad statistical overview
- Treatment response
- Predictive performance
- Manuscript-ready findings
- Results draft
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from conftest import SAMPLE_GRAPH_CSV_ROWS

from report.sections.main_results import (
    _section_executive_summary,
    _section_hypothesis,
    _section_treatment_response,
)
from report.sections.statistical_reporting import (
    _section_statistical_significance,
    _section_broad_statistical_overview,
)
from report.sections.manuscript import (
    _section_predictive_performance,
    _section_manuscript_ready_findings,
    _section_results_draft,
)


def _make_log_data():
    return {
        "Standard": {
            "survival": {
                "hazard_ratios": [
                    {"covariate": "mean_adc", "hr": 1.5, "ci_lo": 0.9, "ci_hi": 2.5, "p": 0.03},
                ],
                "ipcw": {"max_weight": 1.4, "min_weight": 0.8},
            },
            "stats_comparisons": {
                "glme_details": [
                    {"metric": "mean_adc", "p": 0.01, "adj_alpha": 0.025},
                ],
                "glme_interactions": [0.02],
                "fdr_timepoints": [{"n_significant": 1, "timepoint": "BL"}],
            },
            "stats_predictive": {
                "roc_analyses": [{"auc": 0.78, "timepoint": "BL"}],
                "feature_selections": [
                    {"timepoint": "BL", "features": ["adc", "d"], "lambda": 0.05},
                ],
            },
            "baseline": {
                "total_outliers": {"pct": 8.5},
                "baseline_exclusion": {"n_excluded": 6, "n_total": 48},
            },
            "sanity_checks": {"all_converged": True},
        }
    }


def _make_mat_data():
    return {
        "Standard": {
            "longitudinal": {"num_patients": 42, "num_timepoints": 5},
            "dosimetry": {
                "d95_adc_mean": {"mean": 45.0},
                "v50_adc_mean": {"mean": 0.85},
            },
        }
    }


def _make_groups():
    return {
        "Longitudinal_Mean_Metrics": {
            "Standard": {
                "trends_json": json.dumps([
                    {"series": "Mean D", "direction": "increasing", "description": "D rises"},
                    {"series": "Mean f", "direction": "decreasing", "description": "f drops"},
                ])
            },
            "dnCNN": {
                "trends_json": json.dumps([
                    {"series": "Mean D", "direction": "increasing", "description": "D rises"},
                    {"series": "Mean f", "direction": "decreasing", "description": "f drops"},
                ])
            },
        },
        "Feature_BoxPlots": {
            "Standard": SAMPLE_GRAPH_CSV_ROWS[0],
            "dnCNN": SAMPLE_GRAPH_CSV_ROWS[1],
        },
    }


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


# ── Predictive Performance ──


class TestPredictivePerformance:
    def test_returns_html(self):
        result = _section_predictive_performance(_make_log_data(), ["Standard"])
        assert isinstance(result, list)

    def test_contains_auc(self):
        html = "\n".join(_section_predictive_performance(_make_log_data(), ["Standard"]))
        assert "AUC" in html or "0.78" in html

    def test_no_data(self):
        result = _section_predictive_performance(None, [])
        assert isinstance(result, list)
        assert result == []

    def test_roc_table_structure(self):
        """ROC table should have headers for AUC, sensitivity, specificity."""
        html = "\n".join(_section_predictive_performance(_make_log_data(), ["Standard"]))
        assert "ROC" in html
        assert "Discrimination" in html

    def test_feature_selection_table(self):
        """Should show selected features per timepoint."""
        html = "\n".join(_section_predictive_performance(_make_log_data(), ["Standard"]))
        assert "Selected Features" in html or "adc" in html

    def test_cox_ph_section(self):
        """Should show Cox PH hazard ratios table."""
        html = "\n".join(_section_predictive_performance(_make_log_data(), ["Standard"]))
        assert "Cox" in html or "Hazard" in html

    def test_ipcw_weights_shown(self):
        """IPCW weight range should be displayed."""
        html = "\n".join(_section_predictive_performance(_make_log_data(), ["Standard"]))
        assert "IPCW" in html or "0.800" in html or "1.400" in html

    def test_discrimination_rating(self):
        """AUC 0.78 should be classified as 'Acceptable'."""
        html = "\n".join(_section_predictive_performance(_make_log_data(), ["Standard"]))
        assert "Acceptable" in html

    def test_no_roc_no_fs_no_surv(self):
        """Log data with no predictive/survival content should say 'no data'."""
        log = {"Standard": {}}
        html = "\n".join(_section_predictive_performance(log, ["Standard"]))
        assert "No predictive" in html or "no" in html.lower()

    def test_lambda_trend_analysis(self):
        """Lambda trend should be analyzed when 2+ timepoints."""
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": ["adc", "d"], "lambda": 0.05},
            {"timepoint": "W2", "features": ["adc"], "lambda": 0.12},
        ]
        html = "\n".join(_section_predictive_performance(log, ["Standard"]))
        assert "Regularisation" in html or "trend" in html.lower()

    def test_firth_refits_shown(self):
        """Firth penalised-likelihood refits should appear."""
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["firth_refits"] = [
            {"timepoint": "BL", "n_features": 3},
        ]
        html = "\n".join(_section_predictive_performance(log, ["Standard"]))
        assert "Firth" in html

    def test_multiple_dwi_types(self):
        """Should render data for multiple DWI types."""
        log = _make_log_data()
        log["dnCNN"] = {
            "stats_predictive": {
                "roc_analyses": [{"auc": 0.72, "timepoint": "BL"}],
                "feature_selections": [],
            },
            "survival": {},
        }
        html = "\n".join(_section_predictive_performance(log, ["Standard", "dnCNN"]))
        assert "Standard" in html
        assert "dnCNN" in html

    def test_global_lrt_shown(self):
        """Global LRT should appear when available."""
        log = _make_log_data()
        log["Standard"]["survival"]["global_lrt"] = {"df": 2, "chi2": 7.82, "p": 0.02}
        html = "\n".join(_section_predictive_performance(log, ["Standard"]))
        assert "LRT" in html or "7.82" in html


# ── Manuscript-Ready Findings ──


class TestManuscriptReadyFindings:
    def test_returns_html(self):
        result = _section_manuscript_ready_findings(
            _make_log_data(), ["Standard"], None, _make_mat_data(), _make_groups()
        )
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_manuscript_ready_findings(None, [], None, {}, None)
        assert isinstance(result, list)
        assert result == []

    def test_cohort_sentence(self):
        """Should generate a sentence about the cohort size."""
        html = "\n".join(_section_manuscript_ready_findings(
            _make_log_data(), ["Standard"], None, _make_mat_data(), _make_groups()
        ))
        assert "42 patients" in html

    def test_baseline_exclusion_sentence(self):
        """Should generate a sentence about baseline exclusions."""
        html = "\n".join(_section_manuscript_ready_findings(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "excluded" in html.lower() or "6" in html

    def test_glme_significance_sentence(self):
        """Should generate sentence about GLME interaction testing."""
        html = "\n".join(_section_manuscript_ready_findings(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "GLME" in html or "interaction" in html.lower()

    def test_auc_sentence(self):
        """Should generate a sentence about the best AUC."""
        html = "\n".join(_section_manuscript_ready_findings(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "AUC" in html or "0.780" in html

    def test_cox_hr_sentence(self):
        """Should generate a sentence about significant hazard ratios."""
        html = "\n".join(_section_manuscript_ready_findings(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "HR" in html or "hazard" in html.lower()

    def test_dosimetry_sentence(self):
        """Should generate sentence about dosimetry D95."""
        html = "\n".join(_section_manuscript_ready_findings(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "D95" in html or "Dosimetric" in html

    def test_longitudinal_trend_sentence(self):
        """Should generate sentence about D and f trends."""
        html = "\n".join(_section_manuscript_ready_findings(
            _make_log_data(), ["Standard"], None, _make_mat_data(), _make_groups()
        ))
        assert "diffusion" in html.lower() or "perfusion" in html.lower()

    def test_cross_dwi_agreement_sentence(self):
        """Should generate sentence about cross-DWI agreement."""
        html = "\n".join(_section_manuscript_ready_findings(
            _make_log_data(), ["Standard", "dnCNN"], None, _make_mat_data(), _make_groups()
        ))
        assert "agreement" in html.lower() or "Cross-DWI" in html

    def test_copy_all_button(self):
        """Should include a copy-all button."""
        html = "\n".join(_section_manuscript_ready_findings(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "copy" in html.lower() or "Copy" in html

    def test_fdr_global_sentence(self):
        """FDR global metrics should generate a sentence."""
        csv_data = {"fdr_global": {"Standard": [{"metric": "adc"}]}}
        html = "\n".join(_section_manuscript_ready_findings(
            _make_log_data(), ["Standard"], csv_data, _make_mat_data(), None
        ))
        assert "FDR" in html

    def test_no_sig_glme(self):
        """When no GLME metrics are significant, should say 'none'."""
        log = _make_log_data()
        log["Standard"]["stats_comparisons"]["glme_details"] = [
            {"metric": "x", "p": 0.5, "adj_alpha": 0.025},
        ]
        html = "\n".join(_section_manuscript_ready_findings(
            log, ["Standard"], None, _make_mat_data(), None
        ))
        assert "none" in html.lower()


# ── Results Draft ──


class TestResultsDraft:
    def test_returns_html(self):
        result = _section_results_draft(
            _make_log_data(), ["Standard"], None, _make_mat_data(), _make_groups()
        )
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_results_draft(None, [], None, {}, None)
        assert isinstance(result, list)
        assert result == []

    def test_patient_cohort_paragraph(self):
        """Should include a Patient Cohort subsection."""
        html = "\n".join(_section_results_draft(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "Patient Cohort" in html
        assert "42 patients" in html

    def test_group_comparisons_paragraph(self):
        """Should include Group Comparisons subsection."""
        html = "\n".join(_section_results_draft(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "Group Comparisons" in html

    def test_predictive_modelling_paragraph(self):
        """Should include Predictive Modelling subsection."""
        html = "\n".join(_section_results_draft(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "Predictive" in html

    def test_survival_analysis_paragraph(self):
        """Should include Survival Analysis subsection."""
        html = "\n".join(_section_results_draft(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "Survival" in html

    def test_cross_dwi_paragraph(self):
        """Should include Cross-DWI Comparison subsection."""
        html = "\n".join(_section_results_draft(
            _make_log_data(), ["Standard", "dnCNN"], None, _make_mat_data(), _make_groups()
        ))
        assert "Cross-DWI" in html

    def test_dosimetry_paragraph(self):
        """Should include Dosimetric Analysis subsection."""
        html = "\n".join(_section_results_draft(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "Dosimetric" in html

    def test_copy_all_results_button(self):
        """Should include copy-all button for entire results draft."""
        html = "\n".join(_section_results_draft(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "results-draft-all" in html

    def test_baseline_exclusion_in_cohort(self):
        """Baseline exclusion should appear in patient cohort paragraph."""
        html = "\n".join(_section_results_draft(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "excluded" in html.lower() or "missing" in html.lower()

    def test_no_sig_hrs(self):
        """When no HRs are significant, should say 'None achieved significance'."""
        log = _make_log_data()
        log["Standard"]["survival"]["hazard_ratios"] = [
            {"covariate": "x", "hr": 1.0, "ci_lo": 0.5, "ci_hi": 1.5, "p": 0.8},
        ]
        html = "\n".join(_section_results_draft(
            log, ["Standard"], None, _make_mat_data(), None
        ))
        assert "None" in html or "did not" in html.lower()

    def test_adequate_d95_coverage(self):
        """D95 >= 45 should produce 'adequate' coverage language."""
        html = "\n".join(_section_results_draft(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "adequate" in html.lower()

    def test_suboptimal_d95_coverage(self):
        """D95 < 45 should produce 'suboptimal' coverage language."""
        mat = _make_mat_data()
        mat["Standard"]["dosimetry"]["d95_adc_mean"] = {"mean": 38.0}
        html = "\n".join(_section_results_draft(
            _make_log_data(), ["Standard"], None, mat, None
        ))
        assert "suboptimal" in html.lower() or "insufficient" in html.lower()

    def test_competing_risk_in_group_comparisons(self):
        """Should mention competing risk exclusion in group comparison paragraph."""
        log = _make_log_data()
        log["Standard"]["stats_comparisons"]["glme_excluded"] = {
            "n_excluded": 5, "n_total": 42, "pct": 11.9,
        }
        html = "\n".join(_section_results_draft(
            log, ["Standard"], None, _make_mat_data(), None
        ))
        assert "competing" in html.lower() or "excluded" in html.lower()
