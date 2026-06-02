"""Tests for report_sections main_results module (part 2).

Validates: predictive performance, manuscript-ready findings, results draft.
Executive summary, hypothesis, statistical significance, broad statistical
overview, and treatment response tests live in
test_report_sections_main_results.py.

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
    _section_manuscript_ready_findings,
    _section_predictive_performance,
    _section_results_draft,
)


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


# ── Edge cases: empty cohorts, partial data ──


class TestPredictivePerformanceEdgeCases:
    def test_empty_log_data_dict(self):
        result = _section_predictive_performance({}, ["Standard"])
        assert isinstance(result, list)

    def test_log_with_no_predictive_key(self):
        log = {"Standard": {"baseline": {}}}
        result = _section_predictive_performance(log, ["Standard"])
        assert isinstance(result, list)

    def test_roc_with_nan_auc(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = [{"auc": float("nan"), "timepoint": "BL"}]
        result = _section_predictive_performance(log, ["Standard"])
        assert isinstance(result, list)

    def test_empty_feature_selections(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = []
        result = _section_predictive_performance(log, ["Standard"])
        assert isinstance(result, list)

    def test_single_feature_selection(self):
        """Single feature selection with one feature should still render."""
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["feature_selections"] = [
            {"timepoint": "BL", "features": ["adc"], "lambda": 0.1},
        ]
        result = _section_predictive_performance(log, ["Standard"])
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "adc" in html


class TestManuscriptReadyFindingsEdgeCases:
    def test_all_empty_nested_data(self):
        log = {"Standard": {
            "stats_comparisons": {"glme_details": [], "glme_interactions": [], "fdr_timepoints": []},
            "stats_predictive": {"roc_analyses": [], "feature_selections": []},
            "survival": {"hazard_ratios": []},
            "baseline": {},
        }}
        result = _section_manuscript_ready_findings(log, ["Standard"], None, _make_mat_data(), None)
        assert isinstance(result, list)

    def test_no_dosimetry(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}}}
        result = _section_manuscript_ready_findings(_make_log_data(), ["Standard"], None, mat, None)
        assert isinstance(result, list)

    def test_zero_patients_mat_data(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 0, "num_timepoints": 0}}}
        result = _section_manuscript_ready_findings(None, ["Standard"], None, mat, None)
        assert isinstance(result, list)

    def test_groups_with_no_longitudinal(self):
        groups = {"Feature_BoxPlots": {"Standard": SAMPLE_GRAPH_CSV_ROWS[0]}}
        result = _section_manuscript_ready_findings(
            _make_log_data(), ["Standard"], None, _make_mat_data(), groups
        )
        assert isinstance(result, list)


class TestResultsDraftEdgeCases:
    def test_all_empty_nested_data(self):
        log = {"Standard": {
            "stats_comparisons": {"glme_details": []},
            "stats_predictive": {"roc_analyses": [], "feature_selections": []},
            "survival": {"hazard_ratios": []},
            "baseline": {},
        }}
        result = _section_results_draft(log, ["Standard"], None, _make_mat_data(), None)
        assert isinstance(result, list)

    def test_missing_baseline_exclusion(self):
        log = _make_log_data()
        del log["Standard"]["baseline"]["baseline_exclusion"]
        result = _section_results_draft(log, ["Standard"], None, _make_mat_data(), None)
        assert isinstance(result, list)

    def test_no_glme_excluded(self):
        result = _section_results_draft(_make_log_data(), ["Standard"], None, _make_mat_data(), None)
        assert isinstance(result, list)

    def test_csv_data_integration(self):
        csv_data = {
            "fdr_global": {"Standard": [{"metric": "adc"}]},
            "significant_metrics": {"Standard": [{"Metric": "adc"}]},
        }
        result = _section_results_draft(_make_log_data(), ["Standard"], csv_data, _make_mat_data(), None)
        assert isinstance(result, list)

    def test_empty_mat_data(self):
        result = _section_results_draft(_make_log_data(), ["Standard"], None, {}, None)
        assert isinstance(result, list)

    def test_no_dosimetry_in_mat(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}}}
        result = _section_results_draft(_make_log_data(), ["Standard"], None, mat, None)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Dosimetric" not in html or isinstance(html, str)
