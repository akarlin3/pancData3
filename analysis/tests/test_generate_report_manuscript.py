"""Tests for manuscript and publication-related report features.

Covers:
- _section_manuscript_ready_findings: copyable manuscript sentences
- Manuscript findings with effect sizes
- _section_reporting_checklist: STROBE/REMARK compliance
- _section_results_draft: auto-generated Results section
- BibTeX export in references
- Table index (List of Tables)
- _section_journal_guide: journal submission guidance
"""

from __future__ import annotations

from pathlib import Path

import pytest  # type: ignore

from report.generate_report import (  # type: ignore
    _section_journal_guide,
    _section_manuscript_ready_findings,
    _section_reporting_checklist,
    _section_results_draft,
    _section_table_index,
    generate_report,
    get_numbering,
    reset_numbering,
    REPORT_JS,
)


# ---------------------------------------------------------------------------
# Manuscript ready findings
# ---------------------------------------------------------------------------

class TestManuscriptReadyFindings:
    """Verify the Key Findings for Manuscript section."""

    def test_empty_without_data(self):
        """Section returns empty when no data is available."""
        result = _section_manuscript_ready_findings(None, [], None, {}, {})
        assert result == []

    def test_generates_sentences_with_log_data(self, saved_files_with_logs: Path):
        """Section generates copyable sentences from log data."""
        from parsers.parse_log_metrics import parse_all_logs  # type: ignore
        log_data = parse_all_logs(saved_files_with_logs)
        result = _section_manuscript_ready_findings(
            log_data, ["Standard"], None, {}, {}
        )
        html = "\n".join(result)
        assert "manuscript-sentence" in html
        assert "Manuscript" in html
        # Should have copy buttons
        assert "copy-btn" in html

    def test_includes_auc_sentence(self, saved_files_with_logs: Path):
        """Section includes AUC performance sentence when ROC data exists."""
        from parsers.parse_log_metrics import parse_all_logs  # type: ignore
        log_data = parse_all_logs(saved_files_with_logs)
        result = _section_manuscript_ready_findings(
            log_data, ["Standard"], None, {}, {}
        )
        html = "\n".join(result)
        assert "AUC" in html or "elastic-net" in html.lower()

    def test_includes_hazard_ratio_sentence(self, saved_files_with_logs: Path):
        """Section includes HR sentence when survival data exists."""
        from parsers.parse_log_metrics import parse_all_logs  # type: ignore
        log_data = parse_all_logs(saved_files_with_logs)
        result = _section_manuscript_ready_findings(
            log_data, ["Standard"], None, {}, {}
        )
        html = "\n".join(result)
        assert "HR" in html or "hazard" in html.lower() or "Cox" in html

    def test_copy_all_button_present(self, saved_files_with_logs: Path):
        """Section has a copy-all button for the full paragraph."""
        from parsers.parse_log_metrics import parse_all_logs  # type: ignore
        log_data = parse_all_logs(saved_files_with_logs)
        result = _section_manuscript_ready_findings(
            log_data, ["Standard"], None, {}, {}
        )
        html = "\n".join(result)
        assert "all-findings" in html


# ---------------------------------------------------------------------------
# Manuscript Findings with Effect Sizes
# ---------------------------------------------------------------------------

class TestManuscriptFindingsEffectSizes:
    """Verify that manuscript sentences include effect sizes."""

    def test_hr_sentence_includes_effect_size(self):
        """Cox PH manuscript sentences include effect size labels."""
        log_data = {"Standard": {"survival": {
            "hazard_ratios": [
                {"covariate": "mean_adc", "hr": 3.0,
                 "ci_lo": 1.5, "ci_hi": 6.0, "p": 0.002},
            ]
        }}}
        result = _section_manuscript_ready_findings(
            log_data, ["Standard"], None, {}, {})
        html = "\n".join(result)
        assert "effect" in html.lower()
        assert "3.000" in html

    def test_glme_sentence_includes_sample_size(self):
        """GLME sentences include sample size when exclusion data present."""
        log_data = {"Standard": {"stats_comparisons": {
            "glme_details": [
                {"metric": "mean_adc", "p": 0.001, "adj_alpha": 0.01},
            ],
            "glme_excluded": {"n_excluded": 3, "n_total": 25, "pct": 12.0},
        }}}
        result = _section_manuscript_ready_findings(
            log_data, ["Standard"], None, {}, {})
        html = "\n".join(result)
        assert "22 evaluable patients" in html
        assert "excluding 3" in html


# ---------------------------------------------------------------------------
# Reporting checklist
# ---------------------------------------------------------------------------

class TestReportingChecklist:
    """Verify the STROBE/REMARK reporting checklist section."""

    def test_checklist_with_no_data(self):
        """Checklist renders even without any data."""
        result = _section_reporting_checklist(None, [], {}, None, [])
        html = "\n".join(result)
        assert "Reporting Guideline" in html
        assert "STROBE" in html
        assert "REMARK" in html

    def test_checklist_counts_addressed_items(self, saved_files_with_logs: Path):
        """Checklist correctly counts addressed vs partial items."""
        from parsers.parse_log_metrics import parse_all_logs  # type: ignore
        log_data = parse_all_logs(saved_files_with_logs)
        result = _section_reporting_checklist(
            log_data, ["Standard"], {}, None, []
        )
        html = "\n".join(result)
        assert "Addressed" in html
        assert "checklist-done" in html

    def test_checklist_has_dynamic_status(self, saved_files_with_logs: Path):
        """Checklist status changes based on available data."""
        from parsers.parse_log_metrics import parse_all_logs  # type: ignore
        log_data = parse_all_logs(saved_files_with_logs)
        result_with = _section_reporting_checklist(
            log_data, ["Standard"], {}, None, []
        )
        result_without = _section_reporting_checklist(
            None, [], {}, None, []
        )
        # With data should have more 'done' items
        done_with = "\n".join(result_with).count("checklist-done")
        done_without = "\n".join(result_without).count("checklist-done")
        assert done_with >= done_without


# ---------------------------------------------------------------------------
# Table index
# ---------------------------------------------------------------------------

class TestTableIndex:
    """Verify the List of Tables section."""

    def test_empty_without_tables(self):
        """Returns empty when no tables have been numbered."""
        reset_numbering()
        result = _section_table_index()
        assert result == []

    def test_lists_numbered_tables(self):
        """Lists all tables that were numbered during report generation."""
        reset_numbering()
        from report.report_formatters import _table_caption  # type: ignore
        _table_caption("Test Table One", "description")
        _table_caption("Test Table Two")
        result = _section_table_index()
        html = "\n".join(result)
        assert "Test Table One" in html
        assert "Test Table Two" in html
        assert "Table 1" in html
        assert "Table 2" in html
        reset_numbering()


# ---------------------------------------------------------------------------
# BibTeX export
# ---------------------------------------------------------------------------

class TestBibTeXExport:
    """Verify BibTeX export in references section."""

    def test_report_contains_bibtex(self, saved_files_with_graph_csv: Path):
        """Report includes BibTeX export block."""
        report = generate_report(saved_files_with_graph_csv)
        assert "BibTeX" in report
        assert "@article" in report or "@book" in report

    def test_report_contains_js(self, saved_files_with_graph_csv: Path):
        """Report includes JavaScript for copy-to-clipboard."""
        report = generate_report(saved_files_with_graph_csv)
        assert "copyText" in report

    def test_report_js_constant_exists(self):
        """REPORT_JS constant is a non-empty string."""
        assert isinstance(REPORT_JS, str)
        assert len(REPORT_JS) > 50
        assert "copyText" in REPORT_JS


# ---------------------------------------------------------------------------
# Draft Results Section
# ---------------------------------------------------------------------------

class TestResultsDraft:
    """Verify the auto-generated Results section draft."""

    def test_empty_when_no_data(self):
        """No results draft without any data."""
        result = _section_results_draft(None, [], None, {}, {})
        assert result == []

    def test_cohort_paragraph_from_mat_data(self):
        """Cohort paragraph is generated from MAT data."""
        mat_data = {"Standard": {"longitudinal": {
            "num_patients": 25, "num_timepoints": 5}}}
        result = _section_results_draft(None, ["Standard"], None, mat_data, {})
        html = "\n".join(result)
        assert "25 patients" in html
        assert "5 timepoints" in html
        assert "Draft Results" in html

    def test_group_comparisons_paragraph(self):
        """Group comparisons paragraph is generated from GLME data."""
        log_data = {"Standard": {"stats_comparisons": {
            "glme_details": [
                {"metric": "mean_adc", "p": 0.001, "adj_alpha": 0.01},
                {"metric": "mean_d", "p": 0.5, "adj_alpha": 0.025},
            ],
        }}}
        result = _section_results_draft(log_data, ["Standard"], None, {}, {})
        html = "\n".join(result)
        assert "GLME" in html
        assert "mean_adc" in html

    def test_predictive_paragraph(self):
        """Predictive modelling paragraph with AUC."""
        log_data = {"Standard": {"stats_predictive": {
            "roc_analyses": [{"auc": 0.82, "timepoint": "W2",
                              "sensitivity": 85.0, "specificity": 75.0}],
            "feature_selections": [{"timepoint": "W2", "lambda": 0.01,
                                    "features": ["f1", "f2"]}],
        }}}
        result = _section_results_draft(log_data, ["Standard"], None, {}, {})
        html = "\n".join(result)
        assert "0.82" in html or "0.820" in html
        assert "sensitivity" in html.lower()

    def test_survival_paragraph(self):
        """Survival analysis paragraph with hazard ratios."""
        log_data = {"Standard": {"survival": {
            "hazard_ratios": [
                {"covariate": "mean_adc_delta", "hr": 2.5,
                 "ci_lo": 1.2, "ci_hi": 5.1, "p": 0.01},
            ],
            "global_lrt": {"df": 3, "chi2": 12.5, "p": 0.006},
            "ipcw": {"min_weight": 0.8, "max_weight": 1.2},
        }}}
        result = _section_results_draft(log_data, ["Standard"], None, {}, {})
        html = "\n".join(result)
        assert "mean_adc_delta" in html
        assert "HR" in html
        assert "IPCW" in html

    def test_copy_all_button_present(self):
        """Results draft includes a copy-all button."""
        mat_data = {"Standard": {"longitudinal": {
            "num_patients": 10, "num_timepoints": 3}}}
        result = _section_results_draft(None, ["Standard"], None, mat_data, {})
        html = "\n".join(result)
        assert "results-draft-all" in html
        assert "copy" in html.lower()


# ---------------------------------------------------------------------------
# Journal Submission Guidance
# ---------------------------------------------------------------------------

class TestJournalGuide:
    """Verify the journal submission guidance section."""

    def test_basic_output(self):
        """Section includes journal recommendations."""
        result = _section_journal_guide(None, ["Standard"], {})
        html = "\n".join(result)
        assert "Journal Submission Guidance" in html
        assert "Radiotherapy and Oncology" in html
        assert "Manuscript Preparation Checklist" in html

    def test_adds_survival_journal(self):
        """Adds Acta Oncologica when survival data present."""
        log_data = {"Standard": {"survival": {
            "hazard_ratios": [{"covariate": "x", "hr": 1.5, "p": 0.03}]
        }}}
        result = _section_journal_guide(log_data, ["Standard"], {})
        html = "\n".join(result)
        assert "Acta Oncologica" in html

    def test_adds_predictive_journal(self):
        """Adds European Radiology when predictive data present."""
        log_data = {"Standard": {"stats_predictive": {
            "roc_analyses": [{"auc": 0.8, "timepoint": "BL"}]
        }}}
        result = _section_journal_guide(log_data, ["Standard"], {})
        html = "\n".join(result)
        assert "European Radiology" in html

    def test_keywords_copyable(self):
        """Keywords section includes copy button."""
        result = _section_journal_guide(None, ["Standard"], {})
        html = "\n".join(result)
        assert "diffusion-weighted imaging" in html
        assert "copy" in html.lower()

    def test_keywords_include_survival(self):
        """Keywords include survival terms when survival data present."""
        log_data = {"Standard": {"survival": {
            "hazard_ratios": [{"covariate": "x", "hr": 1.5, "p": 0.03}]
        }}}
        result = _section_journal_guide(log_data, ["Standard"], {})
        html = "\n".join(result)
        assert "survival analysis" in html
