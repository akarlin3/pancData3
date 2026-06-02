"""Tests for report_sections.publication module (split from
test_report_sections_discussion.py).

Validates publication-related discussion section builders:
- Reporting checklist
- Journal guide
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from _test_report_sections_discussion_shared import _make_log_data, _make_mat_data

from report.sections.publication import (
    _section_reporting_checklist,
    _section_journal_guide,
)


# ── Reporting Checklist ──


class TestReportingChecklist:
    def test_returns_html(self):
        result = _section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        )
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Checklist" in html or "checklist" in html

    def test_no_data(self):
        result = _section_reporting_checklist(None, [], {}, None, None)
        assert isinstance(result, list)

    def test_contains_strobe_or_remark(self):
        html = "\n".join(_section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        ))
        assert "STROBE" in html or "REMARK" in html or "TRIPOD" in html or "checklist" in html.lower()

    def test_stat_cards_for_counts(self):
        html = "\n".join(_section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        ))
        assert "Addressed" in html
        assert "Partial" in html

    def test_full_data_increases_done_count(self):
        result_full = _section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        )
        result_empty = _section_reporting_checklist(None, [], {}, None, None)
        html_full = "\n".join(result_full)
        html_empty = "\n".join(result_empty)
        # Full data should have more "Addressed" items
        assert html_full.count("checklist-done") >= html_empty.count("checklist-done")

    def test_na_items_shown_when_present(self):
        html = "\n".join(_section_reporting_checklist(None, [], {}, None, None))
        assert "N/A" in html

    def test_checklist_table_structure(self):
        html = "\n".join(_section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        ))
        assert "<thead>" in html
        assert "Item" in html
        assert "Requirement" in html
        assert "Status" in html

    def test_strobe_items_present(self):
        html = "\n".join(_section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        ))
        assert "STROBE 1" in html
        assert "STROBE 12" in html
        assert "STROBE 19" in html

    def test_remark_items_present(self):
        html = "\n".join(_section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        ))
        assert "REMARK 1" in html
        assert "REMARK 5" in html

    def test_with_csv_data(self):
        csv_data = {"fdr_global": {"Standard": [{"metric": "adc"}]}}
        result = _section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), csv_data, None
        )
        assert isinstance(result, list)

    def test_with_graph_rows(self):
        rows = [{"graph_type": "scatter", "dwi_type": "Standard"}]
        result = _section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, rows
        )
        assert isinstance(result, list)

    def test_methods_section_always_done(self):
        html = "\n".join(_section_reporting_checklist(None, [], {}, None, None))
        # STROBE 12 "Statistical methods" should always be "done"
        assert "Full methods section auto-generated" in html


# ── Journal Guide ──


class TestJournalGuide:
    def test_returns_html(self):
        result = _section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data())
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Journal" in html or "journal" in html

    def test_no_data(self):
        result = _section_journal_guide(None, [], {})
        assert isinstance(result, list)

    def test_suggests_journals(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert any(j in html for j in [
            "Radiology", "Physics", "Oncology", "Cancer",
            "Medical", "International Journal", "journal",
        ])

    def test_survival_adds_acta_oncologica(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Acta Oncologica" in html

    def test_predictive_adds_european_radiology(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "European Radiology" in html

    def test_no_survival_no_acta(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": []},
            "stats_predictive": {"roc_analyses": []},
        }}
        html = "\n".join(_section_journal_guide(log, ["Standard"], _make_mat_data()))
        assert "Acta Oncologica" not in html

    def test_no_predictive_no_european_radiology(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": [{"covariate": "x", "hr": 1.5, "p": 0.03}]},
            "stats_predictive": {"roc_analyses": []},
        }}
        html = "\n".join(_section_journal_guide(log, ["Standard"], _make_mat_data()))
        assert "European Radiology" not in html

    def test_manuscript_checklist(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Manuscript Preparation" in html
        assert "Title page" in html
        assert "Abstract" in html
        assert "IRB" in html

    def test_suggested_keywords(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "diffusion-weighted imaging" in html
        assert "IVIM" in html
        assert "pancreatic cancer" in html

    def test_keywords_include_elastic_net_when_predictive(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "elastic net" in html

    def test_keywords_include_survival_when_present(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "survival analysis" in html

    def test_keywords_include_deep_learning_multi_dwi(self):
        html = "\n".join(_section_journal_guide(
            _make_log_data(), ["Standard", "dnCNN"], _make_mat_data()
        ))
        assert "deep learning" in html.lower()

    def test_no_deep_learning_keyword_single_dwi(self):
        html = "\n".join(_section_journal_guide(
            _make_log_data(), ["Standard"], _make_mat_data()
        ))
        assert "deep learning denoising" not in html

    def test_copy_button_for_keywords(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "copy-btn" in html or "Copy" in html

    def test_journal_table_structure(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Word Limit" in html
        assert "Article Type" in html
        assert "Scope Match" in html

    def test_base_journals_always_present(self):
        html = "\n".join(_section_journal_guide(None, [], {}))
        assert "Radiotherapy and Oncology" in html
        assert "Medical Physics" in html
        assert "Physics in Medicine" in html

    def test_empty_dwi_types(self):
        result = _section_journal_guide(None, [], {})
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Journal" in html


# ── Edge cases ──


class TestReportingChecklistEdgeCases:
    def test_all_none_inputs(self):
        result = _section_reporting_checklist(None, None, None, None, None)
        assert isinstance(result, list)

    def test_comprehensive_data(self):
        from conftest import SAMPLE_GRAPH_CSV_ROWS
        csv_data = {"fdr_global": {"Standard": [{"metric": "adc"}]}}
        rows = SAMPLE_GRAPH_CSV_ROWS
        result = _section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), csv_data, rows
        )
        html = "\n".join(result)
        assert "Addressed" in html

    def test_partial_log_data(self):
        log = {"Standard": {"stats_comparisons": {"glme_details": [
            {"metric": "m1", "p": 0.01, "adj_alpha": 0.025},
        ]}}}
        result = _section_reporting_checklist(log, ["Standard"], {}, None, None)
        assert isinstance(result, list)

    def test_empty_dwi_types_and_data(self):
        result = _section_reporting_checklist({}, [], {}, {}, [])
        assert isinstance(result, list)


class TestJournalGuideEdgeCases:
    def test_all_none_inputs(self):
        result = _section_journal_guide(None, [], None)
        assert isinstance(result, list)

    def test_mat_data_zero_patients(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 0, "num_timepoints": 0}}}
        result = _section_journal_guide(None, [], mat)
        assert isinstance(result, list)

    def test_only_survival_data(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": [{"covariate": "x", "hr": 1.5, "p": 0.03}]},
            "stats_predictive": {"roc_analyses": []},
        }}
        html = "\n".join(_section_journal_guide(log, ["Standard"], _make_mat_data()))
        assert "Acta Oncologica" in html
        assert "European Radiology" not in html

    def test_only_predictive_data(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": []},
            "stats_predictive": {"roc_analyses": [{"auc": 0.8}]},
        }}
        html = "\n".join(_section_journal_guide(log, ["Standard"], _make_mat_data()))
        assert "European Radiology" in html
        assert "Acta Oncologica" not in html

    def test_multiple_dwi_types_keywords(self):
        html = "\n".join(_section_journal_guide(
            _make_log_data(), ["Standard", "dnCNN", "IVIMnet"], _make_mat_data()
        ))
        assert "deep learning" in html.lower()

    def test_empty_log_dict(self):
        result = _section_journal_guide({}, ["Standard"], _make_mat_data())
        assert isinstance(result, list)
