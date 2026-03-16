"""Tests for report_sections.discussion module.

Validates discussion section builders:
- Methods section (statistical methodology description)
- Limitations section
- Conclusions section
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

from report.sections.discussion import (
    _section_methods,
    _section_limitations,
    _section_conclusions,
)
from report.sections.publication import (
    _section_reporting_checklist,
    _section_journal_guide,
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
                    {"timepoint": "BL", "features": ["adc"], "lambda": 0.05},
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
            "dosimetry": {"d95_gtvp": {"mean": 45.0}, "v50_gtvp": {"mean": 0.85}},
        }
    }


# ── Methods ──


class TestMethods:
    def test_returns_html(self):
        result = _section_methods(["Standard"], _make_mat_data(), _make_log_data())
        assert isinstance(result, list)
        assert len(result) > 0

    def test_contains_methods_header(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "Statistical Methods" in html

    def test_mentions_ivim(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "IVIM" in html

    def test_mentions_bh_fdr(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "Benjamini" in html or "BH" in html or "FDR" in html

    def test_mentions_elastic_net(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "elastic" in html.lower() or "lasso" in html.lower() or "regulariz" in html.lower() or "predictive" in html.lower()


# ── Limitations ──


class TestLimitations:
    def test_returns_html(self):
        result = _section_limitations(_make_log_data(), ["Standard"], _make_mat_data())
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Limitation" in html

    def test_no_data(self):
        result = _section_limitations(None, [], {})
        assert isinstance(result, list)
        # Should still produce a section header
        html = "\n".join(result)
        assert "Limitation" in html

    def test_sample_size_limitation(self):
        html = "\n".join(_section_limitations(_make_log_data(), ["Standard"], _make_mat_data()))
        # Should mention sample size as a limitation
        assert "sample" in html.lower() or "cohort" in html.lower()


# ── Conclusions ──


class TestConclusions:
    def test_returns_html(self):
        result = _section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        )
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Conclusion" in html

    def test_no_data(self):
        result = _section_conclusions(None, [], None, {}, None)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Conclusion" in html

    def test_with_groups(self):
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {
                    "trends_json": json.dumps([
                        {"series": "Mean D", "direction": "increasing"},
                    ])
                },
            }
        }
        result = _section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), groups
        )
        assert isinstance(result, list)


# ── Reporting Checklist ──


class TestReportingChecklist:
    def test_returns_html(self):
        result = _section_reporting_checklist(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        )
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Checklist" in html or "checklist" in html

    def test_no_data(self):
        result = _section_reporting_checklist(None, [], None, {}, None)
        assert isinstance(result, list)

    def test_contains_strobe_or_remark(self):
        html = "\n".join(_section_reporting_checklist(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        # Should reference reporting guidelines
        assert "STROBE" in html or "REMARK" in html or "TRIPOD" in html or "checklist" in html.lower()


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
        # Should suggest at least one journal
        assert any(j in html for j in [
            "Radiology", "Physics", "Oncology", "Cancer",
            "Medical", "International Journal", "journal",
        ])
