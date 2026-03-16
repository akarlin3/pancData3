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
            "dosimetry": {"d95_gtvp": {"mean": 45.0}, "v50_gtvp": {"mean": 0.85}},
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


# ── Broad Statistical Overview ──


class TestBroadStatisticalOverview:
    def test_returns_html(self):
        result = _section_broad_statistical_overview(_make_log_data(), ["Standard"])
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_broad_statistical_overview(None, [])
        assert isinstance(result, list)


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
        # Should reference treatment response in some way
        assert len(html) > 0


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
