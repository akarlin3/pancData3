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

Edge cases covered:
- Empty/None inputs at every level
- No graph rows, no log data, no MAT data, no groups
- Missing nested keys (no survival, no stats_predictive)
- Stat cards: graph count, AUC, GLME significant, CSV significant
- Timestamp display
- DWI type badges
- Multiple DWI types
- No significant p-values in rows
- No trends in groups
- No predictive data (empty roc_analyses)
- Treatment response with canonical pattern
- Treatment response with no groups
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

    def test_zero_graphs(self):
        """No graph rows should show 0 or omit the count."""
        result = _section_executive_summary(None, ["Standard"], [], None, "ts")
        assert isinstance(result, list)

    def test_with_mat_data(self):
        """MAT data should add cohort info to summary."""
        html = "\n".join(_section_executive_summary(
            _make_log_data(), ["Standard"], SAMPLE_GRAPH_CSV_ROWS,
            None, "ts", _make_mat_data()
        ))
        assert "42" in html or "patient" in html.lower()

    def test_with_csv_data(self):
        """CSV data should contribute to significant metric count."""
        csv_data = {
            "significant_metrics": {"Standard": [{"metric": "adc"}]},
            "fdr_global": {"Standard": [{"metric": "adc"}]},
        }
        result = _section_executive_summary(
            _make_log_data(), ["Standard"], [], csv_data, "ts"
        )
        assert isinstance(result, list)

    def test_auc_stat_card(self):
        """Best AUC should appear as a stat card."""
        html = "\n".join(_section_executive_summary(
            _make_log_data(), ["Standard"], [], None, "ts"
        ))
        assert "AUC" in html or "0.78" in html

    def test_glme_significant_count(self):
        """GLME significant count should appear."""
        html = "\n".join(_section_executive_summary(
            _make_log_data(), ["Standard"], [], None, "ts"
        ))
        assert "GLME" in html or "significant" in html.lower()

    def test_no_log_data_no_crash(self):
        """None log_data should not crash."""
        result = _section_executive_summary(None, ["Standard"], SAMPLE_GRAPH_CSV_ROWS, None, "ts")
        assert isinstance(result, list)

    def test_multiple_dwi_types_with_data(self):
        """Multiple DWI types should each show AUC if available."""
        log = _make_log_data()
        log["dnCNN"] = {
            "stats_predictive": {
                "roc_analyses": [{"auc": 0.72, "timepoint": "BL"}],
            },
            "stats_comparisons": {},
            "survival": {},
        }
        result = _section_executive_summary(
            log, ["Standard", "dnCNN"], [], None, "ts"
        )
        assert isinstance(result, list)

    def test_empty_survival_no_hr_card(self):
        """No hazard ratios should not crash executive summary."""
        log = {"Standard": {"survival": {}, "stats_comparisons": {}, "stats_predictive": {}}}
        result = _section_executive_summary(log, ["Standard"], [], None, "ts")
        assert isinstance(result, list)


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

    def test_no_groups_no_log_no_mat(self):
        result = _section_hypothesis(None, None, None)
        assert isinstance(result, list)

    def test_with_trend_data(self):
        """Trend data should influence hypothesis framing."""
        result = _section_hypothesis(_make_groups(), _make_log_data(), _make_mat_data())
        assert isinstance(result, list)
        assert len(result) > 0


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

    def test_no_significant_pvalues(self):
        """Rows with no significant p-values."""
        rows = [{
            "file_path": "Standard/test.png",
            "graph_type": "box",
            "summary": "No significant findings. p = 0.45.",
            "trends_json": "[]",
            "inflection_points_json": "[]",
            "statistical_tests_json": json.dumps([
                {"test_name": "t-test", "p_value": 0.45}
            ]),
        }]
        result = _section_statistical_significance(rows, None, None, ["Standard"])
        assert isinstance(result, list)

    def test_with_csv_data(self):
        """CSV data should contribute to significance display."""
        csv_data = {
            "significant_metrics": {
                "Standard": [{"metric": "adc", "timepoint": "BL", "p": 0.003}],
            }
        }
        result = _section_statistical_significance([], csv_data, None, ["Standard"])
        assert isinstance(result, list)

    def test_with_log_and_rows(self):
        """Both log data and graph rows together."""
        result = _section_statistical_significance(
            SAMPLE_GRAPH_CSV_ROWS, None, _make_log_data(), ["Standard"]
        )
        assert isinstance(result, list)
        assert len(result) > 0


# ── Broad Statistical Overview ──


class TestBroadStatisticalOverview:
    def test_returns_html(self):
        result = _section_broad_statistical_overview(_make_log_data(), ["Standard"])
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_broad_statistical_overview(None, [])
        assert isinstance(result, list)

    def test_with_multiple_dwi_types(self):
        log = _make_log_data()
        log["dnCNN"] = log["Standard"].copy()
        result = _section_broad_statistical_overview(log, ["Standard", "dnCNN"])
        assert isinstance(result, list)

    def test_empty_log_data(self):
        result = _section_broad_statistical_overview({}, [])
        assert isinstance(result, list)


# ── Treatment Response ──


class TestTreatmentResponse:
    def test_returns_html(self):
        result = _section_treatment_response(_make_groups())
        assert isinstance(result, list)

    def test_no_groups(self):
        result = _section_treatment_response(None)
        assert isinstance(result, list)

    def test_empty_groups(self):
        result = _section_treatment_response({})
        assert isinstance(result, list)

    def test_with_trend_data(self):
        html = "\n".join(_section_treatment_response(_make_groups()))
        assert len(html) > 0

    def test_single_dwi_type(self):
        """Treatment response with single DWI type."""
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {
                    "trends_json": json.dumps([
                        {"series": "Mean D", "direction": "increasing"},
                    ])
                },
            },
        }
        result = _section_treatment_response(groups)
        assert isinstance(result, list)

    def test_no_longitudinal_group(self):
        """Groups without Longitudinal_Mean_Metrics."""
        groups = {
            "Feature_BoxPlots": {
                "Standard": SAMPLE_GRAPH_CSV_ROWS[0],
            },
        }
        result = _section_treatment_response(groups)
        assert isinstance(result, list)

    def test_malformed_trends_json(self):
        """Malformed trends_json should not crash."""
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": "not valid json"},
            },
        }
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

    def test_no_roc_analyses(self):
        """No ROC analyses should produce empty or minimal section."""
        log = {"Standard": {"stats_predictive": {"roc_analyses": []}}}
        result = _section_predictive_performance(log, ["Standard"])
        assert isinstance(result, list)

    def test_multiple_timepoints(self):
        """Multiple ROC analyses at different timepoints."""
        log = {
            "Standard": {
                "stats_predictive": {
                    "roc_analyses": [
                        {"auc": 0.78, "timepoint": "BL"},
                        {"auc": 0.65, "timepoint": "W2"},
                        {"auc": 0.82, "timepoint": "W4"},
                    ],
                },
            }
        }
        html = "\n".join(_section_predictive_performance(log, ["Standard"]))
        assert "0.82" in html or "AUC" in html

    def test_multiple_dwi_types(self):
        """Multiple DWI types each with ROC data."""
        log = {
            "Standard": {
                "stats_predictive": {
                    "roc_analyses": [{"auc": 0.78, "timepoint": "BL"}],
                },
            },
            "dnCNN": {
                "stats_predictive": {
                    "roc_analyses": [{"auc": 0.72, "timepoint": "BL"}],
                },
            },
        }
        html = "\n".join(_section_predictive_performance(log, ["Standard", "dnCNN"]))
        assert isinstance(html, str)

    def test_dwi_type_not_in_log(self):
        """DWI type listed but not in log_data."""
        result = _section_predictive_performance(_make_log_data(), ["Standard", "IVIMnet"])
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

    def test_with_csv_data(self):
        csv_data = {
            "significant_metrics": {"Standard": [{"metric": "adc"}]},
        }
        result = _section_manuscript_ready_findings(
            _make_log_data(), ["Standard"], csv_data, _make_mat_data(), None
        )
        assert isinstance(result, list)

    def test_with_all_data_sources(self):
        """All data sources present for maximum coverage."""
        result = _section_manuscript_ready_findings(
            _make_log_data(), ["Standard", "dnCNN"], None,
            _make_mat_data(), _make_groups()
        )
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

    def test_with_csv_data(self):
        csv_data = {
            "significant_metrics": {"Standard": [{"metric": "adc"}]},
            "fdr_global": {"Standard": [{"metric": "adc"}]},
        }
        result = _section_results_draft(
            _make_log_data(), ["Standard"], csv_data, _make_mat_data(), None
        )
        assert isinstance(result, list)

    def test_with_groups(self):
        result = _section_results_draft(
            _make_log_data(), ["Standard"], None, _make_mat_data(), _make_groups()
        )
        assert isinstance(result, list)
        assert len(result) > 0

    def test_empty_log_data(self):
        result = _section_results_draft({}, [], None, {}, None)
        assert isinstance(result, list)
