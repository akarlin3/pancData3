"""Tests for report_sections.analysis_sections module.

Validates analysis section builders:
- Graph overview (type/DWI counts, signal density)
- Graph issues (severity classification, issue counting)
- Statistics by graph type (aggregation, trend counting)
- Cross-DWI comparison (trend agreement, disagreement records)
- Correlations (extraction, strength classification)
- Feature overlap (cross-DWI feature comparison)
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

from report.sections.analysis_sections import (
    _section_graph_overview,
    _section_graph_issues,
    _section_stats_by_graph_type,
    _section_cross_dwi_comparison,
    _section_correlations,
    _section_feature_overlap,
)


class TestGraphOverview:
    def test_returns_html_list(self):
        result = _section_graph_overview(SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)
        assert len(result) > 0

    def test_empty_rows_returns_empty(self):
        assert _section_graph_overview([]) == []

    def test_contains_type_counts(self):
        html = "\n".join(_section_graph_overview(SAMPLE_GRAPH_CSV_ROWS))
        assert "box" in html.lower()
        assert "line" in html.lower()

    def test_contains_graph_overview_header(self):
        html = "\n".join(_section_graph_overview(SAMPLE_GRAPH_CSV_ROWS))
        assert "Graph Analysis Overview" in html


class TestGraphIssues:
    def test_returns_issues(self):
        result = _section_graph_issues(SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)
        # At least one row has issues
        html = "\n".join(result)
        assert "Graph Issues" in html

    def test_empty_rows(self):
        assert _section_graph_issues([]) == []

    def test_no_issues_returns_empty(self):
        rows = [{"file_path": "Standard/test.png", "graph_type": "line",
                 "issues_json": "[]"}]
        result = _section_graph_issues(rows)
        assert result == []

    def test_error_type_classified_critical(self):
        rows = [{"file_path": "Standard/test.png", "graph_type": "error",
                 "issues_json": "[]"}]
        result = _section_graph_issues(rows)
        html = "\n".join(result)
        assert "Critical" in html

    def test_severity_classification(self):
        rows = [
            {"file_path": "Standard/a.png", "graph_type": "box",
             "issues_json": json.dumps(["Axis cutoff detected"])},
        ]
        result = _section_graph_issues(rows)
        html = "\n".join(result)
        assert "High" in html


class TestStatsByGraphType:
    def test_returns_html(self):
        result = _section_stats_by_graph_type(SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)
        assert len(result) > 0

    def test_empty_rows(self):
        assert _section_stats_by_graph_type([]) == []

    def test_contains_graph_type_rows(self):
        html = "\n".join(_section_stats_by_graph_type(SAMPLE_GRAPH_CSV_ROWS))
        assert "box" in html.lower()


class TestCrossDwiComparison:
    def test_with_groups(self):
        groups = {
            "Feature_BoxPlots": {
                "Standard": SAMPLE_GRAPH_CSV_ROWS[0],
                "dnCNN": SAMPLE_GRAPH_CSV_ROWS[1],
            }
        }
        result = _section_cross_dwi_comparison(groups, None)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Cross-DWI" in html

    def test_empty_groups(self):
        result = _section_cross_dwi_comparison(None, None)
        assert result == []

    def test_agreement_detection(self):
        groups = {
            "TestGraph": {
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
            }
        }
        result = _section_cross_dwi_comparison(groups, None)
        html = "\n".join(result)
        assert "AGREE" in html

    def test_disagreement_detection(self):
        groups = {
            "TestGraph": {
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "down"}])},
            }
        }
        result = _section_cross_dwi_comparison(groups, None)
        html = "\n".join(result)
        assert "DIFFER" in html

    def test_csv_cross_reference(self):
        csv_data = {
            "cross_reference": [
                {"metric": "adc", "timepoint": "BL", "consistent": False,
                 "significant_in": ["Standard"], "not_significant_in": ["dnCNN"]},
            ]
        }
        result = _section_cross_dwi_comparison({"G": {"Standard": {}, "dnCNN": {}}}, csv_data)
        html = "\n".join(result)
        assert "Inconsistencies" in html


class TestCorrelations:
    def test_finds_correlations(self):
        result = _section_correlations(SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)
        html = "\n".join(result)
        # Row 1 summary has r = 0.65
        assert "Correlations" in html

    def test_empty_rows(self):
        assert _section_correlations([]) == []

    def test_strong_vs_moderate_split(self):
        # Row 1 has r=0.65 (strong)
        html = "\n".join(_section_correlations(SAMPLE_GRAPH_CSV_ROWS))
        assert "Strong" in html


class TestFeatureOverlap:
    def test_with_overlap(self):
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["ADC_BL", "D_BL"]},
                    ]
                }
            },
            "dnCNN": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["ADC_BL", "f_BL"]},
                    ]
                }
            },
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Feature Overlap" in html

    def test_single_dwi_type(self):
        result = _section_feature_overlap({"Standard": {}}, ["Standard"])
        assert result == []

    def test_no_log_data(self):
        result = _section_feature_overlap(None, [])
        assert result == []
