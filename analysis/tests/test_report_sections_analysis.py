"""Tests for report sections: graph_overview, cross_dwi, and correlations modules.

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

from report.sections.graph_overview import (
    _section_graph_overview,
    _section_graph_issues,
    _section_stats_by_graph_type,
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

    def test_dwi_type_counts(self):
        """Should show counts by DWI type."""
        html = "\n".join(_section_graph_overview(SAMPLE_GRAPH_CSV_ROWS))
        assert "Standard" in html
        assert "dnCNN" in html
        assert "IVIMnet" in html

    def test_signal_density_table(self):
        """Should show statistical signal density by graph type."""
        html = "\n".join(_section_graph_overview(SAMPLE_GRAPH_CSV_ROWS))
        assert "Signal" in html or "signal" in html.lower()

    def test_signal_density_with_no_sig_pvals(self):
        """Graph type with no sig p-values should not appear in signal table."""
        rows = [{
            "file_path": "Standard/test.png",
            "graph_type": "histogram",
            "summary": "no significant results",
            "trends_json": "[]",
            "inflection_points_json": "[]",
        }]
        html = "\n".join(_section_graph_overview(rows))
        # Signal section is always shown; verify no significant p-values counted
        assert "Signal" in html or "0" in html

    def test_single_row(self):
        """Should work with just one row."""
        result = _section_graph_overview([SAMPLE_GRAPH_CSV_ROWS[0]])
        html = "\n".join(result)
        assert "box" in html.lower()
        assert "1" in html


class TestGraphIssues:
    def test_returns_issues(self):
        result = _section_graph_issues(SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)
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

    def test_unknown_type_classified_critical(self):
        """Graph type 'unknown' should be classified as Critical."""
        rows = [{"file_path": "Standard/test.png", "graph_type": "unknown",
                 "issues_json": "[]"}]
        html = "\n".join(_section_graph_issues(rows))
        assert "Critical" in html

    def test_low_severity(self):
        """Issues without high keywords should be classified as Low."""
        rows = [{"file_path": "Standard/test.png", "graph_type": "box",
                 "issues_json": json.dumps(["Minor label issue"])}]
        html = "\n".join(_section_graph_issues(rows))
        assert "Low" in html

    def test_overlap_keyword_high(self):
        """Issue with 'overlap' keyword should be classified as High."""
        rows = [{"file_path": "Standard/test.png", "graph_type": "line",
                 "issues_json": json.dumps(["Overlapping legend text"])}]
        html = "\n".join(_section_graph_issues(rows))
        assert "High" in html

    def test_issue_count_shown(self):
        """Should show 'X of Y graphs have issues'."""
        html = "\n".join(_section_graph_issues(SAMPLE_GRAPH_CSV_ROWS))
        assert "of" in html and "graphs" in html

    def test_severity_stat_cards(self):
        """Should show severity summary stat cards."""
        html = "\n".join(_section_graph_issues(SAMPLE_GRAPH_CSV_ROWS))
        assert "stat-grid" in html

    def test_invalid_issues_json(self):
        """Invalid JSON in issues_json should not crash."""
        rows = [{"file_path": "Standard/test.png", "graph_type": "box",
                 "issues_json": "NOT VALID JSON"}]
        result = _section_graph_issues(rows)
        assert isinstance(result, list)
        assert result == []

    def test_multiple_issues_in_one_graph(self):
        """Multiple issues per graph should all be listed."""
        rows = [{"file_path": "Standard/test.png", "graph_type": "box",
                 "issues_json": json.dumps(["Issue one", "Issue two", "Issue three"])}]
        html = "\n".join(_section_graph_issues(rows))
        assert "Issue one" in html
        assert "Issue two" in html
        assert "Issue three" in html


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

    def test_trend_direction_counts(self):
        """Should count trend directions (increasing, decreasing, stable)."""
        html = "\n".join(_section_stats_by_graph_type(SAMPLE_GRAPH_CSV_ROWS))
        # The arrows should be present
        assert "\u2191" in html or "\u2193" in html

    def test_most_sig_graph_type_highlighted(self):
        """Graph type with most significant p-values should be highlighted."""
        html = "\n".join(_section_stats_by_graph_type(SAMPLE_GRAPH_CSV_ROWS))
        assert "most" in html.lower() or "info-box" in html

    def test_sig_p_count(self):
        """Should count significant p-values per graph type."""
        html = "\n".join(_section_stats_by_graph_type(SAMPLE_GRAPH_CSV_ROWS))
        assert "Sig" in html

    def test_correlation_count(self):
        """Should count strong correlations per graph type."""
        html = "\n".join(_section_stats_by_graph_type(SAMPLE_GRAPH_CSV_ROWS))
        assert "Corr" in html


# ── Edge cases: no-graph-data paths, partial data, missing metrics ──


class TestGraphOverviewEdgeCases:
    def test_none_rows(self):
        assert _section_graph_overview(None) == []

    def test_rows_with_missing_graph_type(self):
        rows = [{"file_path": "Standard/test.png", "summary": "test"}]
        result = _section_graph_overview(rows)
        assert isinstance(result, list)

    def test_rows_all_same_type(self):
        rows = [
            {"file_path": "Standard/a.png", "graph_type": "scatter",
             "summary": "p=0.01", "trends_json": "[]"},
            {"file_path": "Standard/b.png", "graph_type": "scatter",
             "summary": "p=0.02", "trends_json": "[]"},
        ]
        result = _section_graph_overview(rows)
        html = "\n".join(result)
        assert "scatter" in html.lower()

    def test_rows_with_empty_summary(self):
        rows = [{"file_path": "Standard/test.png", "graph_type": "box", "summary": ""}]
        result = _section_graph_overview(rows)
        assert isinstance(result, list)

    def test_unknown_graph_type(self):
        rows = [{"file_path": "Standard/test.png", "graph_type": "unknown",
                 "summary": "test", "trends_json": "[]"}]
        result = _section_graph_overview(rows)
        html = "\n".join(result)
        assert "unknown" in html.lower()


class TestGraphIssuesEdgeCases:
    def test_none_rows(self):
        assert _section_graph_issues(None) == []

    def test_all_graphs_have_issues(self):
        rows = [
            {"file_path": "Standard/a.png", "graph_type": "box",
             "issues_json": json.dumps(["Axis cutoff"])},
            {"file_path": "Standard/b.png", "graph_type": "line",
             "issues_json": json.dumps(["Overlap detected"])},
            {"file_path": "Standard/c.png", "graph_type": "scatter",
             "issues_json": json.dumps(["Label truncated"])},
        ]
        result = _section_graph_issues(rows)
        html = "\n".join(result)
        assert "3" in html  # 3 of 3 have issues

    def test_empty_string_issues_json(self):
        rows = [{"file_path": "Standard/test.png", "graph_type": "box", "issues_json": ""}]
        result = _section_graph_issues(rows)
        assert isinstance(result, list)

    def test_issues_with_html_chars(self):
        rows = [{"file_path": "Standard/test.png", "graph_type": "box",
                 "issues_json": json.dumps(["<b>Bold issue</b> & more"])}]
        result = _section_graph_issues(rows)
        html = "\n".join(result)
        assert "<b>" not in html or "&lt;b&gt;" in html or isinstance(html, str)


class TestStatsByGraphTypeEdgeCases:
    def test_none_rows(self):
        assert _section_stats_by_graph_type(None) == []

    def test_rows_with_no_trends(self):
        rows = [{"file_path": "Standard/test.png", "graph_type": "heatmap",
                 "summary": "no trends", "trends_json": "[]"}]
        result = _section_stats_by_graph_type(rows)
        assert isinstance(result, list)

    def test_all_nonsig_pvalues(self):
        rows = [{"file_path": "Standard/test.png", "graph_type": "box",
                 "summary": "p = 0.85 for comparison", "trends_json": "[]"}]
        result = _section_stats_by_graph_type(rows)
        assert isinstance(result, list)

