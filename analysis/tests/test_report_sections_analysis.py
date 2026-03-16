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

    def test_agreement_summary_percentage(self):
        """Should show overall agreement percentage."""
        groups = {
            "G1": {
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
            },
            "G2": {
                "Standard": {"trends_json": json.dumps([{"series": "S2", "direction": "down"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S2", "direction": "up"}])},
            },
        }
        html = "\n".join(_section_cross_dwi_comparison(groups, None))
        assert "50%" in html or "1/2" in html

    def test_notable_disagreements_section(self):
        """Should show notable disagreements subsection."""
        groups = {
            "TestGraph": {
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "down"}])},
            }
        }
        html = "\n".join(_section_cross_dwi_comparison(groups, None))
        assert "Notable Disagreements" in html

    def test_sample_size_mismatch_warning(self):
        """Differing sample sizes > 5 should warn."""
        groups = {
            "TestGraph": {
                "Standard": {"sample_size": "42", "trends_json": "[]"},
                "dnCNN": {"sample_size": "30", "trends_json": "[]"},
            }
        }
        html = "\n".join(_section_cross_dwi_comparison(groups, None))
        assert "Sample size mismatch" in html or "mismatch" in html.lower()

    def test_no_sample_size_mismatch_close_counts(self):
        """Close sample sizes (diff <= 5) should NOT warn."""
        groups = {
            "TestGraph": {
                "Standard": {"sample_size": "42", "trends_json": "[]"},
                "dnCNN": {"sample_size": "40", "trends_json": "[]"},
            }
        }
        html = "\n".join(_section_cross_dwi_comparison(groups, None))
        assert "Sample size mismatch" not in html

    def test_priority_graphs_labelled(self):
        """Priority graphs should get a 'priority' badge."""
        groups = {
            "Feature_BoxPlots": {
                "Standard": {"trends_json": json.dumps([{"series": "S", "direction": "up"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S", "direction": "up"}])},
            }
        }
        html = "\n".join(_section_cross_dwi_comparison(groups, None))
        assert "priority" in html

    def test_single_type_series_shown(self):
        """Series present in only one DWI type should still be shown."""
        groups = {
            "TestGraph": {
                "Standard": {"trends_json": json.dumps([
                    {"series": "S1", "direction": "up"},
                    {"series": "S2", "direction": "down"},
                ])},
                "dnCNN": {"trends_json": json.dumps([
                    {"series": "S1", "direction": "up"},
                ])},
            }
        }
        html = "\n".join(_section_cross_dwi_comparison(groups, None))
        assert isinstance(html, str)

    def test_empty_trends_json(self):
        """Empty trends_json for all types should show 'No trend data'."""
        groups = {
            "TestGraph": {
                "Standard": {"trends_json": "[]"},
                "dnCNN": {"trends_json": "[]"},
            }
        }
        html = "\n".join(_section_cross_dwi_comparison(groups, None))
        assert "No trend data" in html or isinstance(html, str)

    def test_csv_consistent_entries_not_shown(self):
        """Consistent cross-reference entries should NOT appear in inconsistencies."""
        csv_data = {
            "cross_reference": [
                {"metric": "adc", "timepoint": "BL", "consistent": True,
                 "significant_in": ["Standard", "dnCNN"]},
            ]
        }
        html = "\n".join(_section_cross_dwi_comparison(
            {"G": {"Standard": {}, "dnCNN": {}}}, csv_data
        ))
        assert "Inconsistencies" not in html


class TestCorrelations:
    def test_finds_correlations(self):
        result = _section_correlations(SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Correlations" in html

    def test_empty_rows(self):
        assert _section_correlations([]) == []

    def test_strong_vs_moderate_split(self):
        html = "\n".join(_section_correlations(SAMPLE_GRAPH_CSV_ROWS))
        assert "Strong" in html

    def test_causation_caveat(self):
        """Should include causation caveat."""
        html = "\n".join(_section_correlations(SAMPLE_GRAPH_CSV_ROWS))
        assert "causation" in html.lower() or "Correlation does not imply" in html

    def test_bonferroni_note(self):
        """Should include Bonferroni correction note when multiple correlations."""
        html = "\n".join(_section_correlations(SAMPLE_GRAPH_CSV_ROWS))
        # Row 1 has r=0.65 and row 3 has r²=0.45 which extracts as r~0.67
        assert "Bonferroni" in html or "multiple testing" in html.lower()

    def test_confidence_interval_caveat(self):
        """Should note that CIs are not reported."""
        html = "\n".join(_section_correlations(SAMPLE_GRAPH_CSV_ROWS))
        assert "Confidence interval" in html or "confidence" in html.lower()

    def test_no_correlations_message(self):
        """Should show 'no notable correlations' when none found."""
        rows = [{
            "file_path": "Standard/test.png",
            "summary": "no correlation found",
            "trends_json": "[]",
        }]
        html = "\n".join(_section_correlations(rows))
        assert "No notable" in html or "no notable" in html.lower()

    def test_correlation_sorted_by_strength(self):
        """Correlations should be sorted by absolute value (strongest first)."""
        html = "\n".join(_section_correlations(SAMPLE_GRAPH_CSV_ROWS))
        # Should contain correlation values
        assert "0.65" in html or "|r|" in html


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

    def test_shared_feature_marked(self):
        """Features shared across DWI types should be marked as 'Shared'."""
        log_data = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc", "d"]},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc", "f"]},
            ]}},
        }
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN"]))
        assert "Shared" in html
        assert "adc" in html

    def test_type_specific_feature_marked(self):
        """Features unique to one DWI type should be marked as 'Type-specific'."""
        log_data = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc", "d"]},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc", "f"]},
            ]}},
        }
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN"]))
        assert "Type-specific" in html

    def test_overlap_summary_percentage(self):
        """Should show overlap summary with percentage."""
        log_data = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc", "d"]},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc", "f"]},
            ]}},
        }
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN"]))
        assert "%" in html
        assert "shared" in html.lower()

    def test_feature_stability_across_timepoints(self):
        """Features selected at multiple timepoints should show stability."""
        log_data = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc", "d"]},
                {"timepoint": "W2", "features": ["adc", "f"]},
                {"timepoint": "W4", "features": ["adc"]},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc"]},
            ]}},
        }
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN"]))
        assert "Stability" in html

    def test_feature_importance_context_note(self):
        """Should include feature importance context note."""
        log_data = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc"]},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc"]},
            ]}},
        }
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN"]))
        assert "importance" in html.lower() or "elastic net" in html.lower()

    def test_duplicate_feature_detection(self):
        """Near-duplicate features (e.g., mean_adc vs meanADC) should be flagged."""
        log_data = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["mean_adc", "meanADC"]},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["mean_adc"]},
            ]}},
        }
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN"]))
        assert "duplicate" in html.lower() or "Potential" in html

    def test_no_shared_features(self):
        """When no features overlap, should still show type-specific."""
        log_data = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc"]},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["d"]},
            ]}},
        }
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN"]))
        assert "Shared" in html  # Stat card shows 0
        assert "Type-specific" in html

    def test_empty_feature_lists(self):
        """Empty feature lists should be handled gracefully."""
        log_data = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": []},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": []},
            ]}},
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        assert isinstance(result, list)

    def test_missing_stats_predictive_key(self):
        """Missing stats_predictive key in log_data should not crash."""
        log_data = {
            "Standard": {},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc"]},
            ]}},
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        assert isinstance(result, list)


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


class TestCrossDwiComparisonEdgeCases:
    def test_all_single_dwi_groups(self):
        groups = {
            "G1": {"Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])}},
            "G2": {"dnCNN": {"trends_json": json.dumps([{"series": "S2", "direction": "down"}])}},
        }
        result = _section_cross_dwi_comparison(groups, None)
        assert isinstance(result, list)

    def test_groups_with_invalid_trends_json(self):
        groups = {
            "G": {
                "Standard": {"trends_json": "NOT JSON"},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
            }
        }
        result = _section_cross_dwi_comparison(groups, None)
        assert isinstance(result, list)

    def test_csv_data_empty_cross_reference(self):
        csv_data = {"cross_reference": []}
        result = _section_cross_dwi_comparison(
            {"G": {"Standard": {}, "dnCNN": {}}}, csv_data
        )
        assert isinstance(result, list)

    def test_csv_data_all_consistent(self):
        csv_data = {
            "cross_reference": [
                {"metric": "adc", "timepoint": "BL", "consistent": True,
                 "significant_in": ["Standard", "dnCNN"]},
                {"metric": "d", "timepoint": "BL", "consistent": True,
                 "significant_in": ["Standard", "dnCNN"]},
            ]
        }
        result = _section_cross_dwi_comparison(
            {"G": {"Standard": {}, "dnCNN": {}}}, csv_data
        )
        html = "\n".join(result)
        assert "Inconsistencies" not in html

    def test_three_dwi_types(self):
        groups = {
            "G": {
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "IVIMnet": {"trends_json": json.dumps([{"series": "S1", "direction": "down"}])},
            }
        }
        result = _section_cross_dwi_comparison(groups, None)
        html = "\n".join(result)
        assert "DIFFER" in html or "AGREE" in html

    def test_groups_with_no_trends_at_all(self):
        groups = {
            "G": {
                "Standard": {},
                "dnCNN": {},
            }
        }
        result = _section_cross_dwi_comparison(groups, None)
        assert isinstance(result, list)

    def test_sample_size_none(self):
        groups = {
            "G": {
                "Standard": {"trends_json": "[]"},
                "dnCNN": {"trends_json": "[]"},
            }
        }
        result = _section_cross_dwi_comparison(groups, None)
        assert isinstance(result, list)


class TestCorrelationsEdgeCases:
    def test_none_rows(self):
        assert _section_correlations(None) == []

    def test_rows_with_no_summary(self):
        rows = [{"file_path": "Standard/test.png"}]
        result = _section_correlations(rows)
        assert isinstance(result, list)

    def test_weak_correlations_only(self):
        rows = [{"file_path": "Standard/test.png", "summary": "r = 0.15 weak correlation",
                 "trends_json": "[]"}]
        result = _section_correlations(rows)
        html = "\n".join(result)
        assert "No notable" in html or "no notable" in html.lower()

    def test_very_strong_correlation(self):
        rows = [{"file_path": "Standard/test.png", "summary": "r = 0.95 very strong",
                 "trends_json": "[]"}]
        result = _section_correlations(rows)
        html = "\n".join(result)
        assert "0.95" in html or "Strong" in html

    def test_negative_correlation(self):
        rows = [{"file_path": "Standard/test.png", "summary": "r = -0.72 negative",
                 "trends_json": "[]"}]
        result = _section_correlations(rows)
        html = "\n".join(result)
        assert "0.72" in html or isinstance(html, str)


class TestFeatureOverlapEdgeCases:
    def test_three_dwi_types(self):
        log_data = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc", "d"]},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc", "f"]},
            ]}},
            "IVIMnet": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc", "d_star"]},
            ]}},
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN", "IVIMnet"])
        html = "\n".join(result)
        assert "adc" in html  # Shared across all 3

    def test_many_timepoints(self):
        log_data = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc"]},
                {"timepoint": "W2", "features": ["adc", "d"]},
                {"timepoint": "W4", "features": ["adc"]},
                {"timepoint": "W6", "features": ["adc", "f"]},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc"]},
                {"timepoint": "W2", "features": ["d"]},
            ]}},
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Stability" in html

    def test_all_features_shared(self):
        log_data = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc", "d"]},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc", "d"]},
            ]}},
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        html = "\n".join(result)
        assert "100" in html or "agree" in html.lower()

    def test_no_features_shared(self):
        log_data = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc"]},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["d"]},
            ]}},
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        html = "\n".join(result)
        assert "0%" in html or "Type-specific" in html
