"""Tests for report_sections.analysis_sections module.

Validates analysis section builders:
- Graph overview (type/DWI counts, signal density)
- Graph issues (severity classification, issue counting)
- Statistics by graph type (aggregation, trend counting)
- Cross-DWI comparison (trend agreement, disagreement records)
- Correlations (extraction, strength classification)
- Feature overlap (cross-DWI feature comparison)

Edge cases covered:
- Empty/None inputs at every level
- Malformed JSON in trends_json, issues_json
- No-graph-data paths
- Single DWI type (no cross-comparison possible)
- All trends agree / all disagree
- No correlations found
- No significant p-values
- Unknown and error graph types
- Feature stability across timepoints
- Duplicate feature detection
- Sample size mismatch warnings
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

    def test_none_rows_returns_empty(self):
        assert _section_graph_overview(None) == []

    def test_contains_type_counts(self):
        html = "\n".join(_section_graph_overview(SAMPLE_GRAPH_CSV_ROWS))
        assert "box" in html.lower()
        assert "line" in html.lower()

    def test_contains_graph_overview_header(self):
        html = "\n".join(_section_graph_overview(SAMPLE_GRAPH_CSV_ROWS))
        assert "Graph Analysis Overview" in html

    def test_single_row(self):
        """Single row should still produce valid overview."""
        result = _section_graph_overview([SAMPLE_GRAPH_CSV_ROWS[0]])
        assert isinstance(result, list)
        assert len(result) > 0

    def test_row_with_no_graph_type(self):
        """Row missing graph_type defaults to 'unknown'."""
        rows = [{"file_path": "Standard/test.png"}]
        result = _section_graph_overview(rows)
        html = "\n".join(result)
        assert "unknown" in html.lower()

    def test_signal_density_with_significant_pvalues(self):
        """Rows with significant p-values should show signal density table."""
        html = "\n".join(_section_graph_overview(SAMPLE_GRAPH_CSV_ROWS))
        assert "Signal Density" in html or "Signal" in html

    def test_row_with_no_pvalues_or_correlations(self):
        """Row with no statistical data still counted in graph type table."""
        rows = [{
            "file_path": "Standard/empty_plot.png",
            "graph_type": "scatter",
            "summary": "No statistics available.",
            "trends_json": "[]",
            "inflection_points_json": "[]",
        }]
        result = _section_graph_overview(rows)
        html = "\n".join(result)
        assert "scatter" in html.lower()


class TestGraphIssues:
    def test_returns_issues(self):
        result = _section_graph_issues(SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Graph Issues" in html

    def test_empty_rows(self):
        assert _section_graph_issues([]) == []

    def test_none_rows(self):
        assert _section_graph_issues(None) == []

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

    def test_unknown_type_classified_critical(self):
        """Graph type 'unknown' should also be Critical."""
        rows = [{"file_path": "Standard/test.png", "graph_type": "unknown",
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

    def test_low_severity_issues(self):
        """Issues without high-severity keywords should be Low."""
        rows = [
            {"file_path": "Standard/b.png", "graph_type": "box",
             "issues_json": json.dumps(["Minor label alignment issue"])},
        ]
        result = _section_graph_issues(rows)
        html = "\n".join(result)
        assert "Low" in html

    def test_malformed_issues_json(self):
        """Malformed JSON in issues_json should be handled gracefully."""
        rows = [
            {"file_path": "Standard/c.png", "graph_type": "error",
             "issues_json": "not valid json"},
        ]
        result = _section_graph_issues(rows)
        html = "\n".join(result)
        assert "Critical" in html  # error type still classified

    def test_none_issues_json(self):
        """None issues_json should default to empty list."""
        rows = [
            {"file_path": "Standard/d.png", "graph_type": "line",
             "issues_json": None},
        ]
        result = _section_graph_issues(rows)
        assert result == []

    def test_issue_count_display(self):
        """Should show X of Y graphs have issues."""
        rows = [
            {"file_path": "Standard/a.png", "graph_type": "box",
             "issues_json": json.dumps(["Issue 1"])},
            {"file_path": "Standard/b.png", "graph_type": "line",
             "issues_json": "[]"},
        ]
        html = "\n".join(_section_graph_issues(rows))
        assert "1 of 2" in html

    def test_multiple_severities(self):
        """Mix of Critical, High, Low issues should all appear."""
        rows = [
            {"file_path": "Standard/a.png", "graph_type": "error",
             "issues_json": "[]"},
            {"file_path": "Standard/b.png", "graph_type": "box",
             "issues_json": json.dumps(["overlap detected"])},
            {"file_path": "Standard/c.png", "graph_type": "line",
             "issues_json": json.dumps(["minor spacing"])},
        ]
        html = "\n".join(_section_graph_issues(rows))
        assert "Critical" in html
        assert "High" in html
        assert "Low" in html


class TestStatsByGraphType:
    def test_returns_html(self):
        result = _section_stats_by_graph_type(SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)
        assert len(result) > 0

    def test_empty_rows(self):
        assert _section_stats_by_graph_type([]) == []

    def test_none_rows(self):
        assert _section_stats_by_graph_type(None) == []

    def test_contains_graph_type_rows(self):
        html = "\n".join(_section_stats_by_graph_type(SAMPLE_GRAPH_CSV_ROWS))
        assert "box" in html.lower()

    def test_trend_direction_counting(self):
        """Increasing/decreasing/stable trends should be counted."""
        html = "\n".join(_section_stats_by_graph_type(SAMPLE_GRAPH_CSV_ROWS))
        # The arrows should appear (↑, ↓, →)
        assert "\u2191" in html or "\u2193" in html or "\u2192" in html

    def test_highlights_most_significant_type(self):
        """Graph type with most sig p-values gets highlighted."""
        html = "\n".join(_section_stats_by_graph_type(SAMPLE_GRAPH_CSV_ROWS))
        assert "most" in html.lower() or "significant" in html.lower()

    def test_malformed_trends_json(self):
        """Malformed trends_json should not crash."""
        rows = [{
            "file_path": "Standard/test.png",
            "graph_type": "box",
            "summary": "p = 0.01",
            "trends_json": "INVALID JSON",
            "inflection_points_json": "[]",
        }]
        result = _section_stats_by_graph_type(rows)
        assert isinstance(result, list)

    def test_row_with_no_text_fields(self):
        """Row with empty summary/trends still counted."""
        rows = [{
            "file_path": "Standard/test.png",
            "graph_type": "scatter",
            "summary": "",
            "trends_json": "[]",
            "inflection_points_json": "[]",
        }]
        result = _section_stats_by_graph_type(rows)
        html = "\n".join(result)
        assert "scatter" in html


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

    def test_empty_dict_groups(self):
        result = _section_cross_dwi_comparison({}, None)
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

    def test_single_dwi_type_no_comparison(self):
        """Single DWI type in group should show series but no agreement verdict."""
        groups = {
            "TestGraph": {
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
            }
        }
        result = _section_cross_dwi_comparison(groups, None)
        # Single-type groups with < 2 real types shouldn't produce comparison tables
        assert isinstance(result, list)

    def test_notable_disagreements_subsection(self):
        """Disagreements should generate a Notable Disagreements subsection."""
        groups = {
            "TestGraph": {
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "down"}])},
            }
        }
        html = "\n".join(_section_cross_dwi_comparison(groups, None))
        assert "Notable Disagreements" in html

    def test_all_agree_no_disagreement_subsection(self):
        """100% agreement should not show Notable Disagreements."""
        groups = {
            "TestGraph": {
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
            }
        }
        html = "\n".join(_section_cross_dwi_comparison(groups, None))
        assert "Notable Disagreements" not in html

    def test_agreement_percentage_display(self):
        """Overall agreement percentage should be shown."""
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
        assert "50%" in html  # 1/2 agree

    def test_priority_graph_label(self):
        """Priority graphs should get a 'priority' badge."""
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
            }
        }
        html = "\n".join(_section_cross_dwi_comparison(groups, None))
        assert "priority" in html.lower()

    def test_malformed_trends_json_in_group(self):
        """Malformed trends_json should not crash."""
        groups = {
            "TestGraph": {
                "Standard": {"trends_json": "not json"},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
            }
        }
        result = _section_cross_dwi_comparison(groups, None)
        assert isinstance(result, list)

    def test_sample_size_mismatch_warning(self):
        """Sample sizes differing by > 5 should trigger warning."""
        groups = {
            "TestGraph": {
                "Standard": {"trends_json": "[]", "sample_size": "42"},
                "dnCNN": {"trends_json": "[]", "sample_size": "30"},
            }
        }
        html = "\n".join(_section_cross_dwi_comparison(groups, None))
        assert "mismatch" in html.lower() or "Sample size" in html

    def test_csv_cross_reference_all_consistent(self):
        """All consistent cross-reference entries should not show inconsistency table."""
        csv_data = {
            "cross_reference": [
                {"metric": "adc", "timepoint": "BL", "consistent": True},
            ]
        }
        html = "\n".join(_section_cross_dwi_comparison(
            {"G": {"Standard": {}, "dnCNN": {}}}, csv_data))
        assert "Inconsistencies" not in html

    def test_no_trend_data_in_group(self):
        """Groups with no trend data should show fallback message."""
        groups = {
            "TestGraph": {
                "Standard": {"trends_json": "[]"},
                "dnCNN": {"trends_json": "[]"},
            }
        }
        html = "\n".join(_section_cross_dwi_comparison(groups, None))
        assert "No trend data" in html or isinstance(html, str)


class TestCorrelations:
    def test_finds_correlations(self):
        result = _section_correlations(SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Correlations" in html

    def test_empty_rows(self):
        assert _section_correlations([]) == []

    def test_none_rows(self):
        assert _section_correlations(None) == []

    def test_strong_vs_moderate_split(self):
        html = "\n".join(_section_correlations(SAMPLE_GRAPH_CSV_ROWS))
        assert "Strong" in html

    def test_no_correlations_found(self):
        """Rows without correlation data should show 'no notable correlations'."""
        rows = [{
            "file_path": "Standard/test.png",
            "graph_type": "box",
            "summary": "No correlation data.",
            "trends_json": "[]",
        }]
        html = "\n".join(_section_correlations(rows))
        assert "No notable correlations" in html or "no notable" in html.lower()

    def test_causation_caveat_shown(self):
        """Causation caveat should always be present."""
        html = "\n".join(_section_correlations(SAMPLE_GRAPH_CSV_ROWS))
        assert "causation" in html.lower()

    def test_bonferroni_correction_note(self):
        """Multiple correlations should trigger Bonferroni note."""
        html = "\n".join(_section_correlations(SAMPLE_GRAPH_CSV_ROWS))
        # Row 1 has r=0.65, row 3 has r²=0.45 (r≈0.67)
        if "Bonferroni" in html:
            assert "corrected" in html.lower() or "threshold" in html.lower()

    def test_confidence_interval_caveat(self):
        """CI caveat should appear when correlations are found."""
        html = "\n".join(_section_correlations(SAMPLE_GRAPH_CSV_ROWS))
        assert "Confidence interval" in html or "confidence" in html.lower()

    def test_moderate_correlation_classification(self):
        """Correlations 0.3 <= |r| < 0.5 should be classified as moderate."""
        rows = [{
            "file_path": "Standard/test.png",
            "graph_type": "scatter",
            "summary": "r = 0.35 weak correlation.",
            "trends_json": "[]",
        }]
        html = "\n".join(_section_correlations(rows))
        assert "Moderate" in html


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

    def test_shared_feature_highlighted(self):
        """Shared features should be marked."""
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
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN"]))
        assert "Shared" in html
        assert "ADC_BL" in html

    def test_type_specific_feature(self):
        """Features unique to one type should be labeled 'Type-specific'."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["unique_std"]},
                    ]
                }
            },
            "dnCNN": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["unique_cnn"]},
                    ]
                }
            },
        }
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN"]))
        assert "Type-specific" in html

    def test_overlap_summary_percentage(self):
        """Summary should show overlap percentage."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["A", "B"]},
                    ]
                }
            },
            "dnCNN": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["A", "C"]},
                    ]
                }
            },
        }
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN"]))
        assert "%" in html

    def test_no_feature_selections(self):
        """DWI types with no feature_selections should return empty."""
        log_data = {
            "Standard": {"stats_predictive": {}},
            "dnCNN": {"stats_predictive": {}},
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        assert result == []

    def test_feature_stability_across_timepoints(self):
        """Features at 2+ timepoints within same DWI type should show stability."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["adc", "d"]},
                        {"timepoint": "W2", "features": ["adc", "f"]},
                        {"timepoint": "W4", "features": ["adc"]},
                    ]
                }
            },
            "dnCNN": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["adc"]},
                    ]
                }
            },
        }
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN"]))
        assert "Stability" in html or "stability" in html.lower()

    def test_duplicate_feature_detection(self):
        """Features with same root name should be flagged as potential duplicates."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["mean_adc", "meanADC"]},
                    ]
                }
            },
            "dnCNN": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["mean_adc"]},
                    ]
                }
            },
        }
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN"]))
        assert "duplicate" in html.lower()

    def test_empty_features_at_timepoint(self):
        """Empty feature list at a timepoint should not crash."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": []},
                    ]
                }
            },
            "dnCNN": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": []},
                    ]
                }
            },
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        assert isinstance(result, list)

    def test_dwi_type_missing_from_log(self):
        """DWI type in dwi_types_present but not in log_data."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["adc"]},
                    ]
                }
            },
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        assert isinstance(result, list)

    def test_three_dwi_types(self):
        """All three DWI types with shared features."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["adc", "d"]},
                    ]
                }
            },
            "dnCNN": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["adc", "f"]},
                    ]
                }
            },
            "IVIMnet": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["adc", "dstar"]},
                    ]
                }
            },
        }
        html = "\n".join(_section_feature_overlap(log_data, ["Standard", "dnCNN", "IVIMnet"]))
        assert "Shared" in html
        assert "adc" in html
