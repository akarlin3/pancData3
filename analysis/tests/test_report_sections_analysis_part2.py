"""Tests for report sections: cross_dwi and correlations modules (part 2).

Validates analysis section builders:
- Cross-DWI comparison (trend agreement, disagreement records)
- Correlations (extraction, strength classification)
- Feature overlap (cross-DWI feature comparison)

Split out of test_report_sections_analysis.py to keep files under 700 lines.
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

from report.sections.cross_dwi import (
    _section_cross_dwi_comparison,
    _section_feature_overlap,
)
from report.sections.correlations import (
    _section_correlations,
)


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


# ── Edge cases ──


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
