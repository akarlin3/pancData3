"""Unit tests for cross_reference_summary.py.

Tests the concise cross-DWI summary output including:
- Priority graph ordering (priority graphs listed first)
- Trend agreement/disagreement detection (AGREE vs DIFFER)
- Series union logic (default to "overall" when series is None)
- Statistical tests comparison section
- Clinical relevance truncation (150 char cap)
- Summary truncation (180 char cap)
- Parameter map counts section
- Inflection points section (longitudinal graphs only)
- Graphs with <2 real DWI types skipped
"""

from __future__ import annotations

import csv
import json
import sys
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from conftest import GRAPH_CSV_COLUMNS, SAMPLE_GRAPH_CSV_ROWS


def _run_main(saved_files_dir: Path) -> str:
    """Run cross_reference_summary.main() and capture stdout."""
    import shared
    shared.reset_config_cache()

    from cross_reference import cross_reference_summary
    buf = StringIO()
    with patch.object(sys, "argv", ["cross_reference_summary.py", str(saved_files_dir)]), \
         patch.object(sys, "stdout", buf):
        cross_reference_summary.main()
    return buf.getvalue()


class TestTrendAgreement:
    """Tests for trend agreement/disagreement detection."""

    def test_differ_when_opposite_directions(self, saved_files_with_graph_csv):
        """Standard=increasing, dnCNN=decreasing → DIFFER."""
        output = _run_main(saved_files_with_graph_csv)
        assert "DIFFER" in output

    def test_agree_when_same_directions(self, saved_files_dir):
        """Two DWI types with same trend direction → AGREE."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["trends_json"] = json.dumps([
            {"series": "ADC", "direction": "increasing", "description": "rises"}
        ])
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["trends_json"] = json.dumps([
            {"series": "ADC", "direction": "increasing", "description": "also rises"}
        ])

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "AGREE" in output

    def test_overall_default_for_none_series(self, saved_files_dir):
        """When series is None, it should default to "overall"."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["trends_json"] = json.dumps([
            {"series": None, "direction": "increasing", "description": "test"}
        ])
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["trends_json"] = json.dumps([
            {"series": None, "direction": "increasing", "description": "test"}
        ])

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "[overall]" in output

    def test_disagree_marker_prefix(self, saved_files_with_graph_csv):
        """Disagreements should have '>>' prefix for visual highlighting."""
        output = _run_main(saved_files_with_graph_csv)
        assert ">>" in output


class TestStatisticalTestsSection:
    """Tests for the statistical tests comparison section."""

    def test_stat_test_shown(self, saved_files_with_graph_csv):
        """Structured test from Standard row should appear with p-value."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Standard: Wilcoxon" in output
        assert "p=0.0030" in output

    def test_stat_test_comparison_groups(self, saved_files_with_graph_csv):
        """Comparison groups should be shown in parentheses."""
        output = _run_main(saved_files_with_graph_csv)
        assert "(LF vs LC)" in output


class TestClinicalRelevanceSection:
    """Tests for clinical relevance comparison."""

    def test_clinical_relevance_shown(self, saved_files_with_graph_csv):
        """Standard has clinical relevance text; it should appear."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Standard clinical:" in output
        assert "Higher ADC in LF patients" in output

    def test_clinical_relevance_truncation_at_150(self, saved_files_dir):
        """Clinical relevance longer than 150 chars should be truncated."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["clinical_relevance"] = "X" * 200
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "X" * 150 + "..." in output
        assert "X" * 151 not in output


class TestSummarySection:
    """Tests for summary comparison."""

    def test_summary_shown_per_dwi_type(self, saved_files_with_graph_csv):
        """Each DWI type's summary should be printed."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Standard: Box plot showing" in output
        assert "dnCNN: Box plot showing" in output

    def test_summary_truncation_at_180(self, saved_files_dir):
        """Summaries longer than 180 chars should be truncated."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["summary"] = "S" * 200
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "S" * 180 + "..." in output
        assert "S" * 181 not in output


class TestParameterMapSection:
    """Tests for parameter map counts section."""

    def test_parameter_map_section_header(self, saved_files_with_graph_csv):
        """Parameter map section header should always appear."""
        output = _run_main(saved_files_with_graph_csv)
        assert "PARAMETER MAP COUNTS BY DWI TYPE" in output

    def test_parameter_maps_listed(self, saved_files_dir):
        """Graphs with 'Parameter_Maps' in name should appear in section."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["file_path"] = "saved_files_20260301_120000/Standard/ADC_Parameter_Maps_Standard.png"
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["file_path"] = "saved_files_20260301_120000/dnCNN/ADC_Parameter_Maps_dnCNN.png"

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "ADC_Parameter_Maps" in output


class TestInflectionPointsSection:
    """Tests for inflection points section."""

    def test_inflection_section_header(self, saved_files_with_graph_csv):
        """Inflection points section header should always appear."""
        output = _run_main(saved_files_with_graph_csv)
        assert "INFLECTION POINTS: Longitudinal Graphs" in output

    def test_inflection_points_shown_for_longitudinal(self, saved_files_dir):
        """Longitudinal graph with inflection points should display them."""
        row_std = dict(SAMPLE_GRAPH_CSV_ROWS[2])
        row_std["file_path"] = "saved_files_20260301_120000/Standard/Longitudinal_Mean_Metrics_Standard.png"
        row_ivim = dict(SAMPLE_GRAPH_CSV_ROWS[2])

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row_std)
            writer.writerow(row_ivim)

        output = _run_main(saved_files_dir)
        assert "Longitudinal_Mean_Metrics:" in output
        assert "x=60.0" in output
        assert "LF diverges from LC at day 60" in output


class TestEdgeCases:
    """Edge case tests."""

    def test_empty_csv_exits(self, saved_files_dir):
        """Empty CSV should cause sys.exit."""
        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()

        with pytest.raises(SystemExit):
            _run_main(saved_files_dir)

    def test_single_dwi_type_skipped(self, saved_files_dir):
        """Graph present in only 1 DWI type should be skipped."""
        row = dict(SAMPLE_GRAPH_CSV_ROWS[2])  # IVIMnet only

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row)

        # No graphs to compare, but should not crash
        import shared
        shared.reset_config_cache()
        from cross_reference import cross_reference_summary
        buf = StringIO()
        with patch.object(sys, "argv", ["test", str(saved_files_dir)]), \
             patch.object(sys, "stdout", buf):
            cross_reference_summary.main()
        output = buf.getvalue()
        assert "Longitudinal_Mean_Metrics" not in output.split("PARAMETER MAP")[0]

    def test_malformed_trends_json(self, saved_files_dir):
        """Malformed trends_json should not crash."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["trends_json"] = "INVALID JSON"
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "Feature_BoxPlots" in output
