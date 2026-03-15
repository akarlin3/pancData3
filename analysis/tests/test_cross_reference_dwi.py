"""Unit tests for cross_reference_dwi.py.

Tests the full cross-DWI comparison output including:
- Graph matching across DWI types (2+ real types required)
- Trend parsing and display (series, direction, magnitude, start/end values)
- Inflection point display
- Statistical test formatting (p-value precision, statistic values)
- Outlier and reference line display
- Summary truncation (250 char cap)
- Clinical relevance truncation (200 char cap)
- Metadata line assembly (sample_size, error_bars, density, etc.)
- JSON parsing robustness (malformed, empty, None)
- Legend and spatial pattern display
"""

from __future__ import annotations

import csv
import json
import sys
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import pytest

# Ensure analysis/ root is on sys.path
ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from conftest import GRAPH_CSV_COLUMNS, SAMPLE_GRAPH_CSV_ROWS


def _run_main(saved_files_dir: Path) -> str:
    """Run cross_reference_dwi.main() and capture stdout."""
    from cross_reference import cross_reference_dwi

    buf = StringIO()
    with patch.object(sys, "argv", ["cross_reference_dwi.py", str(saved_files_dir)]), \
         patch.object(sys, "stdout", buf):
        cross_reference_dwi.main()
    return buf.getvalue()


class TestCrossReferenceDwiOutput:
    """Tests that validate the detailed output content of cross_reference_dwi."""

    def test_matched_graph_count(self, saved_files_with_graph_csv):
        """Only Feature_BoxPlots is present in 2+ DWI types; count should be 1."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Total matched graph sets across DWI types: 1" in output

    def test_graph_header_shows_types_present(self, saved_files_with_graph_csv):
        """Matched graph header should list all DWI types it appears in."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Graph: Feature_BoxPlots" in output
        assert "Standard" in output
        assert "dnCNN" in output

    def test_trend_direction_displayed(self, saved_files_with_graph_csv):
        """Trends should show series name, direction, and description."""
        output = _run_main(saved_files_with_graph_csv)
        assert "[ADC] increasing" in output
        assert "ADC rises over time" in output
        assert "[ADC] decreasing" in output
        assert "ADC drops after denoising" in output

    def test_trend_start_end_values(self, saved_files_with_graph_csv):
        """Start/end values should be formatted to 4 decimal places."""
        output = _run_main(saved_files_with_graph_csv)
        assert "0.0010 -> 0.0015" in output

    def test_trend_magnitude_shown(self, saved_files_with_graph_csv):
        """Trend magnitude should be displayed in parentheses."""
        output = _run_main(saved_files_with_graph_csv)
        assert "(~15% increase)" in output

    def test_statistical_test_pvalue_formatting(self, saved_files_with_graph_csv):
        """Structured test p-values should be formatted to 4 decimal places."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Wilcoxon: p=0.0030" in output
        assert "(LF vs LC)" in output

    def test_outlier_displayed(self, saved_files_with_graph_csv):
        """dnCNN row has an outlier; it should appear in output."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Outliers (1):" in output
        assert "[ADC] extreme high ADC value" in output

    def test_metadata_line(self, saved_files_with_graph_csv):
        """Metadata should include sample size, series count, error bars."""
        output = _run_main(saved_files_with_graph_csv)
        assert "n=42" in output
        assert "2 series" in output
        assert "error bars: IQR" in output

    def test_clinical_relevance_shown(self, saved_files_with_graph_csv):
        """Clinical relevance text should appear for Standard row."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Higher ADC in LF patients" in output

    def test_summary_shown(self, saved_files_with_graph_csv):
        """Summary text should appear in output."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Box plot showing p = 0.003 for ADC comparison" in output

    def test_legend_items_displayed(self, saved_files_with_graph_csv):
        """Legend items should be displayed as comma-separated list."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Legend: LF, LC" in output

    def test_graph_type_shown(self, saved_files_with_graph_csv):
        """Graph type (box) should be displayed per DWI type."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Type: box" in output

    def test_axis_labels_shown(self, saved_files_with_graph_csv):
        """Axis labels with units should be displayed."""
        output = _run_main(saved_files_with_graph_csv)
        assert "X-axis: Feature" in output
        assert "Y-axis: Value (mm²/s" in output

    def test_single_dwi_type_skipped(self, saved_files_with_graph_csv):
        """Longitudinal_Mean_Metrics only has IVIMnet; should NOT appear."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Graph: Longitudinal_Mean_Metrics" not in output


class TestCrossReferenceDwiEdgeCases:
    """Edge case tests for cross_reference_dwi."""

    def test_empty_csv_exits(self, saved_files_dir):
        """Empty CSV (header only) should cause sys.exit."""
        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()

        with pytest.raises(SystemExit):
            _run_main(saved_files_dir)

    def test_malformed_json_handled(self, saved_files_dir):
        """Malformed JSON in trends_json should not crash, just skip trends."""
        row = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row["trends_json"] = "NOT VALID JSON"
        # Create a second row in a different DWI type
        row2 = dict(row)
        row2["file_path"] = "saved_files_20260301_120000/dnCNN/Feature_BoxPlots_dnCNN.png"
        row2["trends_json"] = "[]"

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row)
            writer.writerow(row2)

        # Should not raise
        output = _run_main(saved_files_dir)
        assert "Feature_BoxPlots" in output

    def test_summary_truncation_at_250(self, saved_files_dir):
        """Summaries longer than 250 chars should be truncated with '...'."""
        row = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row["summary"] = "A" * 300
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["summary"] = "B" * 100

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "A" * 250 + "..." in output
        assert "A" * 251 not in output

    def test_clinical_relevance_truncation_at_200(self, saved_files_dir):
        """Clinical relevance longer than 200 chars should be truncated."""
        row = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row["clinical_relevance"] = "C" * 250
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "C" * 200 + "..." in output
        assert "C" * 201 not in output

    def test_reference_lines_displayed(self, saved_files_dir):
        """Reference lines from IVIMnet Longitudinal row should be displayed
        when the graph has 2+ DWI types."""
        # Create Standard + IVIMnet versions of Longitudinal graph
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
        assert "Reference lines (1):" in output
        assert "horizontal" in output
        assert "'ADC threshold'" in output

    def test_inflection_points_displayed(self, saved_files_dir):
        """Inflection points should show coordinates and description."""
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
        assert "Inflection points (1):" in output
        assert "(60.0, 0.0012)" in output
        assert "LF diverges from LC at day 60" in output

    def test_non_dict_trend_skipped(self, saved_files_dir):
        """Non-dict entries in trends_json array should be skipped."""
        row = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row["trends_json"] = json.dumps(["not a dict", {"direction": "up", "description": "valid"}])
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "up: valid" in output

    def test_stat_test_with_non_numeric_pvalue(self, saved_files_dir):
        """Non-numeric p_value in structured test should be displayed as-is."""
        row = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row["statistical_tests_json"] = json.dumps([
            {"test_name": "Fisher", "p_value": "<0.001", "statistic_value": None,
             "comparison_groups": None}
        ])
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "Fisher: p=<0.001" in output
