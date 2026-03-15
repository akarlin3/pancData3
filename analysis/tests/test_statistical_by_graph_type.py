"""Unit tests for statistical_by_graph_type.py.

Tests the graph-type-grouped statistical analysis including:
- Grouping rows by graph_type
- Type ordering by descending row count
- Significant/non-significant p-value separation
- Top-5 non-significant display with "... and N more"
- Correlation extraction and strength classification
- Trend direction keyword matching (increasing/decreasing/stable/other)
- Structured statistical tests section
- Data density and comparison type aggregation
- Graph list with metadata
- Summary table formatting
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
    """Run statistical_by_graph_type.main() and capture stdout."""
    import shared
    shared.reset_config_cache()

    from cross_reference import statistical_by_graph_type
    buf = StringIO()
    with patch.object(sys, "argv", ["statistical_by_graph_type.py", str(saved_files_dir)]), \
         patch.object(sys, "stdout", buf):
        statistical_by_graph_type.main()
    return buf.getvalue()


class TestGraphTypeGrouping:
    """Tests for grouping by graph_type."""

    def test_box_type_section(self, saved_files_with_graph_csv):
        """BOX type section should appear (2 rows are box type)."""
        output = _run_main(saved_files_with_graph_csv)
        assert "GRAPH TYPE: BOX" in output

    def test_line_type_section(self, saved_files_with_graph_csv):
        """LINE type section should appear (1 row is line type)."""
        output = _run_main(saved_files_with_graph_csv)
        assert "GRAPH TYPE: LINE" in output

    def test_type_count_shown(self, saved_files_with_graph_csv):
        """Graph count per type should be shown in header."""
        output = _run_main(saved_files_with_graph_csv)
        assert "(2 graphs)" in output  # box has 2 rows
        assert "(1 graphs)" in output  # line has 1 row

    def test_type_ordering_by_count(self, saved_files_with_graph_csv):
        """Types should be ordered by descending row count."""
        output = _run_main(saved_files_with_graph_csv)
        box_pos = output.find("GRAPH TYPE: BOX")
        line_pos = output.find("GRAPH TYPE: LINE")
        assert box_pos < line_pos, "BOX (2 rows) should appear before LINE (1 row)"


class TestSignificantPValues:
    """Tests for per-type significant p-value extraction."""

    def test_significant_in_box(self, saved_files_with_graph_csv):
        """p=0.003 from Standard Feature_BoxPlots should appear in BOX section."""
        output = _run_main(saved_files_with_graph_csv)
        # Find the BOX section and check for significant finding
        box_section = output[output.find("GRAPH TYPE: BOX"):]
        assert "0.003" in box_section.split("GRAPH TYPE: LINE")[0]

    def test_significance_markers(self, saved_files_with_graph_csv):
        """Significance markers should appear for p-values < threshold."""
        output = _run_main(saved_files_with_graph_csv)
        # p=0.003 should get ** marker (< 0.01 but >= 0.001)
        assert "** " in output


class TestNonSignificantTop5:
    """Tests for non-significant findings display."""

    def test_nonsig_count_shown(self, saved_files_with_graph_csv):
        """Non-significant count should be displayed."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Non-significant:" in output

    def test_top5_overflow_message(self, saved_files_dir):
        """When >5 non-sig findings, show '... and N more' message."""
        rows = []
        for i in range(8):
            row = dict(SAMPLE_GRAPH_CSV_ROWS[0])
            row["file_path"] = f"saved_files_20260301_120000/Standard/Graph{i}_Standard.png"
            row["summary"] = f"Result: p = 0.{50+i} for metric{i}."
            row["statistical_tests_json"] = "[]"
            row["trends_json"] = "[]"
            rows.append(row)

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

        output = _run_main(saved_files_dir)
        assert "... and" in output
        assert "more" in output


class TestCorrelations:
    """Tests for per-type correlation extraction."""

    def test_notable_correlation_in_box(self, saved_files_with_graph_csv):
        """r=0.65 should appear in BOX type's correlations section."""
        output = _run_main(saved_files_with_graph_csv)
        assert "r=0.65" in output

    def test_correlation_strength(self, saved_files_with_graph_csv):
        """r=0.65 >= 0.5 should be STRONG positive."""
        output = _run_main(saved_files_with_graph_csv)
        assert "STRONG" in output
        assert "positive" in output

    def test_correlations_sorted_by_strength(self, saved_files_dir):
        """Correlations should be sorted by absolute strength (strongest first)."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["summary"] = "r = 0.35 and r = 0.82 found."
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["summary"] = "No correlations."

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        idx_82 = output.find("r=0.82")
        idx_35 = output.find("r=0.35")
        assert idx_82 < idx_35, "Strongest correlation should appear first"


class TestTrendDirections:
    """Tests for trend direction keyword matching."""

    def test_increasing_counted(self, saved_files_with_graph_csv):
        """'increasing' direction should be counted."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Increasing:" in output

    def test_decreasing_counted(self, saved_files_with_graph_csv):
        """'decreasing' direction should be counted."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Decreasing:" in output

    def test_stable_trend(self, saved_files_dir):
        """'flat' direction should be counted as stable."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["trends_json"] = json.dumps([
            {"series": "ADC", "direction": "flat", "description": "no change"}
        ])
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["trends_json"] = json.dumps([
            {"series": "ADC", "direction": "stable", "description": "constant"}
        ])

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "Stable: 2" in output

    def test_other_trend_keyword(self, saved_files_dir):
        """Unrecognised direction keywords should count as 'Other'."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["trends_json"] = json.dumps([
            {"series": "ADC", "direction": "oscillating", "description": "irregular"}
        ])
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["trends_json"] = "[]"

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "Other: 1" in output


class TestStructuredTests:
    """Tests for structured statistical tests per graph type."""

    def test_structured_test_shown(self, saved_files_with_graph_csv):
        """Wilcoxon test from Standard row should appear."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Structured statistical tests" in output
        assert "Wilcoxon" in output

    def test_structured_test_pvalue(self, saved_files_with_graph_csv):
        """Structured test p-value should be formatted."""
        output = _run_main(saved_files_with_graph_csv)
        assert "p=0.0030" in output


class TestDataDensityAndComparisonTypes:
    """Tests for density and comparison type aggregation."""

    def test_density_aggregation(self, saved_files_with_graph_csv):
        """Data density counts should appear."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Data density:" in output
        assert "moderate" in output

    def test_comparison_type_aggregation(self, saved_files_with_graph_csv):
        """Comparison type counts should appear."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Comparison types:" in output
        assert "unpaired" in output


class TestGraphList:
    """Tests for graph list per type."""

    def test_graphs_listed(self, saved_files_with_graph_csv):
        """All graphs should be listed with DWI type and name."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Graphs in this type:" in output
        assert "[Standard] Feature_BoxPlots" in output
        assert "[dnCNN] Feature_BoxPlots" in output

    def test_graph_metadata_shown(self, saved_files_with_graph_csv):
        """Sample size and quality should appear in graph list."""
        output = _run_main(saved_files_with_graph_csv)
        assert "n=42" in output
        assert "quality:high" in output


class TestSummaryTable:
    """Tests for the final summary table."""

    def test_summary_table_header(self, saved_files_with_graph_csv):
        """Summary table header should appear."""
        output = _run_main(saved_files_with_graph_csv)
        assert "SUMMARY TABLE: Statistical Findings by Graph Type" in output

    def test_summary_table_columns(self, saved_files_with_graph_csv):
        """Summary table should have Type, Count, Sig p, Non-sig, Corr, Trends columns."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Type" in output
        assert "Count" in output
        assert "Sig p" in output
        assert "Trends" in output

    def test_summary_table_row_for_box(self, saved_files_with_graph_csv):
        """BOX type row should show count=2."""
        output = _run_main(saved_files_with_graph_csv)
        # Find the summary table section
        table_section = output[output.find("SUMMARY TABLE"):]
        assert "box" in table_section


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

    def test_empty_trends_json(self, saved_files_dir):
        """Empty trends_json should not crash or count any trends."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["trends_json"] = "[]"
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["trends_json"] = "[]"

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        # Should not crash; no trend directions section
        output = _run_main(saved_files_dir)
        assert "GRAPH TYPE: BOX" in output

    def test_malformed_trends_json(self, saved_files_dir):
        """Malformed trends_json should not crash."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["trends_json"] = "NOT JSON"
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "GRAPH TYPE: BOX" in output
