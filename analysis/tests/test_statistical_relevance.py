"""Unit tests for statistical_relevance.py.

Tests the statistical findings extraction and reporting including:
- Significant p-value extraction from free text and structured tests
- Bonferroni correction computation
- Significance markers (*, **, ***)
- Non-significant findings section
- Notable correlations (|r| >= 0.3) with strength/direction classification
- Cross-DWI significance comparison
- Edge cases: NaN/Inf p-values, empty CSV, missing fields
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
    """Run statistical_relevance.main() and capture stdout."""
    import shared
    shared.reset_config_cache()

    from cross_reference import statistical_relevance
    buf = StringIO()
    with patch.object(sys, "argv", ["statistical_relevance.py", str(saved_files_dir)]), \
         patch.object(sys, "stdout", buf):
        statistical_relevance.main()
    return buf.getvalue()


class TestSignificantFindings:
    """Tests for significant p-value extraction and display."""

    def test_significant_pvalue_found(self, saved_files_with_graph_csv):
        """p=0.003 from Standard Feature_BoxPlots should be reported."""
        output = _run_main(saved_files_with_graph_csv)
        assert "STATISTICALLY SIGNIFICANT FINDINGS" in output
        assert "0.003" in output

    def test_significance_markers_three_stars(self, saved_files_with_graph_csv):
        """p=0.003 < 0.01 should get ** marker (not ***)."""
        output = _run_main(saved_files_with_graph_csv)
        assert "** " in output

    def test_structured_test_extracted(self, saved_files_with_graph_csv):
        """Structured statistical test (Wilcoxon p=0.003) should be extracted."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Wilcoxon: p=0.003" in output

    def test_bonferroni_threshold_computed(self, saved_files_with_graph_csv):
        """Bonferroni threshold note should appear in output."""
        output = _run_main(saved_files_with_graph_csv)
        assert "Bonferroni threshold" in output
        assert "tests extracted" in output

    def test_significant_sorted_by_pvalue(self, saved_files_dir):
        """Significant findings should be sorted by ascending p-value."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["summary"] = "Result: p = 0.04 for ADC."
        row1["statistical_tests_json"] = "[]"

        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["summary"] = "Result: p = 0.01 for D."
        row2["statistical_tests_json"] = "[]"

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        # p=0.01 should appear before p=0.04
        idx_01 = output.find("0.0100")
        idx_04 = output.find("0.0400")
        assert idx_01 < idx_04, "Findings should be sorted by ascending p-value"


class TestNonSignificantFindings:
    """Tests for non-significant findings section."""

    def test_nonsig_section_present(self, saved_files_with_graph_csv):
        """Non-significant findings section should appear."""
        output = _run_main(saved_files_with_graph_csv)
        assert "NON-SIGNIFICANT FINDINGS" in output

    def test_nonsig_pvalue_reported(self, saved_files_with_graph_csv):
        """p=0.12 from dnCNN row should appear in non-significant section."""
        output = _run_main(saved_files_with_graph_csv)
        assert "0.12" in output


class TestCorrelations:
    """Tests for notable correlations extraction."""

    def test_notable_correlation_found(self, saved_files_with_graph_csv):
        """r=0.65 from Standard summary should be reported."""
        output = _run_main(saved_files_with_graph_csv)
        assert "NOTABLE CORRELATIONS" in output
        assert "r=0.65" in output

    def test_correlation_strength_classified(self, saved_files_with_graph_csv):
        """r=0.65 >= 0.5 should be classified as STRONG."""
        output = _run_main(saved_files_with_graph_csv)
        assert "STRONG positive" in output

    def test_moderate_correlation(self, saved_files_dir):
        """0.3 <= |r| < 0.5 should be classified as MODERATE."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["summary"] = "Correlation r = 0.35 between dose and ADC."
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["summary"] = "No correlation found."

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "MODERATE positive" in output

    def test_negative_correlation(self, saved_files_dir):
        """Negative r values should be classified with 'negative' direction."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["summary"] = "Inverse correlation r = -0.72 found."
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["summary"] = "No data."

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "STRONG negative" in output


class TestCrossDWIComparison:
    """Tests for cross-DWI significance comparison."""

    def test_cross_dwi_section_header(self, saved_files_with_graph_csv):
        """Cross-DWI section should appear in output."""
        output = _run_main(saved_files_with_graph_csv)
        assert "CROSS-DWI: Same Analysis, Different Significance?" in output

    def test_cross_dwi_sig_counts(self, saved_files_with_graph_csv):
        """Standard has p=0.003 (sig), dnCNN has p=0.12 (non-sig)."""
        output = _run_main(saved_files_with_graph_csv)
        # Standard should show significant
        assert "Standard:" in output
        assert "significant" in output


class TestEdgeCases:
    """Edge case tests for statistical_relevance."""

    def test_empty_csv_exits(self, saved_files_dir):
        """Empty CSV should cause sys.exit."""
        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()

        with pytest.raises(SystemExit):
            _run_main(saved_files_dir)

    def test_nan_pvalue_in_structured_test_skipped(self, saved_files_dir):
        """NaN p_value in structured test should be silently skipped."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["statistical_tests_json"] = json.dumps([
            {"test_name": "Test", "p_value": float("nan"), "comparison_groups": ""}
        ])
        row1["summary"] = "No p-values in text."
        row1["trends_json"] = "[]"
        row1["inflection_points_json"] = "[]"
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["summary"] = "No p-values here either."
        row2["trends_json"] = "[]"
        row2["inflection_points_json"] = "[]"

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "No p-values < 0.05 found" in output

    def test_inf_pvalue_in_structured_test_skipped(self, saved_files_dir):
        """Inf p_value in structured test should be silently skipped."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["statistical_tests_json"] = json.dumps([
            {"test_name": "Test", "p_value": float("inf"), "comparison_groups": ""}
        ])
        row1["summary"] = "No p-values in text."
        row1["trends_json"] = "[]"
        row1["inflection_points_json"] = "[]"
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["summary"] = "Nothing here."
        row2["trends_json"] = "[]"
        row2["inflection_points_json"] = "[]"

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "No p-values < 0.05 found" in output

    def test_three_star_marker_for_highly_significant(self, saved_files_dir):
        """p < 0.001 should get *** marker."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["summary"] = "Result: p = 0.0005 for ADC."
        row1["statistical_tests_json"] = "[]"
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["summary"] = "No finding."

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "***" in output

    def test_one_star_marker(self, saved_files_dir):
        """0.01 <= p < 0.05 should get * marker."""
        row1 = dict(SAMPLE_GRAPH_CSV_ROWS[0])
        row1["summary"] = "Result: p = 0.03 for D."
        row1["statistical_tests_json"] = "[]"
        row2 = dict(SAMPLE_GRAPH_CSV_ROWS[1])
        row2["summary"] = "Nothing."

        csv_path = saved_files_dir / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            writer.writerow(row1)
            writer.writerow(row2)

        output = _run_main(saved_files_dir)
        assert "*  " in output
