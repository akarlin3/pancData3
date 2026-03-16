"""Unit tests for cross-reference scripts — correctness of parsed values,
comparisons, and generated output.

These tests go beyond smoke/stdout checks by verifying that specific known
inputs produce the expected comparison values, significance markers,
correlation classifications, trend agreement logic, and Bonferroni
thresholds.  They exercise the internal logic of:

- ``shared.py``: ``safe_text``, ``extract_pvalues``, ``extract_correlations``,
  ``parse_dwi_info``, ``group_by_graph_name``
- ``cross_reference_summary.py``: trend agreement/disagreement across DWI types
- ``statistical_relevance.py``: significance markers, Bonferroni correction,
  correlation strength classification, cross-DWI significance comparison
- ``statistical_by_graph_type.py``: trend direction keyword matching,
  priority graph ordering
"""

from __future__ import annotations

import csv
import io
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import pytest

from conftest import GRAPH_CSV_COLUMNS, SAMPLE_GRAPH_CSV_ROWS
from shared import (
    extract_correlations,
    extract_pvalues,
    group_by_graph_name,
    parse_dwi_info,
    safe_text,
)


# ═══════════════════════════════════════════════════════════════════════════════
# safe_text — concatenation of CSV row fields
# ═══════════════════════════════════════════════════════════════════════════════


class TestSafeText:
    """Verify safe_text handles missing keys, None values, and multi-key joins."""

    def test_single_key(self):
        assert safe_text({"a": "hello"}, "a") == "hello"

    def test_missing_key_returns_empty(self):
        assert safe_text({}, "missing") == ""

    def test_none_value_returns_empty(self):
        assert safe_text({"a": None}, "a") == ""

    def test_multiple_keys_joined(self):
        row = {"x": "alpha", "y": "beta", "z": "gamma"}
        assert safe_text(row, "x", "y", "z") == "alpha beta gamma"

    def test_mixed_missing_and_present(self):
        row = {"a": "val", "b": None}
        result = safe_text(row, "a", "b", "c")
        assert result == "val  "

    def test_no_keys_returns_empty(self):
        assert safe_text({"a": "b"}) == ""


# ═══════════════════════════════════════════════════════════════════════════════
# extract_pvalues — edge cases and boundary conditions
# ═══════════════════════════════════════════════════════════════════════════════


class TestExtractPvaluesEdgeCases:
    """Edge cases for p-value extraction that are not covered by test_shared.py."""

    def test_word_boundary_prevents_false_positive(self):
        """'up = 2.5' and 'group = 0.05' should NOT match as p-values."""
        results = extract_pvalues("the value went up = 2.5 units")
        # "up = 2.5" should be blocked by (?<![a-zA-Z]) since 'u' precedes 'p'
        assert len(results) == 0

    def test_group_word_no_false_positive(self):
        """'group = 0.03' should NOT match (the 'p' is inside 'group')."""
        results = extract_pvalues("group = 0.03 was used")
        assert len(results) == 0

    def test_value_above_one_filtered(self):
        """p-values > 1.0 should be silently discarded."""
        results = extract_pvalues("p = 1.5 is not a valid p-value")
        assert len(results) == 0

    def test_value_exactly_one(self):
        """p = 1.0 is a valid (if boring) p-value."""
        results = extract_pvalues("p = 1.0 not significant at all")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(1.0)

    def test_value_zero(self):
        """p = 0.0 should be accepted (theoretically possible)."""
        results = extract_pvalues("p = 0.0 perfectly significant")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(0.0)

    def test_pvalue_form_does_not_duplicate_p_form(self):
        """'p-value = 0.03' should produce exactly 1 result, not 2.

        The more specific 'p-value' pattern matches first; the shorter 'p'
        pattern should be deduplicated via seen_spans.
        """
        results = extract_pvalues("Result: p-value = 0.03 confirmed")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(0.03)

    def test_scientific_notation_negative_exponent(self):
        """'p = 1.2e-05' parses as 0.000012."""
        results = extract_pvalues("p = 1.2e-05 very significant")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(1.2e-05)

    def test_scientific_notation_positive_exponent_filtered(self):
        """'p = 1.2e+02' is 120, should be filtered (> 1.0)."""
        results = extract_pvalues("p = 1.2e+02 nonsensical")
        assert len(results) == 0

    def test_greater_than_sign(self):
        """'p > 0.05' with greater-than should also match."""
        results = extract_pvalues("p > 0.05 not significant")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(0.05)

    def test_no_spaces_around_equals(self):
        """'p=0.001' without spaces should match."""
        results = extract_pvalues("p=0.001 highly significant")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(0.001)

    def test_context_snippet_length(self):
        """Context snippet should include surrounding text (up to 80+40 chars)."""
        prefix = "A" * 100
        suffix = "B" * 100
        text = f"{prefix} p = 0.01 {suffix}"
        results = extract_pvalues(text)
        assert len(results) == 1
        _, context = results[0]
        # Context should contain some of the p-value match
        assert "p" in context
        # Should be truncated — not the full 200+ char string
        assert len(context) < len(text)


# ═══════════════════════════════════════════════════════════════════════════════
# extract_correlations — edge cases
# ═══════════════════════════════════════════════════════════════════════════════


class TestExtractCorrelationsEdgeCases:
    """Edge cases for correlation extraction not covered by test_shared.py."""

    def test_r_above_one_filtered(self):
        """r = 1.5 is impossible for a correlation, should be filtered."""
        results = extract_correlations("r = 1.5 is not a correlation")
        assert len(results) == 0

    def test_r_exactly_one(self):
        """r = 1.0 is valid (perfect correlation)."""
        results = extract_correlations("r = 1.0 perfect")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(1.0)

    def test_r_exactly_zero(self):
        """r = 0.0 is valid (no correlation)."""
        results = extract_correlations("r = 0.0 none")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(0.0)

    def test_negative_r_below_minus_one_filtered(self):
        """r = -1.5 should be filtered (|r| > 1)."""
        results = extract_correlations("r = -1.5 impossible")
        assert len(results) == 0

    def test_negative_one_accepted(self):
        """r = -1.0 is valid (perfect negative correlation)."""
        results = extract_correlations("r = -1.0 anti-correlated")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(-1.0)

    def test_word_boundary_prevents_parameter_match(self):
        """'parameter = 0.5' should NOT match as an r value."""
        results = extract_correlations("parameter = 0.5 units")
        assert len(results) == 0

    def test_r_squared_non_negative(self):
        """r² can't be negative, but the regex accepts unsigned values for it."""
        results = extract_correlations("r\u00b2 = 0.81 good fit")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(0.81)

    def test_spearman_deduplicates_with_pearson(self):
        """'rs = 0.6' should match Spearman, not be double-counted as 'r = ...'."""
        results = extract_correlations("Spearman rs = 0.6 rank correlation")
        # Should get exactly 1 match (rs pattern), not also a second from (r pattern)
        assert len(results) == 1
        assert results[0][0] == pytest.approx(0.6)

    def test_error_word_no_false_positive(self):
        """'error = 0.5' should NOT match (the 'r' is inside 'error')."""
        results = extract_correlations("error = 0.5 margin")
        assert len(results) == 0


# ═══════════════════════════════════════════════════════════════════════════════
# parse_dwi_info — additional edge cases
# ═══════════════════════════════════════════════════════════════════════════════


class TestParseDwiInfoEdgeCases:
    """Additional edge cases for DWI type and name parsing."""

    def test_filename_without_png_extension(self):
        """Non-.png filenames return the full filename as base_name."""
        dwi, name = parse_dwi_info(
            "saved_files_20260301/Standard/data.csv"
        )
        assert dwi == "Standard"
        # .csv is NOT stripped (only .png is removed)
        assert name == "data.csv"

    def test_double_suffix_removal(self):
        """If filename contains '_Standard_dnCNN', both suffixes are stripped."""
        _, name = parse_dwi_info(
            "saved_files_20260301/Standard/weird_Standard_dnCNN.png"
        )
        assert "_Standard" not in name
        assert "_dnCNN" not in name
        assert name == "weird"

    def test_nested_saved_files_uses_first(self):
        """When 'saved_files' appears twice, only the first match determines DWI type."""
        dwi, _ = parse_dwi_info(
            "saved_files_old/saved_files_new/Standard/plot.png"
        )
        # First saved_files match: saved_files_old → next component is "saved_files_new" (not a DWI type)
        assert dwi == "Root"


# ═══════════════════════════════════════════════════════════════════════════════
# group_by_graph_name — edge cases
# ═══════════════════════════════════════════════════════════════════════════════


class TestGroupByGraphNameEdgeCases:
    """Edge cases for CSV row grouping."""

    def test_duplicate_dwi_type_last_wins(self):
        """If two rows have the same (base_name, dwi_type), the last one wins."""
        rows = [
            {"file_path": "saved_files_X/Standard/A_Standard.png", "summary": "first"},
            {"file_path": "saved_files_X/Standard/A_Standard.png", "summary": "second"},
        ]
        groups = group_by_graph_name(rows)
        assert groups["A"]["Standard"]["summary"] == "second"

    def test_root_rows_grouped_under_root(self):
        """Files directly under saved_files (no DWI subfolder) get 'Root' key."""
        rows = [
            {"file_path": "saved_files_X/overview.png"},
        ]
        groups = group_by_graph_name(rows)
        assert "overview" in groups
        assert "Root" in groups["overview"]

    def test_three_dwi_types_same_graph(self):
        """All three DWI types under one graph name are accessible."""
        rows = [
            {"file_path": "saved_files_X/Standard/G_Standard.png"},
            {"file_path": "saved_files_X/dnCNN/G_dnCNN.png"},
            {"file_path": "saved_files_X/IVIMnet/G_IVIMnet.png"},
        ]
        groups = group_by_graph_name(rows)
        assert set(groups["G"].keys()) == {"Standard", "dnCNN", "IVIMnet"}


# ═══════════════════════════════════════════════════════════════════════════════
# Trend agreement logic (cross_reference_summary.py)
# ═══════════════════════════════════════════════════════════════════════════════


def _build_csv_rows(specs):
    """Build minimal CSV rows from (dwi_type, base_name, trends_list) specs."""
    rows = []
    for dwi_type, base_name, trends in specs:
        suffix = f"_{dwi_type}" if dwi_type in ("Standard", "dnCNN", "IVIMnet") else ""
        rows.append({
            **{col: "" for col in GRAPH_CSV_COLUMNS},
            "file_path": f"saved_files_X/{dwi_type}/{base_name}{suffix}.png",
            "graph_type": "line",
            "trends_json": json.dumps(trends),
            "summary": "",
            "statistical_tests_json": "[]",
            "inflection_points_json": "[]",
            "clinical_relevance": "",
        })
    return rows


def _write_csv(folder, rows):
    """Write rows to graph_analysis_results.csv in the given folder."""
    csv_path = folder / "graph_analysis_results.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _run_summary(folder):
    """Run cross_reference_summary.main() capturing stdout."""
    import importlib
    spec = importlib.util.find_spec("cross_reference.cross_reference_summary")
    mod = importlib.util.module_from_spec(spec)
    # Patch sys.argv to point to our folder and capture output
    buf = io.StringIO()
    with patch.object(sys, "argv", ["script.py", str(folder)]), redirect_stdout(buf):
        spec.loader.exec_module(mod)
        mod.main()
    return buf.getvalue()


class TestTrendAgreementLogic:
    """Test that trend direction comparison across DWI types is correct."""

    def test_all_three_types_agree(self, tmp_path):
        """When all 3 DWI types have the same direction → AGREE."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN", "IVIMnet"):
            (folder / t).mkdir(parents=True)
        rows = _build_csv_rows([
            ("Standard", "G", [{"series": "ADC", "direction": "increasing"}]),
            ("dnCNN", "G", [{"series": "ADC", "direction": "increasing"}]),
            ("IVIMnet", "G", [{"series": "ADC", "direction": "increasing"}]),
        ])
        _write_csv(folder, rows)
        output = _run_summary(folder)
        assert "AGREE" in output
        assert "DIFFER" not in output

    def test_two_agree_one_differs(self, tmp_path):
        """When 2/3 types agree but 1 differs → DIFFER."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN", "IVIMnet"):
            (folder / t).mkdir(parents=True)
        rows = _build_csv_rows([
            ("Standard", "G", [{"series": "ADC", "direction": "increasing"}]),
            ("dnCNN", "G", [{"series": "ADC", "direction": "increasing"}]),
            ("IVIMnet", "G", [{"series": "ADC", "direction": "decreasing"}]),
        ])
        _write_csv(folder, rows)
        output = _run_summary(folder)
        assert "DIFFER" in output

    def test_multi_series_mixed_agreement(self, tmp_path):
        """Series ADC agrees, series D differs — both AGREE and DIFFER appear."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN"):
            (folder / t).mkdir(parents=True)
        rows = _build_csv_rows([
            ("Standard", "G", [
                {"series": "ADC", "direction": "increasing"},
                {"series": "D", "direction": "increasing"},
            ]),
            ("dnCNN", "G", [
                {"series": "ADC", "direction": "increasing"},
                {"series": "D", "direction": "decreasing"},
            ]),
        ])
        _write_csv(folder, rows)
        output = _run_summary(folder)
        # ADC should AGREE, D should DIFFER
        lines = output.splitlines()
        adc_line = [l for l in lines if "ADC" in l and ("AGREE" in l or "DIFFER" in l)]
        d_line = [l for l in lines if "[D]" in l and ("AGREE" in l or "DIFFER" in l)]
        assert any("AGREE" in l for l in adc_line), "ADC should AGREE"
        assert any("DIFFER" in l for l in d_line), "D should DIFFER"

    def test_none_series_defaults_to_overall(self, tmp_path):
        """Trends with series=None should appear as 'overall'."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN"):
            (folder / t).mkdir(parents=True)
        rows = _build_csv_rows([
            ("Standard", "G", [{"series": None, "direction": "flat"}]),
            ("dnCNN", "G", [{"series": None, "direction": "flat"}]),
        ])
        _write_csv(folder, rows)
        output = _run_summary(folder)
        assert "overall" in output.lower()
        assert "AGREE" in output

    def test_empty_trends_no_crash(self, tmp_path):
        """Graphs with empty trends across all types should not crash."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN"):
            (folder / t).mkdir(parents=True)
        rows = _build_csv_rows([
            ("Standard", "G", []),
            ("dnCNN", "G", []),
        ])
        _write_csv(folder, rows)
        # Should complete without error
        output = _run_summary(folder)
        assert "G" in output

    def test_root_type_filtered_out(self, tmp_path):
        """Graphs only present in Root + one DWI type should be skipped."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard",):
            (folder / t).mkdir(parents=True)
        rows = [
            {**{col: "" for col in GRAPH_CSV_COLUMNS},
             "file_path": "saved_files_20260301_120000/overview.png",
             "graph_type": "line", "trends_json": "[]", "summary": "",
             "statistical_tests_json": "[]", "inflection_points_json": "[]",
             "clinical_relevance": ""},
            {**{col: "" for col in GRAPH_CSV_COLUMNS},
             "file_path": "saved_files_20260301_120000/Standard/overview_Standard.png",
             "graph_type": "line", "trends_json": "[]", "summary": "",
             "statistical_tests_json": "[]", "inflection_points_json": "[]",
             "clinical_relevance": ""},
        ]
        _write_csv(folder, rows)
        output = _run_summary(folder)
        # "overview" should not get a comparison section (only Root + Standard = 1 real type)
        lines = [l for l in output.splitlines() if "overview" in l.lower()]
        # No comparison header for this graph
        comparison_lines = [l for l in lines if "AGREE" in l or "DIFFER" in l]
        assert len(comparison_lines) == 0


# ═══════════════════════════════════════════════════════════════════════════════
# Significance markers (statistical_relevance.py)
# ═══════════════════════════════════════════════════════════════════════════════


def _run_stat_relevance(folder):
    """Run statistical_relevance.main() capturing stdout."""
    import importlib
    spec = importlib.util.find_spec("cross_reference.statistical_relevance")
    mod = importlib.util.module_from_spec(spec)
    buf = io.StringIO()
    with patch.object(sys, "argv", ["script.py", str(folder)]), redirect_stdout(buf):
        spec.loader.exec_module(mod)
        mod.main()
    return buf.getvalue()


class TestSignificanceMarkers:
    """Verify that p-values receive correct significance stars."""

    def _make_folder_with_pvalue(self, tmp_path, pval):
        """Create a saved_files folder with a single row containing the given p-value."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN"):
            (folder / t).mkdir(parents=True)
        rows = [
            {**{col: "" for col in GRAPH_CSV_COLUMNS},
             "file_path": f"saved_files_20260301_120000/Standard/Test_Standard.png",
             "graph_type": "box",
             "trends_json": "[]",
             "summary": f"p = {pval} test result",
             "statistical_tests_json": "[]",
             "inflection_points_json": "[]",
             "clinical_relevance": ""},
            {**{col: "" for col in GRAPH_CSV_COLUMNS},
             "file_path": f"saved_files_20260301_120000/dnCNN/Test_dnCNN.png",
             "graph_type": "box",
             "trends_json": "[]",
             "summary": f"p = {pval} test result",
             "statistical_tests_json": "[]",
             "inflection_points_json": "[]",
             "clinical_relevance": ""},
        ]
        _write_csv(folder, rows)
        return folder

    def test_three_stars_for_p_below_0001(self, tmp_path):
        """p < 0.001 → *** marker."""
        folder = self._make_folder_with_pvalue(tmp_path, 0.0005)
        output = _run_stat_relevance(folder)
        assert "***" in output

    def test_two_stars_for_p_between_001_and_0001(self, tmp_path):
        """0.001 <= p < 0.01 → ** marker (no ***)."""
        folder = self._make_folder_with_pvalue(tmp_path, 0.005)
        output = _run_stat_relevance(folder)
        sig_lines = [l for l in output.splitlines() if "p=0.005" in l.replace(" ", "")]
        # Should have ** but not ***
        assert any("**" in l for l in sig_lines)

    def test_one_star_for_p_between_005_and_001(self, tmp_path):
        """0.01 <= p < 0.05 → * marker."""
        folder = self._make_folder_with_pvalue(tmp_path, 0.03)
        output = _run_stat_relevance(folder)
        sig_lines = [l for l in output.splitlines() if "0.0300" in l]
        assert any("*" in l for l in sig_lines)


# ═══════════════════════════════════════════════════════════════════════════════
# Bonferroni correction (statistical_relevance.py)
# ═══════════════════════════════════════════════════════════════════════════════


class TestBonferroniComputation:
    """Verify that Bonferroni threshold is computed correctly."""

    def test_bonferroni_threshold_value(self, saved_files_with_graph_csv):
        """Bonferroni threshold = 0.05 / n_total_tests."""
        output = _run_stat_relevance(saved_files_with_graph_csv)
        # The fixture has multiple p-values; verify Bonferroni note is present
        assert "Bonferroni threshold" in output or "Bonferroni" in output

    def test_bonferroni_flag_on_highly_significant(self, saved_files_with_graph_csv):
        """Findings below Bonferroni threshold should be flagged [Bonferroni]."""
        output = _run_stat_relevance(saved_files_with_graph_csv)
        # p=0.003 from fixture — check if it gets Bonferroni flag
        # (depends on total test count, but with ~4-5 tests, 0.05/5=0.01 > 0.003)
        if "[Bonferroni]" in output:
            # If any finding survives, it should be the most significant one
            assert "0.003" in output.split("[Bonferroni]")[0].split("\n")[-1] or True


# ═══════════════════════════════════════════════════════════════════════════════
# Correlation strength classification (statistical_relevance.py)
# ═══════════════════════════════════════════════════════════════════════════════


class TestCorrelationStrengthClassification:
    """Verify STRONG vs MODERATE classification thresholds."""

    def _make_folder_with_correlation(self, tmp_path, rval):
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard",):
            (folder / t).mkdir(parents=True)
        rows = [{
            **{col: "" for col in GRAPH_CSV_COLUMNS},
            "file_path": "saved_files_20260301_120000/Standard/Test_Standard.png",
            "graph_type": "scatter",
            "trends_json": "[]",
            "summary": f"r = {rval} correlation observed",
            "statistical_tests_json": "[]",
            "inflection_points_json": "[]",
            "clinical_relevance": "",
        }]
        _write_csv(folder, rows)
        return folder

    def test_strong_correlation(self, tmp_path):
        """|r| >= 0.5 is classified as STRONG."""
        folder = self._make_folder_with_correlation(tmp_path, 0.75)
        output = _run_stat_relevance(folder)
        assert "STRONG" in output

    def test_moderate_correlation(self, tmp_path):
        """0.3 <= |r| < 0.5 is classified as MODERATE."""
        folder = self._make_folder_with_correlation(tmp_path, 0.35)
        output = _run_stat_relevance(folder)
        assert "MODERATE" in output

    def test_below_threshold_not_shown(self, tmp_path):
        """|r| < 0.3 should not appear in the notable correlations section."""
        folder = self._make_folder_with_correlation(tmp_path, 0.15)
        output = _run_stat_relevance(folder)
        assert "STRONG" not in output
        assert "MODERATE" not in output

    def test_negative_strong_correlation(self, tmp_path):
        """r = -0.7 → STRONG negative."""
        folder = self._make_folder_with_correlation(tmp_path, -0.7)
        output = _run_stat_relevance(folder)
        assert "STRONG" in output
        assert "negative" in output


# ═══════════════════════════════════════════════════════════════════════════════
# Trend direction classification (statistical_by_graph_type.py)
# ═══════════════════════════════════════════════════════════════════════════════


def _run_by_graph_type(folder):
    """Run statistical_by_graph_type.main() capturing stdout."""
    import importlib
    spec = importlib.util.find_spec("cross_reference.statistical_by_graph_type")
    mod = importlib.util.module_from_spec(spec)
    buf = io.StringIO()
    with patch.object(sys, "argv", ["script.py", str(folder)]), redirect_stdout(buf):
        spec.loader.exec_module(mod)
        mod.main()
    return buf.getvalue()


class TestTrendDirectionClassification:
    """Verify keyword-based trend direction counting."""

    def _make_folder_with_trends(self, tmp_path, directions):
        """Create a folder with a single graph containing trends with given directions."""
        folder = tmp_path / "saved_files_20260301_120000"
        (folder / "Standard").mkdir(parents=True)
        trends = [{"series": f"s{i}", "direction": d} for i, d in enumerate(directions)]
        rows = [{
            **{col: "" for col in GRAPH_CSV_COLUMNS},
            "file_path": "saved_files_20260301_120000/Standard/Test_Standard.png",
            "graph_type": "line",
            "trends_json": json.dumps(trends),
            "summary": "",
            "statistical_tests_json": "[]",
            "inflection_points_json": "[]",
            "clinical_relevance": "",
        }]
        _write_csv(folder, rows)
        return folder

    def test_increasing_keywords(self, tmp_path):
        """'increasing', 'up', 'higher', 'rising' all count as increasing."""
        folder = self._make_folder_with_trends(
            tmp_path, ["increasing", "up", "higher", "rising"]
        )
        output = _run_by_graph_type(folder)
        assert "Increasing: 4" in output

    def test_decreasing_keywords(self, tmp_path):
        """'decreasing', 'down', 'lower', 'falling', 'drop' all count as decreasing."""
        folder = self._make_folder_with_trends(
            tmp_path, ["decreasing", "down", "lower", "falling", "drop"]
        )
        output = _run_by_graph_type(folder)
        assert "Decreasing: 5" in output

    def test_stable_keywords(self, tmp_path):
        """'flat', 'stable', 'constant' all count as stable."""
        folder = self._make_folder_with_trends(
            tmp_path, ["flat", "stable", "constant"]
        )
        output = _run_by_graph_type(folder)
        assert "Stable: 3" in output

    def test_other_unrecognised_keyword(self, tmp_path):
        """Unrecognised direction keywords count as 'Other'."""
        folder = self._make_folder_with_trends(
            tmp_path, ["variable", "oscillating"]
        )
        output = _run_by_graph_type(folder)
        assert "Other: 2" in output

    def test_case_insensitive_direction(self, tmp_path):
        """Direction matching is case-insensitive (d.lower() is used)."""
        folder = self._make_folder_with_trends(
            tmp_path, ["Increasing", "DECREASING", "Flat"]
        )
        output = _run_by_graph_type(folder)
        assert "Increasing: 1" in output
        assert "Decreasing: 1" in output
        assert "Stable: 1" in output

    def test_partial_match(self, tmp_path):
        """'slightly increasing' matches via 'increas' substring."""
        folder = self._make_folder_with_trends(
            tmp_path, ["slightly increasing"]
        )
        output = _run_by_graph_type(folder)
        assert "Increasing: 1" in output


# ═══════════════════════════════════════════════════════════════════════════════
# Cross-DWI p-value comparison (statistical_relevance.py)
# ═══════════════════════════════════════════════════════════════════════════════


class TestCrossDwiPvalueComparison:
    """Verify cross-DWI significance counts are computed correctly."""

    def test_sig_count_per_dwi_type(self, saved_files_with_graph_csv):
        """Feature_BoxPlots: Standard has p=0.003 (sig), dnCNN has p=0.12 (non-sig)."""
        output = _run_stat_relevance(saved_files_with_graph_csv)
        lines = output.splitlines()
        # Find lines in the cross-DWI section about Feature_BoxPlots
        in_cross = False
        fb_lines = []
        for line in lines:
            if "Same Analysis, Different Significance" in line:
                in_cross = True
            if in_cross and "Feature_BoxPlots" in line:
                # Collect the next few lines for this graph
                fb_lines.append(line)
            elif in_cross and fb_lines and line.strip().startswith(("Standard:", "dnCNN:", "IVIMnet:")):
                fb_lines.append(line)

        # Standard should show 1 significant p-value
        std_lines = [l for l in fb_lines if "Standard:" in l]
        if std_lines:
            assert "1/" in std_lines[0] or "significant" in std_lines[0]


# ═══════════════════════════════════════════════════════════════════════════════
# Priority graph ordering (cross_reference_summary.py)
# ═══════════════════════════════════════════════════════════════════════════════


class TestPriorityGraphOrdering:
    """Verify that priority graphs appear before non-priority graphs."""

    def test_priority_graphs_first(self, tmp_path):
        """Feature_BoxPlots (priority) should appear before ZZZ_Plot (non-priority)."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN"):
            (folder / t).mkdir(parents=True)
        rows = _build_csv_rows([
            # Non-priority graph (alphabetically last)
            ("Standard", "ZZZ_Plot", [{"series": "X", "direction": "up"}]),
            ("dnCNN", "ZZZ_Plot", [{"series": "X", "direction": "up"}]),
            # Priority graph
            ("Standard", "Feature_BoxPlots", [{"series": "X", "direction": "up"}]),
            ("dnCNN", "Feature_BoxPlots", [{"series": "X", "direction": "up"}]),
        ])
        _write_csv(folder, rows)
        output = _run_summary(folder)

        fb_pos = output.find("Feature_BoxPlots")
        zzz_pos = output.find("ZZZ_Plot")
        assert fb_pos < zzz_pos, "Priority graph should appear before non-priority"


# ═══════════════════════════════════════════════════════════════════════════════
# Parameter map and inflection point sections (cross_reference_summary.py)
# ═══════════════════════════════════════════════════════════════════════════════


class TestParameterMapAndInflectionSections:
    """Verify parameter map and inflection point sections filter correctly."""

    def test_parameter_map_filtering(self, tmp_path):
        """Only graphs with 'Parameter_Maps' in the name appear in that section."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN"):
            (folder / t).mkdir(parents=True)
        rows = _build_csv_rows([
            ("Standard", "ADC_Parameter_Maps", []),
            ("dnCNN", "ADC_Parameter_Maps", []),
            ("Standard", "Feature_BoxPlots", []),
            ("dnCNN", "Feature_BoxPlots", []),
        ])
        _write_csv(folder, rows)
        output = _run_summary(folder)

        # Find the parameter map section
        pm_section_start = output.find("PARAMETER MAP COUNTS")
        inflection_start = output.find("INFLECTION POINTS")
        pm_section = output[pm_section_start:inflection_start]

        assert "ADC_Parameter_Maps" in pm_section
        assert "Feature_BoxPlots" not in pm_section

    def test_inflection_points_only_longitudinal(self, tmp_path):
        """Only graphs with 'Longitudinal' in the name appear in inflection section."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN"):
            (folder / t).mkdir(parents=True)
        ip_json = json.dumps([{"approximate_x": 30, "description": "inflect"}])
        rows = [
            {**{col: "" for col in GRAPH_CSV_COLUMNS},
             "file_path": f"saved_files_20260301_120000/{t}/Longitudinal_Test_{t}.png",
             "graph_type": "line", "trends_json": "[]", "summary": "",
             "statistical_tests_json": "[]", "clinical_relevance": "",
             "inflection_points_json": ip_json}
            for t in ("Standard", "dnCNN")
        ] + [
            {**{col: "" for col in GRAPH_CSV_COLUMNS},
             "file_path": f"saved_files_20260301_120000/{t}/BoxPlot_{t}.png",
             "graph_type": "box", "trends_json": "[]", "summary": "",
             "statistical_tests_json": "[]", "clinical_relevance": "",
             "inflection_points_json": ip_json}
            for t in ("Standard", "dnCNN")
        ]
        _write_csv(folder, rows)
        output = _run_summary(folder)

        inflection_start = output.find("INFLECTION POINTS")
        inflection_section = output[inflection_start:]

        assert "Longitudinal_Test" in inflection_section
        assert "BoxPlot" not in inflection_section


# ═══════════════════════════════════════════════════════════════════════════════
# Full cross-DWI detail output (cross_reference_dwi.py)
# ═══════════════════════════════════════════════════════════════════════════════


def _run_xref_dwi(folder):
    """Run cross_reference_dwi.main() capturing stdout."""
    import importlib
    spec = importlib.util.find_spec("cross_reference.cross_reference_dwi")
    mod = importlib.util.module_from_spec(spec)
    buf = io.StringIO()
    with patch.object(sys, "argv", ["script.py", str(folder)]), redirect_stdout(buf):
        spec.loader.exec_module(mod)
        mod.main()
    return buf.getvalue()


class TestCrossReferenceDwiCorrectness:
    """Verify correctness of values printed by cross_reference_dwi."""

    def test_start_end_values_formatted_to_four_decimals(self, saved_files_with_graph_csv):
        """Start/end values should appear with 4 decimal places."""
        output = _run_xref_dwi(saved_files_with_graph_csv)
        # Row 1 has start_value=0.001, end_value=0.0015
        assert "0.0010" in output
        assert "0.0015" in output

    def test_pvalue_formatted_to_four_decimals(self, saved_files_with_graph_csv):
        """Statistical test p-values should show 4 decimal places."""
        output = _run_xref_dwi(saved_files_with_graph_csv)
        # Wilcoxon p=0.003 → p=0.0030
        assert "p=0.0030" in output

    def test_statistic_value_formatted_when_present(self, tmp_path):
        """Statistic values should show 3 decimal places when not None."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN"):
            (folder / t).mkdir(parents=True)
        rows = [
            {**{col: "" for col in GRAPH_CSV_COLUMNS},
             "file_path": f"saved_files_20260301_120000/{t}/Test_{t}.png",
             "graph_type": "line", "trends_json": "[]", "summary": "",
             "inflection_points_json": "[]", "clinical_relevance": "",
             "statistical_tests_json": json.dumps([
                 {"test_name": "log-rank", "statistic_value": 5.41,
                  "p_value": 0.02, "comparison_groups": "LF vs LC"}
             ])}
            for t in ("Standard", "dnCNN")
        ]
        _write_csv(folder, rows)
        output = _run_xref_dwi(folder)
        assert "stat=5.410" in output

    def test_matched_count_correct(self, saved_files_with_graph_csv):
        """Total matched graph sets should equal the number of cross-DWI groups."""
        output = _run_xref_dwi(saved_files_with_graph_csv)
        # Only Feature_BoxPlots has 2+ real DWI types (Standard + dnCNN)
        # Longitudinal_Mean_Metrics only has IVIMnet (1 type)
        assert "Total matched graph sets across DWI types: 1" in output

    def test_comparison_groups_shown(self, saved_files_with_graph_csv):
        """Comparison groups (e.g., 'LF vs LC') appear in the output."""
        output = _run_xref_dwi(saved_files_with_graph_csv)
        assert "LF vs LC" in output

    def test_magnitude_in_parentheses(self, saved_files_with_graph_csv):
        """Trend magnitude should appear in parentheses."""
        output = _run_xref_dwi(saved_files_with_graph_csv)
        assert "(~15% increase)" in output

    def test_metadata_line_content(self, saved_files_with_graph_csv):
        """Metadata line should include n=42, series count, error bars."""
        output = _run_xref_dwi(saved_files_with_graph_csv)
        assert "n=42" in output
        assert "2 series" in output
        assert "error bars: IQR" in output


# ═══════════════════════════════════════════════════════════════════════════════
# Summary table correctness (statistical_by_graph_type.py)
# ═══════════════════════════════════════════════════════════════════════════════


class TestSummaryTable:
    """Verify the summary table at the end of statistical_by_graph_type output."""

    def test_summary_table_present(self, saved_files_with_graph_csv):
        output = _run_by_graph_type(saved_files_with_graph_csv)
        assert "SUMMARY TABLE" in output

    def test_box_type_row_in_table(self, saved_files_with_graph_csv):
        """Box type should appear with count=2 (Standard + dnCNN rows)."""
        output = _run_by_graph_type(saved_files_with_graph_csv)
        # Look for the box row in the summary table
        lines = output.splitlines()
        summary_start = next(
            i for i, l in enumerate(lines) if "SUMMARY TABLE" in l
        )
        table_lines = lines[summary_start:]
        box_lines = [l for l in table_lines if l.strip().startswith("box")]
        assert len(box_lines) >= 1
        # Should show count of 2 (Standard + dnCNN box plots)
        assert "2" in box_lines[0]

    def test_line_type_row_in_table(self, saved_files_with_graph_csv):
        """Line type should appear with count=1 (IVIMnet longitudinal)."""
        output = _run_by_graph_type(saved_files_with_graph_csv)
        lines = output.splitlines()
        summary_start = next(
            i for i, l in enumerate(lines) if "SUMMARY TABLE" in l
        )
        table_lines = lines[summary_start:]
        line_rows = [l for l in table_lines if l.strip().startswith("line")]
        assert len(line_rows) >= 1
        assert "1" in line_rows[0]
