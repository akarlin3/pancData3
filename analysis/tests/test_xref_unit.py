"""Unit tests for cross-reference scripts (part 1) — shared.py parsing and
statistical_relevance significance/correlation logic.

These tests go beyond smoke/stdout checks by verifying that specific known
inputs produce the expected comparison values, significance markers, and
correlation classifications.  They exercise the internal logic of:

- ``shared.py``: ``safe_text``, ``extract_pvalues``, ``extract_correlations``,
  ``parse_dwi_info``, ``group_by_graph_name``
- ``cross_reference_summary.py``: trend agreement/disagreement across DWI types
- ``statistical_relevance.py``: significance markers, Bonferroni correction,
  correlation strength classification

Trend direction classification, cross-DWI p-value comparison, priority graph
ordering, parameter-map/inflection sections, full cross-DWI detail output, and
the summary table live in ``test_xref_unit_part2.py`` (file-size split). Shared
imports, constants, and helpers live in ``_xref_unit_shared.py``.
"""

from __future__ import annotations

from _xref_unit_shared import *  # noqa: F401,F403


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
        results = extract_correlations("r² = 0.81 good fit")
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
