"""Tests for report formatting helpers and small utility functions.

Covers:
- _normalize_series_name / _build_normalised_series_map: series name matching
- _sig_tag: p-value to significance star mapping
- _section: Markdown heading generation
- _forest_plot_cell: forest plot HTML generation
- _effect_size_class / _effect_size_label: effect size classification
- _copy_button / _manuscript_sentence: copy-to-clipboard helpers
- _figure_caption / figure index / table index: numbering helpers
"""

from __future__ import annotations

import pytest  # type: ignore

from report.generate_report import (  # type: ignore
    _copy_button,
    _effect_size_class,
    _effect_size_label,
    _figure_caption,
    _forest_plot_cell,
    _manuscript_sentence,
    _section,
    _section_figure_index,
    _section_table_index,
    _sig_tag,
    get_numbering,
    reset_numbering,
)
from report.sections._helpers import (  # type: ignore
    _normalize_series_name,
    _build_normalised_series_map,
    _best_display_name,
)


# ---------------------------------------------------------------------------
# _normalize_series_name / _build_normalised_series_map
# ---------------------------------------------------------------------------


class TestNormaliseSeriesName:
    """Verify that free-form vision-API series names are normalised."""

    def test_identical_names_match(self):
        assert _normalize_series_name("Mean ADC") == _normalize_series_name("Mean ADC")

    def test_delimiter_variants_match(self):
        """Comma, dash, and parenthetical separators should all normalise."""
        a = _normalize_series_name("Mean D - Local Control")
        b = _normalize_series_name("Mean D (Local Control)")
        c = _normalize_series_name("Mean D, Local Control")
        assert a == b == c

    def test_vs_delimiter_matches(self):
        """'vs' separator should be treated as a delimiter."""
        a = _normalize_series_name("LC vs D95")
        b = _normalize_series_name("LC - D95")
        assert a == b

    def test_case_insensitive(self):
        assert _normalize_series_name("mean adc") == _normalize_series_name("Mean ADC")

    def test_delta_unicode_normalised(self):
        a = _normalize_series_name("\u0394 D (%)")
        b = _normalize_series_name("delta D (%)")
        assert a == b


class TestBuildNormalisedSeriesMap:
    """Verify cross-DWI series matching with normalised names."""

    def test_mismatched_names_merged(self):
        """Series with different delimiters should be merged."""
        all_trends = {
            "Standard": [{"series": "Mean D - Local Control", "direction": "increasing", "description": ""}],
            "dnCNN":    [{"series": "Mean D (Local Control)", "direction": "increasing", "description": ""}],
            "IVIMnet":  [{"series": "Mean D, Local Control",  "direction": "decreasing", "description": ""}],
        }
        norm_map = _build_normalised_series_map(all_trends)
        # All three should be in one group
        assert len(norm_map) == 1
        key = list(norm_map.keys())[0]
        assert len(norm_map[key]) == 3

    def test_best_display_name_picks_longest(self):
        all_trends = {
            "Standard": [{"series": "Mean D - Local Control", "direction": "up", "description": ""}],
            "IVIMnet":  [{"series": "Mean D (Local Control)", "direction": "up", "description": ""}],
        }
        norm_map = _build_normalised_series_map(all_trends)
        key = list(norm_map.keys())[0]
        display = _best_display_name(all_trends, key)
        # Should pick the longest raw name
        assert len(display) >= len("Mean D - Local Control")

    def test_disjoint_series_separate(self):
        """Completely different series remain separate entries."""
        all_trends = {
            "Standard": [
                {"series": "Mean ADC", "direction": "increasing", "description": ""},
                {"series": "Mean D",   "direction": "decreasing", "description": ""},
            ],
        }
        nmap = _build_normalised_series_map(all_trends)
        assert len(nmap) == 2


# ---------------------------------------------------------------------------
# _sig_tag
# ---------------------------------------------------------------------------


class TestSigTag:
    """Verify p-value to significance star mapping."""

    def test_three_stars(self):
        """p < 0.001 gets triple stars (highly significant)."""
        assert _sig_tag(0.0001) == "***"
        assert _sig_tag(0.0009) == "***"

    def test_two_stars(self):
        """0.001 <= p < 0.01 gets double stars."""
        assert _sig_tag(0.001) == "**"
        assert _sig_tag(0.005) == "**"
        assert _sig_tag(0.009) == "**"

    def test_one_star(self):
        """0.01 <= p < 0.05 gets a single star."""
        assert _sig_tag(0.01) == "*"
        assert _sig_tag(0.049) == "*"

    def test_not_significant(self):
        """p >= 0.05 gets no star (empty string)."""
        assert _sig_tag(0.05) == ""
        assert _sig_tag(0.5) == ""
        assert _sig_tag(1.0) == ""

    def test_boundary_values(self):
        """Boundary at 0.001, 0.01, 0.05 (exclusive lower bound)."""
        assert _sig_tag(0.001) == "**"   # exactly 0.001 -> **
        assert _sig_tag(0.01) == "*"     # exactly 0.01 -> *
        assert _sig_tag(0.05) == ""      # exactly 0.05 -> not significant


# ---------------------------------------------------------------------------
# _section
# ---------------------------------------------------------------------------


class TestSection:
    """Verify Markdown/HTML heading generation."""

    def test_default_h2(self):
        """Default heading level is 2 (##)."""
        result = _section("My Section")
        assert result == "\n## My Section\n"

    def test_h3(self):
        """Explicit level=3 produces a ### heading."""
        result = _section("Subsection", level=3)
        assert result == "\n### Subsection\n"

    def test_h1(self):
        """level=1 produces a top-level # heading."""
        result = _section("Title", level=1)
        assert result == "\n# Title\n"


# ---------------------------------------------------------------------------
# _forest_plot_cell
# ---------------------------------------------------------------------------


class TestForestPlotCell:
    """Verify forest plot HTML generation."""

    def test_returns_html(self):
        """Forest plot cell should return valid HTML."""
        result = _forest_plot_cell(1.5, 1.1, 2.0, 0.03)
        assert "<div" in result
        assert "forest-row" in result
        assert "forest-point-sig" in result  # p < 0.05

    def test_non_significant_point(self):
        """Non-significant HR should use the ns class."""
        result = _forest_plot_cell(1.1, 0.8, 1.5, 0.12)
        assert "forest-point-ns" in result

    def test_reference_line_present(self):
        """The HR=1.0 reference line should always be present."""
        result = _forest_plot_cell(0.5, 0.3, 0.8, 0.01)
        assert "forest-ref" in result


# ---------------------------------------------------------------------------
# _effect_size_class / _effect_size_label
# ---------------------------------------------------------------------------


class TestEffectSize:
    """Verify effect size classification helpers."""

    def test_large_effect(self):
        assert _effect_size_class(0.9) == "effect-lg"
        assert _effect_size_label(0.9) == "Large"

    def test_medium_effect(self):
        assert _effect_size_class(0.6) == "effect-md"
        assert _effect_size_label(0.6) == "Medium"

    def test_small_effect(self):
        assert _effect_size_class(0.3) == "effect-sm"
        assert _effect_size_label(0.3) == "Small"

    def test_negative_values_use_absolute(self):
        """Negative effect sizes should be classified by absolute value."""
        assert _effect_size_class(-0.9) == "effect-lg"
        assert _effect_size_label(-0.5) == "Medium"


# ---------------------------------------------------------------------------
# Copy button and manuscript sentence helpers
# ---------------------------------------------------------------------------


class TestCopyHelpers:
    """Verify copy-to-clipboard helper functions."""

    def test_copy_button_with_target(self):
        """Copy button generates onclick with target ID."""
        btn = _copy_button("my-target")
        assert "my-target" in btn
        assert "copy-btn" in btn

    def test_copy_button_without_target(self):
        """Copy button without target uses copyText(this)."""
        btn = _copy_button()
        assert "copyText(this)" in btn

    def test_manuscript_sentence_has_copy(self):
        """Manuscript sentence includes copy button and data-copy attr."""
        html = _manuscript_sentence("Test finding here.")
        assert "Test finding here." in html
        assert "data-copy" in html
        assert "copy-btn" in html

    def test_manuscript_sentence_escapes_html(self):
        """HTML characters are escaped in manuscript sentences."""
        html = _manuscript_sentence("p < 0.05 & HR > 1")
        assert "&lt;" in html or "p < 0.05" not in html.split("data-copy")[0]


# ---------------------------------------------------------------------------
# Figure Caption
# ---------------------------------------------------------------------------


class TestFigureCaption:
    """Verify the _figure_caption formatting helper."""

    def test_basic_caption(self):
        """Caption includes figure number and title."""
        reset_numbering()
        cap = _figure_caption("ADC Map")
        assert "Figure 1" in cap
        assert "ADC Map" in cap

    def test_sequential_numbering(self):
        """Multiple captions are sequentially numbered."""
        reset_numbering()
        cap1 = _figure_caption("First")
        cap2 = _figure_caption("Second")
        assert "Figure 1" in cap1
        assert "Figure 2" in cap2

    def test_with_description(self):
        """Caption includes optional description."""
        reset_numbering()
        cap = _figure_caption("Title", "Additional details here.")
        assert "Additional details here" in cap


# ---------------------------------------------------------------------------
# Figure Index (List of Figures)
# ---------------------------------------------------------------------------


class TestFigureIndex:
    """Verify the List of Figures section."""

    def test_empty_when_no_figures(self):
        """No figure index if no figures were captioned."""
        reset_numbering()
        result = _section_figure_index()
        assert result == []

    def test_lists_captioned_figures(self):
        """Lists all figures that received captions."""
        reset_numbering()
        _figure_caption("ADC Parameter Map")
        _figure_caption("Longitudinal Trends")
        result = _section_figure_index()
        html = "\n".join(result)
        assert "List of Figures" in html
        assert "ADC Parameter Map" in html
        assert "Longitudinal Trends" in html


# ---------------------------------------------------------------------------
# Table index
# ---------------------------------------------------------------------------


class TestTableIndex:
    """Verify the List of Tables section."""

    def test_empty_without_tables(self):
        """Returns empty when no tables have been numbered."""
        reset_numbering()
        result = _section_table_index()
        assert result == []

    def test_lists_numbered_tables(self):
        """Lists all tables that were numbered during report generation."""
        reset_numbering()
        from report.report_formatters import _table_caption  # type: ignore
        _table_caption("Test Table One", "description")
        _table_caption("Test Table Two")
        result = _section_table_index()
        html = "\n".join(result)
        assert "Test Table One" in html
        assert "Test Table Two" in html
        assert "Table 1" in html
        assert "Table 2" in html
        reset_numbering()
