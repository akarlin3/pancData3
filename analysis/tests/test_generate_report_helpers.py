"""Tests for report helper and formatting functions.

Covers:
- _normalize_series_name / _build_normalised_series_map: series name matching
- _sig_tag: p-value to significance star mapping
- _section: Markdown heading generation
- _effect_size_class / _effect_size_label: effect size classification
- _copy_button / _manuscript_sentence: clipboard helper HTML
- _figure_caption: numbered figure caption generation
"""

from __future__ import annotations

import pytest  # type: ignore

from report.generate_report import (  # type: ignore
    _copy_button,
    _effect_size_class,
    _effect_size_label,
    _figure_caption,
    _manuscript_sentence,
    _section,
    _sig_tag,
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

    def test_distinct_series_not_merged(self):
        """Genuinely different series should remain separate."""
        all_trends = {
            "Standard": [
                {"series": "Mean ADC", "direction": "increasing", "description": ""},
                {"series": "Mean D",   "direction": "decreasing", "description": ""},
            ],
        }
        norm_map = _build_normalised_series_map(all_trends)
        assert len(norm_map) == 2

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
    """Verify Markdown section header generation."""

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
