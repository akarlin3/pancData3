"""Unit tests for report_formatters.py.

Validates HTML escaping, significance markers, CSS class mapping, DWI badges,
trend tags, navigation bar, stat cards, forest plot cells, effect size helpers,
and consensus determination.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Ensure analysis/ is on sys.path
ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report_formatters import (
    CSS,
    HTML_TEMPLATE,
    NAV_SECTIONS,
    _dwi_badge,
    _effect_size_class,
    _effect_size_label,
    _esc,
    _forest_plot_cell,
    _get_consensus,
    _h2,
    _nav_bar,
    _section,
    _sig_class,
    _sig_tag,
    _stat_card,
    _trend_tag,
)


# ── _esc ──────────────────────────────────────────────────────────────────


class TestEsc:
    def test_escapes_angle_brackets(self):
        assert "&lt;" in _esc("<script>")
        assert "&gt;" in _esc("</script>")

    def test_escapes_ampersand(self):
        assert "&amp;" in _esc("a & b")

    def test_passthrough_plain_text(self):
        assert _esc("hello world") == "hello world"

    def test_non_string_coerced(self):
        assert _esc(42) == "42"


# ── _section ──────────────────────────────────────────────────────────────


class TestSection:
    def test_default_h2(self):
        assert _section("Foo") == "\n## Foo\n"

    def test_h3_level(self):
        assert _section("Bar", 3) == "\n### Bar\n"


# ── _sig_tag ──────────────────────────────────────────────────────────────


class TestSigTag:
    def test_triple_star(self):
        assert _sig_tag(0.0005) == "***"

    def test_double_star(self):
        assert _sig_tag(0.005) == "**"

    def test_single_star(self):
        assert _sig_tag(0.03) == "*"

    def test_no_star(self):
        assert _sig_tag(0.1) == ""

    def test_boundary_001(self):
        # p = 0.001 is NOT < 0.001, so should be "**"
        assert _sig_tag(0.001) == "**"

    def test_boundary_01(self):
        # p = 0.01 is NOT < 0.01, so should be "*"
        assert _sig_tag(0.01) == "*"

    def test_boundary_05(self):
        assert _sig_tag(0.05) == ""


# ── _sig_class ────────────────────────────────────────────────────────────


class TestSigClass:
    def test_sig3(self):
        assert _sig_class(0.0005) == "sig-3"

    def test_sig2(self):
        assert _sig_class(0.005) == "sig-2"

    def test_sig1(self):
        assert _sig_class(0.03) == "sig-1"

    def test_empty(self):
        assert _sig_class(0.1) == ""


# ── _dwi_badge ────────────────────────────────────────────────────────────


class TestDwiBadge:
    def test_standard(self):
        badge = _dwi_badge("Standard")
        assert "badge-standard" in badge
        assert "Standard" in badge

    def test_dncnn(self):
        assert "badge-dncnn" in _dwi_badge("dnCNN")

    def test_ivimnet(self):
        assert "badge-ivimnet" in _dwi_badge("IVIMnet")

    def test_unknown_type(self):
        assert "badge-root" in _dwi_badge("Unknown")

    def test_escapes_html(self):
        badge = _dwi_badge("<script>")
        assert "<script>" not in badge
        assert "&lt;script&gt;" in badge


# ── _trend_tag ────────────────────────────────────────────────────────────


class TestTrendTag:
    def test_increasing(self):
        tag = _trend_tag("increasing")
        assert "trend-incr" in tag
        assert "\u2191" in tag  # up arrow

    def test_decreasing(self):
        tag = _trend_tag("decreasing trend")
        assert "trend-decr" in tag
        assert "\u2193" in tag  # down arrow

    def test_stable(self):
        tag = _trend_tag("flat")
        assert "trend-flat" in tag
        assert "\u2192" in tag  # right arrow

    def test_non_monotonic(self):
        tag = _trend_tag("U-shaped")
        assert "trend-nm" in tag

    def test_case_insensitive(self):
        tag = _trend_tag("INCREASING")
        assert "trend-incr" in tag

    def test_keyword_up(self):
        assert "trend-incr" in _trend_tag("going up")

    def test_keyword_drop(self):
        assert "trend-decr" in _trend_tag("sharp drop")


# ── _nav_bar ──────────────────────────────────────────────────────────────


class TestNavBar:
    def test_returns_nav_element(self):
        nav = _nav_bar()
        assert nav.startswith('<nav class="toc">')
        assert nav.endswith("</nav>")

    def test_contains_all_sections(self):
        nav = _nav_bar()
        for anchor, label in NAV_SECTIONS:
            assert f'href="#{anchor}"' in nav
            assert label in nav


# ── _h2 ──────────────────────────────────────────────────────────────────


class TestH2:
    def test_basic(self):
        result = _h2("Results", "results")
        assert result == '<h2 id="results">Results</h2>'

    def test_escapes_text(self):
        result = _h2("A & B", "ab")
        assert "&amp;" in result


# ── _stat_card ────────────────────────────────────────────────────────────


class TestStatCard:
    def test_basic_card(self):
        card = _stat_card("Best AUC", "0.843")
        assert "stat-card" in card
        assert "Best AUC" in card
        assert "0.843" in card

    def test_with_subtitle(self):
        card = _stat_card("Best AUC", "0.843", "across all DWI")
        assert "across all DWI" in card

    def test_no_subtitle(self):
        card = _stat_card("Count", "42")
        # Should not contain the sub div when no subtitle
        assert 'class="sub"' not in card


# ── _forest_plot_cell ────────────────────────────────────────────────────


class TestForestPlotCell:
    def test_returns_div(self):
        cell = _forest_plot_cell(1.5, 1.1, 2.0, 0.03)
        assert "forest-row" in cell
        assert "forest-bar" in cell
        assert "forest-point" in cell

    def test_significant_coloring(self):
        cell = _forest_plot_cell(1.5, 1.1, 2.0, 0.01)
        assert "forest-point-sig" in cell

    def test_nonsignificant_coloring(self):
        cell = _forest_plot_cell(1.0, 0.8, 1.2, 0.2)
        assert "forest-point-ns" in cell

    def test_reference_line(self):
        cell = _forest_plot_cell(1.5, 1.1, 2.0, 0.03)
        assert "forest-ref" in cell


# ── _effect_size_class ───────────────────────────────────────────────────


class TestEffectSizeClass:
    def test_large(self):
        assert _effect_size_class(0.9) == "effect-lg"

    def test_medium(self):
        assert _effect_size_class(0.6) == "effect-md"

    def test_small(self):
        assert _effect_size_class(0.3) == "effect-sm"

    def test_negative_uses_abs(self):
        assert _effect_size_class(-0.9) == "effect-lg"

    def test_boundary_08(self):
        assert _effect_size_class(0.8) == "effect-lg"

    def test_boundary_05(self):
        assert _effect_size_class(0.5) == "effect-md"


# ── _effect_size_label ───────────────────────────────────────────────────


class TestEffectSizeLabel:
    def test_large(self):
        assert _effect_size_label(1.0) == "Large"

    def test_medium(self):
        assert _effect_size_label(0.6) == "Medium"

    def test_small(self):
        assert _effect_size_label(0.2) == "Small"

    def test_negative(self):
        assert _effect_size_label(-0.8) == "Large"


# ── _get_consensus ───────────────────────────────────────────────────────


class TestGetConsensus:
    def test_empty_returns_unknown(self):
        assert _get_consensus([]) == "unknown"

    def test_all_increasing(self):
        assert _get_consensus(["increasing", "higher", "up"]) == "increasing"

    def test_all_decreasing(self):
        assert _get_consensus(["decreasing", "lower"]) == "decreasing"

    def test_mixed_majority_increasing(self):
        assert _get_consensus(["increasing", "increasing", "decreasing"]) == "increasing"

    def test_tied_returns_stable(self):
        assert _get_consensus(["increasing", "decreasing"]) == "stable"

    def test_no_keywords_stable(self):
        assert _get_consensus(["flat", "constant"]) == "stable"

    def test_rising_counts_as_increasing(self):
        assert _get_consensus(["rising", "rising", "decreasing"]) == "increasing"

    def test_falling_counts_as_decreasing(self):
        assert _get_consensus(["falling", "falling", "increasing"]) == "decreasing"

    def test_drop_counts_as_decreasing(self):
        assert _get_consensus(["drop", "drop", "increasing"]) == "decreasing"


# ── Constants ────────────────────────────────────────────────────────────


class TestConstants:
    def test_css_nonempty(self):
        assert len(CSS) > 100

    def test_html_template_has_placeholders(self):
        assert "{title}" in HTML_TEMPLATE
        assert "{body}" in HTML_TEMPLATE

    def test_nav_sections_nonempty(self):
        assert len(NAV_SECTIONS) > 10
