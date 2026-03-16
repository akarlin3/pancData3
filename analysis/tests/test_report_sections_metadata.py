"""Tests for report_sections.metadata module.

Validates the metadata section builders:
- Cover page HTML generation
- Part break divider
- Table of contents
- Publication header
- Data availability statement
- Table/figure index generation
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.metadata import (
    _section_cover_page,
    _part_break,
    _section_print_toc,
    _section_publication_header,
    _section_data_availability,
    _section_table_index,
    _section_figure_index,
)


class TestCoverPage:
    def test_returns_html_list(self):
        result = _section_cover_page("20260301_120000", ["Standard", "dnCNN"], 42)
        assert isinstance(result, list)
        assert len(result) > 0

    def test_contains_timestamp(self):
        result = _section_cover_page("20260301_120000", ["Standard"], 10)
        html = "\n".join(result)
        assert "20260301_120000" in html

    def test_contains_dwi_types(self):
        result = _section_cover_page("ts", ["Standard", "dnCNN"], 5)
        html = "\n".join(result)
        assert "Standard" in html
        assert "dnCNN" in html

    def test_contains_graph_count(self):
        result = _section_cover_page("ts", ["Standard"], 67)
        html = "\n".join(result)
        assert "67" in html

    def test_zero_graphs_omits_count(self):
        result = _section_cover_page("ts", ["Standard"], 0)
        html = "\n".join(result)
        assert "Graphs analysed" not in html

    def test_empty_dwi_types(self):
        result = _section_cover_page("ts", [], 0)
        assert isinstance(result, list)

    def test_html_escaping(self):
        # Ensure special characters are escaped
        result = _section_cover_page("<script>alert('xss')</script>", [], 0)
        html = "\n".join(result)
        assert "<script>" not in html


class TestPartBreak:
    def test_returns_string(self):
        result = _part_break("Results")
        assert isinstance(result, str)

    def test_contains_label(self):
        result = _part_break("Data")
        assert "Data" in result

    def test_html_class(self):
        result = _part_break("Test")
        assert "part-break" in result

    def test_escapes_html(self):
        result = _part_break("<b>bold</b>")
        assert "<b>" not in result


class TestPrintToc:
    def test_returns_html_list(self):
        result = _section_print_toc()
        assert isinstance(result, list)
        assert len(result) > 0

    def test_contains_toc_id(self):
        html = "\n".join(_section_print_toc())
        assert 'id="toc"' in html

    def test_contains_section_links(self):
        html = "\n".join(_section_print_toc())
        assert "<a href=" in html


class TestPublicationHeader:
    def test_returns_html_list(self):
        result = _section_publication_header()
        assert isinstance(result, list)
        assert len(result) > 0

    def test_contains_placeholders(self):
        html = "\n".join(_section_publication_header())
        assert "Author" in html
        assert "Affiliations" in html
        assert "IRB" in html


class TestDataAvailability:
    def test_returns_html_list(self):
        result = _section_data_availability()
        assert isinstance(result, list)
        assert len(result) > 0

    def test_mentions_hipaa(self):
        html = "\n".join(_section_data_availability())
        assert "HIPAA" in html

    def test_mentions_mit_license(self):
        html = "\n".join(_section_data_availability())
        assert "MIT" in html


class TestTableIndex:
    def test_empty_when_no_tables(self):
        # Reset numbering state
        from report.report_formatters import get_numbering
        numbering = get_numbering()
        numbering.table_titles.clear()
        result = _section_table_index()
        assert result == []


class TestFigureIndex:
    def test_empty_when_no_figures(self):
        from report.report_formatters import get_numbering
        numbering = get_numbering()
        numbering.figure_titles.clear()
        result = _section_figure_index()
        assert result == []


# ── Edge cases: metadata section robustness ──


class TestCoverPageEdgeCases:
    def test_single_dwi_type(self):
        result = _section_cover_page("ts", ["Standard"], 5)
        html = "\n".join(result)
        assert "Standard" in html

    def test_all_three_dwi_types(self):
        result = _section_cover_page("ts", ["Standard", "dnCNN", "IVIMnet"], 100)
        html = "\n".join(result)
        assert "Standard" in html
        assert "dnCNN" in html
        assert "IVIMnet" in html

    def test_large_graph_count(self):
        result = _section_cover_page("ts", ["Standard"], 999999)
        html = "\n".join(result)
        assert "999999" in html

    def test_special_chars_in_timestamp(self):
        result = _section_cover_page("2026/03/16 12:00 & <test>", ["Standard"], 1)
        html = "\n".join(result)
        assert "&lt;" in html or "<test>" not in html

    def test_unicode_in_dwi_type(self):
        result = _section_cover_page("ts", ["\u0394-DWI"], 1)
        html = "\n".join(result)
        assert isinstance(html, str)

    def test_negative_graph_count(self):
        result = _section_cover_page("ts", ["Standard"], -1)
        assert isinstance(result, list)

    def test_empty_timestamp(self):
        result = _section_cover_page("", ["Standard"], 5)
        assert isinstance(result, list)
        assert len(result) > 0


class TestPartBreakEdgeCases:
    def test_empty_label(self):
        result = _part_break("")
        assert "part-break" in result

    def test_unicode_label(self):
        result = _part_break("\u0394 Analysis")
        assert "part-break" in result

    def test_long_label(self):
        result = _part_break("A" * 200)
        assert "part-break" in result

    def test_html_entities_in_label(self):
        result = _part_break("Results & Discussion")
        assert "&amp;" in result or "Results" in result


class TestPrintTocEdgeCases:
    def test_contains_navigation_groups(self):
        html = "\n".join(_section_print_toc())
        # Should have multiple section links
        assert html.count("<a href=") >= 3

    def test_toc_structure(self):
        result = _section_print_toc()
        html = "\n".join(result)
        assert "Table of Contents" in html or "toc" in html.lower()


class TestTableIndexEdgeCases:
    def test_with_table_entries(self):
        from report.report_formatters import get_numbering
        numbering = get_numbering()
        numbering.table_titles.clear()
        numbering.table_titles.append((1, "Test Table 1"))
        numbering.table_titles.append((2, "Test Table 2"))
        result = _section_table_index()
        html = "\n".join(result)
        assert "Test Table 1" in html
        assert "Test Table 2" in html
        numbering.table_titles.clear()

    def test_table_titles_with_html_chars(self):
        from report.report_formatters import get_numbering
        numbering = get_numbering()
        numbering.table_titles.clear()
        numbering.table_titles.append((1, "<b>Bold Table</b>"))
        result = _section_table_index()
        html = "\n".join(result)
        assert "<b>Bold Table</b>" not in html or "&lt;b&gt;" in html
        numbering.table_titles.clear()


class TestFigureIndexEdgeCases:
    def test_with_figure_entries(self):
        from report.report_formatters import get_numbering
        numbering = get_numbering()
        numbering.figure_titles.clear()
        numbering.figure_titles.append((1, "ADC Parameter Map"))
        numbering.figure_titles.append((2, "Survival Curve"))
        result = _section_figure_index()
        html = "\n".join(result)
        assert "ADC Parameter Map" in html
        assert "Survival Curve" in html
        numbering.figure_titles.clear()

    def test_single_figure(self):
        from report.report_formatters import get_numbering
        numbering = get_numbering()
        numbering.figure_titles.clear()
        numbering.figure_titles.append((1, "Only Figure"))
        result = _section_figure_index()
        assert len(result) > 0
        numbering.figure_titles.clear()


class TestDataAvailabilityEdgeCases:
    def test_contains_code_availability(self):
        html = "\n".join(_section_data_availability())
        assert "code" in html.lower() or "pipeline" in html.lower()

    def test_is_nonempty(self):
        result = _section_data_availability()
        assert len(result) > 0
        html = "\n".join(result)
        assert len(html) > 50


class TestPublicationHeaderEdgeCases:
    def test_contains_correspondence(self):
        html = "\n".join(_section_publication_header())
        assert "Corresponding" in html or "correspondence" in html.lower() or "Author" in html
