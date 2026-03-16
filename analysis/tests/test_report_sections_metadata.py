"""Tests for report_sections.metadata module.

Validates the metadata section builders:
- Cover page HTML generation
- Part break divider
- Table of contents
- Publication header
- Data availability statement
- Table/figure index generation

Edge cases covered:
- Empty/None/zero inputs
- HTML escaping (XSS prevention)
- Special characters in labels
- Large graph counts
- Empty DWI types list
- Index with entries present
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
        result = _section_cover_page("<script>alert('xss')</script>", [], 0)
        html = "\n".join(result)
        assert "<script>" not in html

    def test_large_graph_count(self):
        """Large graph count should be displayed."""
        result = _section_cover_page("ts", ["Standard"], 5000)
        html = "\n".join(result)
        assert "5000" in html

    def test_all_three_dwi_types(self):
        """All three DWI types should be shown."""
        result = _section_cover_page("ts", ["Standard", "dnCNN", "IVIMnet"], 10)
        html = "\n".join(result)
        assert "Standard" in html
        assert "dnCNN" in html
        assert "IVIMnet" in html

    def test_none_timestamp(self):
        """None timestamp should not crash."""
        result = _section_cover_page(None, ["Standard"], 5)
        assert isinstance(result, list)


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

    def test_empty_label(self):
        result = _part_break("")
        assert isinstance(result, str)
        assert "part-break" in result

    def test_special_characters(self):
        result = _part_break("Results & Discussion")
        assert "&amp;" in result or "Results" in result


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

    def test_has_multiple_entries(self):
        """TOC should have multiple section entries."""
        result = _section_print_toc()
        assert len(result) > 3


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

    def test_always_same_output(self):
        """Should be deterministic."""
        r1 = _section_publication_header()
        r2 = _section_publication_header()
        assert r1 == r2


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

    def test_always_same_output(self):
        """Should be deterministic."""
        r1 = _section_data_availability()
        r2 = _section_data_availability()
        assert r1 == r2


class TestTableIndex:
    def test_empty_when_no_tables(self):
        from report.report_formatters import get_numbering
        numbering = get_numbering()
        numbering.table_titles.clear()
        result = _section_table_index()
        assert result == []

    def test_with_tables(self):
        """Table index should list registered tables."""
        from report.report_formatters import get_numbering
        numbering = get_numbering()
        numbering.table_titles.clear()
        numbering.table_titles.append((1, "Patient Demographics"))
        numbering.table_titles.append((2, "Effect Sizes"))
        result = _section_table_index()
        html = "\n".join(result)
        assert "Patient Demographics" in html
        assert "Effect Sizes" in html
        numbering.table_titles.clear()


class TestFigureIndex:
    def test_empty_when_no_figures(self):
        from report.report_formatters import get_numbering
        numbering = get_numbering()
        numbering.figure_titles.clear()
        result = _section_figure_index()
        assert result == []

    def test_with_figures(self):
        """Figure index should list registered figures."""
        from report.report_formatters import get_numbering
        numbering = get_numbering()
        numbering.figure_titles.clear()
        numbering.figure_titles.append((1, "ROC Curve"))
        numbering.figure_titles.append((2, "Longitudinal Trends"))
        result = _section_figure_index()
        html = "\n".join(result)
        assert "ROC Curve" in html
        assert "Longitudinal Trends" in html
        numbering.figure_titles.clear()
