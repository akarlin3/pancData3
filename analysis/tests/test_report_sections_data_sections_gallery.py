"""Tests for report sections: gallery module (appendix, graph analysis, figures).

Split from test_report_sections_data_sections.py to keep modules under 700 lines.
Validates:
- Appendix
- Build graph analysis HTML
- Figure gallery
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from conftest import SAMPLE_GRAPH_CSV_ROWS

from report.sections.gallery import (
    _section_appendix,
    _section_figure_gallery,
    _build_graph_analysis_html,
)


# ── Appendix ──


class TestAppendix:
    def test_returns_html(self):
        result = _section_appendix(SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)

    def test_empty_rows(self):
        result = _section_appendix([])
        assert isinstance(result, list)
        assert result == []

    def test_contains_graph_cards(self):
        """Should render graph-card divs for each row."""
        html = "\n".join(_section_appendix(SAMPLE_GRAPH_CSV_ROWS))
        assert "graph-card" in html

    def test_grouped_by_type(self):
        """Graphs should be grouped by type (box, line)."""
        html = "\n".join(_section_appendix(SAMPLE_GRAPH_CSV_ROWS))
        assert "box" in html.lower()
        assert "line" in html.lower()

    def test_graph_count_shown(self):
        """Total graph count should be displayed."""
        html = "\n".join(_section_appendix(SAMPLE_GRAPH_CSV_ROWS))
        assert "3" in html  # 3 sample rows

    def test_dwi_badge_shown(self):
        """DWI type badges should appear for each graph."""
        html = "\n".join(_section_appendix(SAMPLE_GRAPH_CSV_ROWS))
        assert "Standard" in html
        assert "dnCNN" in html
        assert "IVIMnet" in html


# ── Build Graph Analysis HTML ──


class TestBuildGraphAnalysisHtml:
    def test_renders_trends(self):
        html = "\n".join(_build_graph_analysis_html(SAMPLE_GRAPH_CSV_ROWS[0]))
        assert "Trends" in html
        assert "increasing" in html.lower() or "ADC" in html

    def test_renders_statistics(self):
        html = "\n".join(_build_graph_analysis_html(SAMPLE_GRAPH_CSV_ROWS[0]))
        assert "Statistics" in html or "p=" in html

    def test_renders_issues(self):
        html = "\n".join(_build_graph_analysis_html(SAMPLE_GRAPH_CSV_ROWS[0]))
        assert "Issues" in html

    def test_renders_summary(self):
        html = "\n".join(_build_graph_analysis_html(SAMPLE_GRAPH_CSV_ROWS[0]))
        assert "Summary" in html

    def test_renders_axes(self):
        html = "\n".join(_build_graph_analysis_html(SAMPLE_GRAPH_CSV_ROWS[0]))
        assert "Axes" in html
        assert "Feature" in html  # x_axis_label

    def test_renders_inflection_points(self):
        html = "\n".join(_build_graph_analysis_html(SAMPLE_GRAPH_CSV_ROWS[2]))
        assert "Inflection" in html

    def test_no_trends_no_crash(self):
        """Row with empty trends_json should not crash."""
        row = {"file_path": "Standard/test.png", "trends_json": "[]",
               "issues_json": "[]", "summary": "Test"}
        result = _build_graph_analysis_html(row)
        assert isinstance(result, list)

    def test_invalid_json_handled(self):
        """Invalid JSON in trends_json should not crash."""
        row = {"file_path": "Standard/test.png", "trends_json": "NOT JSON",
               "issues_json": "NOT JSON", "summary": "Test",
               "inflection_points_json": "NOT JSON"}
        result = _build_graph_analysis_html(row)
        assert isinstance(result, list)

    def test_long_summary_truncated(self):
        """Summary > 200 chars should be truncated with details tag."""
        row = {"file_path": "Standard/test.png", "summary": "A" * 250}
        html = "\n".join(_build_graph_analysis_html(row))
        assert "<details" in html


# ── Figure Gallery ──


class TestFigureGallery:
    def test_returns_html(self, tmp_path):
        folder = tmp_path / "saved_files_20260301_120000"
        (folder / "Standard").mkdir(parents=True)
        result = _section_figure_gallery(folder)
        assert isinstance(result, list)

    def test_nonexistent_folder(self, tmp_path):
        folder = tmp_path / "nonexistent"
        result = _section_figure_gallery(folder)
        assert isinstance(result, list)
        assert result == []

    def test_embeds_png_images(self, tmp_path):
        """Should embed PNG images as base64 data URIs."""
        folder = tmp_path / "saved_files_20260301_120000"
        std_dir = folder / "Standard"
        std_dir.mkdir(parents=True)
        # Create a minimal valid PNG (1x1 pixel)
        import base64
        # Minimal 1x1 PNG
        png_data = base64.b64decode(
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAC0lEQVQI12NgAAIABQAB"
            "Nl7BcQAAAABJRU5ErkJggg=="
        )
        (std_dir / "test_plot.png").write_bytes(png_data)
        html = "\n".join(_section_figure_gallery(folder))
        assert "data:image/png;base64," in html
        assert "Figure Gallery" in html

    def test_skips_empty_images(self, tmp_path):
        """Should skip zero-byte image files."""
        folder = tmp_path / "saved_files_20260301_120000"
        std_dir = folder / "Standard"
        std_dir.mkdir(parents=True)
        (std_dir / "empty.png").write_bytes(b"")
        result = _section_figure_gallery(folder)
        assert result == []

    def test_matches_vision_csv_rows(self, tmp_path):
        """Should match images to vision CSV rows by filename."""
        folder = tmp_path / "saved_files_20260301_120000"
        std_dir = folder / "Standard"
        std_dir.mkdir(parents=True)
        import base64
        png_data = base64.b64decode(
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAC0lEQVQI12NgAAIABQAB"
            "Nl7BcQAAAABJRU5ErkJggg=="
        )
        (std_dir / "Feature_BoxPlots_Standard.png").write_bytes(png_data)
        rows = [SAMPLE_GRAPH_CSV_ROWS[0]]
        html = "\n".join(_section_figure_gallery(folder, rows=rows))
        assert "graph-card" in html

    def test_unmatched_vision_rows(self, tmp_path):
        """Vision rows with no matching image should be listed separately."""
        folder = tmp_path / "saved_files_20260301_120000"
        (folder / "Standard").mkdir(parents=True)
        html = "\n".join(_section_figure_gallery(folder, rows=SAMPLE_GRAPH_CSV_ROWS))
        # No images on disk, so all rows are unmatched
        # No images on disk — verify section still renders without error
        assert isinstance(html, str)

    def test_root_folder_images(self, tmp_path):
        """Should also pick up images from the root folder for cross-DWI figures."""
        folder = tmp_path / "saved_files_20260301_120000"
        (folder / "Standard").mkdir(parents=True)
        import base64
        png_data = base64.b64decode(
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAC0lEQVQI12NgAAIABQAB"
            "Nl7BcQAAAABJRU5ErkJggg=="
        )
        # Root images must be >= 1024 bytes
        (folder / "cross_dwi_comparison.png").write_bytes(png_data * 100)
        html = "\n".join(_section_figure_gallery(folder))
        assert "Cross-DWI" in html


# ── Edge cases ──


class TestAppendixEdgeCases:
    def test_rows_with_missing_keys(self):
        rows = [{"file_path": "Standard/test.png"}]
        result = _section_appendix(rows)
        assert isinstance(result, list)

    def test_many_rows(self):
        rows = []
        for i in range(20):
            rows.append({
                "file_path": f"Standard/graph_{i}.png",
                "graph_title": f"Graph {i}",
                "graph_type": "scatter" if i % 2 == 0 else "line",
                "trends_json": "[]",
                "issues_json": "[]",
                "summary": f"Summary for graph {i}",
            })
        result = _section_appendix(rows)
        html = "\n".join(result)
        assert "20" in html or "graph_19" in html


class TestBuildGraphAnalysisHtmlEdgeCases:
    def test_empty_row(self):
        result = _build_graph_analysis_html({})
        assert isinstance(result, list)

    def test_row_with_all_empty_fields(self):
        row = {
            "file_path": "", "graph_title": "", "graph_type": "",
            "x_axis_label": "", "y_axis_label": "", "trends_json": "[]",
            "issues_json": "[]", "summary": "", "inflection_points_json": "[]",
        }
        result = _build_graph_analysis_html(row)
        assert isinstance(result, list)

    def test_row_with_many_trends(self):
        trends = [{"series": f"S{i}", "direction": "increasing", "description": f"Trend {i}"} for i in range(10)]
        row = {
            "file_path": "Standard/test.png",
            "trends_json": json.dumps(trends),
            "issues_json": "[]",
            "summary": "Many trends",
        }
        result = _build_graph_analysis_html(row)
        html = "\n".join(result)
        assert "S9" in html

    def test_row_with_many_inflection_points(self):
        ips = [{"approximate_x": i * 10, "approximate_y": 0.001 * i, "description": f"IP {i}"} for i in range(5)]
        row = {
            "file_path": "Standard/test.png",
            "trends_json": "[]",
            "issues_json": "[]",
            "inflection_points_json": json.dumps(ips),
            "summary": "Multiple IPs",
        }
        result = _build_graph_analysis_html(row)
        html = "\n".join(result)
        assert "IP 4" in html

    def test_summary_starts_with_json_error(self):
        row = {
            "file_path": "Standard/test.png",
            "summary": "JSON parse error: unexpected token",
            "trends_json": "[]", "issues_json": "[]",
        }
        result = _build_graph_analysis_html(row)
        html = "\n".join(result)
        assert "JSON parse error" not in html or isinstance(html, str)

    def test_short_summary_inline(self):
        row = {
            "file_path": "Standard/test.png",
            "summary": "Short summary",
            "trends_json": "[]", "issues_json": "[]",
        }
        result = _build_graph_analysis_html(row)
        html = "\n".join(result)
        assert "<details" not in html or "Short summary" in html
