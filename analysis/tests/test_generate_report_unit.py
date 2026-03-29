"""Unit tests for report/generate_report.py — report generation and re-exports.

Tests:
- Re-exported formatter names are importable from generate_report
- Re-exported section builder names are importable from generate_report
- generate_report() handles empty folders gracefully (no crash)
- generate_report() produces HTML output with minimal synthetic data
- Re-exports match their canonical source modules
- markdown_to_html() produces a valid HTML document
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

# Ensure analysis/ is importable.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


# ---------------------------------------------------------------------------
# Import availability
# ---------------------------------------------------------------------------

class TestImportsAvailable:
    """Verify all re-exported names are importable from generate_report."""

    def test_formatter_re_exports(self):
        """All formatter helpers should be importable from generate_report."""
        from report.generate_report import (
            CSS,
            HTML_TEMPLATE,
            NAV_SECTIONS,
            REFERENCES,
            REPORT_JS,
            _cite,
            _copy_button,
            _dwi_badge,
            _effect_size_class,
            _effect_size_label,
            _esc,
            _figure_caption,
            _forest_plot_cell,
            _get_consensus,
            _h2,
            _manuscript_sentence,
            _nav_bar,
            _references_section,
            _section,
            _sig_class,
            _sig_tag,
            _stat_card,
            _table_caption,
            _trend_tag,
            get_numbering,
            reset_numbering,
        )
        # Smoke check: the imported names are callable or non-None.
        assert callable(_esc)
        assert callable(_sig_tag)
        assert callable(_h2)
        assert isinstance(CSS, str)
        assert isinstance(REPORT_JS, str)

    def test_section_builders_importable(self):
        """All section builder functions should be importable from generate_report."""
        from report.generate_report import (
            _part_break,
            _section_cover_page,
            _section_print_toc,
            _section_appendix,
            _section_broad_statistical_overview,
            _section_cohort_overview,
            _section_conclusions,
            _section_correlations,
            _section_cross_dwi_comparison,
            _section_cross_pipeline_dice,
            _section_failure_rates,
            _section_pruning_results,
            _section_core_method_outcomes,
            _section_data_availability,
            _section_data_completeness,
            _section_effect_sizes,
            _section_executive_summary,
            _section_feature_overlap,
            _section_figure_gallery,
            _section_figure_index,
            _section_forest_plot_figure,
            _section_graph_issues,
            _section_graph_overview,
            _section_hypothesis,
            _section_journal_guide,
            _section_limitations,
            _section_manuscript_ready_findings,
            _section_mat_data,
            _section_methods,
            _section_model_diagnostics,
            _section_model_robustness,
            _section_multiple_comparisons,
            _section_patient_flow,
            _section_power_analysis,
            _section_predictive_performance,
            _section_publication_header,
            _section_reporting_checklist,
            _section_results_draft,
            _section_sensitivity_analysis,
            _section_stats_by_graph_type,
            _section_statistical_significance,
            _section_table_index,
            _section_treatment_response,
            build_dca_section,
            build_nri_idi_section,
            build_texture_section,
            build_registration_quality_section,
        )
        # All should be callable.
        assert callable(_section_executive_summary)
        assert callable(_section_appendix)
        assert callable(build_dca_section)
        assert callable(build_registration_quality_section)


# ---------------------------------------------------------------------------
# Re-export identity checks
# ---------------------------------------------------------------------------

class TestReExportsMatchCanonical:
    """Verify re-exported names are the same objects as their canonical sources."""

    def test_esc_matches_formatters(self):
        """_esc from generate_report should be the same as from report_formatters."""
        from report.generate_report import _esc as gen_esc
        from report.report_formatters import _esc as fmt_esc
        assert gen_esc is fmt_esc

    def test_sig_tag_matches_formatters(self):
        """_sig_tag from generate_report should be the same as from report_formatters."""
        from report.generate_report import _sig_tag as gen_tag
        from report.report_formatters import _sig_tag as fmt_tag
        assert gen_tag is fmt_tag

    def test_h2_matches_formatters(self):
        """_h2 from generate_report should be the same as from report_formatters."""
        from report.generate_report import _h2 as gen_h2
        from report.report_formatters import _h2 as fmt_h2
        assert gen_h2 is fmt_h2

    def test_section_builder_matches_sections_package(self):
        """_section_executive_summary should be the same object from both paths."""
        from report.generate_report import _section_executive_summary as gen_fn
        from report.sections import _section_executive_summary as sec_fn
        assert gen_fn is sec_fn


# ---------------------------------------------------------------------------
# generate_report with empty folder
# ---------------------------------------------------------------------------

class TestGenerateReportEmptyFolder:
    """Verify generate_report handles a folder with no data files gracefully."""

    def test_generate_report_with_empty_folder(self, tmp_path):
        """Calling generate_report on an empty temp dir should not crash."""
        from report.generate_report import generate_report

        # Name the folder so the timestamp parser has something to work with.
        folder = tmp_path / "saved_files_20260101_120000"
        folder.mkdir()

        html = generate_report(folder)
        assert isinstance(html, str)
        assert len(html) > 0
        assert "<html" in html
        assert "</html>" in html
        # Should contain the report title with the formatted timestamp.
        assert "Analysis Report" in html


# ---------------------------------------------------------------------------
# generate_report with minimal data
# ---------------------------------------------------------------------------

class TestGenerateReportMinimalData:
    """Verify generate_report produces HTML output with minimal synthetic data."""

    def test_generate_report_with_minimal_data(self, tmp_path):
        """Create minimal CSV/log files and verify HTML output is created."""
        from report.generate_report import generate_report

        folder = tmp_path / "saved_files_20260315_100000"
        folder.mkdir()

        # Create a minimal Standard DWI type subdirectory with a log file.
        std_dir = folder / "Standard"
        std_dir.mkdir()
        log_file = std_dir / "pipeline_log_Standard.txt"
        log_file.write_text(
            "Pipeline started\n"
            "Wilcoxon rank-sum: p = 0.023 (adc_mean)\n"
            "Pipeline complete\n"
        )

        # Create a minimal graph_analysis_results.csv.
        csv_path = folder / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow([
                "file", "dwi_type", "graph_name", "graph_type",
                "x_label", "y_label", "x_type", "y_type",
                "trend_direction", "trend_description",
                "inflection_x", "inflection_description",
                "stat_test_name", "stat_test_value", "stat_test_p",
                "outlier_label", "outlier_value",
                "ref_line_label", "ref_line_value",
                "clinical_relevance",
            ])
            writer.writerow([
                "Feature_BoxPlots_Standard.png", "Standard",
                "Feature_BoxPlots", "box",
                "Group", "ADC", "categorical", "continuous",
                "increasing", "Higher ADC in treatment group",
                "", "",
                "Wilcoxon", "2.31", "0.023",
                "", "",
                "", "",
                "Moderate",
            ])

        html = generate_report(folder)
        assert isinstance(html, str)
        assert "<html" in html
        assert "Analysis Report" in html
        # The DWI type should appear somewhere in the report.
        assert "Standard" in html

    def test_html_output_written_to_disk(self, tmp_path):
        """Verify that generate_report output can be written to disk."""
        from report.generate_report import generate_report

        folder = tmp_path / "saved_files_20260320_080000"
        folder.mkdir()

        html = generate_report(folder)
        out_path = folder / "analysis_report.html"
        out_path.write_text(html, encoding="utf-8")
        assert out_path.exists()
        assert out_path.stat().st_size > 0


# ---------------------------------------------------------------------------
# markdown_to_html
# ---------------------------------------------------------------------------

class TestMarkdownToHtml:
    """Verify Markdown-to-HTML conversion helper."""

    def test_basic_conversion(self):
        """markdown_to_html should produce a complete HTML document."""
        from report.generate_report import markdown_to_html

        md = "# Hello\n\nThis is a **test**."
        html = markdown_to_html(md, "Test Title")
        assert "<html" in html
        assert "Test Title" in html
        assert "<h1>" in html or "<h1" in html
        assert "<strong>test</strong>" in html

    def test_table_extension(self):
        """Markdown tables should be converted to HTML tables."""
        from report.generate_report import markdown_to_html

        md = "| A | B |\n|---|---|\n| 1 | 2 |"
        html = markdown_to_html(md, "Table Test")
        assert "<table>" in html or "<table" in html


# ---------------------------------------------------------------------------
# _wrap_str_builder
# ---------------------------------------------------------------------------

class TestWrapStrBuilder:
    """Verify the internal wrapper that adapts str-returning builders to list[str]."""

    def test_wraps_non_empty_string(self):
        """A function returning a non-empty string should produce a one-element list."""
        from report.generate_report import _wrap_str_builder

        def builder():
            return "<div>hello</div>"

        wrapped = _wrap_str_builder(builder)
        result = wrapped()
        assert result == ["<div>hello</div>"]

    def test_wraps_empty_string(self):
        """A function returning an empty string should produce an empty list."""
        from report.generate_report import _wrap_str_builder

        def builder():
            return ""

        wrapped = _wrap_str_builder(builder)
        result = wrapped()
        assert result == []

    def test_wraps_none(self):
        """A function returning None should produce an empty list."""
        from report.generate_report import _wrap_str_builder

        def builder():
            return None

        wrapped = _wrap_str_builder(builder)
        result = wrapped()
        assert result == []
