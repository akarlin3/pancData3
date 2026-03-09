"""Tests for generate_report.py — HTML report generation.

Covers:
- _sig_tag: p-value → significance star mapping
- _section: Markdown heading generation
- _forest_plot_cell: forest plot HTML generation
- _effect_size_class / _effect_size_label: effect size classification
- generate_report: full report generation with mocked data sources
- New publication-level sections: Methods, Effect Sizes, Multiple
  Comparisons, Model Diagnostics, Limitations, Conclusions
"""

from __future__ import annotations

from pathlib import Path

import pytest

from generate_report import (
    _effect_size_class,
    _effect_size_label,
    _forest_plot_cell,
    _section,
    _section_data_completeness,
    _section_feature_overlap,
    _section_power_analysis,
    _sig_tag,
    generate_report,
)


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
        assert _sig_tag(0.001) == "**"   # exactly 0.001 → **
        assert _sig_tag(0.01) == "*"     # exactly 0.01 → *
        assert _sig_tag(0.05) == ""      # exactly 0.05 → not significant


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
# generate_report (integration)
# ---------------------------------------------------------------------------

class TestGenerateReport:
    """Integration tests for the full Markdown report generation."""

    def test_report_contains_header(self, saved_files_with_graph_csv: Path):
        """Report should contain the folder timestamp in the title."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Analysis Report" in report
        assert "20260301_120000" in report

    def test_report_contains_dwi_types(self, saved_files_with_graph_csv: Path):
        """All three DWI types should be mentioned somewhere in the report."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Standard" in report
        assert "dnCNN" in report
        assert "IVIMnet" in report

    def test_report_contains_graph_overview(self, saved_files_with_graph_csv: Path):
        """The graph type table should appear."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Graph Type" in report
        assert "box" in report
        assert "line" in report

    def test_report_contains_graph_issues(self, saved_files_with_graph_csv: Path):
        """Graph Issues section should show detected issues from the CSV."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Graph Issues" in report
        assert "Y-axis label partially obscured by legend" in report
        assert "Overlapping legend text" in report

    def test_report_extracts_significant_pvalues(self, saved_files_with_graph_csv: Path):
        """P-values < 0.05 from the graph CSV should appear in the report."""
        report = generate_report(saved_files_with_graph_csv)
        # The fixture has p = 0.003 and p-value = 0.02
        assert "0.003" in report or "0.0030" in report

    def test_report_extracts_correlations(self, saved_files_with_graph_csv: Path):
        """Notable correlations should appear in the report."""
        report = generate_report(saved_files_with_graph_csv)
        # The fixture has r = 0.65
        assert "0.65" in report

    def test_report_appendix_lists_all_graphs(self, saved_files_with_graph_csv: Path):
        report = generate_report(saved_files_with_graph_csv)
        assert "Appendix" in report
        assert "Feature_BoxPlots" in report

    def test_report_footer(self, saved_files_with_graph_csv: Path):
        report = generate_report(saved_files_with_graph_csv)
        assert "Report generated by pancData3" in report

    def test_empty_folder_still_generates(self, saved_files_dir: Path):
        """A folder with no graph CSV and no logs should still produce a report."""
        report = generate_report(saved_files_dir)
        assert "Analysis Report" in report
        assert "No significant findings" in report or "No predictive" in report

    def test_report_with_logs(self, saved_files_with_logs: Path):
        """With log data present, predictive performance should appear."""
        report = generate_report(saved_files_with_logs)
        # The fixture has AUC = 0.781
        assert "0.781" in report or "Predictive" in report

    def test_report_cross_dwi_comparison(self, saved_files_with_graph_csv: Path):
        """Cross-DWI comparison section for Feature_BoxPlots (Standard + dnCNN)."""
        report = generate_report(saved_files_with_graph_csv)
        # Feature_BoxPlots has Standard (increasing) and dnCNN (decreasing) → DIFFER
        if "Cross-DWI" in report:
            assert "DIFFER" in report or "AGREE" in report

    def test_report_is_valid_html(self, saved_files_with_graph_csv: Path):
        """Basic structural check: headings and tables."""
        report = generate_report(saved_files_with_graph_csv)
        assert "<h2" in report  # at least a few headings
        assert "<table" in report  # at least one table

    def test_report_contains_methods_section(self, saved_files_with_graph_csv: Path):
        """The Methods section should describe statistical methodology."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Statistical Methods" in report
        assert "Wilcoxon" in report
        assert "Benjamini" in report
        assert "Cox" in report

    def test_report_contains_structured_abstract(self, saved_files_with_graph_csv: Path):
        """Executive summary should have structured abstract subsections."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Objective" in report
        assert "Key Results" in report
        assert "Conclusions" in report

    def test_report_contains_limitations(self, saved_files_with_graph_csv: Path):
        """Limitations section should appear in the report."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Limitations" in report
        assert "Single-institution" in report

    def test_report_contains_model_diagnostics(self, saved_files_with_graph_csv: Path):
        """Model diagnostics section should have assumptions listed."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Model Diagnostics" in report
        assert "Proportional hazards" in report

    def test_report_contains_conclusions_section(self, saved_files_with_graph_csv: Path):
        """Conclusions section should appear."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Conclusions" in report
        assert "Future directions" in report

    def test_report_with_logs_has_effect_sizes(self, saved_files_with_logs: Path):
        """With log data, the Effect Size section should report HR effects."""
        report = generate_report(saved_files_with_logs)
        assert "Effect Size" in report

    def test_report_with_logs_has_multiple_comparisons(self, saved_files_with_logs: Path):
        """With log data, the Multiple Comparisons section should appear."""
        report = generate_report(saved_files_with_logs)
        assert "Multiple Comparison" in report
        assert "BH-FDR" in report or "FDR" in report

    def test_report_with_logs_has_diagnostics_detail(self, saved_files_with_logs: Path):
        """With log data, model diagnostics should show data-driven issues."""
        report = generate_report(saved_files_with_logs)
        assert "IPCW" in report or "weight" in report or "Assumptions" in report

    def test_report_with_logs_has_discrimination_table(self, saved_files_with_logs: Path):
        """AUC discrimination interpretation should appear."""
        report = generate_report(saved_files_with_logs)
        # The fixture has AUC = 0.781, which maps to "Acceptable"
        assert "Acceptable" in report or "Discrimination" in report or "discrimination" in report


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
# _section_data_completeness
# ---------------------------------------------------------------------------

class TestSectionDataCompleteness:
    """Verify the data completeness section builder."""

    def test_empty_when_no_log_data(self):
        """Returns empty list when log_data is None."""
        result = _section_data_completeness(None, ["Standard"])
        assert result == []

    def test_empty_when_no_sanity_data(self):
        """Returns empty list when sanity_checks has no flags."""
        log_data = {
            "Standard": {
                "sanity_checks": {
                    "total_convergence": 0,
                    "all_converged": False,
                    "dim_mismatches": 0,
                    "excessive_nan": [],
                },
            },
        }
        result = _section_data_completeness(log_data, ["Standard"])
        assert result == []

    def test_shows_convergence_passed(self):
        """All-converged status generates the section with a 'Passed' card."""
        log_data = {
            "Standard": {
                "sanity_checks": {
                    "all_converged": True,
                    "total_convergence": 0,
                    "convergence_flags": [],
                    "outliers": [],
                    "total_outliers": 0,
                    "dim_mismatches": 0,
                    "nan_dose_warnings": 0,
                    "excessive_nan": [],
                },
            },
        }
        result = _section_data_completeness(log_data, ["Standard"])
        html = "\n".join(result)
        assert "Data Completeness" in html
        assert "Passed" in html

    def test_shows_convergence_flags(self):
        """Convergence flags generate a warning box."""
        log_data = {
            "Standard": {
                "sanity_checks": {
                    "all_converged": False,
                    "total_convergence": 5,
                    "convergence_flags": [],
                    "outliers": [],
                    "total_outliers": 0,
                    "dim_mismatches": 0,
                    "nan_dose_warnings": 0,
                    "excessive_nan": [],
                },
            },
        }
        result = _section_data_completeness(log_data, ["Standard"])
        html = "\n".join(result)
        assert "5 convergence" in html
        assert "warn-box" in html

    def test_shows_excessive_nan(self):
        """Excessive NaN parameters generate specific warnings."""
        log_data = {
            "Standard": {
                "sanity_checks": {
                    "all_converged": False,
                    "total_convergence": 0,
                    "convergence_flags": [],
                    "outliers": [],
                    "total_outliers": 0,
                    "dim_mismatches": 0,
                    "nan_dose_warnings": 0,
                    "excessive_nan": [
                        {"parameter": "D", "pct_nan": 65.0},
                    ],
                },
            },
        }
        result = _section_data_completeness(log_data, ["Standard"])
        html = "\n".join(result)
        assert "65.0%" in html
        assert "Parameter D" in html


# ---------------------------------------------------------------------------
# _section_feature_overlap
# ---------------------------------------------------------------------------

class TestSectionFeatureOverlap:
    """Verify the cross-DWI feature overlap section builder."""

    def test_empty_when_single_dwi_type(self):
        """Returns empty list with only one DWI type."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "Fx5", "lambda": 0.05, "features": ["ADC", "D"]},
                    ],
                },
            },
        }
        result = _section_feature_overlap(log_data, ["Standard"])
        assert result == []

    def test_empty_when_no_log_data(self):
        """Returns empty list when log_data is None."""
        result = _section_feature_overlap(None, ["Standard", "dnCNN"])
        assert result == []

    def test_shows_shared_features(self):
        """Features selected by multiple DWI types are flagged as 'Shared'."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "Fx5", "lambda": 0.05,
                         "features": ["ADC_mean", "D_mean", "f_mean"]},
                    ],
                },
            },
            "dnCNN": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "Fx5", "lambda": 0.06,
                         "features": ["ADC_mean", "D_mean", "Dstar_vol"]},
                    ],
                },
            },
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        html = "\n".join(result)
        assert "Feature Overlap" in html
        assert "Shared" in html
        assert "ADC_mean" in html
        assert "D_mean" in html
        # f_mean and Dstar_vol are type-specific
        assert "Type-specific" in html

    def test_summary_percentage(self):
        """Summary shows correct shared percentage."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "lambda": 0.05,
                         "features": ["A", "B"]},
                    ],
                },
            },
            "dnCNN": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "lambda": 0.06,
                         "features": ["A", "C"]},
                    ],
                },
            },
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        html = "\n".join(result)
        # A is shared, B and C are unique -> 1/3 = 33%
        assert "1/3" in html
        assert "33%" in html


# ---------------------------------------------------------------------------
# _section_power_analysis
# ---------------------------------------------------------------------------

class TestSectionPowerAnalysis:
    """Verify the statistical power commentary section builder."""

    def test_empty_when_no_patients(self):
        """Returns empty list when no cohort data is available."""
        result = _section_power_analysis(None, ["Standard"], {})
        assert result == []

    def test_shows_power_with_cohort(self):
        """With cohort data, approximate detectable effect sizes are shown."""
        mat_data = {
            "Standard": {
                "longitudinal": {
                    "num_patients": 25,
                    "num_timepoints": 6,
                },
            },
        }
        result = _section_power_analysis(None, ["Standard"], mat_data)
        html = "\n".join(result)
        assert "Statistical Power" in html
        assert "n = 25" in html
        assert "Wilcoxon" in html
        assert "Cox PH" in html

    def test_shows_fdr_penalty(self):
        """When GLME tests exist, FDR correction penalty is shown."""
        log_data = {
            "Standard": {
                "stats_comparisons": {
                    "glme_details": [
                        {"metric": "ADC", "p": 0.01, "adj_alpha": 0.025},
                        {"metric": "D", "p": 0.03, "adj_alpha": 0.05},
                    ],
                },
                "survival": {"hazard_ratios": []},
                "stats_predictive": {"roc_analyses": [], "feature_selections": []},
            },
        }
        mat_data = {
            "Standard": {
                "longitudinal": {"num_patients": 30, "num_timepoints": 6},
            },
        }
        result = _section_power_analysis(log_data, ["Standard"], mat_data)
        html = "\n".join(result)
        assert "FDR correction" in html
        assert "2 tests" in html
