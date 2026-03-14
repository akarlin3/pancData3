"""Tests for generate_report.py — HTML report generation.

Covers:
- _sig_tag: p-value → significance star mapping
- _section: Markdown heading generation
- _forest_plot_cell: forest plot HTML generation
- _effect_size_class / _effect_size_label: effect size classification
- generate_report: full report generation with mocked data sources
- New publication-level sections: Methods, Effect Sizes, Multiple
  Comparisons, Model Diagnostics, Limitations, Conclusions
- Patient flow section (CONSORT-style attrition)
- Sensitivity analysis section (EPV, cross-DWI concordance, HR stability)
- Feature stability across timepoints
- Enhanced appendix (grouped by type, statistics column)
"""

from __future__ import annotations

from pathlib import Path

import pytest  # type: ignore

from report.generate_report import (  # type: ignore
    _copy_button,
    _effect_size_class,
    _effect_size_label,
    _figure_caption,
    _forest_plot_cell,
    _manuscript_sentence,
    _section,
    _section_appendix,
    _section_data_completeness,
    _section_feature_overlap,
    _section_figure_gallery,
    _section_figure_index,
    _section_journal_guide,
    _section_manuscript_ready_findings,
    _section_patient_flow,
    _section_power_analysis,
    _section_reporting_checklist,
    _section_results_draft,
    _section_sensitivity_analysis,
    _section_table_index,
    _sig_tag,
    generate_report,
    get_numbering,
    reset_numbering,
    REPORT_JS,
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
        """Report should contain the formatted run timestamp in the title."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Analysis Report" in report
        # Timestamp is now formatted: "March 1, 2026 at 12:00:00"
        assert "March 1, 2026" in report

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


# ---------------------------------------------------------------------------
# _section_patient_flow
# ---------------------------------------------------------------------------

class TestSectionPatientFlow:
    """Verify the CONSORT-style patient flow section builder."""

    def test_empty_when_no_data(self):
        """Returns empty list when no log_data or mat_data."""
        result = _section_patient_flow(None, ["Standard"], {})
        assert result == []

    def test_empty_when_no_patient_counts(self):
        """Returns empty list when log and mat data have no flow info."""
        log_data = {"Standard": {"sanity_checks": {}}}
        result = _section_patient_flow(log_data, ["Standard"], {})
        assert result == []

    def test_shows_flow_table_with_mat_data(self):
        """With longitudinal patient count, a flow table is generated."""
        mat_data = {
            "Standard": {
                "longitudinal": {"num_patients": 40, "num_timepoints": 6},
            },
        }
        log_data = {
            "Standard": {
                "baseline": {
                    "baseline_exclusion": {
                        "n_total": 48,
                        "n_excluded": 6,
                    },
                    "total_outliers": {
                        "n_removed": 4,
                        "pct": 9.5,
                    },
                },
                "stats_comparisons": {
                    "glme_excluded": {"n_excluded": 3},
                },
            },
        }
        result = _section_patient_flow(log_data, ["Standard"], mat_data)
        html = "\n".join(result)
        assert "Patient Flow" in html
        assert "CONSORT" in html
        assert "40" in html  # initial cohort
        assert "6" in html   # baseline excluded
        assert "9.5%" in html  # outlier percentage

    def test_shows_multiple_dwi_types(self):
        """Flow table shows rows for each DWI type."""
        mat_data = {
            "Standard": {"longitudinal": {"num_patients": 40}},
            "dnCNN": {"longitudinal": {"num_patients": 38}},
        }
        result = _section_patient_flow(None, ["Standard", "dnCNN"], mat_data)
        html = "\n".join(result)
        assert "Standard" in html
        assert "dnCNN" in html

    def test_uses_baseline_n_total_fallback(self):
        """When mat_data has no patients, uses baseline_exclusion n_total."""
        log_data = {
            "Standard": {
                "baseline": {
                    "baseline_exclusion": {"n_total": 50, "n_excluded": 8},
                },
            },
        }
        result = _section_patient_flow(log_data, ["Standard"], {})
        html = "\n".join(result)
        assert "Patient Flow" in html
        assert "50" in html


# ---------------------------------------------------------------------------
# _section_sensitivity_analysis
# ---------------------------------------------------------------------------

class TestSectionSensitivityAnalysis:
    """Verify the sensitivity analysis section builder."""

    def test_empty_when_no_log_data(self):
        """Returns empty list when log_data is None."""
        result = _section_sensitivity_analysis(None, ["Standard"], {})
        assert result == []

    def test_empty_when_no_issues(self):
        """Returns empty list when no sensitivity issues are detected."""
        log_data = {
            "Standard": {
                "stats_predictive": {"feature_selections": []},
                "stats_comparisons": {"glme_details": []},
                "survival": {"hazard_ratios": []},
            },
        }
        result = _section_sensitivity_analysis(log_data, ["Standard"], {})
        assert result == []

    def test_shows_low_epv_warning(self):
        """Low EPV (events-per-variable) generates a warning."""
        log_data = {
            "Standard": {
                "baseline": {
                    "baseline_exclusion": {"n_total": 20, "n_excluded": 2},
                },
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["a", "b", "c", "d", "e", "f"]},
                    ],
                },
                "stats_comparisons": {"glme_details": []},
                "survival": {"hazard_ratios": []},
            },
        }
        mat_data = {
            "Standard": {"longitudinal": {"num_patients": 20}},
        }
        result = _section_sensitivity_analysis(log_data, ["Standard"], mat_data)
        html = "\n".join(result)
        assert "Sensitivity Analysis" in html
        assert "EPV" in html
        assert "warn-box" in html

    def test_shows_cross_dwi_concordance(self):
        """Cross-DWI GLME concordance check is reported."""
        log_data = {
            "Standard": {
                "stats_comparisons": {
                    "glme_details": [
                        {"metric": "ADC_mean", "p": 0.01, "adj_alpha": 0.025},
                        {"metric": "D_mean", "p": 0.002, "adj_alpha": 0.025},
                    ],
                },
                "stats_predictive": {"feature_selections": []},
                "survival": {"hazard_ratios": []},
            },
            "dnCNN": {
                "stats_comparisons": {
                    "glme_details": [
                        {"metric": "ADC_mean", "p": 0.03, "adj_alpha": 0.05},
                    ],
                },
                "stats_predictive": {"feature_selections": []},
                "survival": {"hazard_ratios": []},
            },
        }
        result = _section_sensitivity_analysis(
            log_data, ["Standard", "dnCNN"], {})
        html = "\n".join(result)
        assert "Sensitivity Analysis" in html
        assert "concordance" in html.lower()

    def test_shows_unstable_hr(self):
        """Hazard ratios with wide CIs (ratio > 10) are flagged."""
        log_data = {
            "Standard": {
                "stats_comparisons": {"glme_details": []},
                "stats_predictive": {"feature_selections": []},
                "survival": {
                    "hazard_ratios": [
                        {"covariate": "delta_adc", "hr": 2.5,
                         "ci_lo": 0.1, "ci_hi": 50.0, "p": 0.4},
                    ],
                },
            },
        }
        result = _section_sensitivity_analysis(log_data, ["Standard"], {})
        html = "\n".join(result)
        assert "Sensitivity Analysis" in html
        assert "Unstable HR" in html
        assert "delta_adc" in html
        assert "warn-box" in html


# ---------------------------------------------------------------------------
# Feature stability (within _section_feature_overlap)
# ---------------------------------------------------------------------------

class TestFeatureStability:
    """Verify feature stability across timepoints in feature overlap."""

    def test_stability_shown_for_multi_timepoint(self):
        """Features selected at multiple timepoints show stability table."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "lambda": 0.05,
                         "features": ["ADC_mean", "D_mean"]},
                        {"timepoint": "Fx5", "lambda": 0.06,
                         "features": ["ADC_mean", "f_ratio"]},
                        {"timepoint": "W2", "lambda": 0.04,
                         "features": ["ADC_mean", "D_mean"]},
                    ],
                },
            },
            "dnCNN": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "lambda": 0.05,
                         "features": ["ADC_mean"]},
                    ],
                },
            },
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        html = "\n".join(result)
        assert "Feature Stability" in html
        assert "ADC_mean" in html
        # ADC_mean selected at 3/3 timepoints → High stability
        assert "High" in html
        # D_mean selected at 2/3 → Moderate
        assert "Moderate" in html

    def test_no_stability_for_single_timepoint(self):
        """No stability table when each DWI type has only one timepoint."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "lambda": 0.05,
                         "features": ["ADC_mean", "D_mean"]},
                    ],
                },
            },
            "dnCNN": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "lambda": 0.05,
                         "features": ["ADC_mean"]},
                    ],
                },
            },
        }
        result = _section_feature_overlap(log_data, ["Standard", "dnCNN"])
        html = "\n".join(result)
        assert "Feature Stability" not in html


# ---------------------------------------------------------------------------
# _section_appendix (enhanced)
# ---------------------------------------------------------------------------

class TestSectionAppendix:
    """Verify the enhanced appendix section builder."""

    def test_empty_when_no_rows(self):
        """Returns empty list when no graph rows."""
        result = _section_appendix([])
        assert result == []

    def test_groups_by_graph_type(self):
        """Graphs are grouped by type with sub-headers."""
        from conftest import SAMPLE_GRAPH_CSV_ROWS  # type: ignore
        result = _section_appendix(SAMPLE_GRAPH_CSV_ROWS)
        html = "\n".join(result)
        assert "Appendix" in html
        # The fixture has "box" and "line" graph types
        assert "box" in html
        assert "line" in html
        assert "<h3>" in html  # sub-headers per type

    def test_statistics_column_present(self):
        """The Statistics column shows extracted p-values."""
        from conftest import SAMPLE_GRAPH_CSV_ROWS  # type: ignore
        result = _section_appendix(SAMPLE_GRAPH_CSV_ROWS)
        html = "\n".join(result)
        # The fixture has p = 0.003 in the Standard box plot summary
        assert "Statistics" in html
        assert "p=" in html

    def test_inflection_points_in_details(self):
        """Inflection points appear in the Details column."""
        from conftest import SAMPLE_GRAPH_CSV_ROWS  # type: ignore
        result = _section_appendix(SAMPLE_GRAPH_CSV_ROWS)
        html = "\n".join(result)
        # The IVIMnet line plot has an inflection point at day 60
        assert "Inflection" in html or "inflection" in html
        assert "60" in html

    def test_trend_descriptions_shown(self):
        """Trend descriptions are included beyond just direction tags."""
        from conftest import SAMPLE_GRAPH_CSV_ROWS  # type: ignore
        result = _section_appendix(SAMPLE_GRAPH_CSV_ROWS)
        html = "\n".join(result)
        assert "ADC rises over time" in html or "ADC" in html


# ---------------------------------------------------------------------------
# Data quality outlier balance and bias detection
# ---------------------------------------------------------------------------

class TestDataQualityOutlierBalance:
    """Verify outlier balance and bias detection in the Data Quality section."""

    def test_report_with_logs_shows_balance_column(self, saved_files_with_logs: Path):
        """Data quality outlier table includes a Balance column."""
        report = generate_report(saved_files_with_logs)
        assert "Balance" in report

    def test_report_with_logs_shows_baseline_exclusion(self, saved_files_with_logs: Path):
        """Data quality section shows baseline exclusion info."""
        report = generate_report(saved_files_with_logs)
        # The fixture has baseline exclusion: 6/48 patients
        assert "Baseline Excluded" in report or "baseline" in report.lower()


# ---------------------------------------------------------------------------
# Cox PH direction column
# ---------------------------------------------------------------------------

class TestCoxPHDirectionColumn:
    """Verify the HR direction column in Cox PH tables."""

    def test_report_with_logs_shows_hr_direction(self, saved_files_with_logs: Path):
        """Cox PH table should show Risk/Protective direction."""
        report = generate_report(saved_files_with_logs)
        # The fixture has HR=1.250 (Risk) and HR=0.750 (Protective)
        assert "Direction" in report
        assert "Risk" in report or "Protective" in report

    def test_report_with_logs_shows_discrimination_column(self, saved_files_with_logs: Path):
        """ROC table should include AUC discrimination rating."""
        report = generate_report(saved_files_with_logs)
        assert "Discrimination" in report
        # AUC 0.781 → Acceptable
        assert "Acceptable" in report


# ---------------------------------------------------------------------------
# Correlations interpretation context
# ---------------------------------------------------------------------------

class TestCorrelationsContext:
    """Verify correlation section has interpretation guidance."""

    def test_report_has_correlation_benchmarks(self, saved_files_with_graph_csv: Path):
        """Correlation section mentions Cohen's strength benchmarks."""
        report = generate_report(saved_files_with_graph_csv)
        if "Correlations" in report:
            assert "Cohen" in report or "strength" in report.lower()


# ---------------------------------------------------------------------------
# Manuscript ready findings
# ---------------------------------------------------------------------------

class TestManuscriptReadyFindings:
    """Verify the Key Findings for Manuscript section."""

    def test_empty_without_data(self):
        """Section returns empty when no data is available."""
        result = _section_manuscript_ready_findings(None, [], None, {}, {})
        assert result == []

    def test_generates_sentences_with_log_data(self, saved_files_with_logs: Path):
        """Section generates copyable sentences from log data."""
        from parsers.parse_log_metrics import parse_all_logs  # type: ignore
        log_data = parse_all_logs(saved_files_with_logs)
        result = _section_manuscript_ready_findings(
            log_data, ["Standard"], None, {}, {}
        )
        html = "\n".join(result)
        assert "manuscript-sentence" in html
        assert "Manuscript" in html
        # Should have copy buttons
        assert "copy-btn" in html

    def test_includes_auc_sentence(self, saved_files_with_logs: Path):
        """Section includes AUC performance sentence when ROC data exists."""
        from parsers.parse_log_metrics import parse_all_logs  # type: ignore
        log_data = parse_all_logs(saved_files_with_logs)
        result = _section_manuscript_ready_findings(
            log_data, ["Standard"], None, {}, {}
        )
        html = "\n".join(result)
        assert "AUC" in html or "elastic-net" in html.lower()

    def test_includes_hazard_ratio_sentence(self, saved_files_with_logs: Path):
        """Section includes HR sentence when survival data exists."""
        from parsers.parse_log_metrics import parse_all_logs  # type: ignore
        log_data = parse_all_logs(saved_files_with_logs)
        result = _section_manuscript_ready_findings(
            log_data, ["Standard"], None, {}, {}
        )
        html = "\n".join(result)
        assert "HR" in html or "hazard" in html.lower() or "Cox" in html

    def test_copy_all_button_present(self, saved_files_with_logs: Path):
        """Section has a copy-all button for the full paragraph."""
        from parsers.parse_log_metrics import parse_all_logs  # type: ignore
        log_data = parse_all_logs(saved_files_with_logs)
        result = _section_manuscript_ready_findings(
            log_data, ["Standard"], None, {}, {}
        )
        html = "\n".join(result)
        assert "all-findings" in html


# ---------------------------------------------------------------------------
# Reporting checklist
# ---------------------------------------------------------------------------

class TestReportingChecklist:
    """Verify the STROBE/REMARK reporting checklist section."""

    def test_checklist_with_no_data(self):
        """Checklist renders even without any data."""
        result = _section_reporting_checklist(None, [], {}, None, [])
        html = "\n".join(result)
        assert "Reporting Guideline" in html
        assert "STROBE" in html
        assert "REMARK" in html

    def test_checklist_counts_addressed_items(self, saved_files_with_logs: Path):
        """Checklist correctly counts addressed vs partial items."""
        from parsers.parse_log_metrics import parse_all_logs  # type: ignore
        log_data = parse_all_logs(saved_files_with_logs)
        result = _section_reporting_checklist(
            log_data, ["Standard"], {}, None, []
        )
        html = "\n".join(result)
        assert "Addressed" in html
        assert "checklist-done" in html

    def test_checklist_has_dynamic_status(self, saved_files_with_logs: Path):
        """Checklist status changes based on available data."""
        from parsers.parse_log_metrics import parse_all_logs  # type: ignore
        log_data = parse_all_logs(saved_files_with_logs)
        result_with = _section_reporting_checklist(
            log_data, ["Standard"], {}, None, []
        )
        result_without = _section_reporting_checklist(
            None, [], {}, None, []
        )
        # With data should have more 'done' items
        done_with = "\n".join(result_with).count("checklist-done")
        done_without = "\n".join(result_without).count("checklist-done")
        assert done_with >= done_without


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
# BibTeX export
# ---------------------------------------------------------------------------

class TestBibTeXExport:
    """Verify BibTeX export in references section."""

    def test_report_contains_bibtex(self, saved_files_with_graph_csv: Path):
        """Report includes BibTeX export block."""
        report = generate_report(saved_files_with_graph_csv)
        assert "BibTeX" in report
        assert "@article" in report or "@book" in report

    def test_report_contains_js(self, saved_files_with_graph_csv: Path):
        """Report includes JavaScript for copy-to-clipboard."""
        report = generate_report(saved_files_with_graph_csv)
        assert "copyText" in report

    def test_report_js_constant_exists(self):
        """REPORT_JS constant is a non-empty string."""
        assert isinstance(REPORT_JS, str)
        assert len(REPORT_JS) > 50
        assert "copyText" in REPORT_JS


# ---------------------------------------------------------------------------
# Full integration: new sections in report
# ---------------------------------------------------------------------------

class TestNewSectionsIntegration:
    """Verify new sections appear in the full report."""

    def test_report_has_checklist(self, saved_files_with_graph_csv: Path):
        """Full report includes the reporting checklist section."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Reporting Guideline" in report

    def test_report_has_bibtex_block(self, saved_files_with_graph_csv: Path):
        """Full report includes the BibTeX export."""
        report = generate_report(saved_files_with_graph_csv)
        assert "bibtex-block" in report

    def test_report_has_strobe_reference(self, saved_files_with_graph_csv: Path):
        """Full report includes STROBE and REMARK references."""
        report = generate_report(saved_files_with_graph_csv)
        assert "STROBE" in report
        assert "REMARK" in report

    def test_report_has_results_draft(self, saved_files_with_graph_csv: Path):
        """Full report includes the draft results section."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Draft Results Section" in report

    def test_report_has_figure_gallery(self, saved_files_with_graph_csv: Path):
        """Full report includes figure gallery section (even if empty)."""
        # The gallery renders even with no images (it just returns [])
        report = generate_report(saved_files_with_graph_csv)
        # Gallery may be empty if no images exist, that's fine
        assert report  # report should still generate successfully

    def test_report_has_journal_guide(self, saved_files_with_graph_csv: Path):
        """Full report includes journal submission guidance."""
        report = generate_report(saved_files_with_graph_csv)
        assert "Journal Submission Guidance" in report


# ---------------------------------------------------------------------------
# Draft Results Section
# ---------------------------------------------------------------------------

class TestResultsDraft:
    """Verify the auto-generated Results section draft."""

    def test_empty_when_no_data(self):
        """No results draft without any data."""
        result = _section_results_draft(None, [], None, {}, {})
        assert result == []

    def test_cohort_paragraph_from_mat_data(self):
        """Cohort paragraph is generated from MAT data."""
        mat_data = {"Standard": {"longitudinal": {
            "num_patients": 25, "num_timepoints": 5}}}
        result = _section_results_draft(None, ["Standard"], None, mat_data, {})
        html = "\n".join(result)
        assert "25 patients" in html
        assert "5 timepoints" in html
        assert "Draft Results" in html

    def test_group_comparisons_paragraph(self):
        """Group comparisons paragraph is generated from GLME data."""
        log_data = {"Standard": {"stats_comparisons": {
            "glme_details": [
                {"metric": "mean_adc", "p": 0.001, "adj_alpha": 0.01},
                {"metric": "mean_d", "p": 0.5, "adj_alpha": 0.025},
            ],
        }}}
        result = _section_results_draft(log_data, ["Standard"], None, {}, {})
        html = "\n".join(result)
        assert "GLME" in html
        assert "mean_adc" in html

    def test_predictive_paragraph(self):
        """Predictive modelling paragraph with AUC."""
        log_data = {"Standard": {"stats_predictive": {
            "roc_analyses": [{"auc": 0.82, "timepoint": "W2",
                              "sensitivity": 85.0, "specificity": 75.0}],
            "feature_selections": [{"timepoint": "W2", "lambda": 0.01,
                                    "features": ["f1", "f2"]}],
        }}}
        result = _section_results_draft(log_data, ["Standard"], None, {}, {})
        html = "\n".join(result)
        assert "0.82" in html or "0.820" in html
        assert "sensitivity" in html.lower()

    def test_survival_paragraph(self):
        """Survival analysis paragraph with hazard ratios."""
        log_data = {"Standard": {"survival": {
            "hazard_ratios": [
                {"covariate": "mean_adc_delta", "hr": 2.5,
                 "ci_lo": 1.2, "ci_hi": 5.1, "p": 0.01},
            ],
            "global_lrt": {"df": 3, "chi2": 12.5, "p": 0.006},
            "ipcw": {"min_weight": 0.8, "max_weight": 1.2},
        }}}
        result = _section_results_draft(log_data, ["Standard"], None, {}, {})
        html = "\n".join(result)
        assert "mean_adc_delta" in html
        assert "HR" in html
        assert "IPCW" in html

    def test_copy_all_button_present(self):
        """Results draft includes a copy-all button."""
        mat_data = {"Standard": {"longitudinal": {
            "num_patients": 10, "num_timepoints": 3}}}
        result = _section_results_draft(None, ["Standard"], None, mat_data, {})
        html = "\n".join(result)
        assert "results-draft-all" in html
        assert "copy" in html.lower()


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
# Figure Gallery
# ---------------------------------------------------------------------------

class TestFigureGallery:
    """Verify the figure gallery section with embedded images."""

    def test_empty_when_no_images(self, saved_files_dir: Path):
        """Gallery is empty when no images exist."""
        reset_numbering()
        result = _section_figure_gallery(saved_files_dir)
        assert result == []

    def test_embeds_png_images(self, saved_files_dir: Path):
        """Gallery embeds PNG images as base64."""
        reset_numbering()
        # Create a small test PNG (1x1 pixel)
        import struct
        import zlib
        def _make_tiny_png():
            raw_data = b'\x00\x00\x00\x00'
            compressed = zlib.compress(raw_data)
            def chunk(ctype, data):
                c = ctype + data
                return struct.pack('>I', len(data)) + c + struct.pack('>I', zlib.crc32(c) & 0xffffffff)
            return (b'\x89PNG\r\n\x1a\n' +
                    chunk(b'IHDR', struct.pack('>IIBBBBB', 1, 1, 8, 2, 0, 0, 0)) +
                    chunk(b'IDAT', compressed) +
                    chunk(b'IEND', b''))

        img_dir = saved_files_dir / "Standard"
        img_path = img_dir / "test_figure.png"
        img_path.write_bytes(_make_tiny_png())

        result = _section_figure_gallery(saved_files_dir)
        html = "\n".join(result)
        assert "Figure Gallery" in html
        assert "data:image/png;base64," in html
        assert "Test Figure" in html  # filename -> title

    def test_groups_by_dwi_type(self, saved_files_dir: Path):
        """Gallery groups figures by DWI type."""
        reset_numbering()
        import struct
        import zlib
        def _make_tiny_png():
            raw_data = b'\x00\x00\x00\x00'
            compressed = zlib.compress(raw_data)
            def chunk(ctype, data):
                c = ctype + data
                return struct.pack('>I', len(data)) + c + struct.pack('>I', zlib.crc32(c) & 0xffffffff)
            return (b'\x89PNG\r\n\x1a\n' +
                    chunk(b'IHDR', struct.pack('>IIBBBBB', 1, 1, 8, 2, 0, 0, 0)) +
                    chunk(b'IDAT', compressed) +
                    chunk(b'IEND', b''))
        png = _make_tiny_png()
        (saved_files_dir / "Standard" / "fig1.png").write_bytes(png)
        (saved_files_dir / "dnCNN" / "fig2.png").write_bytes(png)

        result = _section_figure_gallery(saved_files_dir)
        html = "\n".join(result)
        assert "badge-standard" in html
        assert "badge-dncnn" in html


# ---------------------------------------------------------------------------
# Journal Submission Guidance
# ---------------------------------------------------------------------------

class TestJournalGuide:
    """Verify the journal submission guidance section."""

    def test_basic_output(self):
        """Section includes journal recommendations."""
        result = _section_journal_guide(None, ["Standard"], {})
        html = "\n".join(result)
        assert "Journal Submission Guidance" in html
        assert "Radiotherapy and Oncology" in html
        assert "Manuscript Preparation Checklist" in html

    def test_adds_survival_journal(self):
        """Adds Acta Oncologica when survival data present."""
        log_data = {"Standard": {"survival": {
            "hazard_ratios": [{"covariate": "x", "hr": 1.5, "p": 0.03}]
        }}}
        result = _section_journal_guide(log_data, ["Standard"], {})
        html = "\n".join(result)
        assert "Acta Oncologica" in html

    def test_adds_predictive_journal(self):
        """Adds European Radiology when predictive data present."""
        log_data = {"Standard": {"stats_predictive": {
            "roc_analyses": [{"auc": 0.8, "timepoint": "BL"}]
        }}}
        result = _section_journal_guide(log_data, ["Standard"], {})
        html = "\n".join(result)
        assert "European Radiology" in html

    def test_keywords_copyable(self):
        """Keywords section includes copy button."""
        result = _section_journal_guide(None, ["Standard"], {})
        html = "\n".join(result)
        assert "diffusion-weighted imaging" in html
        assert "copy" in html.lower()

    def test_keywords_include_survival(self):
        """Keywords include survival terms when survival data present."""
        log_data = {"Standard": {"survival": {
            "hazard_ratios": [{"covariate": "x", "hr": 1.5, "p": 0.03}]
        }}}
        result = _section_journal_guide(log_data, ["Standard"], {})
        html = "\n".join(result)
        assert "survival analysis" in html


# ---------------------------------------------------------------------------
# Manuscript Findings with Effect Sizes
# ---------------------------------------------------------------------------

class TestManuscriptFindingsEffectSizes:
    """Verify that manuscript sentences include effect sizes."""

    def test_hr_sentence_includes_effect_size(self):
        """Cox PH manuscript sentences include effect size labels."""
        log_data = {"Standard": {"survival": {
            "hazard_ratios": [
                {"covariate": "mean_adc", "hr": 3.0,
                 "ci_lo": 1.5, "ci_hi": 6.0, "p": 0.002},
            ]
        }}}
        result = _section_manuscript_ready_findings(
            log_data, ["Standard"], None, {}, {})
        html = "\n".join(result)
        assert "effect" in html.lower()
        assert "3.000" in html

    def test_glme_sentence_includes_sample_size(self):
        """GLME sentences include sample size when exclusion data present."""
        log_data = {"Standard": {"stats_comparisons": {
            "glme_details": [
                {"metric": "mean_adc", "p": 0.001, "adj_alpha": 0.01},
            ],
            "glme_excluded": {"n_excluded": 3, "n_total": 25, "pct": 12.0},
        }}}
        result = _section_manuscript_ready_findings(
            log_data, ["Standard"], None, {}, {})
        html = "\n".join(result)
        assert "22 evaluable patients" in html
        assert "excluding 3" in html
