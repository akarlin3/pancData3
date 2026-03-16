"""Tests for individual report section builders.

Covers:
- _forest_plot_cell: forest plot HTML generation
- _section_data_completeness: data quality / sanity check section
- _section_feature_overlap: cross-DWI feature overlap section
- Feature stability across timepoints (within feature overlap)
- _section_power_analysis: statistical power commentary
- _section_patient_flow: CONSORT-style patient attrition
- _section_sensitivity_analysis: EPV, concordance, HR stability
- _section_appendix: enhanced appendix with grouping
- Figure index (List of Figures)
- Figure gallery with embedded images
"""

from __future__ import annotations

from pathlib import Path

import pytest  # type: ignore

from conftest import SAMPLE_GRAPH_CSV_ROWS, make_tiny_png  # type: ignore

from report.generate_report import (  # type: ignore
    _forest_plot_cell,
    _section_appendix,
    _section_data_completeness,
    _section_feature_overlap,
    _section_figure_gallery,
    _section_figure_index,
    _section_patient_flow,
    _section_power_analysis,
    _section_sensitivity_analysis,
    _figure_caption,
    reset_numbering,
)


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
        # ADC_mean selected at 3/3 timepoints -> High stability
        assert "High" in html
        # D_mean selected at 2/3 -> Moderate
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
        result = _section_appendix(SAMPLE_GRAPH_CSV_ROWS)
        html = "\n".join(result)
        assert "Appendix" in html
        # The fixture has "box" and "line" graph types
        assert "box" in html
        assert "line" in html
        assert "<h3>" in html  # sub-headers per type

    def test_statistics_column_present(self):
        """The Statistics column shows extracted p-values."""
        result = _section_appendix(SAMPLE_GRAPH_CSV_ROWS)
        html = "\n".join(result)
        # The fixture has p = 0.003 in the Standard box plot summary
        assert "Statistics" in html
        assert "p=" in html

    def test_inflection_points_in_details(self):
        """Inflection points appear in the Details column."""
        result = _section_appendix(SAMPLE_GRAPH_CSV_ROWS)
        html = "\n".join(result)
        # The IVIMnet line plot has an inflection point at day 60
        assert "Inflection" in html or "inflection" in html
        assert "60" in html

    def test_trend_descriptions_shown(self):
        """Trend descriptions are included beyond just direction tags."""
        result = _section_appendix(SAMPLE_GRAPH_CSV_ROWS)
        html = "\n".join(result)
        assert "ADC rises over time" in html or "ADC" in html


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
        img_dir = saved_files_dir / "Standard"
        img_path = img_dir / "test_figure.png"
        img_path.write_bytes(make_tiny_png())

        result = _section_figure_gallery(saved_files_dir)
        html = "\n".join(result)
        assert "Figure Gallery" in html
        assert "data:image/png;base64," in html
        assert "Test Figure" in html  # filename -> title

    def test_groups_by_dwi_type(self, saved_files_dir: Path):
        """Gallery groups figures by DWI type."""
        reset_numbering()
        png = make_tiny_png()
        (saved_files_dir / "Standard" / "fig1.png").write_bytes(png)
        (saved_files_dir / "dnCNN" / "fig2.png").write_bytes(png)

        result = _section_figure_gallery(saved_files_dir)
        html = "\n".join(result)
        assert "badge-standard" in html
        assert "badge-dncnn" in html
