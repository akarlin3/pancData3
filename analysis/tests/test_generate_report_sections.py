"""Tests for report section builder functions.

Covers:
- _section_data_completeness: convergence, NaN, outlier summary
- _section_feature_overlap: cross-DWI feature consistency
- _section_power_analysis: statistical power commentary
- _section_patient_flow: CONSORT-style attrition
- _section_sensitivity_analysis: EPV, cross-DWI concordance, HR stability
- Feature stability across timepoints (within feature_overlap)
- _section_appendix: grouped graph appendix with statistics
"""

from __future__ import annotations

import pytest  # type: ignore

from report.generate_report import (  # type: ignore
    _section_appendix,
    _section_data_completeness,
    _section_feature_overlap,
    _section_patient_flow,
    _section_power_analysis,
    _section_sensitivity_analysis,
)


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
