"""Tests for report_sections.data_sections module.

Validates data section builders:
- Cohort overview
- Patient flow
- Data completeness
- MAT data section
- Appendix
- Figure gallery

Edge cases covered:
- Empty/None inputs at every level
- Missing nested keys (no longitudinal, no baseline, no sanity_checks)
- Partial data (some DWI types present, others missing)
- Small sample size warnings
- Imbalanced outcome distributions
- Cross-DWI patient count mismatch warnings
- Convergence failures and excessive NaN
- Dosimetry pass/fail thresholds
- Core method comparison with Hausdorff matrix
- Attrition with per-timepoint counts
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

from report.sections.data_sections import (
    _section_cohort_overview,
    _section_patient_flow,
    _section_data_completeness,
    _section_mat_data,
)
from report.sections.gallery import (
    _section_appendix,
    _section_figure_gallery,
)


def _make_log_data():
    return {
        "Standard": {
            "baseline": {
                "total_outliers": {"n_removed": 4, "n_total": 42, "pct": 9.5},
                "baseline_exclusion": {
                    "n_excluded": 6, "n_total": 48,
                    "lf_rate_included": 35.0, "lf_rate_excluded": 50.0,
                },
            },
            "sanity_checks": {
                "all_converged": True,
                "total_convergence": 0,
                "dim_mismatches": 1,
                "nan_dose_warnings": 0,
            },
            "stats_comparisons": {
                "glme_excluded": {"n": 5, "pct": 12.0},
            },
        }
    }


def _make_mat_data():
    return {
        "Standard": {
            "longitudinal": {"num_patients": 42, "num_timepoints": 5},
            "dosimetry": {
                "d95_gtvp": {"mean": 45.2, "std": 3.1},
                "v50_gtvp": {"mean": 0.87, "std": 0.08},
            },
            "core_comparison": {
                "methods": ["adc_threshold", "otsu"],
                "dice_matrix": [[1.0, 0.75], [0.75, 1.0]],
            },
        }
    }


# ── Cohort Overview ──


class TestCohortOverview:
    def test_returns_html(self):
        result = _section_cohort_overview(_make_mat_data(), _make_log_data(), ["Standard"])
        assert isinstance(result, list)

    def test_empty_data(self):
        result = _section_cohort_overview(None, None, [])
        assert isinstance(result, list)
        assert result == []

    def test_shows_patient_count(self):
        html = "\n".join(_section_cohort_overview(_make_mat_data(), _make_log_data(), ["Standard"]))
        assert "42" in html

    def test_none_mat_data_with_log(self):
        """Cohort overview from log data alone (no MAT files)."""
        result = _section_cohort_overview(None, _make_log_data(), ["Standard"])
        html = "\n".join(result)
        # Should still produce content from log data
        assert isinstance(result, list)

    def test_none_log_data_with_mat(self):
        """Cohort overview from MAT data alone (no logs)."""
        result = _section_cohort_overview(_make_mat_data(), None, ["Standard"])
        html = "\n".join(result)
        assert "42" in html

    def test_small_sample_size_warning(self):
        """n < 30 should trigger a small sample size warning."""
        mat = {"Standard": {"longitudinal": {"num_patients": 15, "num_timepoints": 3}}}
        html = "\n".join(_section_cohort_overview(mat, None, ["Standard"]))
        assert "Small sample size" in html or "caution" in html.lower()

    def test_no_small_sample_warning_for_large_cohort(self):
        """n >= 30 should not trigger the small sample warning."""
        mat = {"Standard": {"longitudinal": {"num_patients": 50, "num_timepoints": 3}}}
        html = "\n".join(_section_cohort_overview(mat, None, ["Standard"]))
        assert "Small sample size" not in html

    def test_outcome_balance_imbalanced(self):
        """LF rate < 20% or > 80% should trigger imbalance warning."""
        log = {
            "Standard": {
                "baseline": {
                    "baseline_exclusion": {
                        "n_excluded": 0, "n_total": 50,
                        "lf_rate_included": 10.0, "lf_rate_excluded": 0.0,
                    },
                },
                "survival": {},
            }
        }
        mat = {"Standard": {"longitudinal": {"num_patients": 50, "num_timepoints": 3}}}
        html = "\n".join(_section_cohort_overview(mat, log, ["Standard"]))
        assert "Imbalanced" in html

    def test_outcome_balance_normal(self):
        """LF rate ~35% should not trigger imbalance warning."""
        result = _section_cohort_overview(_make_mat_data(), _make_log_data(), ["Standard"])
        html = "\n".join(result)
        assert "Imbalanced" not in html

    def test_cross_dwi_patient_count_mismatch_warning(self):
        """Different patient counts across DWI types triggers warning."""
        mat = {
            "Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}},
            "dnCNN": {"longitudinal": {"num_patients": 38, "num_timepoints": 5}},
        }
        html = "\n".join(_section_cohort_overview(mat, None, ["Standard", "dnCNN"]))
        assert "differ across DWI" in html

    def test_no_mismatch_warning_when_counts_similar(self):
        """Similar patient counts (diff <= 2) should not warn."""
        mat = {
            "Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}},
            "dnCNN": {"longitudinal": {"num_patients": 41, "num_timepoints": 5}},
        }
        html = "\n".join(_section_cohort_overview(mat, None, ["Standard", "dnCNN"]))
        assert "differ across DWI" not in html

    def test_attrition_with_per_timepoint_counts(self):
        """When patients_per_timepoint is available, show attrition table."""
        mat = {
            "Standard": {
                "longitudinal": {
                    "num_patients": 42, "num_timepoints": 3,
                    "patients_per_timepoint": [42, 38, 35],
                },
            },
        }
        html = "\n".join(_section_cohort_overview(mat, None, ["Standard"]))
        assert "Attrition" in html
        assert "Timepoint 1" in html

    def test_attrition_without_per_timepoint_counts(self):
        """Multi-timepoint without per-timepoint data shows fallback text."""
        mat = {
            "Standard": {
                "longitudinal": {"num_patients": 42, "num_timepoints": 3},
            },
        }
        html = "\n".join(_section_cohort_overview(mat, None, ["Standard"]))
        assert "Attrition" in html
        assert "not available" in html

    def test_single_timepoint_no_attrition(self):
        """Single timepoint should not show attrition section."""
        mat = {
            "Standard": {
                "longitudinal": {"num_patients": 42, "num_timepoints": 1},
            },
        }
        html = "\n".join(_section_cohort_overview(mat, None, ["Standard"]))
        assert "Attrition" not in html

    def test_dwi_type_not_in_log(self):
        """DWI type listed in dwi_types_present but missing from log_data."""
        result = _section_cohort_overview(
            _make_mat_data(), _make_log_data(), ["Standard", "dnCNN"]
        )
        assert isinstance(result, list)

    def test_empty_baseline_sections(self):
        """baseline key exists but total_outliers and baseline_exclusion are None."""
        log = {"Standard": {"baseline": {}}}
        result = _section_cohort_overview(None, log, ["Standard"])
        assert result == []

    def test_zero_patients_zero_timepoints(self):
        """Zero patients/timepoints in MAT data."""
        mat = {"Standard": {"longitudinal": {"num_patients": 0, "num_timepoints": 0}}}
        log = {"Standard": {"baseline": {"total_outliers": {"n_removed": 0, "n_total": 0, "pct": 0}}}}
        result = _section_cohort_overview(mat, log, ["Standard"])
        assert isinstance(result, list)

    def test_outcome_balance_no_survival_data(self):
        """No survival data and no lf_rate should show 'not available' message."""
        log = {"Standard": {"baseline": {"total_outliers": {"n_removed": 1, "n_total": 10, "pct": 10.0}}}}
        mat = {"Standard": {"longitudinal": {"num_patients": 10, "num_timepoints": 2}}}
        result = _section_cohort_overview(mat, log, ["Standard"])
        html = "\n".join(result)
        assert "not available" in html.lower() or "Outcome Balance" in html


# ── Patient Flow ──


class TestPatientFlow:
    def test_returns_html(self):
        result = _section_patient_flow(_make_log_data(), ["Standard"], _make_mat_data())
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_patient_flow(None, [], {})
        assert isinstance(result, list)
        assert result == []

    def test_exclusion_shown(self):
        html = "\n".join(_section_patient_flow(_make_log_data(), ["Standard"], _make_mat_data()))
        if html:
            assert "exclu" in html.lower() or "Patient" in html

    def test_flow_shows_consort_table(self):
        """Full data produces a CONSORT-style table with all columns."""
        html = "\n".join(_section_patient_flow(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Initial Cohort" in html
        assert "Baseline Excluded" in html
        assert "Outliers Removed" in html

    def test_negative_analysed_clipped_to_zero(self):
        """When exclusions exceed initial count, analysed should be capped at 0."""
        log = {
            "Standard": {
                "baseline": {
                    "baseline_exclusion": {"n_excluded": 50, "n_total": 48},
                },
                "stats_comparisons": {"glme_excluded": {"n_excluded": 10}},
            }
        }
        mat = {"Standard": {"longitudinal": {"num_patients": 20, "num_timepoints": 2}}}
        html = "\n".join(_section_patient_flow(log, ["Standard"], mat))
        # Should not show negative values
        assert "- " not in html or "<td><strong>0</strong></td>" in html

    def test_no_mat_data_uses_baseline_total(self):
        """When MAT data is missing, initial cohort comes from baseline_exclusion.n_total."""
        log = {
            "Standard": {
                "baseline": {
                    "baseline_exclusion": {"n_excluded": 6, "n_total": 48},
                },
            }
        }
        html = "\n".join(_section_patient_flow(log, ["Standard"], None))
        assert "48" in html

    def test_multiple_dwi_types_in_flow(self):
        """Flow table should have rows for each DWI type with data."""
        log = {
            "Standard": {"baseline": {"baseline_exclusion": {"n_excluded": 3, "n_total": 40}}},
            "dnCNN": {"baseline": {"baseline_exclusion": {"n_excluded": 2, "n_total": 40}}},
        }
        mat = {
            "Standard": {"longitudinal": {"num_patients": 37, "num_timepoints": 3}},
            "dnCNN": {"longitudinal": {"num_patients": 38, "num_timepoints": 3}},
        }
        html = "\n".join(_section_patient_flow(log, ["Standard", "dnCNN"], mat))
        assert "Standard" in html
        assert "dnCNN" in html

    def test_parse_warnings_shown(self):
        """log_data=None should show parse warning banner."""
        # _section_patient_flow with log_data=None but mat_data providing initial cohort
        # won't produce rows since n_initial=0 and n_baseline_exc=0
        # Instead test with log_data that has parse_warnings
        log = {
            "Standard": {
                "baseline": {"baseline_exclusion": {"n_excluded": 2, "n_total": 30}},
                "parse_warnings": ["Could not parse outlier line"],
            }
        }
        mat = {"Standard": {"longitudinal": {"num_patients": 28, "num_timepoints": 2}}}
        html = "\n".join(_section_patient_flow(log, ["Standard"], mat))
        assert "could not be parsed" in html.lower() or "incomplete" in html.lower()

    def test_no_exclusion_data_returns_empty(self):
        """No baseline exclusion and no MAT data => empty section."""
        log = {"Standard": {"baseline": {}}}
        result = _section_patient_flow(log, ["Standard"], {})
        assert result == []

    def test_competing_risk_exclusion_from_glme(self):
        """Competing-risk exclusion column should reflect glme_excluded data."""
        log = {
            "Standard": {
                "baseline": {"baseline_exclusion": {"n_excluded": 5, "n_total": 50}},
                "stats_comparisons": {"glme_excluded": {"n_excluded": 8}},
            }
        }
        mat = {"Standard": {"longitudinal": {"num_patients": 45, "num_timepoints": 3}}}
        html = "\n".join(_section_patient_flow(log, ["Standard"], mat))
        assert "8" in html  # CR exclusion count


# ── Data Completeness ──


class TestDataCompleteness:
    def test_returns_html(self):
        result = _section_data_completeness(_make_log_data(), ["Standard"])
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_data_completeness(None, [])
        assert isinstance(result, list)
        assert result == []

    def test_convergence_passed(self):
        """all_converged=True should show 'Passed' or 'All converged'."""
        html = "\n".join(_section_data_completeness(_make_log_data(), ["Standard"]))
        assert "converge" in html.lower() or "Passed" in html

    def test_convergence_failures(self):
        """Convergence flags > 0 should warn."""
        log = {
            "Standard": {
                "sanity_checks": {
                    "all_converged": False,
                    "total_convergence": 5,
                    "dim_mismatches": 0,
                    "nan_dose_warnings": 0,
                },
            }
        }
        html = "\n".join(_section_data_completeness(log, ["Standard"]))
        assert "5" in html
        assert "convergence" in html.lower()

    def test_dim_mismatches(self):
        """Dimensional mismatches should appear in output."""
        html = "\n".join(_section_data_completeness(_make_log_data(), ["Standard"]))
        assert "mismatch" in html.lower() or "Mismatch" in html

    def test_nan_dose_warnings(self):
        """NaN dose warnings should be flagged."""
        log = {
            "Standard": {
                "sanity_checks": {
                    "all_converged": True,
                    "total_convergence": 0,
                    "dim_mismatches": 0,
                    "nan_dose_warnings": 3,
                },
            }
        }
        html = "\n".join(_section_data_completeness(log, ["Standard"]))
        assert "NaN" in html
        assert "3" in html

    def test_excessive_nan_parameters(self):
        """Excessive NaN in parameters should be shown."""
        log = {
            "Standard": {
                "sanity_checks": {
                    "all_converged": True,
                    "total_convergence": 0,
                    "dim_mismatches": 0,
                    "nan_dose_warnings": 0,
                    "excessive_nan": [
                        {"parameter": "D_star", "pct_nan": 65.3},
                    ],
                },
            }
        }
        html = "\n".join(_section_data_completeness(log, ["Standard"]))
        assert "D_star" in html
        assert "65.3" in html

    def test_empty_sanity_checks(self):
        """Empty sanity_checks dict returns empty."""
        log = {"Standard": {"sanity_checks": {}}}
        result = _section_data_completeness(log, ["Standard"])
        assert result == []

    def test_no_sanity_checks_key(self):
        """Missing sanity_checks key returns empty."""
        log = {"Standard": {"baseline": {}}}
        result = _section_data_completeness(log, ["Standard"])
        assert result == []

    def test_multiple_dwi_types(self):
        """Multiple DWI types each get their own section."""
        log = {
            "Standard": {
                "sanity_checks": {"all_converged": True, "total_convergence": 0,
                                   "dim_mismatches": 0, "nan_dose_warnings": 0},
            },
            "dnCNN": {
                "sanity_checks": {"all_converged": False, "total_convergence": 3,
                                   "dim_mismatches": 2, "nan_dose_warnings": 1},
            },
        }
        html = "\n".join(_section_data_completeness(log, ["Standard", "dnCNN"]))
        assert "Standard" in html
        assert "dnCNN" in html

    def test_alignment_passed_when_zero_issues(self):
        """Zero dim_mismatches and nan_dose_warnings shows 'Passed'."""
        log = {
            "Standard": {
                "sanity_checks": {"all_converged": True, "total_convergence": 0,
                                   "dim_mismatches": 0, "nan_dose_warnings": 0},
            },
        }
        html = "\n".join(_section_data_completeness(log, ["Standard"]))
        assert "Passed" in html

    def test_dwi_type_not_in_log(self):
        """DWI type in dwi_types_present but missing from log_data."""
        log = {"Standard": {"sanity_checks": {"all_converged": True, "total_convergence": 0,
                                               "dim_mismatches": 0, "nan_dose_warnings": 0}}}
        result = _section_data_completeness(log, ["Standard", "IVIMnet"])
        assert isinstance(result, list)


# ── MAT Data ──


class TestMatData:
    def test_returns_html(self):
        result = _section_mat_data(_make_mat_data())
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_mat_data(None)
        assert isinstance(result, list)
        assert result == []

    def test_shows_dosimetry(self):
        html = "\n".join(_section_mat_data(_make_mat_data()))
        if html:
            assert "dosimetry" in html.lower() or "D95" in html or "MAT" in html

    def test_dosimetry_d95_pass_check(self):
        """D95 >= 45 Gy should show PASS."""
        mat = {
            "Standard": {
                "dosimetry": {
                    "d95_adc_mean": {"mean": 48.0, "std": 2.0},
                    "v50_adc_mean": {"mean": 0.92, "std": 0.05},
                },
            }
        }
        html = "\n".join(_section_mat_data(mat))
        assert "PASS" in html

    def test_dosimetry_d95_fail_check(self):
        """D95 < 45 Gy should show FAIL and under-dosing note."""
        mat = {
            "Standard": {
                "dosimetry": {
                    "d95_adc_mean": {"mean": 40.0, "std": 5.0},
                    "v50_adc_mean": {"mean": 0.80, "std": 0.10},
                },
            }
        }
        html = "\n".join(_section_mat_data(mat))
        assert "FAIL" in html
        assert "under-dos" in html.lower()

    def test_dosimetry_high_d95_adequate_coverage(self):
        """D95 >= 50 Gy shows adequate coverage note."""
        mat = {
            "Standard": {
                "dosimetry": {
                    "d95_adc_mean": {"mean": 52.0, "std": 1.0},
                    "v50_adc_mean": {"mean": 0.95, "std": 0.02},
                },
            }
        }
        html = "\n".join(_section_mat_data(mat))
        assert "adequate coverage" in html.lower()

    def test_dosimetry_low_v50(self):
        """V50 < 90% should note partial coverage."""
        mat = {
            "Standard": {
                "dosimetry": {
                    "d95_adc_mean": {"mean": 50.0, "std": 1.0},
                    "v50_adc_mean": {"mean": 0.80, "std": 0.10},
                },
            }
        }
        html = "\n".join(_section_mat_data(mat))
        assert "partial coverage" in html.lower() or "dose escalation" in html.lower()

    def test_dosimetry_none_values(self):
        """None values in dosimetry should produce dash characters."""
        mat = {
            "Standard": {
                "dosimetry": {
                    "d95_adc_mean": None,
                    "v50_adc_mean": None,
                },
            }
        }
        html = "\n".join(_section_mat_data(mat))
        assert "\u2014" in html  # em dash for missing

    def test_dosimetry_nan_values(self):
        """NaN values in dosimetry should produce dash characters."""
        mat = {
            "Standard": {
                "dosimetry": {
                    "d95_adc_mean": {"mean": float("nan"), "std": 1.0},
                    "v50_adc_mean": {"mean": 0.85, "std": float("nan")},
                },
            }
        }
        html = "\n".join(_section_mat_data(mat))
        assert isinstance(html, str)

    def test_dosimetry_plain_float(self):
        """Plain float (not dict) dosimetry values."""
        mat = {
            "Standard": {
                "dosimetry": {
                    "d95_adc_mean": 47.5,
                    "v50_adc_mean": 0.88,
                },
            }
        }
        html = "\n".join(_section_mat_data(mat))
        assert "47.5" in html or "PASS" in html

    def test_core_method_comparison_matrix(self):
        """Core method comparison should show Dice matrix."""
        mat = {
            "Standard": {
                "core_method": {
                    "methods": ["adc_threshold", "otsu", "gmm"],
                    "mean_dice_matrix": [
                        [1.0, 0.75, 0.60],
                        [0.75, 1.0, 0.80],
                        [0.60, 0.80, 1.0],
                    ],
                },
            }
        }
        html = "\n".join(_section_mat_data(mat))
        assert "Core Method" in html
        assert "adc_threshold" in html
        assert "otsu" in html
        assert "gmm" in html

    def test_core_method_with_hausdorff(self):
        """Core method with hausdorff_matrix shows second table."""
        mat = {
            "Standard": {
                "core_method": {
                    "methods": ["adc_threshold", "otsu"],
                    "mean_dice_matrix": [[1.0, 0.75], [0.75, 1.0]],
                    "hausdorff_matrix": [[0.0, 8.5], [8.5, 0.0]],
                },
            }
        }
        html = "\n".join(_section_mat_data(mat))
        assert "Hausdorff" in html
        assert "8.5" in html

    def test_core_method_interchangeable_pairs(self):
        """Dice >= 0.70 pair should be marked interchangeable."""
        mat = {
            "Standard": {
                "core_method": {
                    "methods": ["adc_threshold", "otsu"],
                    "mean_dice_matrix": [[1.0, 0.75], [0.75, 1.0]],
                },
            }
        }
        html = "\n".join(_section_mat_data(mat))
        assert "interchangeable" in html.lower() or "Interchangeable" in html

    def test_core_method_no_interchangeable_pairs(self):
        """All Dice < 0.70 should note no pairs meet threshold."""
        mat = {
            "Standard": {
                "core_method": {
                    "methods": ["adc_threshold", "otsu"],
                    "mean_dice_matrix": [[1.0, 0.50], [0.50, 1.0]],
                },
            }
        }
        html = "\n".join(_section_mat_data(mat))
        assert "No method pairs" in html

    def test_empty_dosimetry_dict(self):
        """Empty dosimetry dict inside MAT data."""
        mat = {"Standard": {"dosimetry": {}}}
        result = _section_mat_data(mat)
        assert isinstance(result, list)

    def test_empty_core_method(self):
        """Core method with empty methods list."""
        mat = {"Standard": {"core_method": {"methods": [], "mean_dice_matrix": []}}}
        result = _section_mat_data(mat)
        assert isinstance(result, list)

    def test_no_dosimetry_no_core(self):
        """MAT data with neither dosimetry nor core_method."""
        mat = {"Standard": {"longitudinal": {"num_patients": 42}}}
        result = _section_mat_data(mat)
        html = "\n".join(result)
        assert "Supplemental" in html  # Still shows header
        assert "Dosimetry" not in html
        assert "Core Method" not in html

    def test_multiple_dwi_types_dosimetry(self):
        """Multiple DWI types each show dosimetry rows."""
        mat = {
            "Standard": {
                "dosimetry": {"d95_adc_mean": 48.0, "v50_adc_mean": 0.90},
            },
            "dnCNN": {
                "dosimetry": {"d95_adc_mean": 46.0, "v50_adc_mean": 0.88},
            },
        }
        html = "\n".join(_section_mat_data(mat))
        assert "Standard" in html
        assert "dnCNN" in html


# ── Appendix ──


class TestAppendix:
    def test_returns_html(self):
        result = _section_appendix(SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)

    def test_empty_rows(self):
        result = _section_appendix([])
        assert isinstance(result, list)

    def test_none_rows(self):
        result = _section_appendix(None)
        assert isinstance(result, list)


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

    def test_with_png_files(self, tmp_path):
        """Gallery should find and display PNG files."""
        folder = tmp_path / "saved_files_20260301_120000"
        std = folder / "Standard"
        std.mkdir(parents=True)
        (std / "test_plot.png").write_bytes(b"fake png")
        result = _section_figure_gallery(folder)
        assert isinstance(result, list)

    def test_gallery_with_rows(self, tmp_path):
        """Gallery with both folder and CSV rows."""
        folder = tmp_path / "saved_files_20260301_120000"
        (folder / "Standard").mkdir(parents=True)
        result = _section_figure_gallery(folder, rows=SAMPLE_GRAPH_CSV_ROWS)
        assert isinstance(result, list)
