"""Tests for report sections: enrollment/supplemental data sections (part 2).

Validates data completeness and MAT data section builders. Cohort overview
and patient flow tests live in test_report_sections_data_sections.py.

Shared helpers (_make_log_data, _make_mat_data) live in
_test_report_sections_data_sections_shared.py.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from _test_report_sections_data_sections_shared import *  # noqa: F401,F403
from _test_report_sections_data_sections_shared import (
    _make_log_data,
    _make_mat_data,
    _section_data_completeness,
    _section_mat_data,
)


# ── Data Completeness ──


class TestDataCompleteness:
    def test_returns_html(self):
        result = _section_data_completeness(_make_log_data(), ["Standard"])
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_data_completeness(None, [])
        assert isinstance(result, list)
        assert result == []

    def test_all_converged_shows_passed(self):
        """When all_converged is True, should show 'Passed'."""
        html = "\n".join(_section_data_completeness(_make_log_data(), ["Standard"]))
        assert "Passed" in html or "All converged" in html

    def test_convergence_flags_shown(self):
        """When total_convergence > 0, should show convergence flags."""
        log = _make_log_data()
        log["Standard"]["sanity_checks"]["all_converged"] = False
        log["Standard"]["sanity_checks"]["total_convergence"] = 7
        html = "\n".join(_section_data_completeness(log, ["Standard"]))
        assert "7" in html
        assert "Conv" in html

    def test_dim_mismatches_shown(self):
        """Dimensional mismatches should appear in the output."""
        html = "\n".join(_section_data_completeness(_make_log_data(), ["Standard"]))
        assert "Mismatch" in html or "mismatch" in html.lower()

    def test_nan_dose_warnings_shown(self):
        """NaN dose warnings should appear when present."""
        log = _make_log_data()
        log["Standard"]["sanity_checks"]["nan_dose_warnings"] = 3
        html = "\n".join(_section_data_completeness(log, ["Standard"]))
        assert "NaN" in html

    def test_excessive_nan_shown(self):
        """Excessive NaN parameters should be listed."""
        log = _make_log_data()
        log["Standard"]["sanity_checks"]["excessive_nan"] = [
            {"parameter": "D_star", "pct_nan": 62.5},
        ]
        html = "\n".join(_section_data_completeness(log, ["Standard"]))
        assert "D_star" in html
        assert "62.5" in html

    def test_alignment_passed_when_no_issues(self):
        """Should show alignment 'Passed' when no mismatches and no NaN warnings."""
        log = _make_log_data()
        log["Standard"]["sanity_checks"]["dim_mismatches"] = 0
        log["Standard"]["sanity_checks"]["nan_dose_warnings"] = 0
        html = "\n".join(_section_data_completeness(log, ["Standard"]))
        # All converged + alignment passed
        assert "Passed" in html

    def test_no_sanity_checks_data(self):
        """No sanity_checks key should gracefully return empty."""
        log = {"Standard": {"baseline": {}}}
        result = _section_data_completeness(log, ["Standard"])
        assert result == []

    def test_empty_sanity_checks(self):
        """Empty sanity_checks dict should be skipped."""
        log = {"Standard": {"sanity_checks": {}}}
        result = _section_data_completeness(log, ["Standard"])
        assert result == []

    def test_stat_cards_present(self):
        """Should render stat cards for convergence and alignment status."""
        html = "\n".join(_section_data_completeness(_make_log_data(), ["Standard"]))
        assert "stat-grid" in html

    def test_multiple_dwi_types(self):
        """Should render sections per DWI type."""
        log = _make_log_data()
        log["dnCNN"] = {
            "sanity_checks": {
                "all_converged": False,
                "total_convergence": 3,
                "dim_mismatches": 0,
                "nan_dose_warnings": 0,
            }
        }
        html = "\n".join(_section_data_completeness(log, ["Standard", "dnCNN"]))
        assert "Standard" in html
        assert "dnCNN" in html

    def test_total_outliers_stat_card(self):
        """Sanity outlier count should appear as a stat card."""
        log = _make_log_data()
        log["Standard"]["sanity_checks"]["total_outliers"] = 5
        html = "\n".join(_section_data_completeness(log, ["Standard"]))
        assert "Sanity Outliers" in html or "5" in html


# ── MAT Data ──


class TestMatData:
    def test_returns_html(self):
        result = _section_mat_data(_make_mat_data())
        assert isinstance(result, list)

    def test_no_data(self):
        result = _section_mat_data(None)
        assert isinstance(result, list)
        assert result == []

    def test_empty_dict(self):
        result = _section_mat_data({})
        assert isinstance(result, list)
        assert result == []

    def test_shows_dosimetry_header(self):
        html = "\n".join(_section_mat_data(_make_mat_data()))
        assert "Dosimetry" in html

    def test_d95_values_formatted(self):
        """D95 values should appear with Gy units."""
        html = "\n".join(_section_mat_data(_make_mat_data()))
        assert "45.2" in html
        assert "Gy" in html

    def test_v50_values_formatted_as_pct(self):
        """V50 values should be converted to percentage."""
        html = "\n".join(_section_mat_data(_make_mat_data()))
        assert "87.0%" in html or "87" in html

    def test_d95_pass_fail_check(self):
        """D95 >= 45 Gy should show PASS."""
        html = "\n".join(_section_mat_data(_make_mat_data()))
        assert "PASS" in html or "✅" in html

    def test_d95_fail_below_threshold(self):
        """D95 < 45 Gy should show FAIL."""
        mat = _make_mat_data()
        mat["Standard"]["dosimetry"]["d95_adc_mean"] = {"mean": 40.0, "std": 2.0}
        html = "\n".join(_section_mat_data(mat))
        assert "FAIL" in html or "❌" in html

    def test_clinical_reference_box(self):
        """Should include clinical reference context for dosimetry."""
        html = "\n".join(_section_mat_data(_make_mat_data()))
        assert "Clinical Reference" in html

    def test_dosimetric_interpretation_notes(self):
        """Should show interpretation notes (e.g., under-dosing, adequate coverage)."""
        html = "\n".join(_section_mat_data(_make_mat_data()))
        assert "adequate" in html.lower() or "coverage" in html.lower()

    def test_core_method_comparison_header(self):
        """Should show Core Method Comparison section."""
        html = "\n".join(_section_mat_data(_make_mat_data()))
        assert "Core Method Comparison" in html

    def test_core_dice_matrix_rendered(self):
        """Mean Dice matrix should be rendered as a table."""
        html = "\n".join(_section_mat_data(_make_mat_data()))
        assert "0.75" in html
        assert "adc_threshold" in html
        assert "otsu" in html

    def test_core_summary_statistics(self):
        """Should compute and show mean pairwise Dice."""
        html = "\n".join(_section_mat_data(_make_mat_data()))
        assert "pairwise Dice" in html

    def test_recommendation_section(self):
        """Should recommend interchangeable methods or reference method."""
        html = "\n".join(_section_mat_data(_make_mat_data()))
        assert "Recommendation" in html

    def test_interchangeable_pairs(self):
        """Methods with Dice >= 0.70 should be listed as interchangeable."""
        html = "\n".join(_section_mat_data(_make_mat_data()))
        # adc_threshold <-> otsu have Dice 0.75 >= 0.70
        assert "interchangeable" in html.lower() or "↔" in html

    def test_hausdorff_matrix_shown_when_present(self):
        """Hausdorff matrix should be rendered when available."""
        mat = _make_mat_data()
        mat["Standard"]["core_method"]["hausdorff_matrix"] = [
            [0.0, 3.5, 8.0],
            [3.5, 0.0, 6.0],
            [8.0, 6.0, 0.0],
        ]
        html = "\n".join(_section_mat_data(mat))
        assert "Hausdorff" in html
        assert "3.5" in html

    def test_no_hausdorff_when_absent(self):
        """No hausdorff section when key is missing."""
        html = "\n".join(_section_mat_data(_make_mat_data()))
        assert "Hausdorff" not in html

    def test_dosimetry_none_values(self):
        """None values in dosimetry should render as dashes."""
        mat = _make_mat_data()
        mat["Standard"]["dosimetry"]["d95_adc_mean"] = None
        mat["Standard"]["dosimetry"]["v50_adc_mean"] = None
        html = "\n".join(_section_mat_data(mat))
        assert "—" in html  # em-dash

    def test_dosimetry_nan_values(self):
        """NaN values should render as dashes."""
        mat = _make_mat_data()
        mat["Standard"]["dosimetry"]["d95_adc_mean"] = float("nan")
        html = "\n".join(_section_mat_data(mat))
        assert "—" in html

    def test_dosimetry_plain_float(self):
        """Plain float values (not dict) should be formatted."""
        mat = _make_mat_data()
        mat["Standard"]["dosimetry"]["d95_adc_mean"] = 46.0
        html = "\n".join(_section_mat_data(mat))
        assert "46.0" in html

    def test_core_method_empty_methods(self):
        """Empty methods list should be skipped gracefully."""
        mat = _make_mat_data()
        mat["Standard"]["core_method"]["methods"] = []
        mat["Standard"]["core_method"]["mean_dice_matrix"] = []
        result = _section_mat_data(mat)
        assert isinstance(result, list)

    def test_multiple_dwi_types(self):
        """Should render dosimetry for multiple DWI types."""
        mat = _make_mat_data()
        mat["dnCNN"] = {
            "dosimetry": {
                "d95_adc_mean": {"mean": 44.0, "std": 2.0},
                "v50_adc_mean": {"mean": 0.82, "std": 0.1},
            },
        }
        html = "\n".join(_section_mat_data(mat))
        assert "Standard" in html
        assert "dnCNN" in html

    def test_dosimetry_only(self):
        """Mat data with only dosimetry (no core_method) should work."""
        mat = {"Standard": {
            "dosimetry": {"d95_adc_mean": {"mean": 50.0}, "v50_adc_mean": {"mean": 0.95}},
        }}
        result = _section_mat_data(mat)
        html = "\n".join(result)
        assert "Dosimetry" in html
        assert "Core Method" not in html

    def test_core_method_only(self):
        """Mat data with only core_method (no dosimetry) should work."""
        mat = {"Standard": {
            "core_method": {
                "methods": ["adc_threshold", "otsu"],
                "mean_dice_matrix": [[1.0, 0.8], [0.8, 1.0]],
            },
        }}
        result = _section_mat_data(mat)
        html = "\n".join(result)
        assert "Dosimetry" not in html
        assert "Core Method" in html

    def test_v50_above_one_treated_as_pct(self):
        """V50 value > 1.0 should NOT be multiplied by 100."""
        mat = _make_mat_data()
        mat["Standard"]["dosimetry"]["v50_adc_mean"] = {"mean": 85.0, "std": 5.0}
        html = "\n".join(_section_mat_data(mat))
        assert "85.0" in html
        # Should NOT show 8500
        assert "8500" not in html


# ── Edge cases: missing metrics, convergence failures ──


class TestDataCompletenessEdgeCases:
    def test_convergence_failure_all_types(self):
        log = {
            "Standard": {"sanity_checks": {
                "all_converged": False, "total_convergence": 15,
                "dim_mismatches": 5, "nan_dose_warnings": 3,
            }},
            "dnCNN": {"sanity_checks": {
                "all_converged": False, "total_convergence": 20,
                "dim_mismatches": 8, "nan_dose_warnings": 0,
            }},
        }
        html = "\n".join(_section_data_completeness(log, ["Standard", "dnCNN"]))
        assert "15" in html
        assert "Standard" in html
        assert "dnCNN" in html

    def test_partial_sanity_data(self):
        log = {
            "Standard": {"sanity_checks": {
                "all_converged": True, "total_convergence": 0,
                "dim_mismatches": 0, "nan_dose_warnings": 0,
            }},
            "dnCNN": {"baseline": {}},  # No sanity_checks key
        }
        result = _section_data_completeness(log, ["Standard", "dnCNN"])
        html = "\n".join(result)
        assert "Standard" in html

    def test_excessive_nan_multiple_params(self):
        log = {"Standard": {"sanity_checks": {
            "all_converged": True, "total_convergence": 0,
            "dim_mismatches": 0, "nan_dose_warnings": 0,
            "excessive_nan": [
                {"parameter": "D_star", "pct_nan": 75.0},
                {"parameter": "f", "pct_nan": 60.0},
                {"parameter": "ADC", "pct_nan": 55.0},
            ],
        }}}
        html = "\n".join(_section_data_completeness(log, ["Standard"]))
        assert "D_star" in html
        assert "f" in html
        assert "ADC" in html

    def test_dwi_type_in_list_but_not_in_log(self):
        log = {"Standard": {"sanity_checks": {
            "all_converged": True, "total_convergence": 0,
            "dim_mismatches": 0, "nan_dose_warnings": 0,
        }}}
        result = _section_data_completeness(log, ["Standard", "IVIMnet"])
        assert isinstance(result, list)

    def test_all_zeros_no_issues(self):
        log = {"Standard": {"sanity_checks": {
            "all_converged": True, "total_convergence": 0,
            "dim_mismatches": 0, "nan_dose_warnings": 0,
        }}}
        html = "\n".join(_section_data_completeness(log, ["Standard"]))
        assert "Passed" in html


class TestMatDataEdgeCases:
    def test_dosimetry_all_none(self):
        mat = {"Standard": {"dosimetry": {
            "d95_adc_mean": None, "v50_adc_mean": None,
            "d95_d_mean": None, "v50_d_mean": None,
        }}}
        result = _section_mat_data(mat)
        html = "\n".join(result)
        assert "—" in html

    def test_dosimetry_partial_data(self):
        mat = {"Standard": {"dosimetry": {
            "d95_adc_mean": {"mean": 46.0, "std": 2.0},
        }}}
        result = _section_mat_data(mat)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "46.0" in html

    def test_core_method_single_method(self):
        mat = {"Standard": {"core_method": {
            "methods": ["adc_threshold"],
            "mean_dice_matrix": [[1.0]],
        }}}
        result = _section_mat_data(mat)
        assert isinstance(result, list)

    def test_core_method_mismatched_matrix(self):
        mat = {"Standard": {"core_method": {
            "methods": ["adc_threshold", "otsu"],
            "mean_dice_matrix": [[1.0]],  # Wrong dimensions
        }}}
        result = _section_mat_data(mat)
        assert isinstance(result, list)

    def test_dosimetry_v50_exactly_one(self):
        mat = {"Standard": {"dosimetry": {
            "v50_adc_mean": {"mean": 1.0, "std": 0.0},
        }}}
        result = _section_mat_data(mat)
        html = "\n".join(result)
        assert "100.0%" in html or "100" in html

    def test_dosimetry_d95_exactly_45(self):
        mat = {"Standard": {"dosimetry": {
            "d95_adc_mean": {"mean": 45.0, "std": 0.0},
        }}}
        result = _section_mat_data(mat)
        html = "\n".join(result)
        assert "PASS" in html or "✅" in html

    def test_no_dosimetry_no_core(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 42}}}
        result = _section_mat_data(mat)
        assert isinstance(result, list)

    def test_core_all_high_dice(self):
        mat = {"Standard": {"core_method": {
            "methods": ["adc_threshold", "otsu", "gmm"],
            "mean_dice_matrix": [
                [1.0, 0.90, 0.85],
                [0.90, 1.0, 0.88],
                [0.85, 0.88, 1.0],
            ],
        }}}
        result = _section_mat_data(mat)
        html = "\n".join(result)
        assert "interchangeable" in html.lower()

    def test_core_all_low_dice(self):
        mat = {"Standard": {"core_method": {
            "methods": ["adc_threshold", "otsu", "gmm"],
            "mean_dice_matrix": [
                [1.0, 0.30, 0.25],
                [0.30, 1.0, 0.28],
                [0.25, 0.28, 1.0],
            ],
        }}}
        result = _section_mat_data(mat)
        html = "\n".join(result)
        assert isinstance(html, str)

    def test_dosimetry_dict_without_std(self):
        mat = {"Standard": {"dosimetry": {
            "d95_adc_mean": {"mean": 47.5},
            "v50_adc_mean": {"mean": 0.9},
        }}}
        result = _section_mat_data(mat)
        html = "\n".join(result)
        assert "47.5" in html
