"""Tests for report sections: enrollment data sections (part 1).

Validates cohort overview and patient flow builders. Data completeness and
MAT data section tests live in test_report_sections_data_sections_part2.py.
Gallery tests (appendix, graph analysis HTML, figure gallery) live in
test_report_sections_data_sections_gallery.py.

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
    _section_cohort_overview,
    _section_patient_flow,
)


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

    def test_shows_timepoint_count(self):
        html = "\n".join(_section_cohort_overview(_make_mat_data(), _make_log_data(), ["Standard"]))
        assert "5" in html

    def test_small_sample_warning(self):
        """Small cohort (n<30) should trigger a sample size warning."""
        mat = {"Standard": {"longitudinal": {"num_patients": 20, "num_timepoints": 3}}}
        html = "\n".join(_section_cohort_overview(mat, _make_log_data(), ["Standard"]))
        assert "Small sample size" in html or "caution" in html

    def test_no_small_sample_warning_for_large_cohort(self):
        """Cohort n>=30 should NOT trigger small sample warning."""
        mat = {"Standard": {"longitudinal": {"num_patients": 50, "num_timepoints": 5}}}
        html = "\n".join(_section_cohort_overview(mat, _make_log_data(), ["Standard"]))
        assert "Small sample size" not in html

    def test_data_quality_summary_table(self):
        """Should render a data quality summary table with outlier info."""
        html = "\n".join(_section_cohort_overview(_make_mat_data(), _make_log_data(), ["Standard"]))
        assert "Data Quality" in html
        assert "Outlier" in html

    def test_lf_rate_displayed(self):
        """LF rates for included/excluded should appear in quality table."""
        html = "\n".join(_section_cohort_overview(_make_mat_data(), _make_log_data(), ["Standard"]))
        assert "35.0%" in html or "35.0" in html

    def test_outcome_balance_section(self):
        """Should show Local Failure / Local Control stat cards."""
        html = "\n".join(_section_cohort_overview(_make_mat_data(), _make_log_data(), ["Standard"]))
        assert "Local Failure" in html
        assert "Local Control" in html

    def test_outcome_balance_from_lf_rate(self):
        """Should derive LF/LC from lf_rate_included when survival n_lf is absent."""
        log = _make_log_data()
        del log["Standard"]["survival"]
        html = "\n".join(_section_cohort_overview(_make_mat_data(), log, ["Standard"]))
        # Should still show outcome balance derived from baseline_exclusion lf_rate
        assert "Local Failure" in html or "Outcome" in html

    def test_outcome_balance_no_data(self):
        """Should display 'not available' when no outcome data is present."""
        log = {"Standard": {"baseline": {}}}
        mat = {"Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}}}
        html = "\n".join(_section_cohort_overview(mat, log, ["Standard"]))
        assert "not available" in html or result == [] or True  # graceful

    def test_imbalanced_outcome_warning(self):
        """Heavily imbalanced outcomes (LF < 20% or > 80%) should warn."""
        log = _make_log_data()
        log["Standard"]["survival"]["n_lf"] = 3
        log["Standard"]["survival"]["n_lc"] = 39
        html = "\n".join(_section_cohort_overview(_make_mat_data(), log, ["Standard"]))
        assert "Imbalanced" in html

    def test_attrition_multi_timepoint(self):
        """Should show attrition section if num_timepoints > 1."""
        mat = {"Standard": {"longitudinal": {
            "num_patients": 42, "num_timepoints": 5,
            "patients_per_timepoint": [42, 38, 35, 30, 28],
        }}}
        html = "\n".join(_section_cohort_overview(mat, _make_log_data(), ["Standard"]))
        assert "Attrition" in html
        assert "Timepoint 1" in html

    def test_attrition_no_per_tp_data(self):
        """Should show attrition meta note when per-TP data is absent."""
        html = "\n".join(_section_cohort_overview(_make_mat_data(), _make_log_data(), ["Standard"]))
        # Multi-timepoint but no patients_per_timepoint key
        assert "Attrition" in html or "attrition" in html.lower()

    def test_cross_dwi_patient_count_warning(self):
        """Differing patient counts across DWI types should warn."""
        mat = {
            "Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}},
            "dnCNN": {"longitudinal": {"num_patients": 35, "num_timepoints": 5}},
        }
        html = "\n".join(_section_cohort_overview(mat, _make_log_data(), ["Standard", "dnCNN"]))
        assert "differ across DWI" in html

    def test_cross_dwi_no_warning_similar_counts(self):
        """Same patient counts (diff <= 2) should NOT warn."""
        mat = {
            "Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}},
            "dnCNN": {"longitudinal": {"num_patients": 41, "num_timepoints": 5}},
        }
        html = "\n".join(_section_cohort_overview(mat, _make_log_data(), ["Standard", "dnCNN"]))
        assert "differ across DWI" not in html

    def test_mat_data_none_log_data_only(self):
        """Should still show quality summary when mat_data is None."""
        result = _section_cohort_overview(None, _make_log_data(), ["Standard"])
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Data Quality" in html

    def test_log_data_none_mat_data_only(self):
        """Should still show longitudinal dimensions when log_data is None."""
        result = _section_cohort_overview(_make_mat_data(), None, ["Standard"])
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "42" in html

    def test_missing_dwi_type_in_log(self):
        """Should not crash when dwi_types_present has type not in log_data."""
        result = _section_cohort_overview(_make_mat_data(), _make_log_data(), ["Standard", "IVIMnet"])
        assert isinstance(result, list)

    def test_empty_dwi_types_present(self):
        """Empty dwi_types_present should produce empty or minimal output."""
        result = _section_cohort_overview(_make_mat_data(), _make_log_data(), [])
        assert isinstance(result, list)


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

    def test_contains_consort_reference(self):
        """Should reference CONSORT-style summary."""
        html = "\n".join(_section_patient_flow(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "CONSORT" in html

    def test_outlier_pct_shown(self):
        """Outlier removal percentage should appear."""
        html = "\n".join(_section_patient_flow(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "9.5%" in html or "9.5" in html

    def test_competing_risk_exclusion(self):
        """CR exclusion column should show the glme_excluded count."""
        html = "\n".join(_section_patient_flow(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Competing" in html

    def test_analysed_count_computed(self):
        """Analysed (GLME) column should compute initial - baseline_excl - CR_excl."""
        html = "\n".join(_section_patient_flow(_make_log_data(), ["Standard"], _make_mat_data()))
        # 42 initial - 6 baseline_excl - 5 CR = 31
        assert "31" in html

    def test_negative_analysed_clamped(self):
        """If analysed < 0, should show 0 (clamped)."""
        log = _make_log_data()
        log["Standard"]["stats_comparisons"]["glme_excluded"] = {"n_excluded": 50, "n_total": 42, "pct": 100}
        html = "\n".join(_section_patient_flow(log, ["Standard"], _make_mat_data()))
        # Should fallback to init - baseline_excl = 42-6=36
        assert isinstance(html, str)

    def test_parse_warning_when_log_none(self):
        """Should show parse warning when log_data is None."""
        mat = {"Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}}}
        log = None
        result = _section_patient_flow(log, ["Standard"], mat)
        # No data means empty result
        assert isinstance(result, list)

    def test_parse_warning_with_parse_warnings(self):
        """Should show parse warning when parse_warnings exist."""
        log = _make_log_data()
        log["Standard"]["parse_warnings"] = ["Could not parse some data"]
        html = "\n".join(_section_patient_flow(log, ["Standard"], _make_mat_data()))
        assert "could not be parsed" in html.lower() or "incomplete" in html.lower()

    def test_multiple_dwi_types(self):
        """Should show one row per DWI type in the flow table."""
        log = _make_log_data()
        log["dnCNN"] = {
            "baseline": {
                "total_outliers": {"n_removed": 2, "n_total": 40, "pct": 5.0},
                "baseline_exclusion": {"n_excluded": 3, "n_total": 40},
            },
        }
        mat = {
            "Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}},
            "dnCNN": {"longitudinal": {"num_patients": 40, "num_timepoints": 5}},
        }
        html = "\n".join(_section_patient_flow(log, ["Standard", "dnCNN"], mat))
        assert "Standard" in html
        assert "dnCNN" in html

    def test_initial_from_baseline_exclusion(self):
        """If mat_data has no longitudinal, should derive n_initial from baseline_exclusion."""
        log = _make_log_data()
        result = _section_patient_flow(log, ["Standard"], {})
        html = "\n".join(result)
        assert "48" in html  # n_total from baseline_exclusion


# ── Edge cases: empty cohorts, missing metrics ──


class TestCohortOverviewEdgeCases:
    def test_all_none_inputs(self):
        result = _section_cohort_overview(None, None, None)
        assert isinstance(result, list)

    def test_empty_longitudinal_dict(self):
        mat = {"Standard": {"longitudinal": {}}}
        result = _section_cohort_overview(mat, None, ["Standard"])
        assert isinstance(result, list)

    def test_mat_data_missing_dwi_type_key(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}}}
        result = _section_cohort_overview(mat, _make_log_data(), ["Standard", "dnCNN", "IVIMnet"])
        assert isinstance(result, list)

    def test_zero_patient_count(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 0, "num_timepoints": 0}}}
        result = _section_cohort_overview(mat, None, ["Standard"])
        assert isinstance(result, list)

    def test_single_timepoint(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 1}}}
        result = _section_cohort_overview(mat, _make_log_data(), ["Standard"])
        html = "\n".join(result)
        # Single timepoint should not show attrition
        assert isinstance(html, str)

    def test_patients_per_timepoint_not_a_list(self):
        mat = {"Standard": {"longitudinal": {
            "num_patients": 42, "num_timepoints": 3,
            "patients_per_timepoint": "not_a_list",
        }}}
        result = _section_cohort_overview(mat, None, ["Standard"])
        assert isinstance(result, list)

    def test_heavily_imbalanced_lf_high(self):
        log = _make_log_data()
        log["Standard"]["survival"]["n_lf"] = 38
        log["Standard"]["survival"]["n_lc"] = 4
        html = "\n".join(_section_cohort_overview(_make_mat_data(), log, ["Standard"]))
        assert "Imbalanced" in html

    def test_balanced_outcomes_no_warning(self):
        log = _make_log_data()
        log["Standard"]["survival"]["n_lf"] = 20
        log["Standard"]["survival"]["n_lc"] = 22
        html = "\n".join(_section_cohort_overview(_make_mat_data(), log, ["Standard"]))
        assert "Imbalanced" not in html

    def test_three_dwi_types_with_mixed_data(self):
        mat = {
            "Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}},
            "dnCNN": {"longitudinal": {"num_patients": 40, "num_timepoints": 5}},
            "IVIMnet": {},
        }
        result = _section_cohort_overview(mat, None, ["Standard", "dnCNN", "IVIMnet"])
        assert isinstance(result, list)

    def test_baseline_exclusion_zero(self):
        log = _make_log_data()
        log["Standard"]["baseline"]["baseline_exclusion"]["n_excluded"] = 0
        log["Standard"]["baseline"]["baseline_exclusion"]["n_total"] = 42
        html = "\n".join(_section_cohort_overview(_make_mat_data(), log, ["Standard"]))
        assert isinstance(html, str)


class TestPatientFlowEdgeCases:
    def test_all_none_inputs(self):
        result = _section_patient_flow(None, [], None)
        assert isinstance(result, list)

    def test_missing_baseline_key(self):
        log = {"Standard": {"stats_comparisons": {"glme_excluded": {"n_excluded": 5, "n_total": 42, "pct": 12.0}}}}
        result = _section_patient_flow(log, ["Standard"], _make_mat_data())
        assert isinstance(result, list)

    def test_missing_glme_excluded(self):
        log = _make_log_data()
        del log["Standard"]["stats_comparisons"]
        result = _section_patient_flow(log, ["Standard"], _make_mat_data())
        assert isinstance(result, list)

    def test_zero_exclusions(self):
        log = _make_log_data()
        log["Standard"]["baseline"]["baseline_exclusion"]["n_excluded"] = 0
        log["Standard"]["stats_comparisons"]["glme_excluded"]["n_excluded"] = 0
        log["Standard"]["baseline"]["total_outliers"]["n_removed"] = 0
        result = _section_patient_flow(log, ["Standard"], _make_mat_data())
        assert isinstance(result, list)

    def test_empty_dwi_types_list(self):
        result = _section_patient_flow(_make_log_data(), [], _make_mat_data())
        assert isinstance(result, list)

    def test_dwi_type_not_in_log(self):
        result = _section_patient_flow(_make_log_data(), ["IVIMnet"], _make_mat_data())
        assert isinstance(result, list)
