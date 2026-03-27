"""Tests for the core method outcomes report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.core_method_outcomes import _section_core_method_outcomes


def _make_cmo_data(significant=True):
    """Build synthetic core_method_outcomes mat_data entry."""
    methods = []

    # Method 1: significant
    m1 = {
        "method_name": "otsu",
        "n_patients": 20,
        "n_events": 8,
        "univariable": [
            {"metric_name": "d95_adc_sub", "hr": 0.72, "hr_ci": [0.54, 0.96],
             "p_value": 0.024 if significant else 0.35, "n": 20, "n_events": 8},
            {"metric_name": "v50_adc_sub", "hr": 0.88, "hr_ci": [0.65, 1.19],
             "p_value": 0.41, "n": 20, "n_events": 8},
        ],
        "km": {
            "best_metric": "d95_adc_sub",
            "logrank_p": 0.018 if significant else 0.42,
            "logrank_chi2": 5.6,
            "median_high": 350,
            "median_low": 90,
        },
    }
    methods.append(m1)

    # Method 2: not significant
    m2 = {
        "method_name": "adc_threshold",
        "n_patients": 20,
        "n_events": 8,
        "univariable": [
            {"metric_name": "d95_adc_sub", "hr": 0.85, "hr_ci": [0.63, 1.15],
             "p_value": 0.291, "n": 20, "n_events": 8},
            {"metric_name": "v50_adc_sub", "hr": 0.91, "hr_ci": [0.68, 1.22],
             "p_value": 0.53, "n": 20, "n_events": 8},
        ],
        "km": {
            "best_metric": "d95_adc_sub",
            "logrank_p": 0.31,
            "logrank_chi2": 1.02,
            "median_high": 280,
            "median_low": 120,
        },
    }
    methods.append(m2)

    ranking = ["otsu", "adc_threshold"] if significant else ["adc_threshold", "otsu"]

    return {
        "method_results": methods,
        "ranking": ranking,
        "active_methods": ["otsu", "adc_threshold"],
    }


def _make_mat_data(with_cmo=True, significant=True):
    """Build mat_data dict with optional core_method_outcomes."""
    data: dict = {
        "Standard": {
            "core_method": {},
            "dosimetry": {},
            "longitudinal": {},
        }
    }
    if with_cmo:
        data["Standard"]["core_method_outcomes"] = _make_cmo_data(significant=significant)
    return data


class TestCoreMethodOutcomesSection:
    """Tests for _section_core_method_outcomes."""

    def test_with_significant_results(self):
        """Section should show ranking table with significant results highlighted."""
        mat_data = _make_mat_data(with_cmo=True, significant=True)
        result = _section_core_method_outcomes(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Core Method Outcome" in html
        assert "<table>" in html
        assert "otsu" in html
        assert "sig-moderate" in html  # significant row highlighted
        assert "Interpretation" in html

    def test_no_significant_results(self):
        """When no method is significant, show appropriate message."""
        mat_data = _make_mat_data(with_cmo=True, significant=False)
        result = _section_core_method_outcomes(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "not significantly predicted" in html.lower() or "No core method" in html

    def test_empty_data(self):
        """Missing data should produce 'not available' message."""
        mat_data = _make_mat_data(with_cmo=False)
        result = _section_core_method_outcomes(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_none_mat_data(self):
        """None mat_data should not crash."""
        result = _section_core_method_outcomes(None, [])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_ranking_table_structure(self):
        """Table should have Rank, Method, Best Metric, HR, p-value, Events columns."""
        mat_data = _make_mat_data(with_cmo=True, significant=True)
        result = _section_core_method_outcomes(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Rank" in html
        assert "Method" in html
        assert "Best Metric" in html
        assert "HR" in html
        assert "p-value" in html
        assert "Events" in html

    def test_hr_ci_shown(self):
        """HR confidence intervals should be shown."""
        mat_data = _make_mat_data(with_cmo=True, significant=True)
        result = _section_core_method_outcomes(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "0.72" in html  # HR value
        assert "0.54" in html  # CI lower

    def test_both_methods_in_table(self):
        """Both methods should appear in the table."""
        mat_data = _make_mat_data(with_cmo=True, significant=True)
        result = _section_core_method_outcomes(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "otsu" in html
        assert "adc_threshold" in html
