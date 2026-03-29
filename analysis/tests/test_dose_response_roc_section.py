"""Tests for the dose-response ROC report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.dose_response_roc import _section_dose_response_roc


def _make_roc_data():
    """Build synthetic dose_response_roc mat_data entry."""
    return {
        "method_results": [
            {
                "method_name": "adc_threshold",
                "best_metric": "d95_adc_sub",
                "best_auc": 0.82,
                "metrics": [
                    {
                        "metric_name": "d95_adc_sub",
                        "auc": 0.82,
                        "auc_ci": [0.71, 0.93],
                        "optimal_threshold": 42.5,
                        "sensitivity": 0.85,
                        "specificity": 0.78,
                    },
                ],
            },
            {
                "method_name": "otsu",
                "best_metric": "d95_adc_sub",
                "best_auc": 0.65,
                "metrics": [
                    {
                        "metric_name": "d95_adc_sub",
                        "auc": 0.65,
                        "auc_ci": [0.50, 0.80],
                        "optimal_threshold": 38.0,
                        "sensitivity": 0.70,
                        "specificity": 0.60,
                    },
                ],
            },
        ],
        "ranking": ["adc_threshold", "otsu"],
    }


def _make_mat_data(with_roc=True):
    data: dict = {"Standard": {"core_method": {}, "dosimetry": {}, "longitudinal": {}}}
    if with_roc:
        data["Standard"]["dose_response_roc"] = _make_roc_data()
    return data


class TestDoseResponseRocSection:

    def test_with_data_returns_html(self):
        mat_data = _make_mat_data(with_roc=True)
        result = _section_dose_response_roc(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "ROC" in html
        assert "<table>" in html
        assert "adc_threshold" in html

    def test_auc_values_shown(self):
        mat_data = _make_mat_data(with_roc=True)
        result = _section_dose_response_roc(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "0.82" in html
        assert "42.5" in html

    def test_ci_shown(self):
        mat_data = _make_mat_data(with_roc=True)
        result = _section_dose_response_roc(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "0.71" in html

    def test_clinical_guidance_for_high_auc(self):
        mat_data = _make_mat_data(with_roc=True)
        result = _section_dose_response_roc(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Clinical guidance" in html

    def test_high_auc_highlighted(self):
        mat_data = _make_mat_data(with_roc=True)
        result = _section_dose_response_roc(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "#d4edda" in html

    def test_empty_data_returns_not_available(self):
        mat_data = _make_mat_data(with_roc=False)
        result = _section_dose_response_roc(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_none_mat_data(self):
        result = _section_dose_response_roc(None, [])
        html = "\n".join(result)
        assert "not available" in html.lower()
