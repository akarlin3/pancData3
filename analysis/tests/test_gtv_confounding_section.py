"""Tests for the GTV confounding report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.gtv_confounding import _section_gtv_confounding


def _make_gtv_data():
    """Build synthetic gtv_confounding mat_data entry."""
    return {
        "method_results": [
            {
                "method_name": "adc_threshold",
                "d95_gtv_correlation": -0.45,
                "d95_gtv_pvalue": 0.02,
                "unadjusted_hr": 0.65,
                "adjusted_hr": 0.78,
                "confounding_flag": True,
            },
            {
                "method_name": "otsu",
                "d95_gtv_correlation": -0.12,
                "d95_gtv_pvalue": 0.55,
                "unadjusted_hr": 0.70,
                "adjusted_hr": 0.72,
                "confounding_flag": False,
            },
        ],
        "summary": "GTV volume is a significant confounder for at least one method.",
    }


def _make_mat_data(with_gtv=True):
    data: dict = {"Standard": {"core_method": {}, "dosimetry": {}, "longitudinal": {}}}
    if with_gtv:
        data["Standard"]["gtv_confounding"] = _make_gtv_data()
    return data


class TestGtvConfoundingSection:

    def test_with_data_returns_html(self):
        mat_data = _make_mat_data(with_gtv=True)
        result = _section_gtv_confounding(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Confounding" in html
        assert "<table>" in html
        assert "adc_threshold" in html

    def test_confounding_flag_shown(self):
        mat_data = _make_mat_data(with_gtv=True)
        result = _section_gtv_confounding(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Yes" in html
        assert "sig-strong" in html

    def test_hr_values_shown(self):
        mat_data = _make_mat_data(with_gtv=True)
        result = _section_gtv_confounding(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "0.65" in html  # unadjusted HR
        assert "0.78" in html  # adjusted HR

    def test_warning_box_when_confounded(self):
        mat_data = _make_mat_data(with_gtv=True)
        result = _section_gtv_confounding(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "warn-box" in html

    def test_empty_data_returns_not_available(self):
        mat_data = _make_mat_data(with_gtv=False)
        result = _section_gtv_confounding(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_none_mat_data(self):
        result = _section_gtv_confounding(None, [])
        html = "\n".join(result)
        assert "not available" in html.lower()
