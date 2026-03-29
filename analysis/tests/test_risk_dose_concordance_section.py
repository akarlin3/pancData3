"""Tests for the risk-dose concordance report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.risk_dose_concordance import _section_risk_dose_concordance


def _make_concordance_data():
    """Build synthetic risk_dose_concordance mat_data entry."""
    return {
        "method_results": [
            {
                "method_name": "adc_threshold",
                "best_dose_metric": "d95_adc_sub",
                "cohen_kappa": 0.55,
                "concordance_pct": 72.0,
                "n_complementary": 6,
                "combined_auc": 0.85,
            },
            {
                "method_name": "otsu",
                "best_dose_metric": "d95_adc_sub",
                "cohen_kappa": 0.15,
                "concordance_pct": 48.0,
                "n_complementary": 12,
                "combined_auc": 0.78,
            },
        ],
        "summary": "Models show moderate agreement (mean kappa=0.35). Combining both may improve prediction.",
    }


def _make_mat_data(with_rdc=True):
    data: dict = {"Standard": {"core_method": {}, "dosimetry": {}, "longitudinal": {}}}
    if with_rdc:
        data["Standard"]["risk_dose_concordance"] = _make_concordance_data()
    return data


class TestRiskDoseConcordanceSection:

    def test_with_data_returns_html(self):
        mat_data = _make_mat_data(with_rdc=True)
        result = _section_risk_dose_concordance(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Concordance" in html
        assert "<table>" in html
        assert "adc_threshold" in html

    def test_kappa_interpretation_shown(self):
        mat_data = _make_mat_data(with_rdc=True)
        result = _section_risk_dose_concordance(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "moderate" in html  # kappa=0.55 -> moderate
        assert "poor" in html  # kappa=0.15 -> poor

    def test_concordance_pct_shown(self):
        mat_data = _make_mat_data(with_rdc=True)
        result = _section_risk_dose_concordance(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "72.0%" in html

    def test_combined_auc_shown(self):
        mat_data = _make_mat_data(with_rdc=True)
        result = _section_risk_dose_concordance(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "0.85" in html

    def test_summary_shown(self):
        mat_data = _make_mat_data(with_rdc=True)
        result = _section_risk_dose_concordance(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "moderate agreement" in html

    def test_empty_data_returns_not_available(self):
        mat_data = _make_mat_data(with_rdc=False)
        result = _section_risk_dose_concordance(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_none_mat_data(self):
        result = _section_risk_dose_concordance(None, [])
        html = "\n".join(result)
        assert "not available" in html.lower()
