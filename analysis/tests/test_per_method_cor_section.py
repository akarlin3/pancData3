"""Tests for the per-method CoR report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.per_method_cor import _section_per_method_cor


def _make_cor_data():
    """Build synthetic per_method_cor mat_data entry."""
    return {
        "method_names": ["adc_threshold", "otsu", "kmeans", "spectral"],
        "median_wcv": [0.05, 0.12, 0.20, 0.35],
        "cor": [13.9, 33.3, 55.4, 97.0],
        "n_patients_with_repeats": 8,
    }


def _make_mat_data(with_cor=True):
    data: dict = {"Standard": {"core_method": {}, "dosimetry": {}, "longitudinal": {}}}
    if with_cor:
        data["Standard"]["per_method_cor"] = _make_cor_data()
    return data


class TestPerMethodCorSection:

    def test_with_data_returns_html(self):
        mat_data = _make_mat_data(with_cor=True)
        result = _section_per_method_cor(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Reproducibility" in html
        assert "<table>" in html
        assert "adc_threshold" in html

    def test_color_coding_green(self):
        mat_data = _make_mat_data(with_cor=True)
        result = _section_per_method_cor(mat_data, ["Standard"])
        html = "\n".join(result)
        # Green for CoR < 15%
        assert "#d4edda" in html

    def test_color_coding_red(self):
        mat_data = _make_mat_data(with_cor=True)
        result = _section_per_method_cor(mat_data, ["Standard"])
        html = "\n".join(result)
        # Red for CoR > 30%
        assert "#f8d7da" in html

    def test_empty_data_returns_not_available(self):
        mat_data = _make_mat_data(with_cor=False)
        result = _section_per_method_cor(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_none_mat_data(self):
        result = _section_per_method_cor(None, [])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_patient_count_shown(self):
        mat_data = _make_mat_data(with_cor=True)
        result = _section_per_method_cor(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "8" in html
