"""Tests for the sub-volume stability report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.subvolume_stability import _section_subvolume_stability


def _make_stab_data():
    """Build synthetic subvolume_stability mat_data entry."""
    return {
        "method_names": ["adc_threshold", "otsu", "kmeans"],
        "n_patients": 5,
        "n_timepoints": 3,
        "dice_vs_baseline": [
            # adc_threshold: stable
            [[1.0, 1.0, 1.0, 1.0, 1.0], [0.85, 0.80, 0.82, 0.78, 0.81], [0.75, 0.70, 0.72, 0.68, 0.71]],
            # otsu: moderate decline
            [[1.0, 1.0, 1.0, 1.0, 1.0], [0.60, 0.55, 0.58, 0.52, 0.57], [0.45, 0.40, 0.42, 0.38, 0.41]],
            # kmeans: unstable
            [[1.0, 1.0, 1.0, 1.0, 1.0], [0.35, 0.30, 0.32, 0.28, 0.31], [0.20, 0.15, 0.18, 0.12, 0.16]],
        ],
    }


def _make_mat_data(with_stab=True):
    data: dict = {"Standard": {"core_method": {}, "dosimetry": {}, "longitudinal": {}}}
    if with_stab:
        data["Standard"]["subvolume_stability"] = _make_stab_data()
    return data


class TestSubvolumeStabilitySection:

    def test_with_data_returns_html(self):
        mat_data = _make_mat_data(with_stab=True)
        result = _section_subvolume_stability(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Stability" in html
        assert "<table>" in html
        assert "adc_threshold" in html

    def test_fx1_column_present(self):
        mat_data = _make_mat_data(with_stab=True)
        result = _section_subvolume_stability(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Fx1" in html

    def test_color_coding(self):
        mat_data = _make_mat_data(with_stab=True)
        result = _section_subvolume_stability(mat_data, ["Standard"])
        html = "\n".join(result)
        # Green for high Dice values
        assert "#d4edda" in html
        # Red for low Dice values
        assert "#f8d7da" in html

    def test_empty_data_returns_not_available(self):
        mat_data = _make_mat_data(with_stab=False)
        result = _section_subvolume_stability(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_none_mat_data(self):
        result = _section_subvolume_stability(None, [])
        html = "\n".join(result)
        assert "not available" in html.lower()

    def test_all_methods_in_table(self):
        mat_data = _make_mat_data(with_stab=True)
        result = _section_subvolume_stability(mat_data, ["Standard"])
        html = "\n".join(result)
        for method in ["adc_threshold", "otsu", "kmeans"]:
            assert method in html
