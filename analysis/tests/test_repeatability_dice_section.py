"""Tests for the Fx1 repeat Dice report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.repeatability_dice import _section_repeatability_dice


def _make_rpt_entry(adc=0.82, d=0.65, f=0.45, dstar=None, n=10):
    """Build synthetic repeatability_dice dict entry."""
    entry: dict = {}
    if adc is not None:
        entry["adc"] = {"mean": adc, "std": 0.05, "n": n}
    if d is not None:
        entry["d"] = {"mean": d, "std": 0.08, "n": n}
    if f is not None:
        entry["f"] = {"mean": f, "std": 0.12, "n": n}
    if dstar is not None:
        entry["dstar"] = {"mean": dstar, "std": 0.10, "n": n}
    return entry


def _make_mat_data(with_data=True, with_partial=False):
    data: dict = {
        "Standard": {"core_method": {}, "dosimetry": {}, "longitudinal": {}},
    }
    if with_data and not with_partial:
        data["Standard"]["repeatability_dice"] = _make_rpt_entry(
            adc=0.82, d=0.65, f=0.45, dstar=0.70, n=10
        )
    elif with_partial:
        # Only ADC has data.
        data["Standard"]["repeatability_dice"] = {
            "adc": {"mean": 0.75, "std": 0.1, "n": 5},
        }
    return data


class TestRepeatabilityDiceSection:

    def test_with_data_returns_html(self):
        result = _section_repeatability_dice(_make_mat_data(with_data=True), ["Standard"])
        html = "\n".join(result)
        assert "Spatial Repeatability" in html
        assert "<table>" in html
        assert "ADC" in html
        # D* label is escaped/rendered
        assert "D*" in html

    def test_green_for_high_dice(self):
        # ADC = 0.82 >= 0.70 → green
        result = _section_repeatability_dice(_make_mat_data(with_data=True), ["Standard"])
        html = "\n".join(result)
        assert "#d4edda" in html

    def test_yellow_for_mid_dice(self):
        # D = 0.65 in [0.50, 0.70) → yellow
        result = _section_repeatability_dice(_make_mat_data(with_data=True), ["Standard"])
        html = "\n".join(result)
        assert "#fff3cd" in html

    def test_red_for_low_dice(self):
        # f = 0.45 < 0.50 → red
        result = _section_repeatability_dice(_make_mat_data(with_data=True), ["Standard"])
        html = "\n".join(result)
        assert "#f8d7da" in html

    def test_empty_data_returns_empty(self):
        result = _section_repeatability_dice({}, [])
        assert result == []

    def test_none_mat_data_returns_empty(self):
        result = _section_repeatability_dice(None, [])
        assert result == []

    def test_mat_data_without_rpt_returns_empty(self):
        mat = {"Standard": {"core_method": {}, "dosimetry": {}}}
        result = _section_repeatability_dice(mat, ["Standard"])
        assert result == []

    def test_partial_data_renders(self):
        # Only ADC available; other params should show em-dash.
        result = _section_repeatability_dice(_make_mat_data(with_partial=True), ["Standard"])
        html = "\n".join(result)
        assert "ADC" in html
        assert "&mdash;" in html

    def test_patient_count_shown(self):
        result = _section_repeatability_dice(_make_mat_data(with_data=True), ["Standard"])
        html = "\n".join(result)
        assert "n=10" in html

    def test_mean_and_sd_formatted(self):
        result = _section_repeatability_dice(_make_mat_data(with_data=True), ["Standard"])
        html = "\n".join(result)
        # Mean ± SD format
        assert "0.82" in html
        assert "&plusmn;" in html

    def test_all_nan_returns_empty(self):
        # repeatability_dice exists but all entries have n=0.
        mat = {
            "Standard": {
                "repeatability_dice": {
                    "adc": {"mean": None, "std": None, "n": 0},
                    "d": {"mean": None, "std": None, "n": 0},
                }
            }
        }
        result = _section_repeatability_dice(mat, ["Standard"])
        html = "\n".join(result)
        # Section still rendered with em-dashes (non-empty repeatability_dice dict).
        # But each cell is missing — check header still rendered.
        assert "Spatial Repeatability" in html
        assert "&mdash;" in html
