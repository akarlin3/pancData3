"""Tests for the sub-volume dose context report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.dose_context import _section_dose_context


def _make_dosi(whole_d95=55.0, whole_dmean=54.0, whole_v50=0.98,
               sub_d95=48.0, sub_v50=0.80):
    """Build a dosimetry dict with whole-GTV and sub-volume fields."""
    return {
        "d95_adc_mean": {"mean": sub_d95, "std": 5.0},
        "v50_adc_mean": {"mean": sub_v50, "std": 0.1},
        "whole_gtv_d95_mean": whole_d95,
        "whole_gtv_d95_std": 3.0,
        "whole_gtv_v50_mean": whole_v50,
        "whole_gtv_v50_std": 0.02,
        "whole_gtv_dmean_mean": whole_dmean,
        "whole_gtv_dmean_std": 2.0,
    }


def _make_mat_data(**kw):
    return {"Standard": {"dosimetry": _make_dosi(**kw)}}


class TestDoseContextSection:

    def test_with_large_deficit_shows_warning(self):
        # Whole D95 = 55, sub D95 = 48 -> deficit = 7 > 5 (flagged)
        result = _section_dose_context(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "Sub-Volume Dose Context" in html
        assert "Clinically significant dose deficit" in html

    def test_no_deficit_no_warning(self):
        # Whole D95 = 50, sub D95 = 49 -> deficit = 1 < 5 (not flagged)
        # Whole V50 98%, sub V50 95% -> deficit = 3 < 10 (not flagged)
        result = _section_dose_context(
            _make_mat_data(whole_d95=50.0, sub_d95=49.0,
                            whole_v50=0.98, sub_v50=0.95),
            ["Standard"],
        )
        html = "\n".join(result)
        assert "Sub-Volume Dose Context" in html
        assert "Clinically significant dose deficit" not in html

    def test_missing_data_returns_empty(self):
        # No dosimetry fields.
        mat = {"Standard": {"dosimetry": {}}}
        result = _section_dose_context(mat, ["Standard"])
        assert result == []

    def test_none_mat_data_returns_empty(self):
        assert _section_dose_context(None, []) == []

    def test_empty_mat_data_returns_empty(self):
        assert _section_dose_context({}, []) == []

    def test_only_whole_gtv_no_sub_returns_empty(self):
        mat = {"Standard": {"dosimetry": {
            "whole_gtv_d95_mean": 55.0,
            "whole_gtv_v50_mean": 0.98,
        }}}
        result = _section_dose_context(mat, ["Standard"])
        assert result == []

    def test_only_sub_no_whole_returns_empty(self):
        mat = {"Standard": {"dosimetry": {
            "d95_adc_mean": {"mean": 50.0, "std": 0.0},
            "v50_adc_mean": {"mean": 0.9, "std": 0.0},
        }}}
        result = _section_dose_context(mat, ["Standard"])
        assert result == []

    def test_table_shows_metric_rows(self):
        result = _section_dose_context(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "D95" in html
        assert "V50" in html
        assert "Dmean" in html

    def test_large_v50_deficit_flagged(self):
        # Whole V50 = 98%, sub V50 = 80% -> deficit 18% > 10%
        result = _section_dose_context(
            _make_mat_data(whole_v50=0.98, sub_v50=0.80),
            ["Standard"],
        )
        html = "\n".join(result)
        assert "Clinically significant dose deficit" in html
