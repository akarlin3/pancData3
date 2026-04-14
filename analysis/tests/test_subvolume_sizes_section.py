"""Tests for the sub-volume sizes report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.subvolume_sizes import _section_subvolume_sizes


def _make_subvol(lf_frac_pct=None, wilcoxon_p=None):
    return {
        "timepoints": [1, 2, 3],
        "lc": {
            "mean_vol_cm3": [5.0, 4.2, 3.5],
            "std_vol_cm3": [1.0, 0.9, 0.8],
            "mean_frac_pct": [0.18, 0.16, 0.12],
        },
        "lf": {
            "mean_vol_cm3": [6.0, 5.5, 5.0],
            "std_vol_cm3": [1.2, 1.1, 1.0],
            "mean_frac_pct": lf_frac_pct if lf_frac_pct is not None else [0.22, 0.20, 0.18],
        },
        "wilcoxon_p": wilcoxon_p if wilcoxon_p is not None else [0.03, 0.12, 0.08],
        "median_voxel_count": [5.2, 4.6, 3.8],
    }


def _make_mat_data(**kw):
    return {"Standard": {"subvolume_sizes": _make_subvol(**kw)}}


class TestSubvolumeSizesSection:

    def test_with_data_returns_html(self):
        result = _section_subvolume_sizes(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "Sub-Volume Sizes" in html
        assert "<table>" in html
        assert "Fx1" in html

    def test_empty_returns_empty(self):
        assert _section_subvolume_sizes({}, []) == []

    def test_none_returns_empty(self):
        assert _section_subvolume_sizes(None, []) == []

    def test_no_subvol_data_returns_empty(self):
        result = _section_subvolume_sizes({"Standard": {}}, ["Standard"])
        assert result == []

    def test_p_value_significant_marked(self):
        result = _section_subvolume_sizes(
            _make_mat_data(wilcoxon_p=[0.001, 0.5, 0.5]),
            ["Standard"],
        )
        html = "\n".join(result)
        # sig class should appear for p=0.001
        assert "sig" in html

    def test_scipy_unavailable_shows_em_dash(self):
        result = _section_subvolume_sizes(
            _make_mat_data(wilcoxon_p=[None, None, None]),
            ["Standard"],
        )
        html = "\n".join(result)
        # em-dash should appear for missing p-values
        assert "&mdash;" in html

    def test_small_fraction_warning(self):
        result = _section_subvolume_sizes(
            _make_mat_data(lf_frac_pct=[0.02, 0.02, 0.02]),  # 2% < 5%
            ["Standard"],
        )
        html = "\n".join(result)
        assert "Sub-volumes" in html and "5%" in html

    def test_no_small_fraction_warning(self):
        result = _section_subvolume_sizes(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "may lack clinical significance" not in html

    def test_percent_scale_fraction(self):
        # Accept fractions stored as 0-1 and convert.
        result = _section_subvolume_sizes(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        # 0.18 fraction -> ~18%
        assert "18" in html or "18.0%" in html
