"""Tests for the baseline vs delta report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.baseline_vs_delta import _section_baseline_vs_delta


def _make_rows():
    return [
        {
            "parameter": "ADC",
            "timepoint": "Fx2",
            "baseline_hr": 1.4,
            "baseline_p": 0.08,
            "baseline_cindex": 0.62,
            "delta_hr": 1.8,
            "delta_p": 0.01,
            "delta_cindex": 0.71,
            "better_predictor": "Delta",
            "n_events": 14,
        },
        {
            "parameter": "D",
            "timepoint": "Fx3",
            "baseline_hr": 1.6,
            "baseline_p": 0.02,
            "baseline_cindex": 0.70,
            "delta_hr": 1.3,
            "delta_p": 0.15,
            "delta_cindex": 0.64,
            "better_predictor": "Baseline",
            "n_events": 14,
        },
    ]


def _make_mat_data(rows=None):
    return {"Standard": {"baseline_vs_delta": rows if rows is not None else _make_rows()}}


class TestBaselineVsDeltaSection:

    def test_with_data_returns_html(self):
        result = _section_baseline_vs_delta(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "Baseline vs Delta" in html
        assert "<table>" in html
        assert "ADC" in html
        assert "Fx2" in html

    def test_winner_highlighted(self):
        result = _section_baseline_vs_delta(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "Delta" in html
        assert "Baseline" in html
        # Green background for winner column
        assert "#d4edda" in html

    def test_better_c_bold(self):
        result = _section_baseline_vs_delta(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        # Winning C-indices bolded
        assert "<strong>0.710</strong>" in html  # delta winner
        assert "<strong>0.700</strong>" in html  # baseline winner

    def test_empty_returns_empty(self):
        assert _section_baseline_vs_delta({}, []) == []

    def test_none_returns_empty(self):
        assert _section_baseline_vs_delta(None, []) == []

    def test_missing_field_returns_empty(self):
        mat = {"Standard": {"dosimetry": {}}}
        assert _section_baseline_vs_delta(mat, ["Standard"]) == []

    def test_empty_rows_returns_empty(self):
        mat = {"Standard": {"baseline_vs_delta": []}}
        assert _section_baseline_vs_delta(mat, ["Standard"]) == []

    def test_p_values_formatted(self):
        result = _section_baseline_vs_delta(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "p=" in html

    def test_sig_class_on_significant_p(self):
        result = _section_baseline_vs_delta(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        # delta_p=0.01 (significant) -> should have sig-related class
        assert "sig" in html

    def test_interpretation_text(self):
        result = _section_baseline_vs_delta(_make_mat_data(), ["Standard"])
        html = "\n".join(result)
        assert "Interpretation" in html
