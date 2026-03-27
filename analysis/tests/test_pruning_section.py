"""Tests for the core method pruning results report section builder."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.pruning_results import _section_pruning_results


def _make_pruning_data(pruned=True):
    """Build synthetic pruning mat_data entry."""
    if pruned:
        return {
            "active_methods": [
                "adc_threshold", "d_threshold", "df_intersection",
                "otsu", "kmeans", "percentile",
            ],
            "pruned_info": [
                {"name": "gmm", "reason": "failure_rate",
                 "failure_rate": 0.32, "pipeline": "IVIMNet"},
                {"name": "active_contours", "reason": "failure_rate",
                 "failure_rate": 0.55, "pipeline": "IVIMNet"},
                {"name": "region_growing", "reason": "failure_rate",
                 "failure_rate": 0.22, "pipeline": "DnCNN"},
                {"name": "spectral", "reason": "insufficient_voxels",
                 "failure_rate": 0.42, "pipeline": "IVIMNet"},
                {"name": "fdm", "reason": "manual_exclusion",
                 "failure_rate": 0.12, "pipeline": "Standard"},
            ],
        }
    else:
        # No pruning — all methods retained
        return {
            "active_methods": [
                "adc_threshold", "d_threshold", "df_intersection",
                "otsu", "gmm", "kmeans", "region_growing", "active_contours",
                "percentile", "spectral", "fdm",
            ],
            "pruned_info": [],
        }


def _make_mat_data(with_pruning=True, pruned=True):
    """Build mat_data dict with optional pruning data."""
    data: dict = {
        "Standard": {
            "core_method": {},
            "dosimetry": {},
            "longitudinal": {},
        }
    }
    if with_pruning:
        data["Standard"]["pruning"] = _make_pruning_data(pruned=pruned)
    return data


class TestPruningSection:
    """Tests for _section_pruning_results."""

    def test_with_pruned_methods(self):
        """Section should show pruned methods table when data is present."""
        mat_data = _make_mat_data(with_pruning=True, pruned=True)
        result = _section_pruning_results(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Pruning" in html
        assert "<table>" in html
        assert "gmm" in html
        assert "active_contours" in html
        assert "Failure Rate" in html

    def test_no_pruning_data(self):
        """Missing pruning data should produce 'not available' message."""
        mat_data = _make_mat_data(with_pruning=False)
        result = _section_pruning_results(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "No pruning data available" in html

    def test_none_mat_data(self):
        """None mat_data should not crash."""
        result = _section_pruning_results(None, [])
        html = "\n".join(result)
        assert "No pruning data available" in html

    def test_no_methods_pruned(self):
        """When no methods are pruned, show 'no methods were pruned'."""
        mat_data = _make_mat_data(with_pruning=True, pruned=False)
        result = _section_pruning_results(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "No methods were pruned" in html

    def test_retained_methods_shown(self):
        """Retained methods should be listed."""
        mat_data = _make_mat_data(with_pruning=True, pruned=True)
        result = _section_pruning_results(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "adc_threshold" in html
        assert "Retained" in html

    def test_pruned_reasons_shown(self):
        """All three reason types should appear in the table."""
        mat_data = _make_mat_data(with_pruning=True, pruned=True)
        result = _section_pruning_results(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "Failure Rate" in html
        assert "Manual Exclusion" in html
        assert "Insufficient Voxels" in html

    def test_all_pruned_except_adc(self):
        """Edge case: only adc_threshold retained."""
        data = _make_pruning_data(pruned=True)
        data["active_methods"] = ["adc_threshold"]
        data["pruned_info"] = [
            {"name": m, "reason": "failure_rate", "failure_rate": 0.99, "pipeline": "Standard"}
            for m in ["d_threshold", "df_intersection", "otsu", "gmm", "kmeans",
                       "region_growing", "active_contours", "percentile", "spectral", "fdm"]
        ]
        mat_data = {"Standard": {"pruning": data}}
        result = _section_pruning_results(mat_data, ["Standard"])
        html = "\n".join(result)
        assert "adc_threshold" in html
        assert len(data["pruned_info"]) == 10
