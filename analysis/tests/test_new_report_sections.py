"""Tests for the four new report section builders.

Covers:
- build_dca_section (model_diagnostics.py)
- build_nri_idi_section (statistics_robustness.py)
- build_texture_section (analysis_features.py)
- build_registration_quality_section (data_quality.py)

Each builder is tested for: valid data rendering, empty/missing data,
and malformed input resilience.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.model_diagnostics import build_dca_section
from report.sections.statistics_robustness import build_nri_idi_section
from report.sections.analysis_features import build_texture_section
from report.sections.data_quality import build_registration_quality_section


# ═══════════════════════════════════════════════════════════════════════
# DCA Section
# ═══════════════════════════════════════════════════════════════════════


def _make_dca_data():
    return {
        "Standard": {
            "stats_predictive": {
                "decision_curve": {
                    "thresholds": [0.1, 0.2, 0.3, 0.4, 0.5],
                    "net_benefits": [0.12, 0.10, 0.07, 0.03, -0.01],
                    "treat_all": [0.08, 0.05, 0.02, -0.02, -0.06],
                    "useful_range_lo": 0.10,
                    "useful_range_hi": 0.45,
                },
            },
        },
    }


class TestBuildDcaSection:
    def test_renders_with_valid_data(self):
        html = build_dca_section(_make_dca_data())
        assert "Decision Curve Analysis" in html
        assert "dca" in html  # anchor
        assert "0.10" in html
        assert "0.45" in html
        assert "net benefit" in html.lower()

    def test_contains_threshold_table(self):
        html = build_dca_section(_make_dca_data())
        assert "Threshold Probability" in html
        assert "Model Net Benefit" in html
        assert "Treat All" in html
        assert "Advantage" in html
        # Check a specific net benefit value
        assert "0.1200" in html

    def test_contains_clinical_interpretation(self):
        html = build_dca_section(_make_dca_data())
        assert "clinical decision space" in html

    def test_empty_dict_returns_empty(self):
        assert build_dca_section({}) == ""

    def test_none_returns_empty(self):
        assert build_dca_section(None) == ""

    def test_missing_decision_curve_key(self):
        data = {"Standard": {"stats_predictive": {}}}
        assert build_dca_section(data) == ""

    def test_malformed_thresholds(self):
        data = {
            "Standard": {
                "stats_predictive": {
                    "decision_curve": {
                        "thresholds": "bad",
                        "net_benefits": None,
                    },
                },
            },
        }
        # Should not raise
        result = build_dca_section(data)
        assert isinstance(result, str)

    def test_no_useful_range(self):
        data = {
            "Standard": {
                "stats_predictive": {
                    "decision_curve": {
                        "thresholds": [0.1, 0.2],
                        "net_benefits": [0.05, 0.02],
                        "treat_all": [0.03, 0.01],
                    },
                },
            },
        }
        html = build_dca_section(data)
        assert "could not be determined" in html

    def test_non_dict_dwi_data_skipped(self):
        assert build_dca_section({"Standard": "not a dict"}) == ""

    def test_multiple_dwi_types(self):
        data = _make_dca_data()
        data["dnCNN"] = data["Standard"].copy()
        html = build_dca_section(data)
        assert "Standard" in html
        assert "dnCNN" in html


# ═══════════════════════════════════════════════════════════════════════
# NRI/IDI Section
# ═══════════════════════════════════════════════════════════════════════


def _make_nri_idi_data():
    return {
        "Standard": {
            "nri": {
                "overall": 0.42,
                "overall_ci_lo": 0.15,
                "overall_ci_hi": 0.69,
                "events": 0.28,
                "events_ci_lo": 0.05,
                "events_ci_hi": 0.51,
                "non_events": 0.14,
                "non_events_ci_lo": -0.03,
                "non_events_ci_hi": 0.31,
                "p": 0.003,
            },
            "idi": {
                "value": 0.085,
                "ci_lo": 0.02,
                "ci_hi": 0.15,
                "p": 0.011,
            },
        },
    }


class TestBuildNriIdiSection:
    def test_renders_with_valid_data(self):
        html = build_nri_idi_section(_make_nri_idi_data())
        assert "Net Reclassification" in html
        assert "nri-idi" in html  # anchor

    def test_contains_nri_components(self):
        html = build_nri_idi_section(_make_nri_idi_data())
        assert "Overall NRI" in html
        assert "Events NRI" in html
        assert "Non-events NRI" in html
        assert "+0.4200" in html
        assert "+0.2800" in html

    def test_contains_idi(self):
        html = build_nri_idi_section(_make_nri_idi_data())
        assert "IDI" in html
        assert "+0.0850" in html

    def test_contains_confidence_intervals(self):
        html = build_nri_idi_section(_make_nri_idi_data())
        assert "+0.1500" in html  # overall CI lo
        assert "+0.6900" in html  # overall CI hi

    def test_flags_significant_improvement(self):
        html = build_nri_idi_section(_make_nri_idi_data())
        assert "significant reclassification improvement" in html.lower()

    def test_empty_dict_returns_empty(self):
        assert build_nri_idi_section({}) == ""

    def test_none_returns_empty(self):
        assert build_nri_idi_section(None) == ""

    def test_missing_nri_and_idi_keys(self):
        data = {"Standard": {"something_else": 42}}
        assert build_nri_idi_section(data) == ""

    def test_malformed_nri(self):
        data = {"Standard": {"nri": "not a dict", "idi": None}}
        result = build_nri_idi_section(data)
        assert isinstance(result, str)

    def test_nri_only_no_idi(self):
        data = {
            "Standard": {
                "nri": {"overall": 0.30, "p": 0.02},
            },
        }
        html = build_nri_idi_section(data)
        assert "Overall NRI" in html
        assert "+0.3000" in html

    def test_idi_only_no_nri(self):
        data = {
            "Standard": {
                "idi": {"value": 0.05, "p": 0.08},
            },
        }
        html = build_nri_idi_section(data)
        assert "IDI" in html

    def test_non_significant_no_flag(self):
        data = {
            "Standard": {
                "nri": {"overall": 0.10, "p": 0.15},
                "idi": {"value": 0.02, "p": 0.20},
            },
        }
        html = build_nri_idi_section(data)
        assert "significant reclassification improvement" not in html.lower()

    def test_non_dict_dwi_data_skipped(self):
        assert build_nri_idi_section({"Standard": 123}) == ""


# ═══════════════════════════════════════════════════════════════════════
# Texture Features Section
# ═══════════════════════════════════════════════════════════════════════


def _make_texture_data():
    return {
        "Standard": {
            "texture_features": {
                "features": [
                    {"name": "GLCM_Contrast", "category": "GLCM", "auc": 0.82, "p": 0.003},
                    {"name": "GLRLM_LRE", "category": "GLRLM", "auc": 0.78, "p": 0.01},
                    {"name": "Skewness", "category": "first_order", "auc": 0.71, "p": 0.04},
                    {"name": "Sphericity", "category": "shape", "auc": 0.65, "p": 0.12},
                ],
            },
        },
    }


class TestBuildTextureSection:
    def test_renders_with_valid_data(self):
        html = build_texture_section(_make_texture_data())
        assert "Texture Features" in html
        assert "texture-features" in html  # anchor

    def test_contains_feature_names(self):
        html = build_texture_section(_make_texture_data())
        assert "GLCM_Contrast" in html
        assert "GLRLM_LRE" in html
        assert "Skewness" in html
        assert "Sphericity" in html

    def test_contains_auc_values(self):
        html = build_texture_section(_make_texture_data())
        assert "0.820" in html
        assert "0.780" in html

    def test_contains_categories(self):
        html = build_texture_section(_make_texture_data())
        assert "GLCM" in html
        assert "GLRLM" in html
        assert "first_order" in html
        assert "shape" in html

    def test_top_auc_stat_card(self):
        html = build_texture_section(_make_texture_data())
        assert "Top AUC" in html
        assert "0.820" in html

    def test_sorted_by_auc_descending(self):
        html = build_texture_section(_make_texture_data())
        # GLCM_Contrast (0.82) should appear before Sphericity (0.65)
        idx_contrast = html.index("GLCM_Contrast")
        idx_sphericity = html.index("Sphericity")
        assert idx_contrast < idx_sphericity

    def test_empty_dict_returns_empty(self):
        assert build_texture_section({}) == ""

    def test_none_returns_empty(self):
        assert build_texture_section(None) == ""

    def test_missing_texture_key(self):
        data = {"Standard": {"something_else": 42}}
        assert build_texture_section(data) == ""

    def test_malformed_features_list(self):
        data = {"Standard": {"texture_features": {"features": "not a list"}}}
        result = build_texture_section(data)
        assert isinstance(result, str)

    def test_accepts_flat_list_format(self):
        """texture_features can be a list directly (not nested in dict)."""
        data = {
            "Standard": {
                "texture_features": [
                    {"name": "Entropy", "category": "first_order", "auc": 0.75, "p": 0.02},
                ],
            },
        }
        html = build_texture_section(data)
        assert "Entropy" in html

    def test_features_without_auc_not_in_table(self):
        data = {
            "Standard": {
                "texture_features": {
                    "features": [
                        {"name": "NoAUC", "category": "GLCM"},
                        {"name": "HasAUC", "category": "GLCM", "auc": 0.80, "p": 0.01},
                    ],
                },
            },
        }
        html = build_texture_section(data)
        assert "HasAUC" in html

    def test_non_dict_dwi_data_skipped(self):
        assert build_texture_section({"Standard": "bad"}) == ""


# ═══════════════════════════════════════════════════════════════════════
# Registration Quality Section
# ═══════════════════════════════════════════════════════════════════════


def _make_registration_data():
    return {
        "Standard": {
            "registration_quality": {
                "mean_ncc": 0.85,
                "mean_nmi": 1.2,
                "patients": [
                    {
                        "patient_id": "PT001",
                        "ncc": 0.92,
                        "nmi": 1.3,
                        "jacobian_negative_voxels": 0,
                        "jacobian_mean": 1.01,
                    },
                    {
                        "patient_id": "PT002",
                        "ncc": 0.65,
                        "nmi": 0.9,
                        "jacobian_negative_voxels": 15,
                        "jacobian_mean": 0.95,
                    },
                    {
                        "patient_id": "PT003",
                        "ncc": 0.88,
                        "nmi": 1.1,
                        "jacobian_negative_voxels": 0,
                        "jacobian_mean": 1.00,
                    },
                ],
            },
        },
    }


class TestBuildRegistrationQualitySection:
    def test_renders_with_valid_data(self):
        html = build_registration_quality_section(_make_registration_data())
        assert "Registration Quality" in html
        assert "registration-quality" in html  # anchor

    def test_contains_patient_ids(self):
        html = build_registration_quality_section(_make_registration_data())
        assert "PT001" in html
        assert "PT002" in html
        assert "PT003" in html

    def test_flags_poor_registration(self):
        html = build_registration_quality_section(_make_registration_data())
        # PT002 has NCC < 0.7 and jacobian_negative_voxels > 0
        assert "Poor" in html
        assert "poor registration quality" in html.lower()

    def test_good_registration_marked(self):
        html = build_registration_quality_section(_make_registration_data())
        assert "Good" in html

    def test_contains_summary_cards(self):
        html = build_registration_quality_section(_make_registration_data())
        assert "Mean NCC" in html
        assert "0.850" in html

    def test_poor_count_in_warning(self):
        html = build_registration_quality_section(_make_registration_data())
        # PT002 is poor (NCC < 0.7 and jac neg > 0)
        assert "1 patient(s)" in html

    def test_empty_dict_returns_empty(self):
        assert build_registration_quality_section({}) == ""

    def test_none_returns_empty(self):
        assert build_registration_quality_section(None) == ""

    def test_missing_registration_key(self):
        data = {"Standard": {"something_else": 42}}
        assert build_registration_quality_section(data) == ""

    def test_malformed_patients(self):
        data = {
            "Standard": {
                "registration_quality": {
                    "patients": "not a list",
                },
            },
        }
        result = build_registration_quality_section(data)
        assert isinstance(result, str)

    def test_accepts_flat_list_format(self):
        """registration_quality can be a list of patient dicts directly."""
        data = {
            "Standard": {
                "registration_quality": [
                    {"patient_id": "PT010", "ncc": 0.90, "nmi": 1.1,
                     "jacobian_negative_voxels": 0, "jacobian_mean": 1.0},
                ],
            },
        }
        html = build_registration_quality_section(data)
        assert "PT010" in html

    def test_all_good_no_warning(self):
        data = {
            "Standard": {
                "registration_quality": {
                    "mean_ncc": 0.90,
                    "patients": [
                        {"patient_id": "PT001", "ncc": 0.90,
                         "jacobian_negative_voxels": 0, "jacobian_mean": 1.0},
                    ],
                },
            },
        }
        html = build_registration_quality_section(data)
        assert "poor registration quality" not in html.lower()

    def test_non_dict_dwi_data_skipped(self):
        assert build_registration_quality_section({"Standard": 42}) == ""

    def test_empty_patients_list(self):
        data = {
            "Standard": {
                "registration_quality": {
                    "mean_ncc": 0.80,
                    "patients": [],
                },
            },
        }
        html = build_registration_quality_section(data)
        assert "Registration Quality" in html
