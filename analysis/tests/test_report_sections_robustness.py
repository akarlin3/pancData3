"""Tests for the model_robustness report section builder.

Validates HTML generation for:
- Imputation sensitivity AUC comparison tables
- Time-varying Cox HR summary tables
- Edge cases (empty data, missing fields, single DWI type)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.model_robustness import _section_model_robustness


# ── Fixtures ──


def _make_log_data_with_imputation():
    """Log data with imputation sensitivity results."""
    return {
        "Standard": {
            "stats_predictive": {
                "imputation_sensitivity": [
                    {"method": "KNN", "auc": 0.843, "n_imputed": 42},
                    {"method": "LOCF", "auc": 0.812, "n_imputed": 38},
                    {"method": "Mean", "auc": 0.798, "n_imputed": 42},
                    {"method": "Linear_Interp", "auc": 0.825, "n_imputed": 40},
                ],
                "roc_analyses": [],
                "feature_selections": [],
            },
            "survival": {},
        },
    }


def _make_log_data_with_tv_cox():
    """Log data with time-varying Cox results."""
    return {
        "Standard": {
            "stats_predictive": {},
            "survival": {
                "time_varying_cox": {
                    "violated_covariates": ["mean_adc", "delta_d"],
                    "stratified_by": "mean_adc",
                    "interaction_models": [
                        {
                            "covariate": "mean_adc",
                            "base_coef": 0.35,
                            "base_p": 0.02,
                            "interaction_coef": -0.08,
                            "interaction_p": 0.015,
                        },
                        {
                            "covariate": "delta_d",
                            "base_coef": -0.21,
                            "base_p": 0.04,
                            "interaction_coef": 0.05,
                            "interaction_p": 0.12,
                        },
                    ],
                },
            },
        },
    }


def _make_log_data_full():
    """Log data with both imputation sensitivity and time-varying Cox."""
    data = _make_log_data_with_imputation()
    data["Standard"]["survival"] = _make_log_data_with_tv_cox()["Standard"]["survival"]
    return data


# ── Tests: Empty / No Data ──


class TestModelRobustnessEmpty:
    def test_none_log_data_returns_empty(self):
        result = _section_model_robustness(None, ["Standard"], {})
        assert result == []

    def test_no_relevant_data_returns_empty(self):
        log_data = {"Standard": {"stats_predictive": {}, "survival": {}}}
        result = _section_model_robustness(log_data, ["Standard"], {})
        assert result == []

    def test_missing_dwi_type_returns_empty(self):
        log_data = _make_log_data_with_imputation()
        result = _section_model_robustness(log_data, ["dnCNN"], {})
        assert result == []


# ── Tests: Imputation Sensitivity ──


class TestImputationSensitivitySection:
    def test_returns_html_list(self):
        result = _section_model_robustness(
            _make_log_data_with_imputation(), ["Standard"], {}
        )
        assert isinstance(result, list)
        assert len(result) > 0

    def test_contains_heading(self):
        result = _section_model_robustness(
            _make_log_data_with_imputation(), ["Standard"], {}
        )
        html = "\n".join(result)
        assert "Model Robustness" in html
        assert "model-robustness" in html

    def test_contains_method_names(self):
        result = _section_model_robustness(
            _make_log_data_with_imputation(), ["Standard"], {}
        )
        html = "\n".join(result)
        assert "KNN" in html
        assert "LOCF" in html
        assert "Mean" in html
        assert "Linear_Interp" in html

    def test_contains_auc_values(self):
        result = _section_model_robustness(
            _make_log_data_with_imputation(), ["Standard"], {}
        )
        html = "\n".join(result)
        assert "0.843" in html
        assert "0.812" in html

    def test_contains_best_auc_card(self):
        result = _section_model_robustness(
            _make_log_data_with_imputation(), ["Standard"], {}
        )
        html = "\n".join(result)
        assert "Best AUC" in html
        assert "KNN" in html

    def test_stable_spread_message(self):
        """AUC spread < 0.05 shows 'robust' info box."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "imputation_sensitivity": [
                        {"method": "KNN", "auc": 0.80, "n_imputed": 10},
                        {"method": "LOCF", "auc": 0.79, "n_imputed": 10},
                        {"method": "Mean", "auc": 0.78, "n_imputed": 10},
                        {"method": "Linear_Interp", "auc": 0.80, "n_imputed": 10},
                    ],
                },
                "survival": {},
            },
        }
        result = _section_model_robustness(log_data, ["Standard"], {})
        html = "\n".join(result)
        assert "robust" in html.lower()

    def test_variable_spread_warning(self):
        """AUC spread >= 0.10 shows warning box."""
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "imputation_sensitivity": [
                        {"method": "KNN", "auc": 0.85, "n_imputed": 10},
                        {"method": "LOCF", "auc": 0.72, "n_imputed": 10},
                        {"method": "Mean", "auc": 0.70, "n_imputed": 10},
                        {"method": "Linear_Interp", "auc": 0.80, "n_imputed": 10},
                    ],
                },
                "survival": {},
            },
        }
        result = _section_model_robustness(log_data, ["Standard"], {})
        html = "\n".join(result)
        assert "warn-box" in html


# ── Tests: Time-Varying Cox ──


class TestTimeVaryingCoxSection:
    def test_returns_html_list(self):
        result = _section_model_robustness(
            _make_log_data_with_tv_cox(), ["Standard"], {}
        )
        assert isinstance(result, list)
        assert len(result) > 0

    def test_contains_heading(self):
        result = _section_model_robustness(
            _make_log_data_with_tv_cox(), ["Standard"], {}
        )
        html = "\n".join(result)
        assert "Time-Varying Cox" in html

    def test_contains_violated_covariates(self):
        result = _section_model_robustness(
            _make_log_data_with_tv_cox(), ["Standard"], {}
        )
        html = "\n".join(result)
        assert "mean_adc" in html
        assert "delta_d" in html

    def test_contains_stratification_info(self):
        result = _section_model_robustness(
            _make_log_data_with_tv_cox(), ["Standard"], {}
        )
        html = "\n".join(result)
        assert "Stratified" in html or "stratified" in html.lower()
        assert "mean_adc" in html

    def test_contains_interaction_coefficients(self):
        result = _section_model_robustness(
            _make_log_data_with_tv_cox(), ["Standard"], {}
        )
        html = "\n".join(result)
        assert "-0.0800" in html
        assert "0.0150" in html

    def test_significant_interaction_interpretation(self):
        """Significant interaction (p<0.05) shows directional interpretation."""
        result = _section_model_robustness(
            _make_log_data_with_tv_cox(), ["Standard"], {}
        )
        html = "\n".join(result)
        # mean_adc has interaction_coef=-0.08, p=0.015 → "HR decreases over time"
        assert "decreases over time" in html

    def test_nonsig_interaction_interpretation(self):
        """Non-significant interaction shows 'No significant time variation'."""
        result = _section_model_robustness(
            _make_log_data_with_tv_cox(), ["Standard"], {}
        )
        html = "\n".join(result)
        # delta_d has interaction_p=0.12 → "No significant time variation"
        assert "No significant time variation" in html


# ── Tests: Combined ──


class TestModelRobustnessCombined:
    def test_both_sections_present(self):
        result = _section_model_robustness(
            _make_log_data_full(), ["Standard"], {}
        )
        html = "\n".join(result)
        assert "Imputation Sensitivity" in html
        assert "Time-Varying Cox" in html

    def test_multiple_dwi_types(self):
        log_data = _make_log_data_full()
        log_data["dnCNN"] = {
            "stats_predictive": {
                "imputation_sensitivity": [
                    {"method": "KNN", "auc": 0.79, "n_imputed": 30},
                    {"method": "LOCF", "auc": 0.77, "n_imputed": 28},
                ],
            },
            "survival": {},
        }
        result = _section_model_robustness(log_data, ["Standard", "dnCNN"], {})
        html = "\n".join(result)
        assert "Standard" in html
        assert "dnCNN" in html

    def test_dwi_badge_present(self):
        result = _section_model_robustness(
            _make_log_data_with_imputation(), ["Standard"], {}
        )
        html = "\n".join(result)
        assert "badge" in html
