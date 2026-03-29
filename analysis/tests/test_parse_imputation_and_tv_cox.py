"""Tests for parsing: imputation sensitivity AUC and time-varying Cox HR.

Validates regex extraction from MATLAB log output for:
- Imputation sensitivity comparison table (KNN/LOCF/Mean/Linear_Interp AUC)
- Time-varying Cox model interaction terms and stratification info
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from parsers.parse_log_metrics import (
    parse_stats_predictive,
    parse_survival,
)


# ── Imputation Sensitivity ──


IMPUTATION_LOG = """\
  --- Imputation Sensitivity Analysis ---
  [1/4] KNN imputation...
  [2/4] LOCF imputation...
  [3/4] Mean imputation...
  [4/4] Linear interpolation...

  Imputation Sensitivity Results:
  Method                      AUC      N_Imputed
  --------------------------------------------
  KNN                     0.843            42
  LOCF                    0.812            38
  Mean                    0.798            42
  Linear_Interp           0.825            40

  Concordance (Spearman rho between risk scores):
                       KNN      LOCF      Mean  Linear_Interp
  KNN                 1.000     0.923     0.890     0.912
  LOCF                0.923     1.000     0.856     0.945
  Mean                0.890     0.856     1.000     0.878
  Linear_Interp       0.912     0.945     0.878     1.000
"""


class TestParseImputationSensitivity:
    """Verify extraction of imputation sensitivity AUC values."""

    def test_extracts_four_methods(self):
        result = parse_stats_predictive(IMPUTATION_LOG)
        imp = result["imputation_sensitivity"]
        assert len(imp) == 4

    def test_method_names(self):
        result = parse_stats_predictive(IMPUTATION_LOG)
        methods = [e["method"] for e in result["imputation_sensitivity"]]
        assert methods == ["KNN", "LOCF", "Mean", "Linear_Interp"]

    def test_auc_values(self):
        result = parse_stats_predictive(IMPUTATION_LOG)
        aucs = [e["auc"] for e in result["imputation_sensitivity"]]
        assert aucs == pytest.approx([0.843, 0.812, 0.798, 0.825])

    def test_n_imputed_values(self):
        result = parse_stats_predictive(IMPUTATION_LOG)
        counts = [e["n_imputed"] for e in result["imputation_sensitivity"]]
        assert counts == [42, 38, 42, 40]

    def test_empty_log_returns_empty_list(self):
        result = parse_stats_predictive("")
        assert result["imputation_sensitivity"] == []

    def test_no_imputation_section_returns_empty_list(self):
        text = "PRIMARY ROC ANALYSIS for BL\nAUC = 0.75\n"
        result = parse_stats_predictive(text)
        assert result["imputation_sensitivity"] == []

    def test_nan_auc_handled(self):
        text = "  KNN                     NaN            10\n"
        result = parse_stats_predictive(text)
        imp = result["imputation_sensitivity"]
        assert len(imp) == 1
        assert imp[0]["method"] == "KNN"
        import math
        assert math.isnan(imp[0]["auc"])

    def test_mixed_with_roc_analysis(self):
        """Imputation sensitivity and ROC blocks coexist in the same log."""
        text = (
            "PRIMARY ROC ANALYSIS for BL\n"
            "AUC = 0.78\n"
            "Sensitivity = 80.0% | Specificity = 70.0%\n\n"
            + IMPUTATION_LOG
        )
        result = parse_stats_predictive(text)
        assert len(result["roc_analyses"]) == 1
        assert result["roc_analyses"][0]["auc"] == pytest.approx(0.78)
        assert len(result["imputation_sensitivity"]) == 4


# ── Time-Varying Cox ──


TV_COX_LOG = """\
  --- Time-Varying Cox Model (PH Violation Follow-Up) ---
  PH violations detected for: mean_adc, delta_d

  Stratified Cox Model (stratified by mean_adc, median=0.0012):
  Covariate        HR        SE         p
  ----------------------------------------
  delta_d       0.720    0.1500    0.0600

  Extended Cox: mean_adc × log(time) interaction
    Base mean_adc: coef=0.3500, p=0.0200
    mean_adc × log(t): coef=-0.0800, p=0.0150

  Extended Cox: delta_d × log(time) interaction
    Base delta_d: coef=-0.2100, p=0.0400
    delta_d × log(t): coef=0.0500, p=0.1200
"""


class TestParseTimeVaryingCox:
    """Verify extraction of time-varying Cox model results."""

    def test_extracts_violated_covariates(self):
        result = parse_survival(TV_COX_LOG)
        tv = result["time_varying_cox"]
        assert tv is not None
        assert tv["violated_covariates"] == ["mean_adc", "delta_d"]

    def test_stratified_by(self):
        result = parse_survival(TV_COX_LOG)
        tv = result["time_varying_cox"]
        assert tv["stratified_by"] == "mean_adc"

    def test_interaction_models_count(self):
        result = parse_survival(TV_COX_LOG)
        tv = result["time_varying_cox"]
        assert len(tv["interaction_models"]) == 2

    def test_first_interaction_model(self):
        result = parse_survival(TV_COX_LOG)
        im = result["time_varying_cox"]["interaction_models"][0]
        assert im["covariate"] == "mean_adc"
        assert im["base_coef"] == pytest.approx(0.35)
        assert im["base_p"] == pytest.approx(0.02)
        assert im["interaction_coef"] == pytest.approx(-0.08)
        assert im["interaction_p"] == pytest.approx(0.015)

    def test_second_interaction_model(self):
        result = parse_survival(TV_COX_LOG)
        im = result["time_varying_cox"]["interaction_models"][1]
        assert im["covariate"] == "delta_d"
        assert im["interaction_coef"] == pytest.approx(0.05)
        assert im["interaction_p"] == pytest.approx(0.12)

    def test_no_tv_cox_returns_none(self):
        text = "Global LRT: chi2(2) = 5.00, p = 0.0800\n"
        result = parse_survival(text)
        assert result["time_varying_cox"] is None

    def test_empty_log_returns_none(self):
        result = parse_survival("")
        assert result["time_varying_cox"] is None

    def test_single_violation(self):
        text = (
            "PH violations detected for: mean_f\n"
            "  Base mean_f: coef=0.1000, p=0.0300\n"
            "  mean_f \u00d7 log(t): coef=-0.0200, p=0.0400\n"
        )
        result = parse_survival(text)
        tv = result["time_varying_cox"]
        assert tv["violated_covariates"] == ["mean_f"]
        assert len(tv["interaction_models"]) == 1
        assert tv["interaction_models"][0]["covariate"] == "mean_f"

    def test_coexists_with_standard_survival(self):
        """Time-varying Cox parsing doesn't interfere with standard HR parsing."""
        text = (
            "  mean_adc    1.500   0.900   2.500   0.0300\n"
            "Global LRT: chi2(2) = 7.82, p = 0.0200\n"
            "IPCW weights applied ... [0.85, 1.42]\n\n"
            + TV_COX_LOG
        )
        result = parse_survival(text)
        assert len(result["hazard_ratios"]) == 1
        assert result["global_lrt"]["p"] == pytest.approx(0.02)
        assert result["ipcw"]["min_weight"] == pytest.approx(0.85)
        assert result["time_varying_cox"] is not None
        assert len(result["time_varying_cox"]["interaction_models"]) == 2
