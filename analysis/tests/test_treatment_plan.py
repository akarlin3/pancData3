"""Tests for the Suggested Treatment Plan in the hypothesis section.

Covers:
- Treatment plan rendering for each response pattern (good/mixed/poor)
- GLME significance qualifying the core recommendation
- GLME detail metrics with p-values and adjusted alpha
- FDR-significant timepoint references
- Survival hazard ratios with numeric HR, CI, and p-values
- Predictive model AUC, sensitivity, specificity, Youden cutoff, features
- Dosimetry D95/V50 integration with dose-concern detection
- Inflection-point timing guidance
- Backward compatibility: _section_hypothesis with no log_data/mat_data
"""

from __future__ import annotations

import json

from report_sections import _section_hypothesis


# ---------------------------------------------------------------------------
# Helper builders
# ---------------------------------------------------------------------------

def _make_groups(d_direction: str, f_direction: str):
    """Return a groups dict with specified D and f trend directions."""
    return {
        "Longitudinal_Mean_Metrics": {
            "Standard": {
                "trends_json": json.dumps([
                    {"series": "Mean D", "direction": d_direction},
                    {"series": "Mean f", "direction": f_direction},
                ]),
                "inflection_points_json": "[]",
            }
        }
    }


def _make_groups_with_inflections(d_direction: str, f_direction: str):
    """Return groups with inflection points for timing guidance."""
    return {
        "Longitudinal_Mean_Metrics": {
            "Standard": {
                "trends_json": json.dumps([
                    {"series": "Mean D", "direction": d_direction},
                    {"series": "Mean f", "direction": f_direction},
                ]),
                "inflection_points_json": json.dumps([
                    {
                        "approximate_x": 5.0,
                        "approximate_y": 0.0015,
                        "description": "significant diffusion increase at Fx5",
                    },
                    {
                        "approximate_x": 3.0,
                        "approximate_y": None,
                        "description": "vascular f drop observed around Fx3",
                    },
                ]),
            }
        }
    }


def _make_log_data(
    *,
    hr_entries=None,
    roc_analyses=None,
    feature_selections=None,
    glme_interactions=None,
    glme_details=None,
    fdr_timepoints=None,
):
    """Build a minimal log_data dict for Standard DWI type."""
    survival = {"hazard_ratios": hr_entries or [], "global_lrt": None, "ipcw": None}
    predictive = {
        "roc_analyses": roc_analyses or [],
        "feature_selections": feature_selections or [],
    }
    stats_comparisons = {
        "glme_interactions": glme_interactions or [],
        "glme_details": glme_details or [],
        "fdr_timepoints": fdr_timepoints or [],
    }
    return {
        "Standard": {
            "survival": survival,
            "stats_predictive": predictive,
            "stats_comparisons": stats_comparisons,
            "baseline": {},
        }
    }


def _make_mat_data(*, d95_adc=None, v50_adc=None, d95_d=None, v50_d=None):
    """Build a minimal mat_data dict with dosimetry for Standard."""
    dosi = {}
    if d95_adc is not None:
        dosi["d95_adc_mean"] = d95_adc
    if v50_adc is not None:
        dosi["v50_adc_mean"] = v50_adc
    if d95_d is not None:
        dosi["d95_d_mean"] = d95_d
    if v50_d is not None:
        dosi["v50_d_mean"] = v50_d
    return {"Standard": {"dosimetry": dosi}} if dosi else {}


# ---------------------------------------------------------------------------
# Tests: Core recommendations
# ---------------------------------------------------------------------------

class TestTreatmentPlanCoreRecommendation:
    """Verify the core treatment recommendation based on D/f trends."""

    def test_good_response_continue_regimen(self):
        groups = _make_groups("increasing", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "Suggested Treatment Plan" in html
        assert "current radiotherapy regimen appears effective" in html

    def test_poor_cellular_good_vascular(self):
        groups = _make_groups("stable", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "dose escalation" in html.lower() or "radiosensitiz" in html.lower()

    def test_good_cellular_poor_vascular(self):
        groups = _make_groups("increasing", "stable")
        html = "\n".join(_section_hypothesis(groups))
        assert "anti-angiogenic" in html.lower()

    def test_poor_both(self):
        groups = _make_groups("stable", "stable")
        html = "\n".join(_section_hypothesis(groups))
        assert "multidisciplinary" in html.lower()

    def test_unknown_trends_fallback(self):
        html = "\n".join(_section_hypothesis({}))
        assert "Suggested Treatment Plan" in html
        assert "multidisciplinary" in html.lower()


# ---------------------------------------------------------------------------
# Tests: GLME significance qualifying the recommendation
# ---------------------------------------------------------------------------

class TestGLMESignificance:
    """Verify GLME interaction p-values qualify the core recommendation."""

    def test_significant_glme_noted(self):
        """Significant GLME interaction should include p-value."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(glme_interactions=[0.023])
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "0.0230" in html
        assert "diverge significantly" in html

    def test_nonsignificant_glme_noted(self):
        """No significant GLME should note lack of significance."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(glme_interactions=[0.12])
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "did not identify" in html
        assert "observed trends only" in html

    def test_no_log_data_no_qualifier(self):
        """Without log_data, no GLME qualifier should appear."""
        groups = _make_groups("increasing", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "GLME" not in html


# ---------------------------------------------------------------------------
# Tests: GLME detail metrics with p-values
# ---------------------------------------------------------------------------

class TestGLMEDetails:
    """Verify per-metric GLME detail rows with p and adjusted alpha."""

    def test_significant_details_shown(self):
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(glme_details=[
            {"metric": "mean_adc", "p": 0.015, "adj_alpha": 0.025},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "mean_adc" in html
        assert "0.0150" in html
        assert "0.0250" in html
        assert "BH-FDR" in html

    def test_nonsignificant_details_omitted(self):
        """Details with p >= adj_alpha should not appear."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(glme_details=[
            {"metric": "mean_d*", "p": 0.042, "adj_alpha": 0.0125},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        plan_idx = html.find("Suggested Treatment Plan")
        plan_section = html[plan_idx:] if plan_idx >= 0 else ""
        assert "mean_d*" not in plan_section

    def test_hypothesis_sig_annotation_cellular(self):
        """Significant cellular metric adds p-value to hypothesis text."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(glme_details=[
            {"metric": "mean_adc", "p": 0.008, "adj_alpha": 0.025},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "statistically significant" in html
        assert "0.0080" in html

    def test_hypothesis_sig_annotation_vascular(self):
        """Significant vascular metric adds p-value to hypothesis text."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(glme_details=[
            {"metric": "mean_d*", "p": 0.003, "adj_alpha": 0.0125},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        # d* with p < adj_alpha → significant → annotated in hypothesis
        assert "0.0030" in html


# ---------------------------------------------------------------------------
# Tests: FDR timepoint references
# ---------------------------------------------------------------------------

class TestFDRTimepoints:
    """Verify FDR-significant timepoints appear in the treatment plan."""

    def test_fdr_timepoints_shown(self):
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(fdr_timepoints=[
            {"timepoint": "BL", "n_significant": 3},
            {"timepoint": "W2", "n_significant": 1},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "BL" in html
        assert "3 significant metrics" in html
        assert "W2" in html
        assert "1 significant metric)" in html  # singular
        assert "critical monitoring windows" in html

    def test_zero_significant_omitted(self):
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(fdr_timepoints=[
            {"timepoint": "W4", "n_significant": 0},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        plan_idx = html.find("Suggested Treatment Plan")
        plan_section = html[plan_idx:] if plan_idx >= 0 else ""
        assert "W4" not in plan_section


# ---------------------------------------------------------------------------
# Tests: Survival hazard ratios with numeric values
# ---------------------------------------------------------------------------

class TestSurvivalHR:
    """Verify HR, CI, and p-values appear in the treatment plan."""

    def test_protective_hr_with_values(self):
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(hr_entries=[
            {"covariate": "delta_d", "hr": 0.750, "ci_lo": 0.55, "ci_hi": 1.02, "p": 0.034},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "delta_d" in html
        assert "0.750" in html
        assert "0.55" in html
        assert "1.02" in html
        assert "0.0340" in html
        assert "Protective factors" in html

    def test_risk_hr_with_values(self):
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(hr_entries=[
            {"covariate": "mean_adc", "hr": 1.850, "ci_lo": 1.10, "ci_hi": 3.12, "p": 0.021},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "mean_adc" in html
        assert "1.850" in html
        assert "Risk factors" in html

    def test_nonsignificant_hr_omitted(self):
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(hr_entries=[
            {"covariate": "mean_adc", "hr": 1.25, "ci_lo": 0.98, "ci_hi": 1.59, "p": 0.068},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        plan_idx = html.find("Suggested Treatment Plan")
        plan_section = html[plan_idx:] if plan_idx >= 0 else ""
        assert "mean_adc" not in plan_section

    def test_mixed_hr_both_sections(self):
        """Both protective and risk HR should appear when present."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(hr_entries=[
            {"covariate": "delta_d", "hr": 0.75, "ci_lo": 0.55, "ci_hi": 1.02, "p": 0.034},
            {"covariate": "mean_adc", "hr": 1.85, "ci_lo": 1.10, "ci_hi": 3.12, "p": 0.021},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "Protective factors" in html
        assert "Risk factors" in html


# ---------------------------------------------------------------------------
# Tests: Predictive model with numeric values
# ---------------------------------------------------------------------------

class TestPredictiveModel:
    """Verify AUC, sensitivity, specificity, Youden, and features appear."""

    def test_high_auc_with_performance(self):
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(
            roc_analyses=[{
                "timepoint": "BL",
                "auc": 0.850,
                "sensitivity": 82.5,
                "specificity": 68.2,
                "youden_cutoff": 0.450,
            }],
            feature_selections=[
                {"timepoint": "BL", "lambda": 0.0532,
                 "features": ["mean_adc", "mean_d", "f_ratio"]},
            ],
        )
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "0.850" in html
        assert "82.5%" in html
        assert "68.2%" in html
        assert "0.450" in html
        assert "high" in html  # high confidence
        assert "mean_adc" in html
        assert "0.0532" in html  # lambda

    def test_moderate_auc(self):
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(
            roc_analyses=[{"timepoint": "W2", "auc": 0.66}],
            feature_selections=[
                {"timepoint": "W2", "lambda": 0.1, "features": ["delta_adc"]},
            ],
        )
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "moderate" in html.lower()
        assert "0.660" in html

    def test_low_auc_omitted(self):
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(
            roc_analyses=[{"timepoint": "BL", "auc": 0.55}],
            feature_selections=[
                {"timepoint": "BL", "lambda": 0.05, "features": ["mean_adc"]},
            ],
        )
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "discriminative" not in html

    def test_multiple_timepoint_features(self):
        """Feature selections from multiple timepoints should appear."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(
            roc_analyses=[{"timepoint": "BL", "auc": 0.78}],
            feature_selections=[
                {"timepoint": "BL", "lambda": 0.05,
                 "features": ["mean_adc", "mean_d"]},
                {"timepoint": "W2", "lambda": 0.12,
                 "features": ["delta_adc", "mean_f"]},
            ],
        )
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "BL:" in html
        assert "W2:" in html


# ---------------------------------------------------------------------------
# Tests: Dosimetry integration
# ---------------------------------------------------------------------------

class TestDosimetry:
    """Verify D95/V50 values and dose-concern detection."""

    def test_adequate_dose(self):
        groups = _make_groups("increasing", "decreasing")
        mat_data = _make_mat_data(d95_adc=48.5, v50_adc=0.95, d95_d=47.2, v50_d=0.92)
        html = "\n".join(_section_hypothesis(groups, mat_data=mat_data))
        assert "48.50" in html
        assert "95.0%" in html
        assert "47.20" in html
        assert "appears adequate" in html
        assert "concern" not in html.lower()

    def test_underdosed_d95(self):
        groups = _make_groups("increasing", "decreasing")
        mat_data = _make_mat_data(d95_adc=38.2, v50_adc=0.95)
        html = "\n".join(_section_hypothesis(groups, mat_data=mat_data))
        assert "38.20" in html
        assert "under-dosed" in html
        assert "below 45" in html
        assert "dose escalation" in html.lower()

    def test_low_v50(self):
        groups = _make_groups("increasing", "decreasing")
        mat_data = _make_mat_data(d95_adc=50.0, v50_adc=0.82)
        html = "\n".join(_section_hypothesis(groups, mat_data=mat_data))
        assert "82.0%" in html
        assert "concern" in html.lower()

    def test_no_mat_data_no_dosimetry(self):
        groups = _make_groups("increasing", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "Dosimetry" not in html
        assert "D95" not in html


# ---------------------------------------------------------------------------
# Tests: Inflection-point timing
# ---------------------------------------------------------------------------

class TestTreatmentPlanTiming:
    """Verify inflection-point timing guidance in the treatment plan."""

    def test_inflection_timing_mentioned(self):
        groups = _make_groups_with_inflections("increasing", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "Fx5" in html
        assert "Fx3" in html
        assert "adaptive replanning" in html.lower()


# ---------------------------------------------------------------------------
# Tests: Backward compatibility
# ---------------------------------------------------------------------------

class TestBackwardCompat:
    """Ensure _section_hypothesis works with old call signatures."""

    def test_no_log_data_still_works(self):
        groups = _make_groups("increasing", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "Suggested Treatment Plan" in html
        assert "Data-Driven Hypothesis" in html

    def test_none_log_data_and_mat_data(self):
        groups = _make_groups("stable", "stable")
        html = "\n".join(_section_hypothesis(groups, None, None))
        assert "Suggested Treatment Plan" in html

    def test_empty_groups_no_crash(self):
        html = "\n".join(_section_hypothesis({}, None, None))
        assert "Data-Driven Hypothesis" in html

    def test_disclaimer_always_present(self):
        groups = _make_groups("increasing", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "treating physician" in html.lower()

    def test_positional_args_backward_compat(self):
        """Old callers passing only groups should still work."""
        groups = _make_groups("increasing", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "Data-Driven Hypothesis" in html
