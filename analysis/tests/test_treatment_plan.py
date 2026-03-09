"""Tests for the Suggested Treatment Plan in the hypothesis section.

Covers:
- Treatment plan rendering for each response pattern (good/mixed/poor)
- Integration of survival hazard ratios into plan recommendations
- Integration of predictive model AUC/features into plan recommendations
- Inflection-point timing guidance
- Backward compatibility: _section_hypothesis with no log_data
"""

from __future__ import annotations

import json

from report_sections import _section_hypothesis


# ---------------------------------------------------------------------------
# Helper: build a minimal groups dict with Longitudinal_Mean_Metrics
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
):
    """Build a minimal log_data dict for Standard DWI type."""
    survival = {"hazard_ratios": hr_entries or [], "global_lrt": None, "ipcw": None}
    predictive = {
        "roc_analyses": roc_analyses or [],
        "feature_selections": feature_selections or [],
    }
    return {
        "Standard": {
            "survival": survival,
            "stats_predictive": predictive,
            "stats_comparisons": {},
            "baseline": {},
        }
    }


# ---------------------------------------------------------------------------
# Tests: Treatment plan core recommendations
# ---------------------------------------------------------------------------

class TestTreatmentPlanCoreRecommendation:
    """Verify the core treatment recommendation based on D/f trends."""

    def test_good_response_continue_regimen(self):
        """D increasing + f decreasing → continue current regimen."""
        groups = _make_groups("increasing", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "Suggested Treatment Plan" in html
        assert "current radiotherapy regimen appears effective" in html

    def test_poor_cellular_good_vascular(self):
        """D stable + f decreasing → dose escalation / radiosensitizer."""
        groups = _make_groups("stable", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "Suggested Treatment Plan" in html
        assert "dose escalation" in html.lower() or "radiosensitiz" in html.lower()

    def test_good_cellular_poor_vascular(self):
        """D increasing + f stable → anti-angiogenic consideration."""
        groups = _make_groups("increasing", "stable")
        html = "\n".join(_section_hypothesis(groups))
        assert "Suggested Treatment Plan" in html
        assert "anti-angiogenic" in html.lower()

    def test_poor_both(self):
        """D stable + f stable → multidisciplinary review."""
        groups = _make_groups("stable", "stable")
        html = "\n".join(_section_hypothesis(groups))
        assert "Suggested Treatment Plan" in html
        assert "multidisciplinary" in html.lower()

    def test_unknown_trends(self):
        """No longitudinal data → still produces a plan (fallback)."""
        html = "\n".join(_section_hypothesis({}))
        assert "Suggested Treatment Plan" in html
        assert "multidisciplinary" in html.lower()


# ---------------------------------------------------------------------------
# Tests: Inflection-point timing guidance
# ---------------------------------------------------------------------------

class TestTreatmentPlanTiming:
    """Verify inflection-point timing guidance in the treatment plan."""

    def test_inflection_timing_mentioned(self):
        """Inflection points should generate timing guidance with Fx labels."""
        groups = _make_groups_with_inflections("increasing", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "Fx5" in html
        assert "Fx3" in html
        assert "adaptive replanning" in html.lower()


# ---------------------------------------------------------------------------
# Tests: Survival-informed guidance
# ---------------------------------------------------------------------------

class TestTreatmentPlanSurvival:
    """Verify survival hazard ratio integration in the treatment plan."""

    def test_protective_hr_mentioned(self):
        """HR < 1 with p < 0.05 should appear as protective."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(hr_entries=[
            {"covariate": "delta_d", "hr": 0.75, "ci_lo": 0.55, "ci_hi": 1.02, "p": 0.034},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "delta_d" in html
        assert "protective" in html.lower()

    def test_risk_hr_mentioned(self):
        """HR > 1 with p < 0.05 should appear as risk factor."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(hr_entries=[
            {"covariate": "mean_adc", "hr": 1.85, "ci_lo": 1.10, "ci_hi": 3.12, "p": 0.021},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "mean_adc" in html
        assert "risk factor" in html.lower()

    def test_nonsignificant_hr_omitted(self):
        """HR with p >= 0.05 should not appear in treatment plan."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(hr_entries=[
            {"covariate": "mean_adc", "hr": 1.25, "ci_lo": 0.98, "ci_hi": 1.59, "p": 0.068},
        ])
        html = "\n".join(_section_hypothesis(groups, log_data))
        # Should not mention the non-significant covariate in treatment plan
        plan_idx = html.find("Suggested Treatment Plan")
        plan_section = html[plan_idx:] if plan_idx >= 0 else ""
        assert "mean_adc" not in plan_section


# ---------------------------------------------------------------------------
# Tests: Predictive model guidance
# ---------------------------------------------------------------------------

class TestTreatmentPlanPredictive:
    """Verify predictive model AUC/feature integration."""

    def test_high_auc_features_mentioned(self):
        """AUC >= 0.8 with features should produce high-confidence guidance."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(
            roc_analyses=[{"timepoint": "BL", "auc": 0.85}],
            feature_selections=[
                {"timepoint": "BL", "lambda": 0.05,
                 "features": ["mean_adc", "delta_d", "f_ratio"]},
            ],
        )
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "0.850" in html
        assert "high" in html.lower()
        assert "mean_adc" in html

    def test_moderate_auc(self):
        """AUC between 0.65 and 0.7 → moderate confidence."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(
            roc_analyses=[{"timepoint": "W2", "auc": 0.66}],
            feature_selections=[
                {"timepoint": "W2", "lambda": 0.1, "features": ["delta_adc"]},
            ],
        )
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "moderate" in html.lower()

    def test_low_auc_omitted(self):
        """AUC < 0.65 should not generate predictive guidance."""
        groups = _make_groups("increasing", "decreasing")
        log_data = _make_log_data(
            roc_analyses=[{"timepoint": "BL", "auc": 0.55}],
            feature_selections=[
                {"timepoint": "BL", "lambda": 0.05, "features": ["mean_adc"]},
            ],
        )
        html = "\n".join(_section_hypothesis(groups, log_data))
        assert "discriminative ability" not in html


# ---------------------------------------------------------------------------
# Tests: Backward compatibility
# ---------------------------------------------------------------------------

class TestTreatmentPlanBackwardCompat:
    """Ensure _section_hypothesis works without log_data (old call signature)."""

    def test_no_log_data_still_works(self):
        """Calling with only groups (old signature) should not error."""
        groups = _make_groups("increasing", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "Suggested Treatment Plan" in html
        assert "Data-Driven Hypothesis" in html

    def test_none_log_data_still_works(self):
        """Explicit None for log_data should not error."""
        groups = _make_groups("stable", "stable")
        html = "\n".join(_section_hypothesis(groups, None))
        assert "Suggested Treatment Plan" in html

    def test_empty_groups_no_crash(self):
        """Empty groups + no log_data should produce a valid section."""
        html = "\n".join(_section_hypothesis({}, None))
        assert "Data-Driven Hypothesis" in html
        assert "Suggested Treatment Plan" in html

    def test_disclaimer_always_present(self):
        """The physician review disclaimer should always appear."""
        groups = _make_groups("increasing", "decreasing")
        html = "\n".join(_section_hypothesis(groups))
        assert "treating physician" in html.lower()
