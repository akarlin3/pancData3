"""Shared imports and fixtures for report sections statistics tests.

NOT a test module (underscore prefix prevents pytest collection). Imported via
`from _test_report_sections_statistics_shared import *` by the split test files.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.effect_sizes import (
    _section_effect_sizes,
    _section_multiple_comparisons,
)
from report.sections.model_diagnostics import (
    _section_model_diagnostics,
    _section_sensitivity_analysis,
)
from report.sections.power_analysis import (
    _section_power_analysis,
)

# Names are underscore-prefixed, so `import *` skips them unless listed here.
__all__ = [
    "pytest",
    "_section_effect_sizes",
    "_section_multiple_comparisons",
    "_section_model_diagnostics",
    "_section_sensitivity_analysis",
    "_section_power_analysis",
    "_make_log_data",
    "_make_mat_data",
]


# ── Shared fixtures ──

def _make_log_data():
    """Create synthetic log data for statistics tests."""
    return {
        "Standard": {
            "survival": {
                "hazard_ratios": [
                    {"covariate": "mean_adc", "hr": 1.5, "ci_lo": 0.9, "ci_hi": 2.5, "p": 0.03},
                    {"covariate": "delta_d", "hr": 0.7, "ci_lo": 0.5, "ci_hi": 1.0, "p": 0.06},
                ],
                "ipcw": {"max_weight": 1.8, "min_weight": 0.6},
            },
            "stats_comparisons": {
                "glme_details": [
                    {"metric": "mean_adc", "p": 0.01, "adj_alpha": 0.025},
                    {"metric": "mean_d", "p": 0.04, "adj_alpha": 0.0125},
                    {"metric": "mean_f", "p": 0.80, "adj_alpha": 0.05},
                ],
                "glme_excluded": {"pct": 12.0},
            },
            "stats_predictive": {
                "roc_analyses": [
                    {"auc": 0.78, "timepoint": "BL"},
                    {"auc": 0.65, "timepoint": "W2"},
                ],
                "feature_selections": [
                    {"timepoint": "BL", "features": ["adc", "d"], "lambda": 0.05},
                    {"timepoint": "W2", "features": ["adc"], "lambda": 0.12},
                ],
            },
            "baseline": {
                "total_outliers": {"pct": 8.5},
                "baseline_exclusion": {
                    "n_excluded": 6, "n_total": 48,
                    "lf_rate_included": 35.0, "lf_rate_excluded": 55.0,
                },
            },
            "sanity_checks": {
                "all_converged": True,
                "total_convergence": 0,
                "dim_mismatches": 0,
                "nan_dose_warnings": 1,
            },
        }
    }


def _make_mat_data():
    return {
        "Standard": {
            "longitudinal": {"num_patients": 42, "num_timepoints": 5},
        }
    }
