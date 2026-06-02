"""Shared fixtures/helpers for test_report_sections_discussion* modules.

Not a test module itself (name does not match ``test_*``); imported via
``from _test_report_sections_discussion_shared import *``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))


def _make_log_data():
    return {
        "Standard": {
            "survival": {
                "hazard_ratios": [
                    {"covariate": "mean_adc", "hr": 1.5, "ci_lo": 0.9, "ci_hi": 2.5, "p": 0.03},
                ],
                "ipcw": {"max_weight": 1.4, "min_weight": 0.8},
            },
            "stats_comparisons": {
                "glme_details": [
                    {"metric": "mean_adc", "p": 0.01, "adj_alpha": 0.025},
                ],
                "glme_interactions": [0.02],
                "fdr_timepoints": [{"n_significant": 1, "timepoint": "BL"}],
            },
            "stats_predictive": {
                "roc_analyses": [{"auc": 0.78, "timepoint": "BL"}],
                "feature_selections": [
                    {"timepoint": "BL", "features": ["adc"], "lambda": 0.05},
                ],
            },
            "baseline": {
                "total_outliers": {"pct": 8.5},
                "baseline_exclusion": {"n_excluded": 6, "n_total": 48},
            },
            "sanity_checks": {"all_converged": True},
        }
    }


def _make_mat_data():
    return {
        "Standard": {
            "longitudinal": {"num_patients": 42, "num_timepoints": 5},
            "dosimetry": {
                "d95_adc_mean": 45.0,
                "v50_adc_mean": 0.85,
                "d95_d_mean": 44.0,
                "v50_d_mean": 0.80,
            },
        }
    }
