"""Shared fixtures/helpers for main_results report section tests.

Imported via `from _test_report_sections_main_results_shared import *`
by the split test modules. This module name does NOT start with `test_`
so pytest does not collect it directly.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from conftest import SAMPLE_GRAPH_CSV_ROWS

from report.sections.main_results import (
    _section_executive_summary,
    _section_hypothesis,
    _section_treatment_response,
)
from report.sections.statistical_reporting import (
    _section_statistical_significance,
    _section_broad_statistical_overview,
)
from report.sections.manuscript import (
    _section_predictive_performance,
    _section_manuscript_ready_findings,
    _section_results_draft,
)


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
                    {"timepoint": "BL", "features": ["adc", "d"], "lambda": 0.05},
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
                "d95_adc_mean": {"mean": 45.0},
                "v50_adc_mean": {"mean": 0.85},
            },
        }
    }


def _make_groups():
    return {
        "Longitudinal_Mean_Metrics": {
            "Standard": {
                "trends_json": json.dumps([
                    {"series": "Mean D", "direction": "increasing", "description": "D rises"},
                    {"series": "Mean f", "direction": "decreasing", "description": "f drops"},
                ])
            },
            "dnCNN": {
                "trends_json": json.dumps([
                    {"series": "Mean D", "direction": "increasing", "description": "D rises"},
                    {"series": "Mean f", "direction": "decreasing", "description": "f drops"},
                ])
            },
        },
        "Feature_BoxPlots": {
            "Standard": SAMPLE_GRAPH_CSV_ROWS[0],
            "dnCNN": SAMPLE_GRAPH_CSV_ROWS[1],
        },
    }
