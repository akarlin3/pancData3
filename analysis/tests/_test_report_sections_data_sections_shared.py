"""Shared fixtures/helpers for data_sections report section tests.

Imported via `from _test_report_sections_data_sections_shared import *`
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

from report.sections.enrollment import (
    _section_cohort_overview,
    _section_patient_flow,
    _section_data_completeness,
)
from report.sections.supplemental import (
    _section_mat_data,
)

# Note: actual signatures from source:
# _section_cohort_overview(mat_data, log_data, dwi_types_present)
# _section_patient_flow(log_data, dwi_types_present, mat_data)
# _section_data_completeness(log_data, dwi_types_present)  -- 2 args only
# _section_mat_data(mat_data)  -- 1 arg only
# _section_appendix(rows)  -- 1 arg only
# _section_figure_gallery(folder, rows=None)  -- folder is a Path


def _make_log_data():
    return {
        "Standard": {
            "baseline": {
                "total_outliers": {"n_removed": 4, "n_total": 42, "pct": 9.5},
                "baseline_exclusion": {
                    "n_excluded": 6, "n_total": 48,
                    "lf_rate_included": 35.0, "lf_rate_excluded": 50.0,
                },
            },
            "sanity_checks": {
                "all_converged": True,
                "total_convergence": 0,
                "dim_mismatches": 1,
                "nan_dose_warnings": 0,
            },
            "stats_comparisons": {
                "glme_excluded": {"n_excluded": 5, "n_total": 42, "pct": 12.0},
            },
            "survival": {
                "n_lf": 15,
                "n_lc": 27,
            },
        }
    }


def _make_mat_data():
    return {
        "Standard": {
            "longitudinal": {"num_patients": 42, "num_timepoints": 5},
            "dosimetry": {
                "d95_adc_mean": {"mean": 45.2, "std": 3.1},
                "v50_adc_mean": {"mean": 0.87, "std": 0.08},
                "d95_d_mean": {"mean": 48.0, "std": 2.5},
                "v50_d_mean": {"mean": 0.92, "std": 0.05},
            },
            "core_method": {
                "methods": ["adc_threshold", "otsu", "kmeans"],
                "mean_dice_matrix": [
                    [1.0, 0.75, 0.60],
                    [0.75, 1.0, 0.68],
                    [0.60, 0.68, 1.0],
                ],
            },
        }
    }
