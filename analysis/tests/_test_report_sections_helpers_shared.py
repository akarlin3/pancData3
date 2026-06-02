"""Shared imports for report_sections._helpers tests.

NOT a test module (underscore prefix prevents pytest collection). Imported via
`from _test_report_sections_helpers_shared import *` by the split test files.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections._helpers import (
    _safe_json_load,
    _normalize_series_name,
    _build_normalised_series_map,
    _best_display_name,
    _get_cohort_size,
    _find_best_auc,
    _aggregate_dwi_statistics,
    _compute_feature_overlap,
    _aggregate_sanity_checks,
    _extract_significant_metrics,
    _extract_dosimetry,
    _scalar_gy,
    _compute_cross_dwi_trend_agreement,
    _compute_all_groups_trend_agreement,
    _extract_longitudinal_trend_consensus,
)

# Names are underscore-prefixed, so `import *` skips them unless listed here.
__all__ = [
    "pytest",
    "_safe_json_load",
    "_normalize_series_name",
    "_build_normalised_series_map",
    "_best_display_name",
    "_get_cohort_size",
    "_find_best_auc",
    "_aggregate_dwi_statistics",
    "_compute_feature_overlap",
    "_aggregate_sanity_checks",
    "_extract_significant_metrics",
    "_extract_dosimetry",
    "_scalar_gy",
    "_compute_cross_dwi_trend_agreement",
    "_compute_all_groups_trend_agreement",
    "_extract_longitudinal_trend_consensus",
]
