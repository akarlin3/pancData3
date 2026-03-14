"""Tests for report_sections._helpers module.

Validates the shared helper functions used by all report section builders:
- JSON safe loading
- Series name normalisation
- Normalised series map building
- Cohort size extraction
- Best AUC finding
- DWI statistics aggregation
- Feature overlap computation
- Sanity check aggregation
- Significant metrics extraction
- Dosimetry extraction
- Scalar Gy parsing
- Cross-DWI trend agreement
- Longitudinal trend consensus
"""

from __future__ import annotations

import pytest
import sys
from pathlib import Path

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


# ── _safe_json_load ──


class TestSafeJsonLoad:
    def test_valid_list(self):
        assert _safe_json_load('[1, 2, 3]') == [1, 2, 3]

    def test_valid_dict(self):
        assert _safe_json_load('{"a": 1}') == {"a": 1}

    def test_invalid_json_returns_default(self):
        assert _safe_json_load('not json') == []

    def test_scalar_rejected(self):
        assert _safe_json_load('"hello"') == []

    def test_number_rejected(self):
        assert _safe_json_load('42') == []

    def test_none_input(self):
        assert _safe_json_load(None) == []

    def test_custom_default(self):
        assert _safe_json_load('bad', default={}) == {}

    def test_empty_string(self):
        assert _safe_json_load('') == []


# ── _normalize_series_name ──


class TestNormalizeSeriesName:
    def test_basic_normalisation(self):
        result = _normalize_series_name("Mean D")
        # Tokens sorted alphabetically: "d" and "mean"
        assert "d" in result
        assert "mean" in result

    def test_delimiters_collapsed(self):
        # Comma, dash, parentheses should all produce same tokens
        a = _normalize_series_name("Mean D - Local Control")
        b = _normalize_series_name("Mean D (Local Control)")
        c = _normalize_series_name("Mean D, Local Control")
        assert a == b == c

    def test_vs_normalised(self):
        a = _normalize_series_name("LC vs D95")
        b = _normalize_series_name("D95, LC")
        assert a == b

    def test_unicode_delta(self):
        result = _normalize_series_name("\u0394ADC")
        assert "delta" in result.lower()

    def test_empty_string(self):
        assert _normalize_series_name("") == ""

    def test_whitespace_only(self):
        assert _normalize_series_name("   ") == ""


# ── _build_normalised_series_map ──


class TestBuildNormalisedSeriesMap:
    def test_basic_merge(self):
        all_trends = {
            "Standard": [{"series": "Mean D", "direction": "increasing", "description": "goes up"}],
            "dnCNN": [{"series": "Mean D", "direction": "decreasing", "description": "goes down"}],
        }
        result = _build_normalised_series_map(all_trends)
        # Both should merge under the same normalised key
        assert len(result) == 1
        key = list(result.keys())[0]
        assert "Standard" in result[key]
        assert "dnCNN" in result[key]

    def test_non_dict_trends_skipped(self):
        all_trends = {"Standard": ["not a dict", 42]}
        result = _build_normalised_series_map(all_trends)
        assert len(result) == 0

    def test_empty_input(self):
        assert _build_normalised_series_map({}) == {}


# ── _best_display_name ──


class TestBestDisplayName:
    def test_picks_longest(self):
        all_trends = {
            "Standard": [{"series": "D"}],
            "dnCNN": [{"series": "Mean D - Local Control"}],
        }
        nk = _normalize_series_name("Mean D - Local Control")
        result = _best_display_name(all_trends, nk)
        assert result == "Mean D - Local Control"


# ── _get_cohort_size ──


class TestGetCohortSize:
    def test_with_data(self):
        mat_data = {
            "Standard": {
                "longitudinal": {"num_patients": 42, "num_timepoints": 5}
            }
        }
        n_pat, n_tp, dt = _get_cohort_size(mat_data)
        assert n_pat == 42
        assert n_tp == 5
        assert dt == "Standard"

    def test_no_data(self):
        assert _get_cohort_size(None) == (0, 0, "")

    def test_empty_dict(self):
        assert _get_cohort_size({}) == (0, 0, "")


# ── _find_best_auc ──


class TestFindBestAuc:
    def test_finds_best(self):
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "roc_analyses": [
                        {"auc": 0.75, "timepoint": "BL"},
                        {"auc": 0.82, "timepoint": "W2"},
                    ]
                }
            }
        }
        roc, dt = _find_best_auc(log_data, ["Standard"])
        assert roc["auc"] == 0.82
        assert dt == "Standard"

    def test_no_data(self):
        roc, dt = _find_best_auc(None, [])
        assert roc is None
        assert dt == ""

    def test_across_dwi_types(self):
        log_data = {
            "Standard": {"stats_predictive": {"roc_analyses": [{"auc": 0.70}]}},
            "dnCNN": {"stats_predictive": {"roc_analyses": [{"auc": 0.85}]}},
        }
        roc, dt = _find_best_auc(log_data, ["Standard", "dnCNN"])
        assert roc["auc"] == 0.85
        assert dt == "dnCNN"


# ── _aggregate_dwi_statistics ──


class TestAggregateDwiStatistics:
    def test_basic_aggregation(self):
        log_data = {
            "Standard": {
                "stats_predictive": {"roc_analyses": [{"auc": 0.78, "timepoint": "BL"}]},
                "stats_comparisons": {
                    "glme_details": [{"p": 0.01, "adj_alpha": 0.025, "metric": "m1"}]
                },
                "survival": {
                    "hazard_ratios": [{"hr": 1.5, "p": 0.03}]
                },
            }
        }
        stats = _aggregate_dwi_statistics(log_data, ["Standard"])
        assert stats["total_sig"] == 1
        assert stats["total_hrs"] == 1
        assert stats["sig_hrs"] == 1
        assert len(stats["auc_cards"]) == 1

    def test_empty_log(self):
        stats = _aggregate_dwi_statistics(None, [])
        assert stats["total_sig"] == 0


# ── _compute_feature_overlap ──


class TestComputeFeatureOverlap:
    def test_shared_features(self):
        log_data = {
            "Standard": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["ADC_BL", "D_BL"]},
                    ]
                }
            },
            "dnCNN": {
                "stats_predictive": {
                    "feature_selections": [
                        {"timepoint": "BL", "features": ["ADC_BL", "f_BL"]},
                    ]
                }
            },
        }
        shared, total = _compute_feature_overlap(log_data, ["Standard", "dnCNN"])
        assert shared == 1  # ADC_BL
        assert total == 3   # ADC_BL, D_BL, f_BL

    def test_single_type(self):
        shared, total = _compute_feature_overlap({}, ["Standard"])
        assert shared == 0
        assert total == 0


# ── _aggregate_sanity_checks ──


class TestAggregateSanityChecks:
    def test_aggregation(self):
        log_data = {
            "Standard": {
                "sanity_checks": {
                    "all_converged": True,
                    "total_convergence": 0,
                    "dim_mismatches": 1,
                    "nan_dose_warnings": 2,
                }
            }
        }
        result = _aggregate_sanity_checks(log_data, ["Standard"])
        assert result["all_converged_count"] == 1
        assert result["total_conv_flags"] == 0
        assert result["total_dim_issues"] == 3
        assert result["sanity_types_checked"] == 1

    def test_no_data(self):
        result = _aggregate_sanity_checks(None, [])
        assert result["all_converged_count"] == 0


# ── _extract_significant_metrics ──


class TestExtractSignificantMetrics:
    def test_extraction(self):
        log_data = {
            "Standard": {
                "stats_comparisons": {
                    "glme_interactions": [0.01, 0.08],
                    "glme_details": [{"p": 0.01, "adj_alpha": 0.025, "metric": "m1"}],
                    "fdr_timepoints": [{"n_significant": 2, "timepoint": "BL"}],
                },
                "survival": {
                    "hazard_ratios": [{"hr": 2.0, "p": 0.02, "covariate": "adc"}],
                },
                "stats_predictive": {
                    "roc_analyses": [{"auc": 0.80, "timepoint": "BL"}],
                    "feature_selections": [{"timepoint": "BL", "features": ["adc"]}],
                },
            }
        }
        result = _extract_significant_metrics(log_data)
        assert len(result["sig_glme"]) == 1
        assert len(result["sig_glme_details"]) == 1
        assert len(result["sig_fdr_timepoints"]) == 1
        assert len(result["sig_hr"]) == 1
        assert result["best_roc"]["auc"] == 0.80
        assert len(result["all_feature_selections"]) == 1

    def test_no_data(self):
        result = _extract_significant_metrics(None)
        assert len(result["sig_glme"]) == 0


# ── _extract_dosimetry ──


class TestExtractDosimetry:
    def test_found(self):
        mat_data = {"Standard": {"dosimetry": {"d95": 45.0}}}
        dosi, dt = _extract_dosimetry(mat_data)
        assert dosi["d95"] == 45.0
        assert dt == "Standard"

    def test_not_found(self):
        dosi, dt = _extract_dosimetry(None)
        assert dosi == {}
        assert dt == ""


# ── _scalar_gy ──


class TestScalarGy:
    def test_plain_float(self):
        assert _scalar_gy(45.0) == 45.0

    def test_dict_with_mean(self):
        assert _scalar_gy({"mean": 42.5}) == 42.5

    def test_none(self):
        assert _scalar_gy(None) is None

    def test_nan_in_dict(self):
        assert _scalar_gy({"mean": float("nan")}) is None

    def test_nan_plain(self):
        assert _scalar_gy(float("nan")) is None


# ── _compute_cross_dwi_trend_agreement ──


class TestComputeCrossDwiTrendAgreement:
    def test_agreement(self):
        import json
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": json.dumps([{"series": "Mean D", "direction": "increasing"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "Mean D", "direction": "increasing"}])},
            }
        }
        n_agree, n_total, pct = _compute_cross_dwi_trend_agreement(groups, ["Standard", "dnCNN"])
        assert n_agree == 1
        assert n_total == 1
        assert pct == 100.0

    def test_disagreement(self):
        import json
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": json.dumps([{"series": "Mean D", "direction": "increasing"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "Mean D", "direction": "decreasing"}])},
            }
        }
        n_agree, n_total, pct = _compute_cross_dwi_trend_agreement(groups, ["Standard", "dnCNN"])
        assert n_agree == 0
        assert n_total == 1
        assert pct == 0.0

    def test_no_groups(self):
        assert _compute_cross_dwi_trend_agreement(None, []) == (0, 0, 0.0)

    def test_single_dwi_type(self):
        assert _compute_cross_dwi_trend_agreement({}, ["Standard"]) == (0, 0, 0.0)


# ── _compute_all_groups_trend_agreement ──


class TestComputeAllGroupsTrendAgreement:
    def test_across_groups(self):
        import json
        groups = {
            "GroupA": {
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
            },
            "GroupB": {
                "Standard": {"trends_json": json.dumps([{"series": "S2", "direction": "down"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S2", "direction": "up"}])},
            },
        }
        n_agree, n_total, pct = _compute_all_groups_trend_agreement(groups, ["Standard", "dnCNN"])
        assert n_total == 2
        assert n_agree == 1
        assert pct == pytest.approx(50.0)


# ── _extract_longitudinal_trend_consensus ──


class TestExtractLongitudinalTrendConsensus:
    def test_no_data(self):
        d, f, dl, fl = _extract_longitudinal_trend_consensus(None)
        assert d == "unknown"
        assert f == "unknown"

    def test_with_trends(self):
        import json
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {
                    "trends_json": json.dumps([
                        {"series": "Mean D", "direction": "increasing"},
                        {"series": "Mean f", "direction": "decreasing"},
                    ])
                },
                "dnCNN": {
                    "trends_json": json.dumps([
                        {"series": "Mean D", "direction": "increasing"},
                        {"series": "Mean f", "direction": "decreasing"},
                    ])
                },
            }
        }
        d_cons, f_cons, d_list, f_list = _extract_longitudinal_trend_consensus(groups)
        assert len(d_list) == 2
        assert len(f_list) == 2
        assert "increasing" in d_list[0]
