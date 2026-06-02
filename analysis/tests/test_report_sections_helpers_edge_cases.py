"""Edge-case tests for report_sections._helpers module.

Covers empty / partial / missing nested structures across every helper:
safe JSON load, series normalisation, normalised series map, best display
name, cohort size, best AUC, DWI statistics aggregation, feature overlap,
sanity-check aggregation, significant metrics, dosimetry, scalar Gy, and
cross-DWI / longitudinal trend agreement.

Base-behaviour tests live in test_report_sections_helpers.py. Shared imports
live in _test_report_sections_helpers_shared.py.
"""

from __future__ import annotations

from _test_report_sections_helpers_shared import *  # noqa: F401,F403


# ── Edge case: empty / partial / missing nested structures ──


class TestSafeJsonLoadEdgeCases:
    def test_null_json_literal(self):
        assert _safe_json_load("null") == []

    def test_boolean_json(self):
        assert _safe_json_load("true") == []

    def test_nested_list(self):
        result = _safe_json_load('[[1, 2], [3, 4]]')
        assert result == [[1, 2], [3, 4]]

    def test_deeply_nested_dict(self):
        result = _safe_json_load('{"a": {"b": {"c": 1}}}')
        assert result == {"a": {"b": {"c": 1}}}

    def test_unicode_json(self):
        result = _safe_json_load('{"key": "\\u0394ADC"}')
        assert result["key"] == "\u0394ADC"


class TestNormalizeSeriesNameEdgeCases:
    def test_multiple_delimiters(self):
        result = _normalize_series_name("A--B,,C")
        assert "a" in result
        assert "b" in result
        assert "c" in result

    def test_times_symbol(self):
        result = _normalize_series_name("\u00d7 factor")
        assert "x" in result

    def test_single_token(self):
        result = _normalize_series_name("ADC")
        assert result == "adc"

    def test_case_insensitivity(self):
        assert _normalize_series_name("MEAN D") == _normalize_series_name("mean d")


class TestBuildNormalisedSeriesMapEdgeCases:
    def test_missing_series_key_uses_overall(self):
        all_trends = {
            "Standard": [{"direction": "increasing", "description": "goes up"}],
        }
        result = _build_normalised_series_map(all_trends)
        assert len(result) == 1
        key = list(result.keys())[0]
        assert "overall" in key

    def test_missing_direction_key(self):
        all_trends = {
            "Standard": [{"series": "Mean D"}],
        }
        result = _build_normalised_series_map(all_trends)
        assert len(result) == 1
        key = list(result.keys())[0]
        dir_val, _ = result[key]["Standard"]
        assert dir_val == ""

    def test_trends_not_a_list(self):
        all_trends = {"Standard": "not a list"}
        result = _build_normalised_series_map(all_trends)
        assert result == {}

    def test_duplicate_normalized_key_uses_first(self):
        all_trends = {
            "Standard": [
                {"series": "Mean D", "direction": "increasing", "description": "A"},
                {"series": "D Mean", "direction": "decreasing", "description": "B"},
            ],
        }
        result = _build_normalised_series_map(all_trends)
        key = list(result.keys())[0]
        dir_val, _ = result[key]["Standard"]
        assert dir_val == "increasing"


class TestBestDisplayNameEdgeCases:
    def test_not_found_returns_key(self):
        result = _best_display_name({"Standard": [{"series": "X"}]}, "nonexistent")
        assert result == "nonexistent"

    def test_empty_trends(self):
        result = _best_display_name({}, "key")
        assert result == "key"


class TestGetCohortSizeEdgeCases:
    def test_missing_longitudinal_key(self):
        mat = {"Standard": {"dosimetry": {}}}
        n_pat, n_tp, dt = _get_cohort_size(mat)
        assert n_pat == 0

    def test_zero_patients_skipped(self):
        mat = {
            "Standard": {"longitudinal": {"num_patients": 0, "num_timepoints": 5}},
            "dnCNN": {"longitudinal": {"num_patients": 30, "num_timepoints": 3}},
        }
        n_pat, n_tp, dt = _get_cohort_size(mat)
        assert n_pat == 30
        assert dt == "dnCNN"

    def test_first_dwi_with_patients_wins(self):
        mat = {
            "Standard": {"longitudinal": {"num_patients": 42, "num_timepoints": 5}},
            "dnCNN": {"longitudinal": {"num_patients": 50, "num_timepoints": 3}},
        }
        n_pat, n_tp, dt = _get_cohort_size(mat)
        assert n_pat == 42


class TestFindBestAucEdgeCases:
    def test_missing_auc_key_uses_zero(self):
        log = {"Standard": {"stats_predictive": {"roc_analyses": [{"timepoint": "BL"}]}}}
        roc, dt = _find_best_auc(log, ["Standard"])
        assert roc is not None
        assert dt == "Standard"

    def test_empty_roc_list(self):
        log = {"Standard": {"stats_predictive": {"roc_analyses": []}}}
        roc, dt = _find_best_auc(log, ["Standard"])
        assert roc is None

    def test_missing_stats_predictive(self):
        log = {"Standard": {"baseline": {}}}
        roc, dt = _find_best_auc(log, ["Standard"])
        assert roc is None

    def test_dwi_type_not_in_log(self):
        log = {"Standard": {"stats_predictive": {"roc_analyses": [{"auc": 0.8}]}}}
        roc, dt = _find_best_auc(log, ["dnCNN"])
        assert roc is None


class TestAggregateDwiStatisticsEdgeCases:
    def test_multiple_dwi_types(self):
        log = {
            "Standard": {
                "stats_predictive": {"roc_analyses": [{"auc": 0.78, "timepoint": "BL"}]},
                "stats_comparisons": {"glme_details": [{"p": 0.01, "adj_alpha": 0.025, "metric": "m1"}]},
                "survival": {"hazard_ratios": [{"hr": 1.5, "p": 0.03}]},
            },
            "dnCNN": {
                "stats_predictive": {"roc_analyses": [{"auc": 0.72, "timepoint": "BL"}]},
                "stats_comparisons": {"glme_details": [{"p": 0.5, "adj_alpha": 0.05, "metric": "m2"}]},
                "survival": {"hazard_ratios": [{"hr": 1.1, "p": 0.5}]},
            },
        }
        stats = _aggregate_dwi_statistics(log, ["Standard", "dnCNN"])
        assert stats["total_sig"] >= 1
        assert stats["total_hrs"] == 2

    def test_missing_survival_key(self):
        log = {"Standard": {
            "stats_comparisons": {"glme_details": [{"p": 0.01, "adj_alpha": 0.025, "metric": "m1"}]},
            "stats_predictive": {"roc_analyses": []},
        }}
        stats = _aggregate_dwi_statistics(log, ["Standard"])
        assert stats["total_hrs"] == 0
        assert stats["sig_hrs"] == 0

    def test_missing_glme_details(self):
        log = {"Standard": {
            "stats_comparisons": {},
            "stats_predictive": {"roc_analyses": []},
            "survival": {"hazard_ratios": []},
        }}
        stats = _aggregate_dwi_statistics(log, ["Standard"])
        assert stats["total_sig"] == 0

    def test_csv_data_passthrough(self):
        csv_data = {"significant_metrics": {"Standard": [{"Metric": "x"}]}}
        stats = _aggregate_dwi_statistics(None, [], csv_data)
        assert stats["total_sig"] == 0


class TestComputeFeatureOverlapEdgeCases:
    def test_no_log_data(self):
        shared, total = _compute_feature_overlap(None, ["Standard", "dnCNN"])
        assert shared == 0 and total == 0

    def test_empty_feature_lists(self):
        log = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": []},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": []},
            ]}},
        }
        shared, total = _compute_feature_overlap(log, ["Standard", "dnCNN"])
        assert shared == 0 and total == 0

    def test_no_overlapping_timepoints(self):
        log = {
            "Standard": {"stats_predictive": {"feature_selections": [
                {"timepoint": "BL", "features": ["adc"]},
            ]}},
            "dnCNN": {"stats_predictive": {"feature_selections": [
                {"timepoint": "W2", "features": ["adc"]},
            ]}},
        }
        shared, total = _compute_feature_overlap(log, ["Standard", "dnCNN"])
        assert shared == 0


class TestAggregateSanityChecksEdgeCases:
    def test_not_converged(self):
        log = {"Standard": {"sanity_checks": {
            "all_converged": False,
            "total_convergence": 7,
            "dim_mismatches": 2,
            "nan_dose_warnings": 3,
        }}}
        result = _aggregate_sanity_checks(log, ["Standard"])
        assert result["all_converged_count"] == 0
        assert result["total_conv_flags"] == 7
        assert result["total_dim_issues"] == 5

    def test_missing_sanity_checks_key(self):
        log = {"Standard": {"baseline": {}}}
        result = _aggregate_sanity_checks(log, ["Standard"])
        assert result["sanity_types_checked"] == 0

    def test_multiple_dwi_types(self):
        log = {
            "Standard": {"sanity_checks": {
                "all_converged": True, "total_convergence": 0,
                "dim_mismatches": 0, "nan_dose_warnings": 0,
            }},
            "dnCNN": {"sanity_checks": {
                "all_converged": False, "total_convergence": 3,
                "dim_mismatches": 1, "nan_dose_warnings": 1,
            }},
        }
        result = _aggregate_sanity_checks(log, ["Standard", "dnCNN"])
        assert result["sanity_types_checked"] == 2
        assert result["all_converged_count"] == 1
        assert result["total_conv_flags"] == 3


class TestExtractSignificantMetricsEdgeCases:
    def test_empty_glme_interactions(self):
        log = {"Standard": {
            "stats_comparisons": {
                "glme_interactions": [],
                "glme_details": [],
                "fdr_timepoints": [],
            },
            "survival": {"hazard_ratios": []},
            "stats_predictive": {"roc_analyses": [], "feature_selections": []},
        }}
        result = _extract_significant_metrics(log)
        assert len(result["sig_glme"]) == 0
        assert result["best_roc"] is None

    def test_all_nonsig_glme(self):
        log = {"Standard": {
            "stats_comparisons": {
                "glme_interactions": [0.5, 0.8],
                "glme_details": [{"p": 0.5, "adj_alpha": 0.05, "metric": "m1"}],
                "fdr_timepoints": [],
            },
            "survival": {"hazard_ratios": []},
            "stats_predictive": {"roc_analyses": [], "feature_selections": []},
        }}
        result = _extract_significant_metrics(log)
        assert len(result["sig_glme"]) == 0
        assert len(result["sig_glme_details"]) == 0

    def test_all_nonsig_hr(self):
        log = {"Standard": {
            "stats_comparisons": {"glme_interactions": [], "glme_details": [], "fdr_timepoints": []},
            "survival": {"hazard_ratios": [
                {"hr": 1.0, "p": 0.8, "covariate": "x"},
            ]},
            "stats_predictive": {"roc_analyses": [], "feature_selections": []},
        }}
        result = _extract_significant_metrics(log)
        assert len(result["sig_hr"]) == 0


class TestExtractDosimetryEdgeCases:
    def test_empty_mat_data(self):
        dosi, dt = _extract_dosimetry({})
        assert dosi == {} and dt == ""

    def test_dosimetry_is_none(self):
        mat = {"Standard": {"dosimetry": None}}
        dosi, dt = _extract_dosimetry(mat)
        assert dosi == {} or dt == ""

    def test_multiple_dwi_types_first_with_dosi(self):
        mat = {
            "Standard": {},
            "dnCNN": {"dosimetry": {"d95": 48.0}},
        }
        dosi, dt = _extract_dosimetry(mat)
        assert dt == "dnCNN"
        assert dosi["d95"] == 48.0


class TestScalarGyEdgeCases:
    def test_string_value(self):
        assert _scalar_gy("not a number") is None

    def test_dict_with_none_mean(self):
        assert _scalar_gy({"mean": None}) is None

    def test_zero_value(self):
        assert _scalar_gy(0.0) == 0.0

    def test_negative_value(self):
        assert _scalar_gy(-5.0) == -5.0


class TestCrossDwiTrendAgreementEdgeCases:
    def test_missing_longitudinal_key(self):
        import json
        groups = {
            "Feature_BoxPlots": {
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
            }
        }
        n_agree, n_total, pct = _compute_cross_dwi_trend_agreement(groups, ["Standard", "dnCNN"])
        assert n_total == 0

    def test_empty_trends_in_groups(self):
        import json
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": "[]"},
                "dnCNN": {"trends_json": "[]"},
            }
        }
        n_agree, n_total, pct = _compute_cross_dwi_trend_agreement(groups, ["Standard", "dnCNN"])
        assert n_total == 0
        assert pct == 0.0


class TestAllGroupsTrendAgreementEdgeCases:
    def test_single_dwi_type(self):
        import json
        groups = {
            "G": {"Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])}},
        }
        n_agree, n_total, pct = _compute_all_groups_trend_agreement(groups, ["Standard"])
        assert n_total == 0

    def test_no_groups(self):
        n_agree, n_total, pct = _compute_all_groups_trend_agreement(None, ["Standard", "dnCNN"])
        assert n_total == 0

    def test_root_key_filtered(self):
        import json
        groups = {
            "G": {
                "Root": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
            }
        }
        n_agree, n_total, pct = _compute_all_groups_trend_agreement(groups, ["Standard"])
        assert n_total == 0

    def test_invalid_trends_json(self):
        import json
        groups = {
            "G": {
                "Standard": {"trends_json": "NOT JSON"},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
            }
        }
        n_agree, n_total, pct = _compute_all_groups_trend_agreement(groups, ["Standard", "dnCNN"])
        assert isinstance(n_total, int)


class TestLongitudinalTrendConsensusEdgeCases:
    def test_only_d_trends(self):
        import json
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": json.dumps([
                    {"series": "Mean D", "direction": "increasing"},
                ])},
            }
        }
        d_cons, f_cons, d_list, f_list = _extract_longitudinal_trend_consensus(groups)
        assert d_cons != "unknown"
        assert f_cons == "unknown"
        assert len(f_list) == 0

    def test_only_f_trends(self):
        import json
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": json.dumps([
                    {"series": "Mean f", "direction": "decreasing"},
                ])},
            }
        }
        d_cons, f_cons, d_list, f_list = _extract_longitudinal_trend_consensus(groups)
        assert d_cons == "unknown"
        assert f_cons != "unknown"

    def test_mixed_d_directions(self):
        import json
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": json.dumps([
                    {"series": "Mean D", "direction": "increasing"},
                ])},
                "dnCNN": {"trends_json": json.dumps([
                    {"series": "Mean D", "direction": "decreasing"},
                ])},
            }
        }
        d_cons, f_cons, d_list, f_list = _extract_longitudinal_trend_consensus(groups)
        assert len(d_list) == 2

    def test_invalid_json_in_trends(self):
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {"trends_json": "INVALID"},
            }
        }
        d_cons, f_cons, d_list, f_list = _extract_longitudinal_trend_consensus(groups)
        assert d_cons == "unknown"
