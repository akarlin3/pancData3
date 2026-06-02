#!/usr/bin/env python3
"""Outcome-analysis MAT-file parsers for the analysis suite.

These per-block parsers were extracted from
:mod:`parsers.parse_mat_metrics` to keep that module under the project
file-length limit.  Each function reads one ``.mat`` file for a single DWI
type and mutates the shared ``out_data`` dict in place.  They are
re-imported back into ``parse_mat_metrics`` so its public API is unchanged.

Covers: threshold optimization, baseline-vs-delta predictor comparison,
per-method coefficient of repeatability, sub-volume stability, dose-response
ROC, GTV confounding, and risk-dose concordance.
"""

import sys
import typing
from pathlib import Path

# Ensure analysis/ root is on sys.path so 'parsers' is importable when this
# module is imported from a subprocess-executed sibling.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from parsers.parse_mat_helpers import (  # type: ignore
    METHOD_DESC,
    _array_to_list,
    _nested_safe,
    _safe_float,
    numpy_np,
    scipy_io,
)


def _parse_threshold_optimization(folder, dwi, out_data):
    # ── Threshold Optimization ──
    # The adc_threshold_optimization MAT file contains an opt_results struct
    # with arrays of candidate thresholds and per-threshold Dice statistics.
    opt_mat = folder / f"{dwi}" / f"adc_threshold_optimization_{dwi}.mat"
    if opt_mat.exists() and scipy_io and numpy_np:
        try:
            mat: typing.Any = scipy_io.loadmat(str(opt_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            results = mat.get("opt_results")
            if results:
                opt_entry: dict = {}
                for field in (
                    "thresholds",
                    "median_dice",
                    "mean_dice",
                    "std_dice",
                    "median_vol_frac",
                    "vol_frac_curvature",     # Tactic 2
                    "significance_pvalues",   # Tactic 3
                ):
                    if hasattr(results, field):
                        opt_entry[field] = _array_to_list(getattr(results, field))
                for scalar_field in (
                    "optimal_thresh",
                    "optimal_dice",
                    "optimal_vol_frac",
                    "inflection_thresh",      # Tactic 2
                    "inflection_curvature",   # Tactic 2
                    "significance_thresh",    # Tactic 3
                    "significance_pvalue",    # Tactic 3
                ):
                    if hasattr(results, scalar_field):
                        opt_entry[scalar_field] = _safe_float(getattr(results, scalar_field))
                for int_field in (
                    "n_patients",
                    "inflection_idx",
                    "significance_n_lc",
                    "significance_n_lf",
                ):
                    if hasattr(results, int_field):
                        try:
                            opt_entry[int_field] = int(getattr(results, int_field))
                        except Exception:
                            pass
                if hasattr(results, "significance_metric"):
                    try:
                        opt_entry["significance_metric"] = str(getattr(results, "significance_metric"))
                    except Exception:
                        pass
                # Honest Tactic-3 numbers (nested structs from CP0/CP1):
                # permutation selection-adjusted p and nested-CV out-of-fold.
                for sub_name, sub_fields in (
                    ("permutation", ("perm_adjusted_min_p", "observed_min_p",
                                     "n_notable_naive", "n_notable_adjusted")),
                    ("nested_cv", ("oof_pvalue", "oof_auc", "recommended_thresh",
                                   "modal_fraction", "n_unique_selected")),
                    ("nested_cv_repeated", ("oof_pvalue", "oof_auc",
                                            "recommended_thresh", "modal_fraction",
                                            "n_unique_selected")),
                ):
                    if hasattr(results, sub_name):
                        sub_obj = getattr(results, sub_name)
                        sub_d: dict = {}
                        for sf in sub_fields:
                            if hasattr(sub_obj, sf):
                                sub_d[sf] = _safe_float(getattr(sub_obj, sf))
                        if sub_d:
                            opt_entry[sub_name] = sub_d
                out_data["threshold_optimization"] = opt_entry
        except Exception as e:
            print(f"Error parsing {opt_mat.name}: {e}")



def _parse_baseline_vs_delta(folder, dwi, out_data):
    # ── Baseline vs Delta ──
    # The baseline_vs_delta MAT file contains a comparison struct with a
    # nested results struct array (one entry per parameter/timepoint).
    bvd_mat = folder / f"{dwi}" / f"baseline_vs_delta_{dwi}.mat"
    if bvd_mat.exists() and scipy_io and numpy_np:
        try:
            mat: typing.Any = scipy_io.loadmat(str(bvd_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            comparison = mat.get("comparison")
            bvd_list: list[dict] = []
            if comparison is not None and hasattr(comparison, "results"):
                results_arr = comparison.results
                if results_arr is not None:
                    if not hasattr(results_arr, "__len__") or isinstance(results_arr, str):
                        results_arr = [results_arr]
                    for item in results_arr:
                        entry: dict = {}
                        for field in (
                            "parameter",
                            "timepoint",
                            "better_predictor",
                        ):
                            if hasattr(item, field):
                                entry[field] = str(getattr(item, field))
                        for num_field in (
                            "baseline_hr",
                            "baseline_p",
                            "baseline_cindex",
                            "delta_hr",
                            "delta_p",
                            "delta_cindex",
                        ):
                            if hasattr(item, num_field):
                                entry[num_field] = _safe_float(getattr(item, num_field))
                        if hasattr(item, "n_events"):
                            try:
                                entry["n_events"] = int(getattr(item, "n_events"))
                            except Exception:
                                entry["n_events"] = None
                        if entry:
                            bvd_list.append(entry)
            if bvd_list:
                out_data["baseline_vs_delta"] = bvd_list
        except Exception as e:
            print(f"Error parsing {bvd_mat.name}: {e}")



def _parse_per_method_cor(folder, dwi, out_data):
    # ── Per-Method CoR ──
    cor_mat = folder / f"{dwi}" / f"per_method_cor_{dwi}.mat"
    if cor_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(cor_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            results = mat.get("cor_results")
            if results:
                methods = results.method_names.tolist()
                if isinstance(methods, str):
                    methods = [methods]
                cor_entry: dict = {
                    "method_names": [str(m) for m in methods],
                    "median_wcv": _array_to_list(results.median_wcv) if hasattr(results, "median_wcv") else [],
                    "cor": _array_to_list(results.cor) if hasattr(results, "cor") else [],
                    "n_patients_with_repeats": int(results.n_patients_with_repeats) if hasattr(results, "n_patients_with_repeats") else 0,
                }
                out_data["per_method_cor"] = cor_entry
        except Exception as e:
            print(f"Error parsing {cor_mat.name}: {e}")



def _parse_subvolume_stability(folder, dwi, out_data):
    # ── Sub-Volume Stability ──
    stab_mat = folder / f"{dwi}" / f"subvolume_stability_{dwi}.mat"
    if stab_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(stab_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            results = mat.get("stability")
            if results:
                methods = results.method_names.tolist()
                if isinstance(methods, str):
                    methods = [methods]
                stab_entry: dict = {
                    "method_names": [str(m) for m in methods],
                    "n_patients": int(results.n_patients) if hasattr(results, "n_patients") else 0,
                    "n_timepoints": int(results.n_timepoints) if hasattr(results, "n_timepoints") else 0,
                    "dice_vs_baseline": _array_to_list(results.dice_vs_baseline) if hasattr(results, "dice_vs_baseline") else [],
                }
                out_data["subvolume_stability"] = _nested_safe(stab_entry)
        except Exception as e:
            print(f"Error parsing {stab_mat.name}: {e}")



def _parse_dose_response_roc(folder, dwi, out_data):
    # ── Dose-Response ROC ──
    roc_mat = folder / f"{dwi}" / f"dose_response_roc_{dwi}.mat"
    if roc_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(roc_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            results = mat.get("roc_results")
            if results:
                roc_entry: dict = {"method_results": [], "ranking": []}
                if hasattr(results, "ranking"):
                    ranking = results.ranking
                    if isinstance(ranking, str):
                        ranking = [ranking]
                    else:
                        ranking = [str(r) for r in ranking.tolist()] if hasattr(ranking, "tolist") else list(ranking)
                    roc_entry["ranking"] = ranking
                if hasattr(results, "method_results"):
                    mrs = results.method_results
                    if not hasattr(mrs, "__len__"):
                        mrs = [mrs]
                    for mr in mrs:
                        mr_dict: dict = {
                            "method_name": str(mr.method_name) if hasattr(mr, "method_name") else "",
                            "best_metric": str(mr.best_metric) if hasattr(mr, "best_metric") else "",
                            "best_auc": _safe_float(mr.best_auc) if hasattr(mr, "best_auc") else None,
                            "metrics": [],
                        }
                        if hasattr(mr, "metrics"):
                            mets = mr.metrics
                            if not hasattr(mets, "__len__"):
                                mets = [mets]
                            for met in mets:
                                met_dict: dict = {
                                    "metric_name": str(met.metric_name) if hasattr(met, "metric_name") else "",
                                    "auc": _safe_float(met.auc) if hasattr(met, "auc") else None,
                                    "auc_ci": _array_to_list(met.auc_ci) if hasattr(met, "auc_ci") else [],
                                    "optimal_threshold": _safe_float(met.optimal_threshold) if hasattr(met, "optimal_threshold") else None,
                                    "sensitivity": _safe_float(met.sensitivity) if hasattr(met, "sensitivity") else None,
                                    "specificity": _safe_float(met.specificity) if hasattr(met, "specificity") else None,
                                }
                                mr_dict["metrics"].append(met_dict)
                        roc_entry["method_results"].append(mr_dict)
                out_data["dose_response_roc"] = _nested_safe(roc_entry)
        except Exception as e:
            print(f"Error parsing {roc_mat.name}: {e}")



def _parse_gtv_confounding(folder, dwi, out_data):
    # ── GTV Confounding ──
    gtv_mat = folder / f"{dwi}" / f"gtv_confounding_{dwi}.mat"
    if gtv_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(gtv_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            results = mat.get("confound")
            if results:
                gtv_entry: dict = {"method_results": [], "summary": ""}
                if hasattr(results, "summary"):
                    gtv_entry["summary"] = str(results.summary)
                if hasattr(results, "method_results"):
                    mrs = results.method_results
                    if not hasattr(mrs, "__len__"):
                        mrs = [mrs]
                    for mr in mrs:
                        mr_dict = {
                            "method_name": str(mr.method_name) if hasattr(mr, "method_name") else "",
                            "d95_gtv_correlation": _safe_float(mr.d95_gtv_correlation) if hasattr(mr, "d95_gtv_correlation") else None,
                            "d95_gtv_pvalue": _safe_float(mr.d95_gtv_pvalue) if hasattr(mr, "d95_gtv_pvalue") else None,
                            "d95_gtv_correlation_full": _safe_float(mr.d95_gtv_correlation_full) if hasattr(mr, "d95_gtv_correlation_full") else None,
                            "d95_gtv_pvalue_full": _safe_float(mr.d95_gtv_pvalue_full) if hasattr(mr, "d95_gtv_pvalue_full") else None,
                            "unadjusted_hr": _safe_float(mr.unadjusted_hr) if hasattr(mr, "unadjusted_hr") else None,
                            "adjusted_hr": _safe_float(mr.adjusted_hr) if hasattr(mr, "adjusted_hr") else None,
                            "confounding_flag": bool(mr.confounding_flag) if hasattr(mr, "confounding_flag") else False,
                        }
                        gtv_entry["method_results"].append(mr_dict)
                out_data["gtv_confounding"] = _nested_safe(gtv_entry)
        except Exception as e:
            print(f"Error parsing {gtv_mat.name}: {e}")



def _parse_risk_dose_concordance(folder, dwi, out_data):
    # ── Risk-Dose Concordance ──
    rdc_mat = folder / f"{dwi}" / f"risk_dose_concordance_{dwi}.mat"
    if rdc_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(rdc_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            results = mat.get("concordance")
            if results:
                rdc_entry: dict = {"method_results": [], "summary": ""}
                if hasattr(results, "summary"):
                    rdc_entry["summary"] = str(results.summary)
                if hasattr(results, "method_results"):
                    mrs = results.method_results
                    if not hasattr(mrs, "__len__"):
                        mrs = [mrs]
                    for mr in mrs:
                        mr_dict = {
                            "method_name": str(mr.method_name) if hasattr(mr, "method_name") else "",
                            "best_dose_metric": str(mr.best_dose_metric) if hasattr(mr, "best_dose_metric") else "",
                            "cohen_kappa": _safe_float(mr.cohen_kappa) if hasattr(mr, "cohen_kappa") else None,
                            "concordance_pct": _safe_float(mr.concordance_pct) if hasattr(mr, "concordance_pct") else None,
                            "n_complementary": _safe_float(mr.n_complementary) if hasattr(mr, "n_complementary") else None,
                            "combined_auc": _safe_float(mr.combined_auc) if hasattr(mr, "combined_auc") else None,
                        }
                        rdc_entry["method_results"].append(mr_dict)
                out_data["risk_dose_concordance"] = _nested_safe(rdc_entry)
        except Exception as e:
            print(f"Error parsing {rdc_mat.name}: {e}")

