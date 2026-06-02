#!/usr/bin/env python3
"""Core-method MAT-file parsers for the analysis suite.

These per-block parsers were extracted from
:mod:`parsers.parse_mat_metrics` to keep that module under the project
file-length limit.  Each function reads one ``.mat`` file for a single DWI
type and mutates the shared ``out_data`` dict in place.  They are
re-imported back into ``parse_mat_metrics`` so its public API is unchanged.

Covers: core-method Dice/Hausdorff comparison, per-method outcomes,
core-method pruning, failure rates, and cross-pipeline Dice.
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


def _parse_core_method(folder, dwi, out_data):
    # ── Core Method Comparison ──
    # The compare_core_results MAT file contains a struct with:
    #   - method_names: cell array of method label strings
    #   - mean_dice_matrix: NxN matrix of pairwise mean Dice coefficients
    #   - hausdorff_matrix / mean_hausdorff_matrix: optional NxN Hausdorff matrix
    core_mat = folder / f"{dwi}" / f"compare_core_results_{dwi}.mat"
    if core_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(core_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            results = mat.get("compare_results")
            if results:
                # Extract method names; squeeze_me may collapse single-element
                # cell arrays to a bare string.
                methods = results.method_names.tolist()
                if isinstance(methods, str):
                    methods = [methods]
                method_strs = [str(m) for m in methods]
                n_methods = len(method_strs)

                # Build method descriptions from the MAT file (if present) or
                # fall back to the built-in METHOD_DESC lookup table.
                descriptions: list[str] = []
                if hasattr(results, "method_descriptions"):
                    raw_descs = results.method_descriptions
                    try:
                        desc_list = raw_descs.tolist()
                        if isinstance(desc_list, str):
                            desc_list = [desc_list]
                        descriptions = [str(d) for d in desc_list]
                    except Exception:
                        descriptions = []
                if not descriptions:
                    descriptions = [METHOD_DESC.get(m, m) for m in method_strs]

                mean_dice = results.mean_dice_matrix

                # Validate the Dice matrix is square and matches method count.
                dice_list: list = []
                if hasattr(mean_dice, "shape"):
                    rows_d, cols_d = (mean_dice.shape[0], mean_dice.shape[1]) if mean_dice.ndim >= 2 else (0, 0)
                    if rows_d != n_methods or cols_d != n_methods:
                        print(
                            f"Warning: mean_dice_matrix shape {mean_dice.shape} does not match "
                            f"{n_methods} methods. Skipping dice matrix."
                        )
                    else:
                        dice_list = _array_to_list(mean_dice)

                core_entry: dict = {
                    "methods": method_strs,
                    "method_descriptions": descriptions,
                    "mean_dice_matrix": dice_list,
                }

                # Extract Hausdorff matrix if present (either field name).
                for hd_field in ("hausdorff_matrix", "mean_hausdorff_matrix"):
                    if hasattr(results, hd_field):
                        hd = getattr(results, hd_field)
                        hd_list: list = []
                        if hasattr(hd, "shape"):
                            rows_h, cols_h = (hd.shape[0], hd.shape[1]) if hd.ndim >= 2 else (0, 0)
                            if rows_h != n_methods or cols_h != n_methods:
                                print(
                                    f"Warning: {hd_field} shape {hd.shape} does not match "
                                    f"{n_methods} methods. Skipping."
                                )
                            else:
                                hd_list = _array_to_list(hd)
                        core_entry["mean_hausdorff_matrix"] = hd_list
                        break  # Use the first matching field only.

                out_data["core_method"] = core_entry
        except Exception as e:
            print(f"Error parsing {core_mat.name}: {e}")



def _parse_core_method_outcomes(folder, dwi, out_data):
    # ── Core Method Outcomes ──
    # The core_method_outcomes MAT file contains:
    #   - outcome_results.method_results: struct array with per-method Cox/KM results
    #   - outcome_results.ranking: cell array of method names sorted by p-value
    #   - outcome_results.active_methods: cell array of analyzed methods
    cmo_mat = folder / f"{dwi}" / f"core_method_outcomes_{dwi}.mat"
    if cmo_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(cmo_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            results = mat.get("outcome_results")
            if results:
                cmo_entry: dict = {}

                # Parse active_methods
                am = getattr(results, "active_methods", None)
                if am is not None:
                    am_list = am.tolist() if hasattr(am, "tolist") else am
                    if isinstance(am_list, str):
                        am_list = [am_list]
                    cmo_entry["active_methods"] = [str(m) for m in am_list]

                # Parse ranking
                rk = getattr(results, "ranking", None)
                if rk is not None:
                    rk_list = rk.tolist() if hasattr(rk, "tolist") else rk
                    if isinstance(rk_list, str):
                        rk_list = [rk_list]
                    cmo_entry["ranking"] = [str(r) for r in rk_list]
                else:
                    cmo_entry["ranking"] = []

                # Parse method_results
                mr_raw = getattr(results, "method_results", None)
                method_list: list[dict] = []
                if mr_raw is not None:
                    items = mr_raw if hasattr(mr_raw, "__len__") and not isinstance(mr_raw, str) else [mr_raw]
                    for item in items:
                        me: dict = {}
                        me["method_name"] = str(getattr(item, "method_name", ""))
                        me["n_patients"] = int(getattr(item, "n_patients", 0))
                        me["n_events"] = int(getattr(item, "n_events", 0))

                        # Parse univariable results
                        uv_raw = getattr(item, "univariable", None)
                        uv_list: list[dict] = []
                        if uv_raw is not None and uv_raw is not None:
                            uv_items = uv_raw if hasattr(uv_raw, "__len__") and not isinstance(uv_raw, str) else [uv_raw]
                            for uv in uv_items:
                                uv_entry: dict = {}
                                uv_entry["metric_name"] = str(getattr(uv, "metric_name", ""))
                                uv_entry["hr"] = _safe_float(getattr(uv, "hr", None))
                                ci = getattr(uv, "hr_ci", None)
                                if ci is not None and hasattr(ci, "tolist"):
                                    ci_list = ci.tolist()
                                    uv_entry["hr_ci"] = [_safe_float(c) for c in ci_list] if isinstance(ci_list, list) else []
                                else:
                                    uv_entry["hr_ci"] = []
                                uv_entry["p_value"] = _safe_float(getattr(uv, "p_value", None))
                                uv_entry["n"] = int(getattr(uv, "n", 0)) if hasattr(uv, "n") else 0
                                uv_entry["n_events"] = int(getattr(uv, "n_events", 0)) if hasattr(uv, "n_events") else 0
                                uv_list.append(uv_entry)
                        me["univariable"] = uv_list

                        # Parse KM results
                        km_raw = getattr(item, "km", None)
                        km_entry: dict = {}
                        if km_raw is not None:
                            km_entry["best_metric"] = str(getattr(km_raw, "best_metric", ""))
                            km_entry["logrank_p"] = _safe_float(getattr(km_raw, "logrank_p", None))
                            km_entry["logrank_chi2"] = _safe_float(getattr(km_raw, "logrank_chi2", None))
                            km_entry["median_high"] = _safe_float(getattr(km_raw, "median_high", None))
                            km_entry["median_low"] = _safe_float(getattr(km_raw, "median_low", None))
                        me["km"] = km_entry

                        method_list.append(me)
                cmo_entry["method_results"] = method_list

                out_data["core_method_outcomes"] = cmo_entry
        except Exception as e:
            print(f"Error parsing {cmo_mat.name}: {e}")



def _parse_pruning(folder, dwi, out_data):
    # ── Core Method Pruning ──
    # The core_pruning MAT file contains:
    #   - active_methods: cell array of retained method name strings
    #   - pruned_info: struct array with .name, .reason, .failure_rate, .pipeline
    pruning_mat = folder / f"{dwi}" / f"core_pruning_{dwi}.mat"
    if pruning_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(pruning_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            pruning_entry: dict = {}

            # Parse active_methods
            am = mat.get("active_methods")
            if am is not None:
                am_list = am.tolist() if hasattr(am, "tolist") else am
                if isinstance(am_list, str):
                    am_list = [am_list]
                pruning_entry["active_methods"] = [str(m) for m in am_list]
            else:
                pruning_entry["active_methods"] = []

            # Parse pruned_info struct array
            pi_raw = mat.get("pruned_info")
            pruned_list: list[dict] = []
            if pi_raw is not None:
                # Handle both single struct and array of structs
                if hasattr(pi_raw, "__len__") and not isinstance(pi_raw, str):
                    items = pi_raw if hasattr(pi_raw, "__iter__") else [pi_raw]
                else:
                    items = [pi_raw]
                for item in items:
                    entry: dict = {}
                    if hasattr(item, "name"):
                        entry["name"] = str(item.name)
                    if hasattr(item, "reason"):
                        entry["reason"] = str(item.reason)
                    if hasattr(item, "failure_rate"):
                        entry["failure_rate"] = _safe_float(item.failure_rate)
                    if hasattr(item, "pipeline"):
                        entry["pipeline"] = str(item.pipeline)
                    if entry:
                        pruned_list.append(entry)
            pruning_entry["pruned_info"] = pruned_list

            # Parse retained_with_warning struct array
            rw_raw = mat.get("retained_with_warning")
            rw_list: list[dict] = []
            if rw_raw is not None:
                if hasattr(rw_raw, "__len__") and not isinstance(rw_raw, str):
                    rw_items = rw_raw if hasattr(rw_raw, "__iter__") else [rw_raw]
                else:
                    rw_items = [rw_raw]
                for item in rw_items:
                    rw_entry: dict = {}
                    if hasattr(item, "name"):
                        rw_entry["name"] = str(item.name)
                    if hasattr(item, "failure_rate"):
                        rw_entry["failure_rate"] = _safe_float(item.failure_rate)
                    if hasattr(item, "reason"):
                        rw_entry["reason"] = str(item.reason)
                    if rw_entry:
                        rw_list.append(rw_entry)
            pruning_entry["retained_with_warning"] = rw_list

            # Parse min_core_voxels_used
            mcv = mat.get("min_core_voxels_used")
            if mcv is not None:
                pruning_entry["min_core_voxels"] = int(_safe_float(mcv))

            out_data["pruning"] = pruning_entry
        except Exception as e:
            print(f"Error parsing {pruning_mat.name}: {e}")



def _parse_failure_rates(folder, dwi, out_data):
    # ── Core Method Failure Rates ──
    # The core_failure_rates MAT file contains a struct with:
    #   - method_names: cell array of 11 method label strings
    #   - pipeline_names: cell array of 3 pipeline label strings
    #   - fallback_rate, empty_rate, insufficient_rate, all_nan_rate,
    #     any_failure_rate: [11 x 3] rate matrices
    #   - median_core_voxels: [11 x 3] median voxel counts
    cfr_mat = folder / f"{dwi}" / f"core_failure_rates_{dwi}.mat"
    if cfr_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(cfr_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            results = mat.get("failure_table")
            if results:
                methods = results.method_names.tolist()
                if isinstance(methods, str):
                    methods = [methods]
                method_strs = [str(m) for m in methods]

                pipelines = results.pipeline_names.tolist()
                if isinstance(pipelines, str):
                    pipelines = [pipelines]
                pipeline_strs = [str(p) for p in pipelines]

                cfr_entry: dict = {
                    "method_names": method_strs,
                    "pipeline_names": pipeline_strs,
                }
                for rate_field in ("fallback_rate", "empty_rate", "insufficient_rate",
                                   "all_nan_rate", "any_failure_rate", "median_core_voxels"):
                    if hasattr(results, rate_field):
                        cfr_entry[rate_field] = _array_to_list(getattr(results, rate_field))

                out_data["failure_rates"] = cfr_entry
        except Exception as e:
            print(f"Error parsing {cfr_mat.name}: {e}")



def _parse_cross_pipeline_dice(folder, dwi, out_data):
    # ── Cross-Pipeline Dice ──
    # The cross_pipeline_dice MAT file contains a struct with:
    #   - dice: [11 x 3 x N] array of Dice coefficients
    #   - method_names: cell array of 11 method label strings
    #   - pipeline_pair_labels: cell array of 3 pair label strings
    #   - n_voxels: [11 x 3 x N] voxel counts
    #   - fallback_flags: [11 x 3 x N] logical fallback indicators
    cpd_mat = folder / f"{dwi}" / f"cross_pipeline_dice_{dwi}.mat"
    if cpd_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(cpd_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            results = mat.get("dice_results")
            if results:
                # Extract method names
                methods = results.method_names.tolist()
                if isinstance(methods, str):
                    methods = [methods]
                method_strs = [str(m) for m in methods]

                # Extract pair labels
                pair_labels = results.pipeline_pair_labels.tolist()
                if isinstance(pair_labels, str):
                    pair_labels = [pair_labels]
                pair_strs = [str(p) for p in pair_labels]

                # Extract dice array and compute mean/std across patients (axis 2)
                dice_arr = results.dice
                if hasattr(dice_arr, "shape") and dice_arr.ndim == 3:
                    mean_dice = _array_to_list(numpy_np.nanmean(dice_arr, axis=2))  # type: ignore
                    std_dice = _array_to_list(numpy_np.nanstd(dice_arr, axis=2))  # type: ignore
                elif hasattr(dice_arr, "shape") and dice_arr.ndim == 2:
                    # Single patient: no std
                    mean_dice = _array_to_list(dice_arr)
                    std_dice = None
                else:
                    mean_dice = []
                    std_dice = None

                cpd_entry: dict = {
                    "methods": method_strs,
                    "pair_labels": pair_strs,
                    "mean_dice": mean_dice,
                    "std_dice": std_dice,
                }

                # Extract n_voxels and fallback_flags if present
                if hasattr(results, "n_voxels"):
                    cpd_entry["n_voxels"] = _array_to_list(results.n_voxels)
                if hasattr(results, "fallback_flags"):
                    cpd_entry["fallback_flags"] = _array_to_list(results.fallback_flags)

                out_data["cross_pipeline_dice"] = cpd_entry
        except Exception as e:
            print(f"Error parsing {cpd_mat.name}: {e}")

