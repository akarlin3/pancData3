#!/usr/bin/env python3
"""Parse metrics from MATLAB ``.mat`` files for the analysis suite.

This script reads binary MATLAB ``.mat`` files produced by the pipeline
and converts them to JSON for consumption by :mod:`generate_report`.

**Input files** (per DWI type subfolder):

- ``compare_core_results_{dwi}.mat`` -- pairwise Dice/Hausdorff matrix
  from the 11-method core comparison step.
- ``metrics_dosimetry_results_{dwi}.mat`` -- D95 and V50 dose coverage
  arrays for ADC and D sub-volumes.
- ``summary_metrics_{dwi}.mat`` -- per-patient summary metric struct
  with ADC_abs, D_abs, etc.  Used here only to extract cohort dimensions
  (number of patients and timepoints).

**Output files** (written to the output folder root):

- ``parsed_mat_metrics_{dwi}.json`` -- one per DWI type, containing
  ``core_method``, ``dosimetry``, and ``longitudinal`` keys.

**Dependencies**: Requires ``scipy`` and ``numpy`` (``pip install scipy numpy``).
If not available, the script returns empty dicts without crashing.

Usage:
    python parse_mat_metrics.py /path/to/saved_files_*
"""

import argparse
import json
import math
import os
from pathlib import Path

import typing

from tqdm import tqdm  # type: ignore

# Attempt to import scipy and numpy.  These are optional dependencies --
# if not installed, the script degrades gracefully (returns empty data).
try:
    import scipy.io  # type: ignore
    import numpy as np  # type: ignore
    scipy_io: typing.Any = scipy.io
    numpy_np: typing.Any = np
except ImportError:
    print("Warning: scipy or numpy not found. Please pip install scipy numpy.")
    scipy_io = None
    numpy_np = None


# Known method name → human-readable description mapping.
METHOD_DESC: dict[str, str] = {
    "adc_threshold": "ADC threshold (< threshold mm\u00b2/s)",
    "d_threshold": "D threshold (< threshold mm\u00b2/s)",
    "df_intersection": "D\u2013f intersection (D low AND f high)",
    "otsu": "Otsu automatic thresholding",
    "gmm": "Gaussian mixture model (2-component)",
    "kmeans": "K-means clustering (k=2)",
    "region_growing": "Region growing from seed voxels",
    "active_contours": "Active contour (snake) segmentation",
    "percentile": "Percentile-based threshold",
    "spectral": "Spectral clustering",
    "fdm": "Functional diffusion map",
}


def _safe_float(v: typing.Any) -> typing.Optional[float]:
    """Convert a value to a JSON-safe float, returning None for NaN/Inf.

    Parameters
    ----------
    v : any
        Input value to convert.

    Returns
    -------
    float or None
        Rounded float, or ``None`` if the value is NaN, Inf, or unconvertible.
    """
    if v is None:
        return None
    try:
        f = float(v)
        if math.isnan(f) or math.isinf(f):
            return None
        return float(f"{f:.4g}")
    except (TypeError, ValueError):
        return None


def _array_to_list(arr: typing.Any) -> typing.Any:
    """Recursively convert a numpy array to a nested list of JSON-safe floats.

    NaN and Inf values are replaced with ``None`` (JSON ``null``).

    Parameters
    ----------
    arr : numpy array or scalar
        The value to convert.

    Returns
    -------
    list or None
        Nested Python lists with ``None`` in place of NaN/Inf.
    """
    if not hasattr(arr, "tolist"):
        return _safe_float(arr)
    raw = arr.tolist()
    return _nested_safe(raw)


def _nested_safe(obj: typing.Any) -> typing.Any:
    """Recursively replace float NaN/Inf with None in nested lists/scalars."""
    if isinstance(obj, list):
        return [_nested_safe(x) for x in obj]
    if isinstance(obj, float):
        return _safe_float(obj)
    return obj


def parse_mat_files_for_dwi(folder: Path, dwi: str):
    """Parse all ``.mat`` files for a single DWI type and extract key metrics.

    Parameters
    ----------
    folder : Path
        Path to the ``saved_files_*`` output folder.
    dwi : str
        DWI type name (``"Standard"``, ``"dnCNN"``, or ``"IVIMnet"``).

    Returns
    -------
    dict
        Dictionary with three keys:

        - ``core_method`` -- method names, descriptions, mean Dice matrix,
          and (if present) mean Hausdorff matrix.
        - ``dosimetry`` -- mean±std D95/V50 values for ADC and D sub-volumes.
        - ``longitudinal`` -- cohort dimensions (num_patients, num_timepoints)
          extracted from the summary metrics struct.
    """
    out_data = {
        "core_method": {},
        "core_method_outcomes": {},
        "cross_pipeline_dice": {},
        "failure_rates": {},
        "pruning": {},
        "dosimetry": {},
        "longitudinal": {},
    }
    if not scipy_io or not numpy_np:
        return out_data

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

    # ── Dosimetry ──
    # The dosimetry MAT file contains 1-D arrays of per-patient dose
    # coverage metrics.  We compute the cohort mean and std for each.
    dosimetry_mat = folder / f"{dwi}" / f"metrics_dosimetry_results_{dwi}.mat"
    if dosimetry_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(dosimetry_mat), squeeze_me=True)  # type: ignore

            def _safe_nanmean_std(arr: typing.Any) -> dict:
                """Compute nanmean and nanstd; return dict with mean and std.

                NaN results are returned as None so JSON serialisation does
                not produce invalid output.
                """
                mean_val = _safe_float(numpy_np.nanmean(arr))  # type: ignore
                std_val = _safe_float(numpy_np.nanstd(arr))  # type: ignore
                return {"mean": mean_val, "std": std_val}

            nan_scalar = numpy_np.array([float("nan")])  # type: ignore
            out_data["dosimetry"] = {  # type: ignore
                "d95_adc_mean": _safe_nanmean_std(mat.get("d95_adc_sub", nan_scalar)),  # type: ignore
                "v50_adc_mean": _safe_nanmean_std(mat.get("v50_adc_sub", nan_scalar)),  # type: ignore
                "d95_d_mean": _safe_nanmean_std(mat.get("d95_d_sub", nan_scalar)),  # type: ignore
                "v50_d_mean": _safe_nanmean_std(mat.get("v50_d_sub", nan_scalar)),  # type: ignore
            }
        except Exception as e:
            print(f"Error parsing {dosimetry_mat.name}: {e}")

    # ── Longitudinal scope (from summary_metrics) ──
    # We only extract the array shape to determine how many patients and
    # timepoints were processed.
    # Expected shape: (N_patients, N_params, N_timepoints) or (N_patients, N_timepoints).
    # num_patients = shape[0]; num_timepoints = shape[-1] (last axis).
    summary_mat = folder / f"{dwi}" / f"summary_metrics_{dwi}.mat"
    if summary_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(summary_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            summary = mat.get("summary_metrics")
            if summary and hasattr(summary, "ADC_abs"):
                adc = summary.ADC_abs
                if hasattr(adc, "shape") and len(adc.shape) >= 2:
                    out_data["longitudinal"]["num_patients"] = int(adc.shape[0])
                    out_data["longitudinal"]["num_timepoints"] = int(adc.shape[-1])

            # ── Fx1 Repeat Dice (spatial repeatability) ──
            # summary_metrics stores dice_rpt_{adc,d,f,dstar} as
            # [nPatients x 3] arrays (3 DWI types: Standard=0, DnCNN=1, IVIMnet=2).
            # Extract the column for the current DWI type and compute
            # mean/std across patients with non-NaN entries.
            if summary:
                dwi_col_map = {"Standard": 0, "dnCNN": 1, "DnCNN": 1, "IVIMnet": 2}
                dwi_col = dwi_col_map.get(dwi)
                if dwi_col is not None:
                    rpt_fields = {
                        "adc": "dice_rpt_adc",
                        "d": "dice_rpt_d",
                        "f": "dice_rpt_f",
                        "dstar": "dice_rpt_dstar",
                    }
                    rpt_entry: dict = {}
                    for param, field_name in rpt_fields.items():
                        if not hasattr(summary, field_name):
                            continue
                        arr = getattr(summary, field_name)
                        if not hasattr(arr, "shape"):
                            continue
                        # Select column for this DWI type (may be 1-D if single DWI type).
                        try:
                            if arr.ndim >= 2 and arr.shape[-1] > dwi_col:
                                col = arr[:, dwi_col]
                            elif arr.ndim == 1:
                                col = arr
                            else:
                                continue
                        except Exception:
                            continue
                        col = numpy_np.asarray(col, dtype=float)  # type: ignore
                        valid = col[~numpy_np.isnan(col)]  # type: ignore
                        n_valid = int(valid.size)
                        if n_valid == 0:
                            rpt_entry[param] = {"mean": None, "std": None, "n": 0}
                        else:
                            rpt_entry[param] = {
                                "mean": _safe_float(numpy_np.mean(valid)),  # type: ignore
                                "std": _safe_float(numpy_np.std(valid)),  # type: ignore
                                "n": n_valid,
                            }
                    if rpt_entry:
                        out_data["repeatability_dice"] = rpt_entry
        except Exception as e:
            print(f"Error parsing {summary_mat.name}: {e}")

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
                            "unadjusted_hr": _safe_float(mr.unadjusted_hr) if hasattr(mr, "unadjusted_hr") else None,
                            "adjusted_hr": _safe_float(mr.adjusted_hr) if hasattr(mr, "adjusted_hr") else None,
                            "confounding_flag": bool(mr.confounding_flag) if hasattr(mr, "confounding_flag") else False,
                        }
                        gtv_entry["method_results"].append(mr_dict)
                out_data["gtv_confounding"] = _nested_safe(gtv_entry)
        except Exception as e:
            print(f"Error parsing {gtv_mat.name}: {e}")

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

    return out_data


def main():
    """CLI entry point: parse MAT files and write per-DWI-type JSON."""
    parser = argparse.ArgumentParser(description="Parse MATLAB .mat output files to JSON")
    parser.add_argument("folder", type=str, help="Path to saved_files_* folder")
    args = parser.parse_args()

    folder = Path(args.folder)
    if not folder.is_dir():
        print(f"Folder not found: {folder}")
        return

    # Discover which DWI-type subdirectories exist.
    dwi_types = [d.name for d in folder.iterdir() if d.is_dir() and d.name in ["Standard", "dnCNN", "IVIMnet"]]

    if not dwi_types:
        print("No DWI type folders found.")
        return

    pbar = tqdm(
        dwi_types,
        desc="Parsing MAT files",
        unit="type",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}] {postfix}",
    )
    for dwi in pbar:
        pbar.set_postfix_str(dwi, refresh=True)
        data = parse_mat_files_for_dwi(folder, dwi)

        # Write JSON output alongside the source folder.
        out_json = folder / f"parsed_mat_metrics_{dwi}.json"
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)


if __name__ == "__main__":
    main()
