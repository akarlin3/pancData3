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
                a = numpy_np.asarray(arr, dtype=float).ravel()  # type: ignore
                if a.size == 0 or numpy_np.all(numpy_np.isnan(a)):  # type: ignore
                    return {"mean": None, "std": None}
                mean_val = _safe_float(numpy_np.nanmean(a))  # type: ignore
                std_val = _safe_float(numpy_np.nanstd(a))  # type: ignore
                return {"mean": mean_val, "std": std_val}

            nan_scalar = numpy_np.array([float("nan")])  # type: ignore
            out_data["dosimetry"] = {  # type: ignore
                "d95_adc_mean": _safe_nanmean_std(mat.get("d95_adc_sub", nan_scalar)),  # type: ignore
                "v50_adc_mean": _safe_nanmean_std(mat.get("v50_adc_sub", nan_scalar)),  # type: ignore
                "dmean_adc_mean": _safe_nanmean_std(mat.get("dmean_adc_sub", nan_scalar)),  # type: ignore
                "d95_d_mean": _safe_nanmean_std(mat.get("d95_d_sub", nan_scalar)),  # type: ignore
                "v50_d_mean": _safe_nanmean_std(mat.get("v50_d_sub", nan_scalar)),  # type: ignore
                "dmean_d_mean": _safe_nanmean_std(mat.get("dmean_d_sub", nan_scalar)),  # type: ignore
            }
        except Exception as e:
            print(f"Error parsing {dosimetry_mat.name}: {e}")

    # ── Whole-GTV dosimetry (from baseline_results) ──
    # The baseline results MAT file stores m_d95_gtvp and m_v50gy_gtvp
    # (per-patient whole-GTV D95 / V50 arrays).  We compute cohort
    # nanmean / nanstd and add them to out_data["dosimetry"] alongside
    # the sub-volume metrics.  dmean_gtvp lives on summary_metrics.
    baseline_mat = folder / f"{dwi}" / f"metrics_baseline_results_{dwi}.mat"
    if baseline_mat.exists() and scipy_io and numpy_np:
        try:
            mat: typing.Any = scipy_io.loadmat(str(baseline_mat), squeeze_me=True, struct_as_record=False)  # type: ignore

            def _whole_mean_std(arr_in: typing.Any) -> typing.Tuple[typing.Optional[float], typing.Optional[float]]:
                if arr_in is None:
                    return None, None
                try:
                    a = numpy_np.asarray(arr_in, dtype=float)  # type: ignore
                    # Use baseline (column 0) if 2-D; else flatten.
                    if a.ndim >= 2:
                        col = a[:, 0]
                    else:
                        col = a
                    return (
                        _safe_float(numpy_np.nanmean(col)),  # type: ignore
                        _safe_float(numpy_np.nanstd(col)),   # type: ignore
                    )
                except Exception:
                    return None, None

            dosi_entry = out_data.get("dosimetry") or {}
            d95_m, d95_s = _whole_mean_std(mat.get("m_d95_gtvp"))
            v50_m, v50_s = _whole_mean_std(mat.get("m_v50gy_gtvp"))
            if d95_m is not None:
                dosi_entry["whole_gtv_d95_mean"] = d95_m
                dosi_entry["whole_gtv_d95_std"] = d95_s
            if v50_m is not None:
                dosi_entry["whole_gtv_v50_mean"] = v50_m
                dosi_entry["whole_gtv_v50_std"] = v50_s

            # m_lf may live in baseline file (used by subvolume_sizes).
            m_lf_arr = mat.get("m_lf")
            if m_lf_arr is not None:
                try:
                    out_data["_baseline_m_lf"] = _array_to_list(numpy_np.asarray(m_lf_arr).ravel())  # type: ignore
                except Exception:
                    pass

            out_data["dosimetry"] = dosi_entry
        except Exception as e:
            print(f"Error parsing {baseline_mat.name}: {e}")

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

            # ── Whole-GTV mean dose (dmean_gtvp on summary_metrics) ──
            if summary and hasattr(summary, "dmean_gtvp"):
                try:
                    dm = numpy_np.asarray(summary.dmean_gtvp, dtype=float)  # type: ignore
                    if dm.ndim >= 2:
                        dm_col = dm[:, 0]
                    else:
                        dm_col = dm
                    dosi_entry = out_data.get("dosimetry") or {}
                    dosi_entry["whole_gtv_dmean_mean"] = _safe_float(numpy_np.nanmean(dm_col))  # type: ignore
                    dosi_entry["whole_gtv_dmean_std"] = _safe_float(numpy_np.nanstd(dm_col))  # type: ignore
                    out_data["dosimetry"] = dosi_entry
                except Exception:
                    pass

            # ── Sub-Volume Sizes (adc_sub_vol / adc_sub_vol_pc) ──
            # Shapes: [nPatients x nTimepoints x 3].  Extract the column for
            # this DWI type and stratify by outcome (m_lf from baseline).
            if summary:
                dwi_col_map = {"Standard": 0, "dnCNN": 1, "DnCNN": 1, "IVIMnet": 2}
                dwi_col = dwi_col_map.get(dwi)
                if dwi_col is not None and hasattr(summary, "adc_sub_vol"):
                    try:
                        vol = numpy_np.asarray(summary.adc_sub_vol, dtype=float)  # type: ignore
                        if vol.ndim == 3 and vol.shape[-1] > dwi_col:
                            vol2d = vol[:, :, dwi_col]  # [nPat x nTp]
                        elif vol.ndim == 2:
                            vol2d = vol
                        else:
                            vol2d = None

                        pc2d = None
                        if hasattr(summary, "adc_sub_vol_pc"):
                            pc = numpy_np.asarray(summary.adc_sub_vol_pc, dtype=float)  # type: ignore
                            if pc.ndim == 3 and pc.shape[-1] > dwi_col:
                                pc2d = pc[:, :, dwi_col]
                            elif pc.ndim == 2:
                                pc2d = pc

                        if vol2d is not None and vol2d.ndim == 2 and vol2d.size > 0:
                            n_tp = int(vol2d.shape[1])
                            # Outcome vector — prefer the raw clinical
                            # encoding from summary_metrics.lf (binary
                            # 0 = LC / 1 = LF, read straight from the
                            # clinical spreadsheet) over the baseline
                            # m_lf variant, which may re-encode competing-
                            # risk deaths as 2 and therefore empty the LF
                            # group for cohorts with few pure local-failure
                            # events.  Fall back to baseline m_lf only when
                            # summary.lf is unavailable or sized wrong.
                            n_pat_summary = int(vol2d.shape[0])
                            m_lf = None
                            if hasattr(summary, "lf"):
                                try:
                                    cand = numpy_np.asarray(summary.lf, dtype=float).ravel()  # type: ignore
                                    if cand.size == n_pat_summary:
                                        m_lf = cand
                                except Exception:
                                    m_lf = None
                            if m_lf is None:
                                m_lf_list = out_data.get("_baseline_m_lf") or []
                                if m_lf_list:
                                    cand = numpy_np.asarray(m_lf_list, dtype=float)  # type: ignore
                                    if cand.size == n_pat_summary:
                                        m_lf = cand

                            # Attempt scipy Mann-Whitney (optional dependency).
                            try:
                                from scipy.stats import mannwhitneyu  # type: ignore
                                have_mwu = True
                            except Exception:
                                mannwhitneyu = None  # type: ignore
                                have_mwu = False

                            timepoints: list[int] = []
                            lc_mean_vol: list[typing.Optional[float]] = []
                            lc_std_vol: list[typing.Optional[float]] = []
                            lf_mean_vol: list[typing.Optional[float]] = []
                            lf_std_vol: list[typing.Optional[float]] = []
                            lc_mean_pct: list[typing.Optional[float]] = []
                            lf_mean_pct: list[typing.Optional[float]] = []
                            wilcoxon_p: list[typing.Optional[float]] = []
                            median_vox: list[typing.Optional[float]] = []
                            overall_mean_vol: list[typing.Optional[float]] = []
                            overall_std_vol: list[typing.Optional[float]] = []
                            overall_mean_pct: list[typing.Optional[float]] = []

                            for k in range(n_tp):
                                timepoints.append(k + 1)
                                col_vol = vol2d[:, k]
                                col_pc = pc2d[:, k] if pc2d is not None else None

                                # Overall median voxel count (vol is cm^3; use as-is).
                                try:
                                    median_vox.append(_safe_float(numpy_np.nanmedian(col_vol)))  # type: ignore
                                except Exception:
                                    median_vox.append(None)

                                # Overall mean/std across all patients
                                # (not stratified by outcome).
                                try:
                                    finite_vol = col_vol[~numpy_np.isnan(col_vol)]  # type: ignore
                                    overall_mean_vol.append(
                                        _safe_float(numpy_np.nanmean(finite_vol)) if finite_vol.size else None)  # type: ignore
                                    overall_std_vol.append(
                                        _safe_float(numpy_np.nanstd(finite_vol)) if finite_vol.size else None)  # type: ignore
                                except Exception:
                                    overall_mean_vol.append(None)
                                    overall_std_vol.append(None)
                                try:
                                    if col_pc is not None:
                                        finite_pc = col_pc[~numpy_np.isnan(col_pc)]  # type: ignore
                                        overall_mean_pct.append(
                                            _safe_float(numpy_np.nanmean(finite_pc)) if finite_pc.size else None)  # type: ignore
                                    else:
                                        overall_mean_pct.append(None)
                                except Exception:
                                    overall_mean_pct.append(None)

                                if m_lf is not None and m_lf.size == col_vol.size:
                                    lf_mask = (m_lf == 1) & ~numpy_np.isnan(col_vol)  # type: ignore
                                    lc_mask = (m_lf == 0) & ~numpy_np.isnan(col_vol)  # type: ignore
                                    lf_vals = col_vol[lf_mask]
                                    lc_vals = col_vol[lc_mask]
                                else:
                                    lf_vals = numpy_np.array([])  # type: ignore
                                    lc_vals = numpy_np.array([])  # type: ignore

                                lc_mean_vol.append(_safe_float(numpy_np.nanmean(lc_vals)) if lc_vals.size else None)  # type: ignore
                                lc_std_vol.append(_safe_float(numpy_np.nanstd(lc_vals)) if lc_vals.size else None)  # type: ignore
                                lf_mean_vol.append(_safe_float(numpy_np.nanmean(lf_vals)) if lf_vals.size else None)  # type: ignore
                                lf_std_vol.append(_safe_float(numpy_np.nanstd(lf_vals)) if lf_vals.size else None)  # type: ignore

                                if col_pc is not None and m_lf is not None and m_lf.size == col_pc.size:
                                    lf_pc_mask = (m_lf == 1) & ~numpy_np.isnan(col_pc)  # type: ignore
                                    lc_pc_mask = (m_lf == 0) & ~numpy_np.isnan(col_pc)  # type: ignore
                                    lf_pc_vals = col_pc[lf_pc_mask]
                                    lc_pc_vals = col_pc[lc_pc_mask]
                                    lc_mean_pct.append(_safe_float(numpy_np.nanmean(lc_pc_vals)) if lc_pc_vals.size else None)  # type: ignore
                                    lf_mean_pct.append(_safe_float(numpy_np.nanmean(lf_pc_vals)) if lf_pc_vals.size else None)  # type: ignore
                                else:
                                    lc_mean_pct.append(None)
                                    lf_mean_pct.append(None)

                                p_val = None
                                if have_mwu and lf_vals.size > 0 and lc_vals.size > 0:
                                    try:
                                        _, p_val = mannwhitneyu(lf_vals, lc_vals, alternative="two-sided")  # type: ignore
                                        p_val = _safe_float(p_val)
                                    except Exception:
                                        p_val = None
                                wilcoxon_p.append(p_val)

                            out_data["subvolume_sizes"] = {
                                "timepoints": timepoints,
                                "overall": {
                                    "mean_vol_cm3": overall_mean_vol,
                                    "std_vol_cm3": overall_std_vol,
                                    "mean_frac_pct": overall_mean_pct,
                                },
                                "lc": {
                                    "mean_vol_cm3": lc_mean_vol,
                                    "std_vol_cm3": lc_std_vol,
                                    "mean_frac_pct": lc_mean_pct,
                                },
                                "lf": {
                                    "mean_vol_cm3": lf_mean_vol,
                                    "std_vol_cm3": lf_std_vol,
                                    "mean_frac_pct": lf_mean_pct,
                                },
                                "wilcoxon_p": wilcoxon_p,
                                "median_voxel_count": median_vox,
                            }
                    except Exception:
                        pass
        except Exception as e:
            print(f"Error parsing {summary_mat.name}: {e}")

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
                ):
                    if hasattr(results, field):
                        opt_entry[field] = _array_to_list(getattr(results, field))
                for scalar_field in (
                    "optimal_thresh",
                    "optimal_dice",
                    "optimal_vol_frac",
                ):
                    if hasattr(results, scalar_field):
                        opt_entry[scalar_field] = _safe_float(getattr(results, scalar_field))
                if hasattr(results, "n_patients"):
                    try:
                        opt_entry["n_patients"] = int(getattr(results, "n_patients"))
                    except Exception:
                        pass
                out_data["threshold_optimization"] = opt_entry
        except Exception as e:
            print(f"Error parsing {opt_mat.name}: {e}")

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
