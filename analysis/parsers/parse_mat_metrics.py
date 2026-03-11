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
        except Exception as e:
            print(f"Error parsing {summary_mat.name}: {e}")

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
