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
import os
from pathlib import Path

import typing

from tqdm import tqdm

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

        - ``core_method`` -- method names and mean Dice matrix (or empty).
        - ``dosimetry`` -- mean D95/V50 values for ADC and D (or empty).
        - ``longitudinal`` -- cohort dimensions (num_patients, num_timepoints)
          extracted from the summary metrics struct (or empty).
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
    core_mat = folder / f"{dwi}" / f"compare_core_results_{dwi}.mat"
    if core_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(core_mat), squeeze_me=True, struct_as_record=False)
            results = mat.get("compare_results")
            if results:
                # Extract method names; squeeze_me may collapse single-element
                # cell arrays to a bare string.
                methods = results.method_names.tolist()
                if isinstance(methods, str):
                    methods = [methods]
                mean_dice = results.mean_dice_matrix
                # Replace NaN with 0 for JSON serialisation.
                out_data["core_method"] = {
                    "methods": [str(m) for m in methods],
                    "mean_dice_matrix": numpy_np.nan_to_num(mean_dice).tolist() if hasattr(mean_dice, "tolist") else [],
                }
        except Exception as e:
            print(f"Error parsing {core_mat.name}: {e}")

    # ── Dosimetry ──
    # The dosimetry MAT file contains 1-D arrays of per-patient dose
    # coverage metrics.  We compute the cohort mean for each.
    dosimetry_mat = folder / f"{dwi}" / f"metrics_dosimetry_results_{dwi}.mat"
    if dosimetry_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(dosimetry_mat), squeeze_me=True)
            # d95_adc_sub: D95 (minimum dose to 95% of sub-volume) for ADC
            # v50_adc_sub: V50 (fraction of sub-volume receiving >=50 Gy) for ADC
            out_data["dosimetry"] = {
                "d95_adc_mean": numpy_np.nanmean(mat.get("d95_adc_sub", numpy_np.nan)),
                "v50_adc_mean": numpy_np.nanmean(mat.get("v50_adc_sub", numpy_np.nan)),
                "d95_d_mean": numpy_np.nanmean(mat.get("d95_d_sub", numpy_np.nan)),
                "v50_d_mean": numpy_np.nanmean(mat.get("v50_d_sub", numpy_np.nan)),
            }
        except Exception as e:
            print(f"Error parsing {dosimetry_mat.name}: {e}")

    # ── Longitudinal scope (from summary_metrics) ──
    # We only extract the array shape to determine how many patients and
    # timepoints were processed.
    summary_mat = folder / f"{dwi}" / f"summary_metrics_{dwi}.mat"
    if summary_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(summary_mat), squeeze_me=True, struct_as_record=False)
            summary = mat.get("summary_metrics")
            if summary and hasattr(summary, "ADC_abs"):
                adc = summary.ADC_abs
                if hasattr(adc, "shape"):
                    out_data["longitudinal"]["num_patients"] = adc.shape[0]
                    out_data["longitudinal"]["num_timepoints"] = adc.shape[1] if len(adc.shape) > 1 else 1
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
