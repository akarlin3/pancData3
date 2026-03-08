#!/usr/bin/env python3
"""Parse metrics from MATLAB .mat files to extend the analysis suite.

Reads:
- compare_core_results_{dwi}.mat
- metrics_dosimetry_results_{dwi}.mat
- summary_metrics_{dwi}.mat

Outputs:
- parsed_mat_metrics_{dwi}.json
"""

import argparse
import json
import os
from pathlib import Path

import typing

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
    """Parse .mat files and extract requested data."""
    out_data = {
        "core_method": {},
        "dosimetry": {},
        "longitudinal": {},
    }
    if not scipy_io or not numpy_np:
        return out_data

    # Parse Core Method Comparison
    core_mat = folder / f"{dwi}" / f"compare_core_results_{dwi}.mat"
    if core_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(core_mat), squeeze_me=True, struct_as_record=False)
            results = mat.get("compare_results")
            if results:
                # Extract matrix and labels
                methods = results.method_names.tolist()
                if isinstance(methods, str):
                    methods = [methods]
                mean_dice = results.mean_dice_matrix
                out_data["core_method"] = {
                    "methods": [str(m) for m in methods],
                    "mean_dice_matrix": numpy_np.nan_to_num(mean_dice).tolist() if hasattr(mean_dice, "tolist") else [],
                }
        except Exception as e:
            print(f"Error parsing {core_mat.name}: {e}")

    # Parse Dosimetry
    dosimetry_mat = folder / f"{dwi}" / f"metrics_dosimetry_results_{dwi}.mat"
    if dosimetry_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(dosimetry_mat), squeeze_me=True)
            # Find arrays like 'd95_adc_sub', 'v50_adc_sub'
            out_data["dosimetry"] = {
                "d95_adc_mean": numpy_np.nanmean(mat.get("d95_adc_sub", numpy_np.nan)),
                "v50_adc_mean": numpy_np.nanmean(mat.get("v50_adc_sub", numpy_np.nan)),
                "d95_d_mean": numpy_np.nanmean(mat.get("d95_d_sub", numpy_np.nan)),
                "v50_d_mean": numpy_np.nanmean(mat.get("v50_d_sub", numpy_np.nan)),
            }
        except Exception as e:
            print(f"Error parsing {dosimetry_mat.name}: {e}")

    # Parse Longitudinal (from summary_metrics)
    summary_mat = folder / f"{dwi}" / f"summary_metrics_{dwi}.mat"
    if summary_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(summary_mat), squeeze_me=True, struct_as_record=False)
            # Check for fields like 'ADC_abs', 'D_abs'
            summary = mat.get("summary_metrics")
            # For simplicity, we just confirm its presence to show we processed it,
            # or extract the shape of the data to give basic longitudinal scope.
            if summary and hasattr(summary, "ADC_abs"):
                adc = summary.ADC_abs
                if hasattr(adc, "shape"):
                    out_data["longitudinal"]["num_patients"] = adc.shape[0]
                    out_data["longitudinal"]["num_timepoints"] = adc.shape[1] if len(adc.shape) > 1 else 1
        except Exception as e:
            print(f"Error parsing {summary_mat.name}: {e}")

    return out_data

def main():
    parser = argparse.ArgumentParser(description="Parse MATLAB output logs")
    parser.add_argument("folder", type=str, help="Path to saved_files_* folder")
    args = parser.parse_args()

    folder = Path(args.folder)
    if not folder.is_dir():
        print(f"Folder not found: {folder}")
        return

    # Find DWI folders
    dwi_types = [d.name for d in folder.iterdir() if d.is_dir() and d.name in ["Standard", "dnCNN", "IVIMnet"]]

    if not dwi_types:
        print("No DWI type folders found.")
        return

    for dwi in dwi_types:
        data = parse_mat_files_for_dwi(folder, dwi)
        
        # Save JSON
        out_json = folder / f"parsed_mat_metrics_{dwi}.json"
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)

        print(f"Extracted MAT data for {dwi}")

if __name__ == "__main__":
    main()
