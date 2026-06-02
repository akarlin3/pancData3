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
import sys
from pathlib import Path

from tqdm import tqdm  # type: ignore

# Ensure analysis/ root is on sys.path so 'parsers' is importable when run as
# a subprocess (the per-block parsers live in sibling modules).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Shared helpers and optional scipy/numpy shims live in a sibling module so
# this file stays under the project length limit.  They are re-imported here
# (and re-exported) so the public API of ``parsers.parse_mat_metrics`` is
# unchanged: callers and tests can still reference
# ``parse_mat_metrics.scipy_io``, ``_safe_float``, ``METHOD_DESC`` etc.
from parsers.parse_mat_helpers import (  # type: ignore  # noqa: F401
    METHOD_DESC,
    _array_to_list,
    _nested_safe,
    _safe_float,
    numpy_np,
    scipy_io,
)

# The individual ``.mat``-file parsing blocks were extracted into sibling
# modules to keep this file short.  Each ``_parse_*`` helper reads one file
# for a single DWI type and mutates the shared ``out_data`` dict in place.
from parsers.parse_mat_core import (  # type: ignore  # noqa: F401
    _parse_core_method,
    _parse_core_method_outcomes,
    _parse_cross_pipeline_dice,
    _parse_failure_rates,
    _parse_pruning,
)
from parsers.parse_mat_dosimetry import (  # type: ignore  # noqa: F401
    _parse_baseline_dosimetry,
    _parse_dosimetry,
    _parse_summary,
)
from parsers.parse_mat_outcomes import (  # type: ignore  # noqa: F401
    _parse_baseline_vs_delta,
    _parse_dose_response_roc,
    _parse_gtv_confounding,
    _parse_per_method_cor,
    _parse_risk_dose_concordance,
    _parse_subvolume_stability,
    _parse_threshold_optimization,
)


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

    # Each parsing block reads one ``.mat`` file (if present) and mutates
    # ``out_data`` in place.  They run in the original order because some
    # blocks read keys written by earlier ones (e.g. dosimetry / _baseline_m_lf).
    _parse_core_method(folder, dwi, out_data)
    _parse_core_method_outcomes(folder, dwi, out_data)
    _parse_pruning(folder, dwi, out_data)
    _parse_failure_rates(folder, dwi, out_data)
    _parse_cross_pipeline_dice(folder, dwi, out_data)

    _parse_dosimetry(folder, dwi, out_data)
    _parse_baseline_dosimetry(folder, dwi, out_data)
    _parse_summary(folder, dwi, out_data)

    _parse_threshold_optimization(folder, dwi, out_data)
    _parse_baseline_vs_delta(folder, dwi, out_data)
    _parse_per_method_cor(folder, dwi, out_data)
    _parse_subvolume_stability(folder, dwi, out_data)
    _parse_dose_response_roc(folder, dwi, out_data)
    _parse_gtv_confounding(folder, dwi, out_data)
    _parse_risk_dose_concordance(folder, dwi, out_data)
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
