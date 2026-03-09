#!/usr/bin/env python3
"""Concise cross-DWI summary focusing on clinically priority graphs.

Unlike :mod:`cross_reference_dwi` (which prints every matched graph in full
detail), this script focuses on a curated list of clinically important graph
types and provides a compact comparison:

1. **Trend agreement/disagreement** -- per data series, shows whether
   Standard, dnCNN, and IVIMnet agree on the direction ("AGREE" vs "DIFFER").
2. **Truncated summaries** -- side-by-side plain-English summaries (capped
   at 180 characters each).
3. **Parameter map counts** -- which DWI types produced parameter maps.
4. **Inflection points** -- longitudinal inflection points compared across
   DWI types.

Usage:
    python cross_reference_summary.py [saved_files_path]
"""

from __future__ import annotations

import json
import sys

from shared import (
    DWI_TYPES,
    group_by_graph_name,
    load_graph_csv,
    resolve_folder,
    setup_utf8_stdout,
)

setup_utf8_stdout()


def main():
    """CLI entry point: load vision CSV and print concise cross-DWI summary."""
    folder = resolve_folder(sys.argv)
    rows = load_graph_csv(folder)
    if not rows:
        sys.exit(f"ERROR: No graph_analysis_results.csv found in {folder}")

    groups = group_by_graph_name(rows)

    sep = "=" * 80
    print(sep)
    print("  CROSS-DWI TYPE SUMMARY: Key Differences")
    print(sep)

    # Curated list of clinically interesting graphs to compare across DWI types.
    # These cover dose-response, longitudinal trajectories, feature distributions,
    # and tumor core method agreement.
    priority_graphs = [
        "Dose_vs_Diffusion",
        "Longitudinal_Mean_Metrics",
        "Longitudinal_Mean_Metrics_ByOutcome",
        "Longitudinal_Mean_Metrics_LC",
        "Longitudinal_Mean_Metrics_LF",
        "Feature_BoxPlots",
        "Feature_Histograms",
        "core_method_dice_heatmap",
        "core_method_volume_comparison",
    ]

    for base_name in priority_graphs:
        if base_name not in groups:
            continue
        dwi_dict = groups[base_name]
        # Need at least 2 real DWI types (not Root) for comparison.
        real = [t for t in dwi_dict if t != "Root"]
        if len(real) < 2:
            continue

        print(f"\n{'-' * 80}")
        print(f"  {base_name}")
        print(f"{'-' * 80}")

        # ── Collect trends per DWI type ──
        all_trends = {}
        for dwi_type in DWI_TYPES:
            if dwi_type not in dwi_dict:
                continue
            r = dwi_dict[dwi_type]
            trends = json.loads(r["trends_json"])
            all_trends[dwi_type] = trends

        # ── Compare trend directions across DWI types ──
        if all_trends:
            # Build union of all series names across DWI types.
            all_series: set[str] = set()
            for dwi_type, trends in all_trends.items():
                for t in trends:
                    s = t.get("series") or "overall"
                    all_series.add(s)

            print()
            for series in sorted(all_series):
                # Collect the direction for this series from each DWI type.
                directions: dict[str, str] = {}
                for dwi_type in DWI_TYPES:
                    if dwi_type not in all_trends:
                        continue
                    for t in all_trends[dwi_type]:
                        s = t.get("series") or "overall"
                        if s == series:
                            directions[dwi_type] = t["direction"]

                if len(directions) >= 2:
                    vals = list(directions.values())
                    # All directions identical = AGREE; otherwise DIFFER.
                    match = "AGREE" if len(set(vals)) == 1 else "DIFFER"
                    # ">>" prefix highlights disagreements visually.
                    tag = "  " if match == "AGREE" else ">>"
                    compact = ", ".join(f"{k}={v}" for k, v in directions.items())
                    print(f"  {tag} [{series}]: {compact}  ({match})")

        # ── Print truncated summaries side by side ──
        print()
        for dwi_type in DWI_TYPES:
            if dwi_type not in dwi_dict:
                continue
            summary = dwi_dict[dwi_type]["summary"]
            if len(summary) > 180:
                summary = summary[:180] + "..."
            print(f"    {dwi_type}: {summary}")

    # ── Parameter map counts by DWI type ──
    print(f"\n{sep}")
    print("  PARAMETER MAP COUNTS BY DWI TYPE")
    print(sep)
    for base_name in sorted(groups.keys()):
        if "Parameter_Maps" not in base_name:
            continue
        dwi_dict = groups[base_name]
        types = [t for t in dwi_dict if t != "Root"]
        if types:
            print(f"  {base_name}: {', '.join(sorted(types))}")

    # ── Inflection points comparison for longitudinal graphs ──
    print(f"\n{sep}")
    print("  INFLECTION POINTS: Longitudinal Graphs")
    print(sep)
    for base_name in sorted(groups.keys()):
        if "Longitudinal" not in base_name:
            continue
        dwi_dict = groups[base_name]
        has_ip = False
        for dwi_type in DWI_TYPES:
            if dwi_type not in dwi_dict:
                continue
            ips = json.loads(dwi_dict[dwi_type]["inflection_points_json"])
            if ips:
                if not has_ip:
                    print(f"\n  {base_name}:")
                    has_ip = True
                for ip in ips:
                    x = ip.get("approximate_x", "?")
                    print(f"    {dwi_type}: x={x} - {ip['description']}")

    print()


if __name__ == "__main__":
    main()
