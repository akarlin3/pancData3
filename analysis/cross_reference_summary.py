#!/usr/bin/env python3
"""Concise cross-DWI summary: key differences between Standard, dnCNN, IVIMnet."""

import csv
import io
import json
import os
import sys
from collections import defaultdict

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


def parse_dwi_info(fp):
    fp = fp.replace("\\", "/")
    parts = fp.split("/")
    dwi_type = "Root"
    base_name = parts[-1]
    for i, p in enumerate(parts):
        if "saved_files" in p and i + 1 < len(parts):
            if parts[i + 1] in ("Standard", "dnCNN", "IVIMnet"):
                dwi_type = parts[i + 1]
            break
    for t in ["_Standard", "_dnCNN", "_IVIMnet"]:
        base_name = base_name.replace(t, "")
    return dwi_type, base_name.replace(".png", "")


def main():
    rows = []
    csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "saved_files_20260308_010713", "graph_analysis_results.csv")
    with open(csv_path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append(r)

    groups = defaultdict(dict)
    for r in rows:
        dwi_type, base_name = parse_dwi_info(r["file_path"])
        groups[base_name][dwi_type] = r

    sep = "=" * 80
    print(sep)
    print("  CROSS-DWI TYPE SUMMARY: Key Differences")
    print(sep)

    # Focus on the most clinically interesting graphs
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
        real = [t for t in dwi_dict if t != "Root"]
        if len(real) < 2:
            continue

        print(f"\n{'-' * 80}")
        print(f"  {base_name}")
        print(f"{'-' * 80}")

        # Collect trends per DWI type
        all_trends = {}
        for dwi_type in ["Standard", "dnCNN", "IVIMnet"]:
            if dwi_type not in dwi_dict:
                continue
            r = dwi_dict[dwi_type]
            trends = json.loads(r["trends_json"])
            all_trends[dwi_type] = trends

        # Compare trend directions
        if all_trends:
            # Get union of series names
            all_series = set()
            for dwi_type, trends in all_trends.items():
                for t in trends:
                    s = t.get("series") or "overall"
                    all_series.add(s)

            print()
            for series in sorted(all_series):
                directions = {}
                for dwi_type in ["Standard", "dnCNN", "IVIMnet"]:
                    if dwi_type not in all_trends:
                        continue
                    for t in all_trends[dwi_type]:
                        s = t.get("series") or "overall"
                        if s == series:
                            directions[dwi_type] = t["direction"]

                if len(directions) >= 2:
                    vals = list(directions.values())
                    match = "AGREE" if len(set(vals)) == 1 else "DIFFER"
                    tag = "  " if match == "AGREE" else ">>"
                    compact = ", ".join(f"{k}={v}" for k, v in directions.items())
                    print(f"  {tag} [{series}]: {compact}  ({match})")

        # Print summaries side by side
        print()
        for dwi_type in ["Standard", "dnCNN", "IVIMnet"]:
            if dwi_type not in dwi_dict:
                continue
            summary = dwi_dict[dwi_type]["summary"]
            if len(summary) > 180:
                summary = summary[:180] + "..."
            print(f"    {dwi_type}: {summary}")

    # Overall stats
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

    # Inflection points comparison for longitudinal graphs
    print(f"\n{sep}")
    print("  INFLECTION POINTS: Longitudinal Graphs")
    print(sep)
    for base_name in sorted(groups.keys()):
        if "Longitudinal" not in base_name:
            continue
        dwi_dict = groups[base_name]
        has_ip = False
        for dwi_type in ["Standard", "dnCNN", "IVIMnet"]:
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
