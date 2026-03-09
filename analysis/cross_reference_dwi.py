#!/usr/bin/env python3
"""Full cross-reference of graph analysis results across DWI types.

For every graph that appears in two or more DWI-type subfolders (Standard,
dnCNN, IVIMnet), this script prints a detailed side-by-side comparison
including:

- Graph type, axis labels, and units
- All extracted trends (series name, direction, description)
- All inflection points (approximate coordinates and description)
- Plain-English summary (truncated to 250 characters)

This is the verbose counterpart to :mod:`cross_reference_summary`, which
focuses only on clinically priority graphs and trend agreement.

Usage:
    python cross_reference_dwi.py [saved_files_path]
"""

from __future__ import annotations

import json
import sys

from tqdm import tqdm

from shared import (
    DWI_TYPES,
    group_by_graph_name,
    load_graph_csv,
    resolve_folder,
    setup_utf8_stdout,
)

setup_utf8_stdout()


def main():
    """CLI entry point: load vision CSV and print full cross-DWI comparison."""
    folder = resolve_folder(sys.argv)
    rows = load_graph_csv(folder)
    if not rows:
        sys.exit(f"ERROR: No graph_analysis_results.csv found in {folder}")

    # Group rows by normalised graph name for cross-DWI matching.
    groups = group_by_graph_name(rows)

    sep = "=" * 90
    print(sep)
    print("CROSS-DWI TYPE COMPARISON: Standard vs dnCNN vs IVIMnet")
    print(sep)

    matched = 0
    for base_name in tqdm(sorted(groups.keys()), desc="Cross-referencing graphs", unit="graph",
                          bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}]"):
        dwi_dict = groups[base_name]
        types_present = sorted(dwi_dict.keys())
        if len(types_present) < 2:
            continue
        # Skip graphs that only exist in Root (not inside a DWI-type subfolder).
        real_types = [t for t in types_present if t != "Root"]
        if len(real_types) < 2:
            continue

        matched += 1
        print(f"\n{sep}")
        print(f"  Graph: {base_name}")
        print(f"  Present in: {', '.join(types_present)}")
        print(sep)

        # Print detailed info for each DWI type.
        for dwi_type in DWI_TYPES:
            if dwi_type not in dwi_dict:
                continue
            r = dwi_dict[dwi_type]
            print(f"\n  [{dwi_type}]")
            print(f"    Type: {r.get('graph_type', 'unknown')}")
            if r.get("x_axis_label"):
                print(f"    X-axis: {r['x_axis_label']} ({r.get('x_axis_units') or 'no units'})")
            if r.get("y_axis_label"):
                print(f"    Y-axis: {r['y_axis_label']} ({r.get('y_axis_units') or 'no units'})")

            # Parse and display trend list from JSON string.
            try:
                trends = json.loads(r.get("trends_json", "[]") or "[]")
            except Exception:
                trends = []
            if trends:
                print(f"    Trends ({len(trends)}):")
                for t in trends:
                    if not isinstance(t, dict):
                        continue
                    series = t.get("series") or ""
                    pfx = f"[{series}] " if series else ""
                    print(f"      - {pfx}{t.get('direction', '')}: {t.get('description', '')}")

            # Parse and display inflection points from JSON string.
            try:
                inflections = json.loads(r.get("inflection_points_json", "[]") or "[]")
            except Exception:
                inflections = []
            if inflections:
                print(f"    Inflection points ({len(inflections)}):")
                for ip in inflections:
                    if not isinstance(ip, dict):
                        continue
                    x = ip.get("approximate_x", "?")
                    y = ip.get("approximate_y", "?")
                    print(f"      - ({x}, {y}): {ip.get('description', '')}")

            # Truncate long summaries for readability.
            summary = r.get("summary", "") or ""
            if len(summary) > 250:
                summary = summary[:250] + "..."
            print(f"    Summary: {summary}")

    print(f"\n{sep}")
    print(f"Total matched graph sets across DWI types: {matched}")
    print(sep)


if __name__ == "__main__":
    main()
