#!/usr/bin/env python3
"""Cross-reference graph analysis results across DWI types."""

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
            next_part = parts[i + 1]
            if next_part in ("Standard", "dnCNN", "IVIMnet"):
                dwi_type = next_part
            break
    # Normalize: remove DWI type suffix from filename
    for t in ["_Standard", "_dnCNN", "_IVIMnet"]:
        base_name = base_name.replace(t, "")
    base_name = base_name.replace(".png", "")
    return dwi_type, base_name


def main():
    csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "saved_files_20260308_010713", "graph_analysis_results.csv")
    rows = []
    with open(csv_path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append(r)

    # Group by base graph name
    groups = defaultdict(dict)
    for r in rows:
        dwi_type, base_name = parse_dwi_info(r["file_path"])
        groups[base_name][dwi_type] = r

    sep = "=" * 90
    print(sep)
    print("CROSS-DWI TYPE COMPARISON: Standard vs dnCNN vs IVIMnet")
    print(sep)

    matched = 0
    for base_name in sorted(groups.keys()):
        dwi_dict = groups[base_name]
        types_present = sorted(dwi_dict.keys())
        if len(types_present) < 2:
            continue
        # Skip Root-only
        real_types = [t for t in types_present if t != "Root"]
        if len(real_types) < 2:
            continue

        matched += 1
        print(f"\n{sep}")
        print(f"  Graph: {base_name}")
        print(f"  Present in: {', '.join(types_present)}")
        print(sep)

        for dwi_type in ["Standard", "dnCNN", "IVIMnet"]:
            if dwi_type not in dwi_dict:
                continue
            r = dwi_dict[dwi_type]
            print(f"\n  [{dwi_type}]")
            print(f"    Type: {r['graph_type']}")
            if r["x_axis_label"]:
                print(f"    X-axis: {r['x_axis_label']} ({r['x_axis_units'] or 'no units'})")
            if r["y_axis_label"]:
                print(f"    Y-axis: {r['y_axis_label']} ({r['y_axis_units'] or 'no units'})")

            trends = json.loads(r["trends_json"])
            if trends:
                print(f"    Trends ({len(trends)}):")
                for t in trends:
                    series = t.get("series") or ""
                    pfx = f"[{series}] " if series else ""
                    print(f"      - {pfx}{t['direction']}: {t['description']}")

            inflections = json.loads(r["inflection_points_json"])
            if inflections:
                print(f"    Inflection points ({len(inflections)}):")
                for ip in inflections:
                    x = ip.get("approximate_x", "?")
                    y = ip.get("approximate_y", "?")
                    print(f"      - ({x}, {y}): {ip['description']}")

            summary = r["summary"]
            if len(summary) > 250:
                summary = summary[:250] + "..."
            print(f"    Summary: {summary}")

    print(f"\n{sep}")
    print(f"Total matched graph sets across DWI types: {matched}")
    print(sep)


if __name__ == "__main__":
    main()
