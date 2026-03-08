#!/usr/bin/env python3
"""Cross-reference graph analysis results across DWI types."""

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
    folder = resolve_folder(sys.argv)
    rows = load_graph_csv(folder)
    if not rows:
        sys.exit(f"ERROR: No graph_analysis_results.csv found in {folder}")

    groups = group_by_graph_name(rows)

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

        for dwi_type in DWI_TYPES:
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
