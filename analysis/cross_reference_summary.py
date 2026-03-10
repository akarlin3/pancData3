#!/usr/bin/env python3
"""Concise cross-DWI summary comparing all graphs across DWI types.

Unlike :mod:`cross_reference_dwi` (which prints every matched graph in full
detail), this script provides a compact comparison of all graphs, with
priority graphs listed first:

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

from tqdm import tqdm  # type: ignore

from shared import (  # type: ignore
    DWI_TYPES,
    get_config,
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

    # Compare all graph groups that have at least 2 DWI types.
    # Priority graphs are listed first, then remaining graphs in sorted order.
    priority_graphs = get_config()["priority_graphs"]
    priority_set = set(priority_graphs)
    ordered_names = [g for g in priority_graphs if g in groups]
    ordered_names += sorted(g for g in groups if g not in priority_set)

    for base_name in tqdm(ordered_names, desc="Comparing graphs", unit="graph",
                          bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}]"):
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
            try:
                trends = json.loads(r.get("trends_json", "[]") or "[]")
            except Exception:
                trends = []
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
                    if dwi_type not in all_trends:  # type: ignore
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

        # ── Statistical tests comparison ──
        stat_tests_by_dwi: dict[str, list] = {}
        for dwi_type in DWI_TYPES:
            if dwi_type not in dwi_dict:
                continue
            try:
                tests = json.loads(dwi_dict[dwi_type].get("statistical_tests_json", "[]") or "[]")
            except Exception:
                tests = []
            if tests:
                stat_tests_by_dwi[dwi_type] = tests

        if stat_tests_by_dwi:
            print()
            for dwi_type, tests in stat_tests_by_dwi.items():
                for st in tests:
                    if not isinstance(st, dict):
                        continue
                    name = st.get("test_name", "?")
                    pval = f"p={st['p_value']:.4f}" if st.get("p_value") is not None else ""
                    cmp_groups = f" ({st['comparison_groups']})" if st.get("comparison_groups") else ""
                    print(f"    {dwi_type}: {name} {pval}{cmp_groups}")

        # ── Clinical relevance comparison ──
        clin_by_dwi = {}
        for dwi_type in DWI_TYPES:
            if dwi_type not in dwi_dict:
                continue
            clin = dwi_dict[dwi_type].get("clinical_relevance", "") or ""
            if clin:
                clin_by_dwi[dwi_type] = clin

        if clin_by_dwi:
            print()
            for dwi_type, clin in clin_by_dwi.items():
                if len(clin) > 150:
                    clin = clin[:150] + "..."
                print(f"    {dwi_type} clinical: {clin}")

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
