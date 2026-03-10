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

from tqdm import tqdm  # type: ignore

from shared import (  # type: ignore
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

        matched += 1  # type: ignore
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
                scale = f", {r['x_axis_scale_type']}" if r.get("x_axis_scale_type") else ""
                print(f"    X-axis: {r['x_axis_label']} ({r.get('x_axis_units') or 'no units'}{scale})")
            if r.get("y_axis_label"):
                scale = f", {r['y_axis_scale_type']}" if r.get("y_axis_scale_type") else ""
                print(f"    Y-axis: {r['y_axis_label']} ({r.get('y_axis_units') or 'no units'}{scale})")

            # Metadata line: sample size, series count, error bars, density, quality.
            meta_parts = []
            if r.get("sample_size"):
                meta_parts.append(f"n={r['sample_size']}")
            if r.get("data_series_count"):
                meta_parts.append(f"{r['data_series_count']} series")
            if r.get("error_bars"):
                meta_parts.append(f"error bars: {r['error_bars']}")
            if r.get("data_density"):
                meta_parts.append(f"density: {r['data_density']}")
            if r.get("comparison_type"):
                meta_parts.append(f"comparison: {r['comparison_type']}")
            if r.get("figure_quality"):
                meta_parts.append(f"quality: {r['figure_quality']}")
            if r.get("subpanel_count"):
                meta_parts.append(f"{r['subpanel_count']} subpanels")
            if meta_parts:
                print(f"    Metadata: {', '.join(meta_parts)}")

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
                    mag = f" ({t['magnitude']})" if t.get("magnitude") else ""
                    sig = f" [{t['statistical_significance']}]" if t.get("statistical_significance") else ""
                    vals = ""
                    if t.get("start_value") is not None and t.get("end_value") is not None:
                        vals = f" [{t['start_value']:.4f} -> {t['end_value']:.4f}]"
                    print(f"      - {pfx}{t.get('direction', '')}: {t.get('description', '')}{mag}{sig}{vals}")

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
                    mag = f" (magnitude: {ip['magnitude']})" if ip.get("magnitude") is not None else ""
                    print(f"      - ({x}, {y}): {ip.get('description', '')}{mag}")

            # Parse and display statistical tests.
            try:
                stat_tests = json.loads(r.get("statistical_tests_json", "[]") or "[]")
            except Exception:
                stat_tests = []
            if stat_tests:
                print(f"    Statistical tests ({len(stat_tests)}):")
                for st in stat_tests:
                    if not isinstance(st, dict):
                        continue
                    name = st.get("test_name", "unknown")
                    pval = f"p={st['p_value']:.4f}" if st.get("p_value") is not None else ""
                    cmp_groups = f" ({st['comparison_groups']})" if st.get("comparison_groups") else ""
                    stat_val = f", stat={st['statistic_value']:.3f}" if st.get("statistic_value") is not None else ""
                    print(f"      - {name}: {pval}{stat_val}{cmp_groups}")

            # Parse and display outliers.
            try:
                outliers = json.loads(r.get("outliers_json", "[]") or "[]")
            except Exception:
                outliers = []
            if outliers:
                print(f"    Outliers ({len(outliers)}):")
                for o in outliers:
                    if not isinstance(o, dict):
                        continue
                    series = f"[{o['series']}] " if o.get("series") else ""
                    print(f"      - {series}{o.get('description', '')}")

            # Parse and display reference lines.
            try:
                ref_lines = json.loads(r.get("reference_lines_json", "[]") or "[]")
            except Exception:
                ref_lines = []
            if ref_lines:
                print(f"    Reference lines ({len(ref_lines)}):")
                for rl in ref_lines:
                    if not isinstance(rl, dict):
                        continue
                    lbl = f" '{rl['label']}'" if rl.get("label") else ""
                    val = f" at {rl['value']}" if rl.get("value") is not None else ""
                    style = f" ({rl['style']})" if rl.get("style") else ""
                    print(f"      - {rl.get('orientation', 'unknown')}{val}{lbl}{style}")

            # Legend items
            try:
                legend = json.loads(r.get("legend_items_json", "[]") or "[]")
            except Exception:
                legend = []
            if legend:
                print(f"    Legend: {', '.join(legend)}")

            # Spatial pattern (for parameter maps)
            if r.get("spatial_pattern"):
                print(f"    Spatial pattern: {r['spatial_pattern']}")

            # Clinical relevance
            if r.get("clinical_relevance"):
                clin = r["clinical_relevance"]
                if len(clin) > 200:
                    clin = clin[:200] + "..."
                print(f"    Clinical relevance: {clin}")

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
