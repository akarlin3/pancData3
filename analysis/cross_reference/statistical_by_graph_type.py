#!/usr/bin/env python3
"""Filter statistical relevance findings by graph type.

Groups the vision analysis results by ``graph_type`` (line, scatter, box,
heatmap, histogram, parameter_map, etc.) and for each type reports:

- **Significant p-values** (p < 0.05) with significance markers.
- **Non-significant p-values** (top 5 closest to significance shown).
- **Notable correlations** (|r| >= 0.3) classified by strength/direction.
- **Trend direction summary** -- counts of increasing, decreasing, stable,
  and other trends.
- **Graph list** -- all graphs belonging to this type.

A final summary table aggregates counts across all graph types.

Usage:
    python statistical_by_graph_type.py [saved_files_path]
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

from tqdm import tqdm  # type: ignore

# Ensure analysis/ root is on sys.path so 'shared' is importable when run as subprocess.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import (  # type: ignore
    extract_correlations,
    extract_pvalues,
    get_config,
    load_graph_csv,
    parse_dwi_info,
    resolve_folder,
    safe_text,
    setup_utf8_stdout,
)

setup_utf8_stdout()


def main():
    """CLI entry point: group by graph type and print statistical findings."""
    folder = resolve_folder(sys.argv)
    rows = load_graph_csv(folder)
    if not rows:
        sys.exit(f"ERROR: No graph_analysis_results.csv found in {folder}")

    stats_cfg = get_config()["statistics"]
    p_threshold = stats_cfg["p_noteworthy"]
    corr_threshold = stats_cfg["correlation_threshold"]

    # Group rows by their graph_type field (line, scatter, box, etc.).
    by_type: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_type[r["graph_type"]].append(r)

    sep = "=" * 80
    thin = "-" * 80

    # Sort graph types by count (most common first).
    type_order = sorted(by_type.keys(), key=lambda t: -len(by_type[t]))

    for graph_type in tqdm(type_order, desc="Analyzing graph types", unit="type",
                           bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}]"):
        type_rows = by_type[graph_type]
        print(f"\n{sep}")
        print(f"  GRAPH TYPE: {graph_type.upper()} ({len(type_rows)} graphs)")
        print(sep)

        # ── Significant p-values ──
        sig: list[dict] = []
        nonsig: list[dict] = []
        for r in type_rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            # Concatenate all text fields to search for p-value patterns.
            all_text = safe_text(r, "summary", "trends_json", "inflection_points_json")
            pvals = extract_pvalues(all_text)
            for pval, context in pvals:
                entry = {
                    "dwi": dwi_type,
                    "graph": base_name,
                    "p": pval,
                    "context": context,
                }
                if pval < p_threshold:
                    sig.append(entry)
                else:
                    nonsig.append(entry)

        print(f"\n  {thin}")
        print(f"  Significant findings (p < {p_threshold}): {len(sig)}")
        print(f"  {thin}")

        if sig:
            sig.sort(key=lambda x: x["p"])
            for f in sig:
                p_hi = stats_cfg["p_highly_significant"]  # type: ignore
                p_sig = stats_cfg["p_significant"]  # type: ignore
                tag = "***" if f["p"] < p_hi else "** " if f["p"] < p_sig else "*  "
                print(f"\n  {tag} p={f['p']:.4f}  [{f['dwi']}] {f['graph']}")
                print(f"      {f['context']}")
        else:
            print("\n  (none)")

        # Show top 5 non-significant findings closest to significance threshold.
        print(f"\n  Non-significant: {len(nonsig)}")
        if nonsig:
            nonsig.sort(key=lambda x: x["p"])
            for f in nonsig[:5]:  # type: ignore
                print(f"     p={f['p']:.4f}  [{f['dwi']}] {f['graph']}")
            if len(nonsig) > 5:
                print(f"     ... and {len(nonsig) - 5} more")

        # ── Correlations ──
        corrs_found: list[dict] = []
        for r in type_rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            all_text = safe_text(r, "summary", "trends_json")
            corrs = extract_correlations(all_text)
            for rval, context in corrs:
                if abs(rval) >= corr_threshold:
                    corrs_found.append({
                        "dwi": dwi_type,
                        "graph": base_name,
                        "r": rval,
                        "context": context,
                    })

        if corrs_found:
            print(f"\n  {thin}")
            print(f"  Notable correlations (|r| >= {corr_threshold}): {len(corrs_found)}")
            print(f"  {thin}")
            # Sort by absolute correlation strength (strongest first).
            corrs_found.sort(key=lambda x: -abs(x["r"]))  # type: ignore
            for c in corrs_found:
                strength = "STRONG" if abs(c["r"]) >= 0.5 else "MODERATE"
                direction = "positive" if c["r"] > 0 else "negative"
                print(f"\n    r={c['r']:.2f} ({strength} {direction})  [{c['dwi']}] {c['graph']}")
                print(f"      {c['context']}")

        # ── Trend direction summary ──
        # Classify each trend by keyword matching on its direction string.
        trends_up = 0
        trends_down = 0
        trends_stable = 0
        trends_other = 0
        for r in type_rows:
            try:
                trends = json.loads(r.get("trends_json", "[]") or "[]")
            except Exception:
                trends = []
            for t in trends:
                if not isinstance(t, dict):
                    continue
                d = (t.get("direction") or "").lower()
                if "increas" in d or "up" in d or "higher" in d or "rising" in d:
                    trends_up += 1  # type: ignore
                elif "decreas" in d or "down" in d or "lower" in d or "falling" in d or "drop" in d:
                    trends_down += 1  # type: ignore
                elif "flat" in d or "stable" in d or "constant" in d:
                    trends_stable += 1  # type: ignore
                else:
                    trends_other += 1  # type: ignore

        total_trends = trends_up + trends_down + trends_stable + trends_other  # type: ignore
        if total_trends > 0:
            print(f"\n  {thin}")
            print(f"  Trend directions ({total_trends} total):")
            print(f"  {thin}")
            print(f"    Increasing: {trends_up}  |  Decreasing: {trends_down}  |  Stable: {trends_stable}  |  Other: {trends_other}")

        # ── Structured statistical tests ──
        structured_tests: list[dict] = []
        for r in type_rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            try:
                tests = json.loads(r.get("statistical_tests_json", "[]") or "[]")
            except Exception:
                tests = []
            for st in tests:
                if isinstance(st, dict) and st.get("test_name"):
                    structured_tests.append({**st, "dwi": dwi_type, "graph": base_name})

        if structured_tests:
            print(f"\n  {thin}")
            print(f"  Structured statistical tests ({len(structured_tests)}):")
            print(f"  {thin}")
            for st in structured_tests:
                pval_str = f"p={st['p_value']:.4f}" if st.get("p_value") is not None else "p=N/A"
                groups = f" ({st['comparison_groups']})" if st.get("comparison_groups") else ""
                print(f"    [{st['dwi']}] {st['graph']}: {st['test_name']} {pval_str}{groups}")

        # ── Data density / comparison type breakdown ──
        densities: dict[str, int] = {}
        comp_types: dict[str, int] = {}
        for r in type_rows:
            d = r.get("data_density", "") or ""
            if d:
                densities[d] = densities.get(d, 0) + 1
            c = r.get("comparison_type", "") or ""
            if c:
                comp_types[c] = comp_types.get(c, 0) + 1

        if densities:
            density_str = ", ".join(f"{k}: {v}" for k, v in sorted(densities.items()))
            print(f"\n  Data density: {density_str}")
        if comp_types:
            comp_str = ", ".join(f"{k}: {v}" for k, v in sorted(comp_types.items()))
            print(f"  Comparison types: {comp_str}")

        # ── Graph list ──
        print(f"\n  Graphs in this type:")
        for r in type_rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            meta = []
            if r.get("sample_size"):
                meta.append(r["sample_size"])
            if r.get("figure_quality"):
                meta.append(f"quality:{r['figure_quality']}")
            suffix = f"  ({', '.join(meta)})" if meta else ""
            print(f"    [{dwi_type}] {base_name}{suffix}")

    # ── Overall summary table ──
    # Compact tabular view aggregating counts across all graph types.
    print(f"\n\n{sep}")
    print("  SUMMARY TABLE: Statistical Findings by Graph Type")
    print(sep)
    print(f"\n  {'Type':<16} {'Count':>5} {'Sig p':>6} {'Non-sig':>8} {'Corr':>5} {'Trends':>7}")
    print(f"  {'-'*14:<16} {'-----':>5} {'------':>6} {'--------':>8} {'-----':>5} {'-------':>7}")

    for graph_type in type_order:
        type_rows = by_type[graph_type]
        n_sig = 0
        n_nonsig = 0
        n_corr = 0
        n_trends = 0

        for r in type_rows:
            all_text = safe_text(r, "summary", "trends_json", "inflection_points_json")
            pvals = extract_pvalues(all_text)
            for pval, _ in pvals:
                if pval < p_threshold:
                    n_sig += 1  # type: ignore
                else:
                    n_nonsig += 1  # type: ignore

            corrs = extract_correlations(safe_text(r, "summary", "trends_json"))
            n_corr += sum(1 for rval, _ in corrs if abs(rval) >= corr_threshold)  # type: ignore

            try:
                trends = json.loads(r.get("trends_json", "[]") or "[]")
            except Exception:
                trends = []
            n_trends += len(trends)  # type: ignore

        print(f"  {graph_type:<16} {len(type_rows):>5} {n_sig:>6} {n_nonsig:>8} {n_corr:>5} {n_trends:>7}")

    print()


if __name__ == "__main__":
    main()
