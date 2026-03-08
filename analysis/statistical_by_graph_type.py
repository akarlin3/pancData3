#!/usr/bin/env python3
"""Filter statistical relevance findings by graph type."""

import csv
import io
import json
import os
import re
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


def extract_pvalues(text):
    results = []
    patterns = [
        r'p\s*[=<>]\s*([\d.]+(?:e[+-]?\d+)?)',
        r'p-value\s*[=<>]\s*([\d.]+(?:e[+-]?\d+)?)',
    ]
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            try:
                val = float(m.group(1))
                start = max(0, m.start() - 80)
                end = min(len(text), m.end() + 40)
                context = text[start:end].strip()
                results.append((val, context))
            except ValueError:
                pass
    return results


def extract_correlations(text):
    results = []
    patterns = [
        r'r\s*=\s*(-?[\d.]+)',
        r'rs\s*=\s*(-?[\d.]+)',
        r'r\xb2\s*=\s*([\d.]+)',
    ]
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            try:
                val = float(m.group(1))
                start = max(0, m.start() - 60)
                end = min(len(text), m.end() + 60)
                context = text[start:end].strip()
                results.append((val, context))
            except ValueError:
                pass
    return results


def main():
    rows = []
    csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "saved_files_20260308_010713", "graph_analysis_results.csv")
    with open(csv_path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append(r)

    # Group rows by graph_type
    by_type = defaultdict(list)
    for r in rows:
        by_type[r["graph_type"]].append(r)

    sep = "=" * 80
    thin = "-" * 80

    # Type ordering by count (descending)
    type_order = sorted(by_type.keys(), key=lambda t: -len(by_type[t]))

    for graph_type in type_order:
        type_rows = by_type[graph_type]
        print(f"\n{sep}")
        print(f"  GRAPH TYPE: {graph_type.upper()} ({len(type_rows)} graphs)")
        print(sep)

        # ── Significant p-values ──
        sig = []
        nonsig = []
        for r in type_rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            all_text = r["summary"] + " " + r["trends_json"] + " " + r["inflection_points_json"]
            pvals = extract_pvalues(all_text)
            for pval, context in pvals:
                entry = {
                    "dwi": dwi_type,
                    "graph": base_name,
                    "p": pval,
                    "context": context,
                }
                if pval < 0.05:
                    sig.append(entry)
                else:
                    nonsig.append(entry)

        print(f"\n  {thin}")
        print(f"  Significant findings (p < 0.05): {len(sig)}")
        print(f"  {thin}")

        if sig:
            sig.sort(key=lambda x: x["p"])
            for f in sig:
                tag = "***" if f["p"] < 0.001 else "** " if f["p"] < 0.01 else "*  "
                print(f"\n  {tag} p={f['p']:.4f}  [{f['dwi']}] {f['graph']}")
                print(f"      {f['context']}")
        else:
            print("\n  (none)")

        print(f"\n  Non-significant: {len(nonsig)}")
        if nonsig:
            nonsig.sort(key=lambda x: x["p"])
            for f in nonsig[:5]:  # Show top 5 closest to significance
                print(f"     p={f['p']:.4f}  [{f['dwi']}] {f['graph']}")
            if len(nonsig) > 5:
                print(f"     ... and {len(nonsig) - 5} more")

        # ── Correlations ──
        corrs_found = []
        for r in type_rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            all_text = r["summary"] + " " + r["trends_json"]
            corrs = extract_correlations(all_text)
            for rval, context in corrs:
                if abs(rval) >= 0.3:
                    corrs_found.append({
                        "dwi": dwi_type,
                        "graph": base_name,
                        "r": rval,
                        "context": context,
                    })

        if corrs_found:
            print(f"\n  {thin}")
            print(f"  Notable correlations (|r| >= 0.3): {len(corrs_found)}")
            print(f"  {thin}")
            corrs_found.sort(key=lambda x: -abs(x["r"]))
            for c in corrs_found:
                strength = "STRONG" if abs(c["r"]) >= 0.5 else "MODERATE"
                direction = "positive" if c["r"] > 0 else "negative"
                print(f"\n    r={c['r']:.2f} ({strength} {direction})  [{c['dwi']}] {c['graph']}")
                print(f"      {c['context']}")

        # ── Trend summary ──
        trends_up = 0
        trends_down = 0
        trends_stable = 0
        trends_other = 0
        for r in type_rows:
            trends = json.loads(r["trends_json"])
            for t in trends:
                d = t["direction"].lower()
                if "increas" in d or "upward" in d or "rising" in d:
                    trends_up += 1
                elif "decreas" in d or "downward" in d or "declining" in d or "drop" in d:
                    trends_down += 1
                elif "stable" in d or "flat" in d or "no " in d or "constant" in d:
                    trends_stable += 1
                else:
                    trends_other += 1

        total_trends = trends_up + trends_down + trends_stable + trends_other
        if total_trends > 0:
            print(f"\n  {thin}")
            print(f"  Trend directions ({total_trends} total):")
            print(f"  {thin}")
            print(f"    Increasing: {trends_up}  |  Decreasing: {trends_down}  |  Stable: {trends_stable}  |  Other: {trends_other}")

        # ── Graph list ──
        print(f"\n  Graphs in this type:")
        for r in type_rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            print(f"    [{dwi_type}] {base_name}")

    # ── Overall summary table ──
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
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            all_text = r["summary"] + " " + r["trends_json"] + " " + r["inflection_points_json"]
            pvals = extract_pvalues(all_text)
            for pval, _ in pvals:
                if pval < 0.05:
                    n_sig += 1
                else:
                    n_nonsig += 1

            corrs = extract_correlations(r["summary"] + " " + r["trends_json"])
            n_corr += sum(1 for rval, _ in corrs if abs(rval) >= 0.3)

            trends = json.loads(r["trends_json"])
            n_trends += len(trends)

        print(f"  {graph_type:<16} {len(type_rows):>5} {n_sig:>6} {n_nonsig:>8} {n_corr:>5} {n_trends:>7}")

    print()


if __name__ == "__main__":
    main()
