#!/usr/bin/env python3
"""Extract statistically significant findings from graph analysis CSV."""

from __future__ import annotations

import json
import sys
from collections import defaultdict

from shared import (
    DWI_TYPES,
    extract_correlations,
    extract_pvalues,
    load_graph_csv,
    parse_dwi_info,
    resolve_folder,
    setup_utf8_stdout,
)

setup_utf8_stdout()


def main():
    folder = resolve_folder(sys.argv)
    rows = load_graph_csv(folder)
    if not rows:
        sys.exit(f"ERROR: No graph_analysis_results.csv found in {folder}")

    groups = defaultdict(dict)
    for r in rows:
        dwi_type, base_name = parse_dwi_info(r["file_path"])
        groups[base_name][dwi_type] = r

    sep = "=" * 80

    # ── 1. Significant p-values ──────────────────────────────────────────────
    print(sep)
    print("  STATISTICALLY SIGNIFICANT FINDINGS (p < 0.05)")
    print(sep)

    sig_findings = []
    nonsig_findings = []

    for r in rows:
        dwi_type, base_name = parse_dwi_info(r["file_path"])
        # Search in trends, summary, and inflection points
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
                sig_findings.append(entry)
            else:
                nonsig_findings.append(entry)

    # Print significant
    if sig_findings:
        sig_findings.sort(key=lambda x: x["p"])
        for f in sig_findings:
            tag = "***" if f["p"] < 0.001 else "** " if f["p"] < 0.01 else "*  "
            print(f"\n  {tag} p={f['p']:.4f}  [{f['dwi']}] {f['graph']}")
            print(f"      {f['context']}")
    else:
        print("\n  No p-values < 0.05 found in extracted data.")

    # ── 2. Non-significant (for comparison) ──────────────────────────────────
    print(f"\n{sep}")
    print("  NON-SIGNIFICANT FINDINGS (p >= 0.05)")
    print(sep)

    if nonsig_findings:
        nonsig_findings.sort(key=lambda x: x["p"])
        for f in nonsig_findings:
            print(f"\n     p={f['p']:.4f}  [{f['dwi']}] {f['graph']}")
            print(f"      {f['context']}")

    # ── 3. Strong correlations ───────────────────────────────────────────────
    print(f"\n{sep}")
    print("  NOTABLE CORRELATIONS (|r| >= 0.3)")
    print(sep)

    for r in rows:
        dwi_type, base_name = parse_dwi_info(r["file_path"])
        all_text = r["summary"] + " " + r["trends_json"]
        corrs = extract_correlations(all_text)

        for rval, context in corrs:
            if abs(rval) >= 0.3:
                strength = "STRONG" if abs(rval) >= 0.5 else "MODERATE"
                direction = "positive" if rval > 0 else "negative"
                print(f"\n  [{dwi_type}] {base_name}")
                print(f"    r={rval:.2f} ({strength} {direction})")
                print(f"    {context}")

    # ── 4. Cross-DWI significance comparison ─────────────────────────────────
    print(f"\n{sep}")
    print("  CROSS-DWI: Same Analysis, Different Significance?")
    print(sep)

    for base_name in sorted(groups.keys()):
        dwi_dict = groups[base_name]
        real = [t for t in dwi_dict if t != "Root"]
        if len(real) < 2:
            continue

        # Collect p-values per DWI type for this graph
        pvals_by_dwi = {}
        for dwi_type in DWI_TYPES:
            if dwi_type not in dwi_dict:
                continue
            r = dwi_dict[dwi_type]
            all_text = r["summary"] + " " + r["trends_json"] + " " + r["inflection_points_json"]
            pvals = extract_pvalues(all_text)
            if pvals:
                pvals_by_dwi[dwi_type] = pvals

        if len(pvals_by_dwi) >= 2:
            print(f"\n  {base_name}:")
            for dwi_type in DWI_TYPES:
                if dwi_type not in pvals_by_dwi:
                    continue
                pvals = pvals_by_dwi[dwi_type]
                sig_count = sum(1 for p, _ in pvals if p < 0.05)
                total = len(pvals)
                pval_strs = [f"p={p:.3f}{'*' if p < 0.05 else ''}" for p, _ in pvals]
                print(f"    {dwi_type}: {sig_count}/{total} significant  [{', '.join(pval_strs)}]")

    print()


if __name__ == "__main__":
    main()
