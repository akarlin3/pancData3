#!/usr/bin/env python3
"""Extract statistically significant findings from the vision analysis CSV.

This script performs a comprehensive statistical survey of the
``graph_analysis_results.csv`` produced by :mod:`batch_graph_analysis`.
It uses regex-based extraction (:func:`shared.extract_pvalues`,
:func:`shared.extract_correlations`) to mine p-values and correlation
coefficients from the free-text summaries, trend descriptions, and
inflection-point descriptions.

Output sections:

1. **Significant findings** (p < 0.05) -- sorted by p-value with
   significance markers (``*``, ``**``, ``***``).
2. **Non-significant findings** (p >= 0.05) -- for comparison.
3. **Notable correlations** (|r| >= 0.3) -- classified as moderate or
   strong, positive or negative.
4. **Cross-DWI significance comparison** -- for graphs present in multiple
   DWI types, shows whether the same analysis yields significance in
   different DWI-type runs.

Usage:
    python statistical_relevance.py [saved_files_path]
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict

from shared import (
    DWI_TYPES,
    extract_correlations,
    extract_pvalues,
    get_config,
    load_graph_csv,
    parse_dwi_info,
    resolve_folder,
    setup_utf8_stdout,
)

setup_utf8_stdout()


def main():
    """CLI entry point: extract and print all statistical findings."""
    folder = resolve_folder(sys.argv)
    rows = load_graph_csv(folder)
    if not rows:
        sys.exit(f"ERROR: No graph_analysis_results.csv found in {folder}")

    # Group rows by normalised graph name for cross-DWI comparison (section 4).
    groups = defaultdict(dict)
    for r in rows:
        dwi_type, base_name = parse_dwi_info(r["file_path"])
        groups[base_name][dwi_type] = r

    stats_cfg = get_config()["statistics"]
    p_threshold = stats_cfg["p_noteworthy"]
    corr_threshold = stats_cfg["correlation_threshold"]

    sep = "=" * 80

    # ── 1. Significant p-values ──────────────────────────────────────────────
    print(sep)
    print(f"  STATISTICALLY SIGNIFICANT FINDINGS (p < {p_threshold})")
    print(sep)

    sig_findings = []
    nonsig_findings = []

    for r in rows:
        dwi_type, base_name = parse_dwi_info(r["file_path"])
        # Concatenate all text fields to search for p-value patterns.
        all_text = r["summary"] + " " + r["trends_json"] + " " + r["inflection_points_json"]
        pvals = extract_pvalues(all_text)

        for pval, context in pvals:
            entry = {
                "dwi": dwi_type,
                "graph": base_name,
                "p": pval,
                "context": context,
            }
            if pval < p_threshold:
                sig_findings.append(entry)
            else:
                nonsig_findings.append(entry)

    # Print significant findings sorted by ascending p-value.
    if sig_findings:
        sig_findings.sort(key=lambda x: x["p"])
        for f in sig_findings:
            p_hi = stats_cfg["p_highly_significant"]
            p_sig = stats_cfg["p_significant"]
            tag = "***" if f["p"] < p_hi else "** " if f["p"] < p_sig else "*  "
            print(f"\n  {tag} p={f['p']:.4f}  [{f['dwi']}] {f['graph']}")
            print(f"      {f['context']}")
    else:
        print("\n  No p-values < 0.05 found in extracted data.")

    # ── 2. Non-significant findings (for comparison) ─────────────────────────
    print(f"\n{sep}")
    print(f"  NON-SIGNIFICANT FINDINGS (p >= {p_threshold})")
    print(sep)

    if nonsig_findings:
        nonsig_findings.sort(key=lambda x: x["p"])
        for f in nonsig_findings:
            print(f"\n     p={f['p']:.4f}  [{f['dwi']}] {f['graph']}")
            print(f"      {f['context']}")

    # ── 3. Strong correlations ───────────────────────────────────────────────
    print(f"\n{sep}")
    print(f"  NOTABLE CORRELATIONS (|r| >= {corr_threshold})")
    print(sep)

    for r in rows:
        dwi_type, base_name = parse_dwi_info(r["file_path"])
        all_text = r["summary"] + " " + r["trends_json"]
        corrs = extract_correlations(all_text)

        for rval, context in corrs:
            if abs(rval) >= corr_threshold:
                # Classify correlation strength and direction.
                strength = "STRONG" if abs(rval) >= 0.5 else "MODERATE"
                direction = "positive" if rval > 0 else "negative"
                print(f"\n  [{dwi_type}] {base_name}")
                print(f"    r={rval:.2f} ({strength} {direction})")
                print(f"    {context}")

    # ── 4. Cross-DWI significance comparison ─────────────────────────────────
    # For graphs present in multiple DWI types, compare how many p-values
    # reach significance in each type -- highlights processing-method effects.
    print(f"\n{sep}")
    print("  CROSS-DWI: Same Analysis, Different Significance?")
    print(sep)

    for base_name in sorted(groups.keys()):
        dwi_dict = groups[base_name]
        real = [t for t in dwi_dict if t != "Root"]
        if len(real) < 2:
            continue

        # Collect p-values per DWI type for this graph.
        pvals_by_dwi: dict[str, list] = {}
        for dwi_type in DWI_TYPES:
            if dwi_type not in dwi_dict:
                continue
            r = dwi_dict[dwi_type]
            all_text = r["summary"] + " " + r["trends_json"] + " " + r["inflection_points_json"]
            pvals = extract_pvalues(all_text)
            if pvals:
                pvals_by_dwi[dwi_type] = pvals

        # Only print if we have p-values from at least 2 DWI types.
        if len(pvals_by_dwi) >= 2:
            print(f"\n  {base_name}:")
            for dwi_type in DWI_TYPES:
                if dwi_type not in pvals_by_dwi:
                    continue
                pvals = pvals_by_dwi[dwi_type]
                sig_count = sum(1 for p, _ in pvals if p < p_threshold)
                total = len(pvals)
                pval_strs = [f"p={p:.3f}{'*' if p < p_threshold else ''}" for p, _ in pvals]
                print(f"    {dwi_type}: {sig_count}/{total} significant  [{', '.join(pval_strs)}]")

    print()


if __name__ == "__main__":
    main()
