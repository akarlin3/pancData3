#!/usr/bin/env python3
"""Extract statistically significant findings from graph analysis CSV."""

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
    """Extract all p-values from a text string."""
    results = []
    # Match p=0.XXX, p<0.XXX, p = 0.XXX, p-value=0.XXX, etc.
    patterns = [
        r'p\s*[=<>]\s*([\d.]+(?:e[+-]?\d+)?)',
        r'p-value\s*[=<>]\s*([\d.]+(?:e[+-]?\d+)?)',
    ]
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            try:
                val = float(m.group(1))
                # Get context around the match
                start = max(0, m.start() - 80)
                end = min(len(text), m.end() + 40)
                context = text[start:end].strip()
                results.append((val, context))
            except ValueError:
                pass
    return results


def extract_correlations(text):
    """Extract correlation coefficients."""
    results = []
    patterns = [
        r'r\s*=\s*(-?[\d.]+)',
        r'rs\s*=\s*(-?[\d.]+)',
        r'r²\s*=\s*([\d.]+)',
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
        for dwi_type in ["Standard", "dnCNN", "IVIMnet"]:
            if dwi_type not in dwi_dict:
                continue
            r = dwi_dict[dwi_type]
            all_text = r["summary"] + " " + r["trends_json"] + " " + r["inflection_points_json"]
            pvals = extract_pvalues(all_text)
            if pvals:
                pvals_by_dwi[dwi_type] = pvals

        if len(pvals_by_dwi) >= 2:
            print(f"\n  {base_name}:")
            for dwi_type in ["Standard", "dnCNN", "IVIMnet"]:
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
