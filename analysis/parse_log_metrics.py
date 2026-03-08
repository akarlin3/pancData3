#!/usr/bin/env python3
"""Parse MATLAB pipeline log files to extract structured metrics.

Reads the diary-captured log files from each DWI type subfolder and extracts:
- Wilcoxon rank-sum p-values and effect sizes
- GLME interaction p-values
- AUC values, selected features, hazard ratios
- Cox PH results and sensitivity analysis
- Baseline outlier summaries
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

from shared import DWI_TYPES, resolve_folder, setup_utf8_stdout

setup_utf8_stdout()


# ── Regex patterns ────────────────────────────────────────────────────────────

# metrics_stats_comparisons
RE_GLME_INTERACTION = re.compile(
    r"Interaction P-Value.*?:\s*([0-9.]+(?:e[+-]?\d+)?)", re.IGNORECASE
)
RE_GLME_DETAIL = re.compile(
    r"(\w[\w* ]*?):\s*p=([0-9.]+(?:e[+-]?\d+)?),\s*adj_alpha=([0-9.]+)"
)
RE_FDR_TIMEPOINT = re.compile(
    r"Timepoint:\s*(\S+)\s*—\s*(\d+)\s+significant"
)
RE_GLME_EXCLUDED = re.compile(
    r"Excluded\s+(\d+)/(\d+)\s+\(([0-9.]+)%\)\s+competing-risk"
)

# metrics_stats_predictive
RE_ELASTIC_NET = re.compile(
    r"Elastic Net Selected Features for (\S+)\s+\(Opt Lambda=([0-9.]+)\):\s+(.*)"
)
RE_FIRTH = re.compile(r"Firth refit successful for (\S+)\s+\((\d+) features\)")
RE_ROC_HEADER = re.compile(r"PRIMARY ROC ANALYSIS.*for (\S+)")
RE_AUC = re.compile(r"AUC\s*=\s*([0-9.]+)")
RE_YOUDEN = re.compile(r"Youden Optimal Score Cutoff\s*=\s*([0-9.]+)")
RE_SENS_SPEC = re.compile(
    r"Sensitivity\s*=\s*([0-9.]+)%\s*\|\s*Specificity\s*=\s*([0-9.]+)%"
)

# metrics_survival
RE_HR_ROW = re.compile(
    r"^\s*([\w*]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$",
    re.MULTILINE,
)
RE_GLOBAL_LRT = re.compile(
    r"Global LRT:\s*chi2\((\d+)\)\s*=\s*([0-9.]+),\s*p\s*=\s*([0-9.]+)"
)
RE_IPCW_WEIGHTS = re.compile(
    r"IPCW weights applied.*\[([0-9.]+),\s*([0-9.]+)\]"
)

# metrics_baseline
RE_OUTLIER_FLAG = re.compile(
    r"Outlier flag \((\w+)\):\s*(\d+) flagged \(LF=(\d+),\s*LC=(\d+),\s*CR=(\d+)\)"
)
RE_TOTAL_OUTLIERS = re.compile(
    r"Total outliers removed:\s*(\d+)\s*/\s*(\d+)\s+\(([0-9.]+)%\)"
)
RE_BASELINE_EXCLUDED = re.compile(
    r"Excluded\s+(\d+)/(\d+) patients due to missing baseline"
)
RE_LF_RATE = re.compile(
    r"LF rate:\s*included=([0-9.]+)%,\s*excluded=([0-9.]+)%"
)


def _read_log(path: Path) -> str:
    """Read a log file, return empty string if missing."""
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def parse_stats_comparisons(text: str) -> dict:
    """Parse metrics_stats_comparisons output."""
    result = {
        "glme_interactions": [],
        "glme_details": [],
        "fdr_timepoints": [],
        "glme_excluded": None,
    }

    for m in RE_GLME_INTERACTION.finditer(text):
        result["glme_interactions"].append(float(m.group(1)))

    for m in RE_GLME_DETAIL.finditer(text):
        result["glme_details"].append({
            "metric": m.group(1).strip(),
            "p": float(m.group(2)),
            "adj_alpha": float(m.group(3)),
        })

    for m in RE_FDR_TIMEPOINT.finditer(text):
        result["fdr_timepoints"].append({
            "timepoint": m.group(1),
            "n_significant": int(m.group(2)),
        })

    m = RE_GLME_EXCLUDED.search(text)
    if m:
        result["glme_excluded"] = {
            "n_excluded": int(m.group(1)),
            "n_total": int(m.group(2)),
            "pct": float(m.group(3)),
        }

    return result


def parse_stats_predictive(text: str) -> dict:
    """Parse metrics_stats_predictive output."""
    result = {"feature_selections": [], "roc_analyses": []}

    for m in RE_ELASTIC_NET.finditer(text):
        result["feature_selections"].append({
            "timepoint": m.group(1),
            "lambda": float(m.group(2)),
            "features": [f.strip() for f in m.group(3).split(",") if f.strip()],
        })

    # Parse ROC blocks: each starts with a header line
    roc_headers = list(RE_ROC_HEADER.finditer(text))
    for i, hdr in enumerate(roc_headers):
        start = hdr.start()
        end = roc_headers[i + 1].start() if i + 1 < len(roc_headers) else len(text)
        block = text[start:end]

        entry = {"timepoint": hdr.group(1)}
        m = RE_AUC.search(block)
        if m:
            entry["auc"] = float(m.group(1))
        m = RE_YOUDEN.search(block)
        if m:
            entry["youden_cutoff"] = float(m.group(1))
        m = RE_SENS_SPEC.search(block)
        if m:
            entry["sensitivity"] = float(m.group(1))
            entry["specificity"] = float(m.group(2))
        result["roc_analyses"].append(entry)

    return result


def parse_survival(text: str) -> dict:
    """Parse metrics_survival output."""
    result = {"hazard_ratios": [], "global_lrt": None, "ipcw": None}

    for m in RE_HR_ROW.finditer(text):
        result["hazard_ratios"].append({
            "covariate": m.group(1),
            "hr": float(m.group(2)),
            "ci_lo": float(m.group(3)),
            "ci_hi": float(m.group(4)),
            "p": float(m.group(5)),
        })

    m = RE_GLOBAL_LRT.search(text)
    if m:
        result["global_lrt"] = {
            "df": int(m.group(1)),
            "chi2": float(m.group(2)),
            "p": float(m.group(3)),
        }

    m = RE_IPCW_WEIGHTS.search(text)
    if m:
        result["ipcw"] = {
            "min_weight": float(m.group(1)),
            "max_weight": float(m.group(2)),
        }

    return result


def parse_baseline(text: str) -> dict:
    """Parse metrics_baseline output."""
    result = {"outlier_flags": [], "total_outliers": None, "baseline_exclusion": None}

    for m in RE_OUTLIER_FLAG.finditer(text):
        result["outlier_flags"].append({
            "metric": m.group(1),
            "n_flagged": int(m.group(2)),
            "n_lf": int(m.group(3)),
            "n_lc": int(m.group(4)),
            "n_cr": int(m.group(5)),
        })

    m = RE_TOTAL_OUTLIERS.search(text)
    if m:
        result["total_outliers"] = {
            "n_removed": int(m.group(1)),
            "n_total": int(m.group(2)),
            "pct": float(m.group(3)),
        }

    m = RE_BASELINE_EXCLUDED.search(text)
    if m:
        result["baseline_exclusion"] = {
            "n_excluded": int(m.group(1)),
            "n_total": int(m.group(2)),
        }
        m2 = RE_LF_RATE.search(text)
        if m2:
            result["baseline_exclusion"]["lf_rate_included"] = float(m2.group(1))
            result["baseline_exclusion"]["lf_rate_excluded"] = float(m2.group(2))

    return result


def parse_all_logs(folder: Path) -> dict:
    """Parse all log files from a saved_files folder.

    Returns dict keyed by DWI type, each containing parsed sub-dicts.
    """
    results = {}
    for dwi_type in DWI_TYPES:
        dwi_dir = folder / dwi_type
        if not dwi_dir.is_dir():
            continue

        results[dwi_type] = {
            "stats_comparisons": parse_stats_comparisons(
                _read_log(dwi_dir / f"metrics_stats_comparisons_output_{dwi_type}.txt")
            ),
            "stats_predictive": parse_stats_predictive(
                _read_log(dwi_dir / f"metrics_stats_predictive_output_{dwi_type}.txt")
            ),
            "survival": parse_survival(
                _read_log(dwi_dir / f"metrics_survival_output_{dwi_type}.txt")
            ),
            "baseline": parse_baseline(
                _read_log(dwi_dir / f"metrics_baseline_output_{dwi_type}.txt")
            ),
        }

    return results


def main():
    folder = resolve_folder(sys.argv)
    results = parse_all_logs(folder)

    if not results:
        print(f"No DWI type subfolders found in {folder}")
        return

    sep = "=" * 80

    for dwi_type, data in results.items():
        print(f"\n{sep}")
        print(f"  {dwi_type} — Parsed Log Metrics")
        print(sep)

        # Baseline
        bl = data["baseline"]
        if bl["outlier_flags"]:
            print("\n  Outlier Flags:")
            for o in bl["outlier_flags"]:
                print(f"    {o['metric']}: {o['n_flagged']} flagged (LF={o['n_lf']}, LC={o['n_lc']}, CR={o['n_cr']})")
        if bl["total_outliers"]:
            t = bl["total_outliers"]
            print(f"  Total outliers: {t['n_removed']}/{t['n_total']} ({t['pct']}%)")

        # Stats comparisons
        sc = data["stats_comparisons"]
        if sc["glme_details"]:
            print("\n  GLME Interaction Tests:")
            for g in sc["glme_details"]:
                sig = "*" if g["p"] < g["adj_alpha"] else ""
                print(f"    {g['metric']}: p={g['p']:.4f} (adj_alpha={g['adj_alpha']:.4f}) {sig}")
        if sc["fdr_timepoints"]:
            print("\n  FDR-Significant Timepoints:")
            for f in sc["fdr_timepoints"]:
                print(f"    {f['timepoint']}: {f['n_significant']} significant")

        # Predictive
        sp = data["stats_predictive"]
        if sp["roc_analyses"]:
            print("\n  ROC/AUC Performance:")
            for r in sp["roc_analyses"]:
                auc = r.get("auc", "N/A")
                sens = r.get("sensitivity", "N/A")
                spec = r.get("specificity", "N/A")
                print(f"    {r['timepoint']}: AUC={auc}, Sens={sens}%, Spec={spec}%")
        if sp["feature_selections"]:
            print("\n  Selected Features:")
            for fs in sp["feature_selections"]:
                print(f"    {fs['timepoint']}: {', '.join(fs['features'])}")

        # Survival
        sv = data["survival"]
        if sv["hazard_ratios"]:
            print("\n  Hazard Ratios (Cox PH):")
            print(f"    {'Covariate':<12} {'HR':>6} {'CI_lo':>6} {'CI_hi':>6} {'p':>8}")
            for hr in sv["hazard_ratios"]:
                print(f"    {hr['covariate']:<12} {hr['hr']:>6.3f} {hr['ci_lo']:>6.3f} {hr['ci_hi']:>6.3f} {hr['p']:>8.4f}")
        if sv["global_lrt"]:
            g = sv["global_lrt"]
            print(f"\n  Global LRT: chi2({g['df']}) = {g['chi2']:.2f}, p = {g['p']:.4f}")

    print()


if __name__ == "__main__":
    main()
