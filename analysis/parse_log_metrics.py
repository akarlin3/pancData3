#!/usr/bin/env python3
"""Parse MATLAB pipeline log files to extract structured metrics.

The MATLAB pipeline captures console output (via ``diary``) into per-module
log files inside each DWI-type subfolder.  This script reads those text files
and applies targeted regular expressions to extract quantitative results:

- **metrics_stats_comparisons**: GLME interaction p-values, BH-FDR timepoints,
  competing-risk exclusion counts.
- **metrics_stats_predictive**: Elastic-net selected features, ROC/AUC
  performance, Youden optimal cutoffs, sensitivity/specificity.
- **metrics_survival**: Cox PH hazard ratios with confidence intervals,
  global likelihood-ratio test, IPCW weight ranges.
- **metrics_baseline**: Per-metric outlier flags (LF/LC/CR breakdown),
  total outlier removal, baseline exclusion counts and LF-rate comparison.

The parsed data is returned as nested dicts (keyed by DWI type) and printed
as a human-readable summary when the script is run directly.

Usage:
    python parse_log_metrics.py [saved_files_path]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

from shared import DWI_TYPES, resolve_folder, setup_utf8_stdout

setup_utf8_stdout()


# ── Regex patterns ────────────────────────────────────────────────────────────
# Each compiled regex targets a specific line format produced by the
# corresponding MATLAB core module.  Capture groups extract the numeric
# values (and labels) that are assembled into structured dicts.

# ----- metrics_stats_comparisons -----

# Matches: "Interaction P-Value (LF vs LC): 0.0234"
# Captures the numeric p-value (group 1), supporting scientific notation.
RE_GLME_INTERACTION = re.compile(
    r"Interaction P-Value.*?:\s*([0-9.]+(?:e[+-]?\d+)?)", re.IGNORECASE
)

# Matches: "Mean_ADC: p=0.0012, adj_alpha=0.0250"
# Captures metric name (group 1), raw p-value (group 2), and BH-adjusted
# significance threshold (group 3).
RE_GLME_DETAIL = re.compile(
    r"(\w[\w* ]*?):\s*p=([0-9.]+(?:e[+-]?\d+)?),\s*adj_alpha=([0-9.]+)"
)

# Matches: "Timepoint: Fx5 — 3 significant"
# Captures timepoint label (group 1) and count of significant metrics (group 2).
RE_FDR_TIMEPOINT = re.compile(
    r"Timepoint:\s*(\S+)\s*\u2014\s*(\d+)\s+significant"
)

# Matches: "Excluded 4/25 (16.0%) competing-risk"
# Captures counts and percentage of patients excluded from GLME analysis.
RE_GLME_EXCLUDED = re.compile(
    r"Excluded\s+(\d+)/(\d+)\s+\(([0-9.]+)%\)\s+competing-risk"
)

# ----- metrics_stats_predictive -----

# Matches: "Elastic Net Selected Features for Fx10 (Opt Lambda=0.0523): feat1, feat2"
# Captures timepoint (group 1), optimal lambda (group 2), and comma-separated
# feature names (group 3).
RE_ELASTIC_NET = re.compile(
    r"Elastic Net Selected Features for (\S+)\s+\(Opt Lambda=([0-9.]+)\):\s+(.*)"
)

# Matches: "Firth refit successful for Fx5 (3 features)"
# Captures timepoint (group 1) and number of features (group 2).
RE_FIRTH = re.compile(r"Firth refit successful for (\S+)\s+\((\d+) features\)")

# Matches: "PRIMARY ROC ANALYSIS ... for Fx10"
# Captures the timepoint label (group 1); used to delimit ROC blocks.
RE_ROC_HEADER = re.compile(r"PRIMARY ROC ANALYSIS.*for (\S+)")

# Matches: "AUC = 0.843"
RE_AUC = re.compile(r"AUC\s*=\s*([0-9.]+)")

# Matches: "Youden Optimal Score Cutoff = 0.512"
RE_YOUDEN = re.compile(r"Youden Optimal Score Cutoff\s*=\s*([0-9.]+)")

# Matches: "Sensitivity = 85.7% | Specificity = 72.3%"
RE_SENS_SPEC = re.compile(
    r"Sensitivity\s*=\s*([0-9.]+)%\s*\|\s*Specificity\s*=\s*([0-9.]+)%"
)

# ----- metrics_survival -----

# Matches tabular Cox PH rows, e.g.:
#   "Mean_D*    2.340   1.120   4.890   0.0231"
# Captures covariate name (group 1), HR (2), CI lower (3), CI upper (4), p (5).
RE_HR_ROW = re.compile(
    r"^\s*([\w*]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$",
    re.MULTILINE,
)

# Matches: "Global LRT: chi2(3) = 12.45, p = 0.0061"
# Captures degrees of freedom (1), chi-squared statistic (2), p-value (3).
RE_GLOBAL_LRT = re.compile(
    r"Global LRT:\s*chi2\((\d+)\)\s*=\s*([0-9.]+),\s*p\s*=\s*([0-9.]+)"
)

# Matches: "IPCW weights applied ... [0.850, 1.230]"
# Captures the minimum (1) and maximum (2) IPCW weight values.
RE_IPCW_WEIGHTS = re.compile(
    r"IPCW weights applied.*\[([0-9.]+),\s*([0-9.]+)\]"
)

# ----- metrics_baseline -----

# Matches: "Outlier flag (ADC): 3 flagged (LF=1, LC=2, CR=0)"
# Captures metric name (1), total flagged (2), and per-outcome counts (3-5).
RE_OUTLIER_FLAG = re.compile(
    r"Outlier flag \((\w+)\):\s*(\d+) flagged \(LF=(\d+),\s*LC=(\d+),\s*CR=(\d+)\)"
)

# Matches: "Total outliers removed: 5 / 200 (2.50%)"
RE_TOTAL_OUTLIERS = re.compile(
    r"Total outliers removed:\s*(\d+)\s*/\s*(\d+)\s+\(([0-9.]+)%\)"
)

# Matches: "Excluded 3/25 patients due to missing baseline"
RE_BASELINE_EXCLUDED = re.compile(
    r"Excluded\s+(\d+)/(\d+) patients due to missing baseline"
)

# Matches: "LF rate: included=40.0%, excluded=20.0%"
# Used to compare local-failure rates between included and excluded patients.
RE_LF_RATE = re.compile(
    r"LF rate:\s*included=([0-9.]+)%,\s*excluded=([0-9.]+)%"
)


def _read_log(path: Path) -> str:
    """Read a log file and return its contents as a string.

    Parameters
    ----------
    path : Path
        Path to the log file.

    Returns
    -------
    str
        File contents, or an empty string if the file does not exist.
        Encoding errors are replaced to avoid crashes on malformed logs.
    """
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


# ── Per-module parsers ───────────────────────────────────────────────────────


def parse_stats_comparisons(text: str) -> dict:
    """Parse ``metrics_stats_comparisons`` log output.

    Extracts GLME interaction p-values, per-metric detail rows (with
    BH-adjusted alpha), FDR-significant timepoint counts, and
    competing-risk exclusion statistics.

    Parameters
    ----------
    text : str
        Full text of the ``metrics_stats_comparisons_output_*.txt`` log.

    Returns
    -------
    dict
        Keys: ``glme_interactions``, ``glme_details``, ``fdr_timepoints``,
        ``glme_excluded``.
    """
    result = {
        "glme_interactions": [],
        "glme_details": [],
        "fdr_timepoints": [],
        "glme_excluded": None,
    }

    # Collect raw interaction p-values (one per GLME model).
    for m in RE_GLME_INTERACTION.finditer(text):
        result["glme_interactions"].append(float(m.group(1)))

    # Collect per-metric detail rows with raw p and adjusted alpha.
    for m in RE_GLME_DETAIL.finditer(text):
        result["glme_details"].append({
            "metric": m.group(1).strip(),
            "p": float(m.group(2)),
            "adj_alpha": float(m.group(3)),
        })

    # Collect per-timepoint FDR significance counts.
    for m in RE_FDR_TIMEPOINT.finditer(text):
        result["fdr_timepoints"].append({
            "timepoint": m.group(1),
            "n_significant": int(m.group(2)),
        })

    # Competing-risk exclusion line (appears at most once).
    m = RE_GLME_EXCLUDED.search(text)
    if m:
        result["glme_excluded"] = {
            "n_excluded": int(m.group(1)),
            "n_total": int(m.group(2)),
            "pct": float(m.group(3)),
        }

    return result


def parse_stats_predictive(text: str) -> dict:
    """Parse ``metrics_stats_predictive`` log output.

    Extracts elastic-net feature selections (with optimal lambda) and
    ROC analysis blocks (AUC, Youden cutoff, sensitivity, specificity).

    Parameters
    ----------
    text : str
        Full text of the ``metrics_stats_predictive_output_*.txt`` log.

    Returns
    -------
    dict
        Keys: ``feature_selections``, ``roc_analyses``.
    """
    result = {"feature_selections": [], "roc_analyses": []}

    # Extract elastic-net feature selections per timepoint.
    for m in RE_ELASTIC_NET.finditer(text):
        result["feature_selections"].append({
            "timepoint": m.group(1),
            "lambda": float(m.group(2)),
            "features": [f.strip() for f in m.group(3).split(",") if f.strip()],
        })

    # Parse ROC blocks.  Each block starts with a "PRIMARY ROC ANALYSIS"
    # header line and extends until the next header (or end of text).
    roc_headers = list(RE_ROC_HEADER.finditer(text))
    for i, hdr in enumerate(roc_headers):
        start = hdr.start()
        # The block ends at the start of the next ROC header, or EOF.
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
    """Parse ``metrics_survival`` log output.

    Extracts Cox proportional-hazards table rows (covariate, HR, 95% CI,
    p-value), the global likelihood-ratio test, and IPCW weight ranges.

    Parameters
    ----------
    text : str
        Full text of the ``metrics_survival_output_*.txt`` log.

    Returns
    -------
    dict
        Keys: ``hazard_ratios``, ``global_lrt``, ``ipcw``.
    """
    result = {"hazard_ratios": [], "global_lrt": None, "ipcw": None}

    # Parse whitespace-delimited Cox PH table rows.
    for m in RE_HR_ROW.finditer(text):
        result["hazard_ratios"].append({
            "covariate": m.group(1),
            "hr": float(m.group(2)),
            "ci_lo": float(m.group(3)),
            "ci_hi": float(m.group(4)),
            "p": float(m.group(5)),
        })

    # Global likelihood-ratio test (single occurrence).
    m = RE_GLOBAL_LRT.search(text)
    if m:
        result["global_lrt"] = {
            "df": int(m.group(1)),
            "chi2": float(m.group(2)),
            "p": float(m.group(3)),
        }

    # IPCW (Inverse Probability of Censoring Weighting) range.
    m = RE_IPCW_WEIGHTS.search(text)
    if m:
        result["ipcw"] = {
            "min_weight": float(m.group(1)),
            "max_weight": float(m.group(2)),
        }

    return result


def parse_baseline(text: str) -> dict:
    """Parse ``metrics_baseline`` log output.

    Extracts per-metric outlier flags with outcome breakdown (LF/LC/CR),
    total outlier removal statistics, and baseline exclusion counts with
    associated LF-rate comparison.

    Parameters
    ----------
    text : str
        Full text of the ``metrics_baseline_output_*.txt`` log.

    Returns
    -------
    dict
        Keys: ``outlier_flags``, ``total_outliers``, ``baseline_exclusion``.
    """
    result = {"outlier_flags": [], "total_outliers": None, "baseline_exclusion": None}

    # Per-metric outlier flags with LF/LC/CR breakdown.
    for m in RE_OUTLIER_FLAG.finditer(text):
        result["outlier_flags"].append({
            "metric": m.group(1),
            "n_flagged": int(m.group(2)),
            "n_lf": int(m.group(3)),
            "n_lc": int(m.group(4)),
            "n_cr": int(m.group(5)),
        })

    # Aggregate outlier removal line.
    m = RE_TOTAL_OUTLIERS.search(text)
    if m:
        result["total_outliers"] = {
            "n_removed": int(m.group(1)),
            "n_total": int(m.group(2)),
            "pct": float(m.group(3)),
        }

    # Baseline exclusion count and optional LF-rate comparison.
    m = RE_BASELINE_EXCLUDED.search(text)
    if m:
        result["baseline_exclusion"] = {
            "n_excluded": int(m.group(1)),
            "n_total": int(m.group(2)),
        }
        # Look for the LF-rate comparison line that follows exclusion info.
        m2 = RE_LF_RATE.search(text)
        if m2:
            result["baseline_exclusion"]["lf_rate_included"] = float(m2.group(1))
            result["baseline_exclusion"]["lf_rate_excluded"] = float(m2.group(2))

    return result


def parse_all_logs(folder: Path) -> dict:
    """Parse all log files from a ``saved_files_*`` folder.

    Iterates over each known DWI type subfolder and applies all four
    module-specific parsers to their respective log files.

    Parameters
    ----------
    folder : Path
        Path to the ``saved_files_*`` output folder.

    Returns
    -------
    dict
        Mapping of ``dwi_type -> {"stats_comparisons": ...,
        "stats_predictive": ..., "survival": ..., "baseline": ...}``.
        Only DWI types whose subfolder exists are included.
    """
    results = {}
    for dwi_type in DWI_TYPES:
        dwi_dir = folder / dwi_type
        if not dwi_dir.is_dir():
            continue

        # Each log file follows the naming convention:
        #   <module_name>_output_<DWI_type>.txt
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
    """CLI entry point: parse all logs and print a human-readable summary."""
    folder = resolve_folder(sys.argv)
    results = parse_all_logs(folder)

    if not results:
        print(f"No DWI type subfolders found in {folder}")
        return

    sep = "=" * 80

    for dwi_type, data in results.items():
        print(f"\n{sep}")
        print(f"  {dwi_type} \u2014 Parsed Log Metrics")
        print(sep)

        # ── Baseline quality ──
        bl = data["baseline"]
        if bl["outlier_flags"]:
            print("\n  Outlier Flags:")
            for o in bl["outlier_flags"]:
                print(f"    {o['metric']}: {o['n_flagged']} flagged (LF={o['n_lf']}, LC={o['n_lc']}, CR={o['n_cr']})")
        if bl["total_outliers"]:
            t = bl["total_outliers"]
            print(f"  Total outliers: {t['n_removed']}/{t['n_total']} ({t['pct']}%)")

        # ── GLME interaction tests ──
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

        # ── Predictive performance ──
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

        # ── Survival analysis ──
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
