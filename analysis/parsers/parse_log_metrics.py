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
- **sanity_checks**: Convergence flags (Inf/NaN/Neg counts), outlier
  detection, dimensional alignment, and excessive NaN warnings.

The parsed data is returned as nested dicts (keyed by DWI type) and printed
as a human-readable summary when the script is run directly.

Usage:
    python parse_log_metrics.py [saved_files_path]
"""

from __future__ import annotations

import sys
from pathlib import Path

from tqdm import tqdm  # type: ignore

# Ensure analysis/ root is on sys.path so 'shared' is importable when run as subprocess.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import DWI_TYPES, resolve_folder, setup_utf8_stdout  # type: ignore

# Compiled regexes and float/log helpers live in a sibling module so this
# file stays focused on the per-module parsing logic.  They are re-imported
# here (and re-exported) so the public API of ``parsers.parse_log_metrics``
# is unchanged: callers and tests can still ``from parsers.parse_log_metrics
# import _read_log`` (and any RE_* pattern) exactly as before.
from parsers.parse_log_patterns import (  # type: ignore  # noqa: F401
    RE_ALL_CONVERGED,
    RE_AUC,
    RE_BASELINE_EXCLUDED,
    RE_CONVERGENCE_ISSUE,
    RE_DIM_MISMATCH,
    RE_DIM_MISMATCH_TOTAL,
    RE_ELASTIC_NET,
    RE_EXCESSIVE_NAN,
    RE_FDR_TIMEPOINT,
    RE_FG_COMPARISON_ROW,
    RE_FG_COMPLETED,
    RE_FIRTH,
    RE_GLME_DETAIL,
    RE_GLME_EXCLUDED,
    RE_GLME_INTERACTION,
    RE_GLOBAL_LRT,
    RE_HR_ROW,
    RE_HR_ROW_BRACKET,
    RE_IMPUTATION_AUC,
    RE_IMPUTATION_CONCORDANCE,
    RE_IPCW_WEIGHTS,
    RE_LF_RATE,
    RE_NAN_DOSE_WARNINGS,
    RE_OUTLIER_FLAG,
    RE_ROC_HEADER,
    RE_SANITY_OUTLIER,
    RE_SANITY_OUTLIER_TOTAL,
    RE_SCHOENFELD_HEADER,
    RE_SCHOENFELD_ROW,
    RE_SENS_SPEC,
    RE_TD_PANEL,
    RE_TOTAL_CONVERGENCE,
    RE_TOTAL_OUTLIERS,
    RE_TV_BASE_COEF,
    RE_TV_INTERACTION_COEF,
    RE_TV_STRATIFIED,
    RE_TV_VIOLATED,
    RE_YOUDEN,
    _parse_float,
    _read_log,
)

# ``parse_survival`` and ``parse_sanity_checks`` live in a sibling module to
# keep this file under the project length limit.  Re-imported here so they
# remain part of the public ``parsers.parse_log_metrics`` API.
from parsers.parse_log_survival_sanity import (  # type: ignore  # noqa: F401
    parse_sanity_checks,
    parse_survival,
)

setup_utf8_stdout()


# ── Per-module parsers ───────────────────────────────────────────────────────


def parse_stats_comparisons(text: str, log_path: str = "") -> dict:
    """Parse ``metrics_stats_comparisons`` log output.

    Extracts GLME interaction p-values, per-metric detail rows (with
    BH-adjusted alpha), FDR-significant timepoint counts, and
    competing-risk exclusion statistics.

    Parameters
    ----------
    text : str
        Full text of the ``metrics_stats_comparisons_output_*.txt`` log.
    log_path : str, optional
        Path of the log file being parsed, used in parse-failure warning
        messages so callers can identify which file had no matches.

    Returns
    -------
    dict
        Keys: ``glme_interactions``, ``glme_details``, ``fdr_timepoints``,
        ``glme_excluded``, ``parse_warnings``.
    """
    result: dict = {
        "glme_interactions": [],
        "glme_details": [],
        "fdr_timepoints": [],
        "glme_excluded": None,
        "parse_warnings": [],
    }

    # Collect raw interaction p-values (one per GLME model).
    for m in RE_GLME_INTERACTION.finditer(text):
        result["glme_interactions"].append(_parse_float(m.group(1)))  # type: ignore

    # Warn when a log file exists but yields no GLME interaction matches —
    # this typically indicates the log format has changed.
    if text and not result["glme_interactions"]:
        label = log_path or "metrics_stats_comparisons log"
        result["parse_warnings"].append(
            f"GLME interaction p-values not found in {label} (format may have changed)"
        )

    # Collect per-metric detail rows with raw p and adjusted alpha.
    for m in RE_GLME_DETAIL.finditer(text):
        result["glme_details"].append({  # type: ignore
            "metric": m.group(1).strip(),
            "p": _parse_float(m.group(2)),
            "adj_alpha": _parse_float(m.group(3)),
        })

    # Collect per-timepoint FDR significance counts.
    for m in RE_FDR_TIMEPOINT.finditer(text):
        result["fdr_timepoints"].append({  # type: ignore
            "timepoint": m.group(1),
            "n_significant": int(m.group(2)),
        })

    if text and not result["fdr_timepoints"]:
        label = log_path or "metrics_stats_comparisons log"
        result["parse_warnings"].append(
            f"FDR timepoint counts not found in {label} (format may have changed)"
        )

    # Competing-risk exclusion line (appears at most once).
    m = RE_GLME_EXCLUDED.search(text)
    if m:
        n_excluded = int(m.group(1))
        n_total = int(m.group(2))
        pct = float(m.group(3))
        # Validate that exclusion count does not exceed total.
        if n_excluded > n_total:
            label = log_path or "metrics_stats_comparisons log"
            result["parse_warnings"].append(
                f"Competing-risk exclusion count {n_excluded} > total {n_total} in {label}"
            )
        result["glme_excluded"] = {  # type: ignore
            "n_excluded": n_excluded,
            "n_total": n_total,
            "pct": pct,
        }

    return result


def parse_stats_predictive(text: str) -> dict:
    """Parse ``metrics_stats_predictive`` log output.

    Extracts elastic-net feature selections (with optimal lambda),
    ROC analysis blocks (AUC, Youden cutoff, sensitivity, specificity),
    and imputation sensitivity results.

    Parameters
    ----------
    text : str
        Full text of the ``metrics_stats_predictive_output_*.txt`` log.

    Returns
    -------
    dict
        Keys: ``feature_selections``, ``roc_analyses``, ``firth_refits``,
        ``imputation_sensitivity``.
    """
    result = {"feature_selections": [], "roc_analyses": [], "firth_refits": []}

    # Extract Firth penalised-likelihood refits per timepoint.
    for m in RE_FIRTH.finditer(text):
        result["firth_refits"].append({
            "timepoint": m.group(1),
            "n_features": int(m.group(2)),
        })

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
        block = text[start:end]  # type: ignore

        entry = {"timepoint": hdr.group(1)}
        m = RE_AUC.search(block)
        if m:
            entry["auc"] = _parse_float(m.group(1))
        m = RE_YOUDEN.search(block)
        if m:
            entry["youden_cutoff"] = _parse_float(m.group(1))
        m = RE_SENS_SPEC.search(block)
        if m:
            entry["sensitivity"] = float(m.group(1))
            entry["specificity"] = float(m.group(2))
        result["roc_analyses"].append(entry)

    # ── Imputation Sensitivity ──
    # Parses the comparison table printed by imputation_sensitivity.m.
    imp_methods = list(RE_IMPUTATION_AUC.finditer(text))
    if imp_methods:
        imp_results: list[dict] = []
        for m in imp_methods:
            imp_results.append({
                "method": m.group(1),
                "auc": _parse_float(m.group(2)),
                "n_imputed": int(m.group(3)),
            })
        result["imputation_sensitivity"] = imp_results
    else:
        result["imputation_sensitivity"] = []

    return result


def parse_baseline(text: str, log_path: str = "") -> dict:
    """Parse ``metrics_baseline`` log output.

    Extracts per-metric outlier flags with outcome breakdown (LF/LC/CR),
    total outlier removal statistics, and baseline exclusion counts with
    associated LF-rate comparison.

    Parameters
    ----------
    text : str
        Full text of the ``metrics_baseline_output_*.txt`` log.
    log_path : str, optional
        Path of the log file, used in parse-failure warning messages.

    Returns
    -------
    dict
        Keys: ``outlier_flags``, ``total_outliers``, ``baseline_exclusion``,
        ``parse_warnings``.
    """
    result: dict = {
        "outlier_flags": [],
        "total_outliers": None,
        "baseline_exclusion": None,
        "parse_warnings": [],
    }

    # Per-metric outlier flags with LF/LC/CR breakdown.
    for m in RE_OUTLIER_FLAG.finditer(text):
        n_flagged = int(m.group(2))
        n_lf = int(m.group(3))
        n_lc = int(m.group(4))
        n_cr = int(m.group(5))
        flag_entry: dict = {
            "metric": m.group(1),
            "n_flagged": n_flagged,
            "n_lf": n_lf,
            "n_lc": n_lc,
            "n_cr": n_cr,
        }
        # Validate that per-group breakdown doesn't exceed the reported total.
        group_sum = n_lf + n_lc + n_cr
        if group_sum > n_flagged:
            label = log_path or "metrics_baseline log"
            result["parse_warnings"].append(
                f"Outlier group sum ({n_lf}+{n_lc}+{n_cr}={group_sum}) "
                f"> n_flagged ({n_flagged}) for metric '{m.group(1)}' in {label}"
            )
        result["outlier_flags"].append(flag_entry)  # type: ignore

    # Aggregate outlier removal line.
    m = RE_TOTAL_OUTLIERS.search(text)
    if m:
        result["total_outliers"] = {  # type: ignore
            "n_removed": int(m.group(1)),
            "n_total": int(m.group(2)),
            "pct": float(m.group(3)),
        }

    # Baseline exclusion count and optional LF-rate comparison.
    m = RE_BASELINE_EXCLUDED.search(text)
    if m:
        result["baseline_exclusion"] = {  # type: ignore
            "n_excluded": int(m.group(1)),
            "n_total": int(m.group(2)),
        }
        # Look for the LF-rate comparison line that follows exclusion info.
        m2 = RE_LF_RATE.search(text)
        if m2:
            result["baseline_exclusion"]["lf_rate_included"] = float(m2.group(1))  # type: ignore
            result["baseline_exclusion"]["lf_rate_excluded"] = float(m2.group(2))  # type: ignore

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
    # Determine which DWI type subfolders exist for accurate progress.
    active_types = [d for d in DWI_TYPES if (folder / d).is_dir()]
    pbar = tqdm(
        active_types,
        desc="Parsing logs",
        unit="type",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}] {postfix}",
    )
    for dwi_type in pbar:
        dwi_dir = folder / dwi_type
        pbar.set_postfix_str(dwi_type, refresh=True)

        # Each log file follows the naming convention:
        #   <module_name>_output_<DWI_type>.txt
        # Pass the file path string to parsers so parse-failure warnings
        # include a useful location for debugging.
        sc_path = dwi_dir / f"metrics_stats_comparisons_output_{dwi_type}.txt"
        sp_path = dwi_dir / f"metrics_stats_predictive_output_{dwi_type}.txt"
        sv_path = dwi_dir / f"metrics_survival_output_{dwi_type}.txt"
        bl_path = dwi_dir / f"metrics_baseline_output_{dwi_type}.txt"

        results[dwi_type] = {
            "stats_comparisons": parse_stats_comparisons(
                _read_log(sc_path), log_path=str(sc_path)
            ),
            "stats_predictive": parse_stats_predictive(
                _read_log(sp_path)
            ),
            "survival": parse_survival(
                _read_log(sv_path), log_path=str(sv_path)
            ),
            "baseline": parse_baseline(
                _read_log(bl_path), log_path=str(bl_path)
            ),
            "sanity_checks": parse_sanity_checks(
                _read_log(dwi_dir / "sanity_checks_output.txt")
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

        # ── Time-Varying Cox ──
        tv = sv.get("time_varying_cox")
        if tv:
            print(f"\n  Time-Varying Cox (PH violations: {', '.join(tv['violated_covariates'])}):")
            if tv.get("stratified_by"):
                print(f"    Stratified by: {tv['stratified_by']}")
            for im in tv.get("interaction_models", []):
                base_p = im.get("base_p", float("nan"))
                print(f"    {im['covariate']}: base_coef={im.get('base_coef', 'N/A')}, "
                      f"interaction_coef={im['interaction_coef']:.4f}, "
                      f"interaction_p={im['interaction_p']:.4f}")

        # ── Imputation Sensitivity ──
        imp = sp.get("imputation_sensitivity", [])
        if imp:
            print("\n  Imputation Sensitivity:")
            for entry in imp:
                print(f"    {entry['method']}: AUC={entry['auc']:.3f}, N_Imputed={entry['n_imputed']}")

        # ── Sanity checks ──
        san = data.get("sanity_checks", {})
        if san:
            if san.get("all_converged"):
                print("\n  Convergence: All values converged")
            elif san.get("total_convergence", 0) > 0:
                print(f"\n  Convergence: {san['total_convergence']} flags raised")
            if san.get("dim_mismatches", 0) > 0:
                print(f"  Dimensional mismatches: {san['dim_mismatches']}")
            if san.get("nan_dose_warnings", 0) > 0:
                print(f"  NaN dose warnings: {san['nan_dose_warnings']}")
            if san.get("excessive_nan"):
                print("  Excessive NaN parameters:")
                for en in san["excessive_nan"]:
                    print(f"    {en['parameter']}: {en['pct_nan']:.1f}%")

    print()


if __name__ == "__main__":
    main()
