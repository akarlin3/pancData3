#!/usr/bin/env python3
"""Survival and sanity-check log parsers for the analysis suite.

These two per-module parsers were extracted from
:mod:`parsers.parse_log_metrics` to keep that module under the project
file-length limit.  They are re-imported back into ``parse_log_metrics``
so its public API is unchanged.

- :func:`parse_survival` -- Cox PH hazard ratios with confidence intervals,
  global likelihood-ratio test, IPCW weight ranges, Fine-Gray competing
  risks, and time-varying Cox model results.
- :func:`parse_sanity_checks` -- convergence flags (Inf/NaN/Neg counts),
  outlier detection, dimensional alignment, and excessive NaN warnings.

This module is import-only (no side effects at import time) so it can be
imported safely from the main parser.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure analysis/ root is on sys.path so 'parsers' is importable when this
# module is imported from a subprocess-executed sibling.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from parsers.parse_log_patterns import (  # type: ignore
    RE_ALL_CONVERGED,
    RE_CONVERGENCE_ISSUE,
    RE_DIM_MISMATCH_TOTAL,
    RE_EXCESSIVE_NAN,
    RE_FG_COMPARISON_ROW,
    RE_FG_COMPLETED,
    RE_GLOBAL_LRT,
    RE_HR_ROW,
    RE_HR_ROW_BRACKET,
    RE_IPCW_WEIGHTS,
    RE_NAN_DOSE_WARNINGS,
    RE_SANITY_OUTLIER,
    RE_SANITY_OUTLIER_TOTAL,
    RE_SCHOENFELD_HEADER,
    RE_SCHOENFELD_ROW,
    RE_TD_PANEL,
    RE_TOTAL_CONVERGENCE,
    RE_TV_BASE_COEF,
    RE_TV_INTERACTION_COEF,
    RE_TV_STRATIFIED,
    RE_TV_VIOLATED,
)


def parse_survival(text: str, log_path: str = "") -> dict:
    """Parse ``metrics_survival`` log output.

    Extracts Cox proportional-hazards table rows (covariate, HR, 95% CI,
    p-value), the global likelihood-ratio test, IPCW weight ranges, and
    time-varying Cox model results.

    Parameters
    ----------
    text : str
        Full text of the ``metrics_survival_output_*.txt`` log.
    log_path : str, optional
        Path of the log file, used in parse-failure warning messages.

    Returns
    -------
    dict
        Keys: ``hazard_ratios``, ``cox_covariates``, ``global_lrt``,
        ``ipcw``, ``fine_gray``, ``time_varying_cox``, ``ph_tests``,
        ``schoenfeld_tested``, ``schoenfeld_violated``,
        ``time_varying_cox_fitted``, ``n_intervals``, ``n_events``,
        ``n_patients``, ``n_competing``, ``parse_warnings``.
    """
    result: dict = {
        "hazard_ratios": [],
        "cox_covariates": [],
        "global_lrt": None,
        "ipcw": None,
        "fine_gray": None,
        "time_varying_cox": None,
        "ph_tests": [],
        "schoenfeld_tested": False,
        "schoenfeld_violated": [],
        "time_varying_cox_fitted": False,
        "n_intervals": None,
        "n_events": None,
        "n_patients": None,
        "n_competing": None,
        "parse_warnings": [],
    }

    # Parse whitespace-delimited Cox PH table rows (legacy format).
    for m in RE_HR_ROW.finditer(text):
        result["hazard_ratios"].append({  # type: ignore
            "covariate": m.group(1),
            "hr": float(m.group(2)),
            "ci_lo": float(m.group(3)),
            "ci_hi": float(m.group(4)),
            "p": float(m.group(5)),
        })

    # Parse bracket-format Cox PH table rows (current MATLAB output).
    for m in RE_HR_ROW_BRACKET.finditer(text):
        entry = {
            "name": m.group(1),
            "coeff": float(m.group(2)),
            "hr": float(m.group(3)),
            "ci_lo": float(m.group(4)),
            "ci_hi": float(m.group(5)),
            "p": float(m.group(6)),
        }
        result["cox_covariates"].append(entry)  # type: ignore
        # Also populate hazard_ratios for backward compatibility.
        result["hazard_ratios"].append({  # type: ignore
            "covariate": entry["name"],
            "hr": entry["hr"],
            "ci_lo": entry["ci_lo"],
            "ci_hi": entry["ci_hi"],
            "p": entry["p"],
        })

    # Parse Schoenfeld residuals PH test table.
    if RE_SCHOENFELD_HEADER.search(text):
        result["schoenfeld_tested"] = True
        for m in RE_SCHOENFELD_ROW.finditer(text):
            violated = m.group(5) is not None
            entry = {
                "name": m.group(1),
                "rho": float(m.group(2)),
                "chi2": float(m.group(3)),
                "p": float(m.group(4)),
                "violated": violated,
            }
            result["ph_tests"].append(entry)  # type: ignore
            if violated:
                result["schoenfeld_violated"].append(entry["name"])  # type: ignore

    # Parse TD Panel summary line.
    m_td = RE_TD_PANEL.search(text)
    if m_td:
        result["n_patients"] = int(m_td.group(1))
        result["n_intervals"] = int(m_td.group(2))
        result["n_events"] = int(m_td.group(3))
        result["n_competing"] = int(m_td.group(4))

    # Global likelihood-ratio test (single occurrence).
    m = RE_GLOBAL_LRT.search(text)
    if m:
        result["global_lrt"] = {  # type: ignore
            "df": int(m.group(1)),
            "chi2": float(m.group(2)),
            "p": float(m.group(3)),
        }

    # IPCW (Inverse Probability of Censoring Weighting) range.
    m = RE_IPCW_WEIGHTS.search(text)
    if m:
        w_min = float(m.group(1))
        w_max = float(m.group(2))
        # Compute range ratio: how much weights vary.  A ratio > 5 suggests
        # extreme censoring imbalance and warrants a warning.
        range_ratio: float | None = None
        if w_min > 0:
            range_ratio = float(f"{w_max / w_min:.4f}")  # type: ignore
            if range_ratio > 5:
                label = log_path or "metrics_survival log"
                result["parse_warnings"].append(
                    f"IPCW weight range ratio {range_ratio:.2f} > 5 in {label}; "
                    "extreme censoring imbalance suspected"
                )
        result["ipcw"] = {  # type: ignore
            "min_weight": w_min,
            "max_weight": w_max,
            "weight_range_ratio": range_ratio,
        }

    # ── Fine-Gray Competing Risks (v2.2) ──
    m_fg = RE_FG_COMPLETED.search(text)
    if m_fg:
        fg_data: dict = {
            "n_competing": int(m_fg.group(1)),
            "n_primary": int(m_fg.group(2)),
            "comparison_table": [],
        }
        for m_row in RE_FG_COMPARISON_ROW.finditer(text):
            fg_data["comparison_table"].append({
                "covariate": m_row.group(1),
                "csh_hr": float(m_row.group(2)),
                "csh_ci_lo": float(m_row.group(3)),
                "csh_ci_hi": float(m_row.group(4)),
                "csh_p": float(m_row.group(5)),
                "shr": float(m_row.group(6)),
                "shr_ci_lo": float(m_row.group(7)),
                "shr_ci_hi": float(m_row.group(8)),
                "shr_p": float(m_row.group(9)),
            })
        result["fine_gray"] = fg_data  # type: ignore

    # ── Time-Varying Cox ──
    # Parses output from fit_time_varying_cox.m in the metrics_survival log.
    m = RE_TV_VIOLATED.search(text)
    if m:
        violated_names = [n.strip() for n in m.group(1).split(",") if n.strip()]
        tv_data: dict = {
            "violated_covariates": violated_names,
            "interaction_models": [],
            "stratified_by": None,
        }

        # Stratified model info
        m_strat = RE_TV_STRATIFIED.search(text)
        if m_strat:
            tv_data["stratified_by"] = m_strat.group(1)

        # Extended Cox interaction models: pair base + interaction lines
        base_matches = {m.group(1): m for m in RE_TV_BASE_COEF.finditer(text)}
        for m_int in RE_TV_INTERACTION_COEF.finditer(text):
            cov_name = m_int.group(1)
            entry: dict = {
                "covariate": cov_name,
                "interaction_coef": float(m_int.group(2)),
                "interaction_p": float(m_int.group(3)),
            }
            m_base = base_matches.get(cov_name)
            if m_base:
                entry["base_coef"] = float(m_base.group(2))
                entry["base_p"] = float(m_base.group(3))
            tv_data["interaction_models"].append(entry)

        result["time_varying_cox"] = tv_data  # type: ignore
        result["time_varying_cox_fitted"] = True

    return result


def parse_sanity_checks(text: str) -> dict:
    """Parse ``sanity_checks`` log output.

    Extracts convergence issue counts, outlier flags, dimensional
    alignment results, and excessive NaN warnings.

    Parameters
    ----------
    text : str
        Full text of the ``sanity_checks_output.txt`` log.

    Returns
    -------
    dict
        Keys: ``convergence_flags``, ``total_convergence``,
        ``all_converged``, ``outliers``, ``total_outliers``,
        ``dim_mismatches``, ``nan_dose_warnings``, ``excessive_nan``.
    """
    result: dict = {
        "convergence_flags": [],
        "total_convergence": 0,
        "all_converged": False,
        "outliers": [],
        "total_outliers": 0,
        "dim_mismatches": 0,
        "nan_dose_warnings": 0,
        "excessive_nan": [],
    }

    # Check if all values converged (no issues).
    if RE_ALL_CONVERGED.search(text):
        result["all_converged"] = True

    # Total convergence flags count.
    m = RE_TOTAL_CONVERGENCE.search(text)
    if m:
        result["total_convergence"] = int(m.group(1))

    # Individual convergence issue lines.
    for m in RE_CONVERGENCE_ISSUE.finditer(text):
        entry: dict = {
            "patient": m.group(1),
            "timepoint": m.group(2),
            "parameter": m.group(3).strip(),
        }
        if m.group(4) is not None:
            entry["n_inf"] = int(m.group(4))
            entry["n_total_inf"] = int(m.group(5))
        if m.group(6) is not None:
            entry["n_nan"] = int(m.group(6))
            entry["n_total_nan"] = int(m.group(7))
        if m.group(8) is not None:
            entry["n_neg"] = int(m.group(8))
            entry["n_total_neg"] = int(m.group(9))
        result["convergence_flags"].append(entry)

    # Outlier detections.
    for m in RE_SANITY_OUTLIER.finditer(text):
        result["outliers"].append({
            "patient": m.group(1),
            "timepoint": m.group(2),
            "metric": m.group(3),
            "value": float(m.group(4)),
        })

    m = RE_SANITY_OUTLIER_TOTAL.search(text)
    if m:
        result["total_outliers"] = int(m.group(1))

    # Dimensional mismatches.
    m = RE_DIM_MISMATCH_TOTAL.search(text)
    if m:
        result["dim_mismatches"] = int(m.group(1))

    # NaN dose warnings.
    m = RE_NAN_DOSE_WARNINGS.search(text)
    if m:
        result["nan_dose_warnings"] = int(m.group(1))

    # Excessive NaN parameters.
    for m in RE_EXCESSIVE_NAN.finditer(text):
        result["excessive_nan"].append({
            "parameter": m.group(1),
            "pct_nan": float(m.group(2)),
        })

    return result
