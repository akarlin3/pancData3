#!/usr/bin/env python3
"""Compiled regular expressions and float helpers for log-metric parsing.

This module holds the targeted regular expressions used by
:mod:`parsers.parse_log_metrics` to extract quantitative results from
MATLAB pipeline log files, together with the small string/float helper
utilities shared across the per-module parsers.  It is a pure data /
helper module with no side effects so it can be imported safely from the
main parser as well as from sibling parser modules.
"""

from __future__ import annotations

import re
from pathlib import Path

# ── Regex patterns ────────────────────────────────────────────────────────────
# Each compiled regex targets a specific line format produced by the
# corresponding MATLAB core module.  Capture groups extract the numeric
# values (and labels) that are assembled into structured dicts.

# ----- metrics_stats_comparisons -----

# Matches: "Interaction P-Value (LF vs LC): 0.0234" or "... Inf" / "... NaN"
# Captures the numeric p-value or special float word (group 1).
RE_GLME_INTERACTION = re.compile(
    r"Interaction P-Value.*?:\s*([0-9.]+(?:e[+-]?\d+)?|Inf|NaN|inf|nan)",
    re.IGNORECASE,
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
    r"Timepoint:\s*(\S+)\s*—\s*(\d+)\s+significant"
)

# Matches: "Excluded 4/25 (16.0%) competing-risk"
# Captures counts and percentage of patients excluded from GLME analysis.
RE_GLME_EXCLUDED = re.compile(
    r"Excluded\s+(\d+)/(\d+)\s+\(([0-9.]+)%\)\s+competing-risk"
)

# ----- metrics_stats_predictive -----

# Matches: "Elastic Net Selected Features for Fx10 (Opt Lambda=0.0523): feat1, feat2"
# Lambda also accepts Inf/NaN for degenerate cases.
# Captures timepoint (group 1), optimal lambda (group 2), and comma-separated
# feature names (group 3).
RE_ELASTIC_NET = re.compile(
    r"Elastic Net Selected Features for (\S+)\s+\(Opt Lambda=([0-9.]+(?:e[+-]?\d+)?|Inf|NaN)\):\s+(.*)",
    re.IGNORECASE,
)

# Matches: "Firth refit successful for Fx5 (3 features)"
# Captures timepoint (group 1) and number of features (group 2).
RE_FIRTH = re.compile(r"Firth refit successful for (\S+)\s+\((\d+) features\)")

# Matches: "PRIMARY ROC ANALYSIS ... for Fx10"
# Captures the timepoint label (group 1); used to delimit ROC blocks.
RE_ROC_HEADER = re.compile(r"PRIMARY ROC ANALYSIS.*for (\S+)")

# Matches: "AUC = 0.843" (also Inf/NaN for edge cases)
RE_AUC = re.compile(r"AUC\s*=\s*([0-9.]+(?:e[+-]?\d+)?|Inf|NaN)", re.IGNORECASE)

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

# Matches bracket-format Cox PH rows produced by metrics_survival.m, e.g.:
#   "  ADC           0.768    2.156 [ 0.28-16.70]   0.4621"
# Format: %-10s %8.3f %8.3f [%5.2f-%5.2f] %8.4f
# Captures: feature (1), coeff (2), HR (3), CI_lo (4), CI_hi (5), p (6).
RE_HR_ROW_BRACKET = re.compile(
    r"^\s{2,}(\S+)\s+([-\d.]+)\s+([\d.]+)\s+\[\s*([\d.]+)-\s*([\d.]+)\]\s+([\d.]+)",
    re.MULTILINE,
)

# Matches Schoenfeld residuals PH test table rows, e.g.:
#   "  ADC         0.0909    0.0826    0.7737             "
#   "  D           0.7091    5.0281    0.0249          ***"
# Captures: covariate (1), rho (2), chi2 (3), p-value (4), violation marker (5).
RE_SCHOENFELD_ROW = re.compile(
    r"^\s{2,}(\S+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)\s*(\*{3})?\s*$",
    re.MULTILINE,
)

# Matches the Schoenfeld section header to confirm its presence.
RE_SCHOENFELD_HEADER = re.compile(
    r"Schoenfeld Residuals.*PH Assumption Test",
)

# Matches TD Panel summary line from build_td_panel.m, e.g.:
#   "  [TD Panel] 35 patients → 210 intervals (10 events of interest, 0 competing)"
RE_TD_PANEL = re.compile(
    r"\[TD Panel\]\s*(\d+)\s*patients?\s*→\s*(\d+)\s*intervals?\s*"
    r"\((\d+)\s*events?\s*of interest,\s*(\d+)\s*competing\)"
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

# ----- sanity_checks -----

# Matches convergence issue lines, e.g.:
#   "Patient Pat001  Fx1  ADC : Inf=5/100 (5.0%)  NaN=10/100 (10.0%)"
# Captures patient ID (1), timepoint (2), parameter (3).
RE_CONVERGENCE_ISSUE = re.compile(
    r"Patient\s+(\S+)\s+(\S+)\s+(.+?)\s*:"
    r"(?:\s+Inf=(\d+)/(\d+)(?:\s+\([0-9.]+%\))?)?"
    r"(?:\s+NaN=(\d+)/(\d+)(?:\s+\([0-9.]+%\))?)?"
    r"(?:\s+Neg=(\d+)/(\d+)(?:\s+\([0-9.]+%\))?)?"
)

# Matches: "Total convergence flags raised: 15"
RE_TOTAL_CONVERGENCE = re.compile(
    r"Total convergence flags raised:\s*(\d+)"
)

# Matches: "All voxel-level fit values are finite, non-NaN, and non-negative."
RE_ALL_CONVERGED = re.compile(
    r"All voxel-level fit values are finite"
)

# Matches outlier lines: "OUTLIER: Patient Pat015  Fx3  ADC_mean = 2.5e-3 ..."
RE_SANITY_OUTLIER = re.compile(
    r"OUTLIER:\s+Patient\s+(\S+)\s+(\S+)\s+(\w+)\s*=\s*([0-9.eE+-]+)"
)

# Matches: "Total outlier flags: 5"
RE_SANITY_OUTLIER_TOTAL = re.compile(
    r"Total outlier flags:\s*(\d+)"
)

# Matches dimensional mismatch lines
RE_DIM_MISMATCH = re.compile(
    r"MISMATCH:\s+Patient\s+(\S+)\s+(\S+)"
)

# Matches: "Dimensional mismatches: 3"
RE_DIM_MISMATCH_TOTAL = re.compile(
    r"Dimensional mismatches:\s*(\d+)"
)

# Matches: "NaN dose warnings: 2"
RE_NAN_DOSE_WARNINGS = re.compile(
    r"NaN dose warnings:\s*(\d+)"
)

# Matches excessive NaN lines:
#   "Excessive NaN fraction in D: 65.0% (threshold: 50%)"
RE_EXCESSIVE_NAN = re.compile(
    r"Excessive NaN fraction in (\w[\w*]*?):\s*([0-9.]+)%"
)

# ----- imputation sensitivity -----

# Matches imputation sensitivity table rows, e.g.:
#   "KNN                     0.843            42"
# The table is printed by imputation_sensitivity.m inside the
# metrics_stats_predictive log file.  Captures method name (1), AUC (2),
# and number of imputed values (3).
RE_IMPUTATION_AUC = re.compile(
    r"^\s{2}(KNN|LOCF|Mean|Linear_Interp)\s+([0-9.]+(?:e[+-]?\d+)?|NaN|Inf)\s+(\d+)\s*$",
    re.MULTILINE | re.IGNORECASE,
)

# Matches the concordance matrix header/rows printed after the AUC table:
#   "KNN              1.000     0.923     0.890     0.912"
# Captures method name (1) and the rest of the row as a single string (2).
RE_IMPUTATION_CONCORDANCE = re.compile(
    r"^\s{2}(KNN|LOCF|Mean|Linear_Interp)\s+((?:[0-9.]+\s*)+)$",
    re.MULTILINE | re.IGNORECASE,
)

# ----- time-varying Cox -----

# Matches PH violation covariate list:
#   "PH violations detected for: mean_adc, delta_d"
RE_TV_VIOLATED = re.compile(
    r"PH violations detected for:\s*(.+)",
)

# Matches extended Cox base coefficient lines:
#   "    Base mean_adc: coef=0.1234, p=0.0456"
RE_TV_BASE_COEF = re.compile(
    r"Base\s+(\S+):\s*coef=([0-9.eE+-]+),\s*p=([0-9.eE+-]+)"
)

# Matches extended Cox interaction coefficient lines:
#   "    mean_adc × log(t): coef=-0.0567, p=0.0123"
RE_TV_INTERACTION_COEF = re.compile(
    r"(\S+)\s*×\s*log\(t\):\s*coef=([0-9.eE+-]+),\s*p=([0-9.eE+-]+)"
)

# Matches stratified Cox model header:
#   "Stratified Cox Model (stratified by mean_adc, median=0.0012):"
RE_TV_STRATIFIED = re.compile(
    r"Stratified Cox Model \(stratified by (\S+),\s*median=([0-9.eE+-]+)\)"
)

# ----- Fine-Gray competing risks -----

# Matches Fine-Gray comparison table rows (CSH HR | sHR format):
#   "  ADC         1.234    0.890    1.456    0.0231  |     1.198    0.856    1.423    0.0345"
RE_FG_COMPARISON_ROW = re.compile(
    r"^\s*([\w*]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*\|\s*"
    r"([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$",
    re.MULTILINE,
)

# Matches: "Fine-Gray model completed: 5 competing events, 8 primary events."
RE_FG_COMPLETED = re.compile(
    r"Fine-Gray model completed:\s*(\d+)\s*competing events?,\s*(\d+)\s*primary events?"
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


def _parse_float(s: str) -> float:
    """Convert a string to float, handling ``Inf`` and ``NaN`` spellings.

    Parameters
    ----------
    s : str
        Numeric string, optionally ``"Inf"``, ``"inf"``, ``"NaN"``, or
        ``"nan"`` (case-insensitive).

    Returns
    -------
    float
        Parsed value including ``float('inf')`` or ``float('nan')``.
    """
    lower = s.lower()
    if lower in ("inf", "+inf"):
        return float("inf")
    if lower in ("-inf",):
        return float("-inf")
    if lower == "nan":
        return float("nan")
    return float(s)
