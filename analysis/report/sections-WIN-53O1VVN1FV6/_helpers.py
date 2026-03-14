"""Shared helper functions for report section builders.

These are cross-cutting utilities used by multiple section functions
to avoid code duplication within the report_sections package.
"""

from __future__ import annotations

import json
import re

import sys
from pathlib import Path

# Ensure analysis/ root is on sys.path so sibling packages are importable.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from shared import DWI_TYPES  # type: ignore


def _safe_json_load(json_str: str, default=None):
    """Safely parse a JSON string, returning *default* on failure.

    Also rejects results that are not lists or dicts.
    """
    if default is None:
        default = []
    try:
        result = json.loads(str(json_str))
        if isinstance(result, (list, dict)):
            return result
    except (json.JSONDecodeError, TypeError, ValueError):
        pass
    return default


def _normalize_series_name(name: str) -> str:
    """Reduce a free-form vision-API series name to a canonical key.

    The Gemini vision API returns inconsistent series names across DWI
    types (e.g. ``"Mean D - Local Control"`` vs ``"Mean D (Local Control)"``
    vs ``"Mean D, Local Control"``).  This function strips separators,
    parenthetical qualifiers, and whitespace differences so that
    semantically identical series can be matched across types.

    Canonical form: lower-case tokens joined by single spaces.
    Known parameter aliases (``adc``, ``d*``/``dstar``, ``d``, ``f``,
    ``gtvvol``) are kept as-is; everything between common delimiters
    (``,``, ``-``, ``(``, ``)``  , ``vs``) is treated as separate tokens
    and reassembled in sorted order so that ``"LC vs D95 (Mean ADC)"``
    and ``"Mean ADC, LC, D95"`` yield the same key.
    """
    s = name.strip()
    # Normalise unicode delta/times symbols before lowering
    s = s.replace("\u0394", "delta").replace("\u00d7", "x")
    s = s.lower()
    # Also handle lowercase delta (from pre-lowered input)
    s = s.replace("\u03b4", "delta")
    # Replace common delimiters with a pipe for uniform splitting
    s = re.sub(r"\s*[,\-\(\)]\s*", " | ", s)
    s = re.sub(r"\s+vs\.?\s+", " | ", s)
    # Split into tokens, drop empties, sort for order-independence
    tokens = sorted(t.strip() for t in s.split("|") if t.strip())
    return " ".join(tokens)


def _build_normalised_series_map(
    all_trends: dict[str, list],
) -> dict[str, dict[str, tuple[str, str]]]:
    """Group trend entries by normalised series name.

    Parameters
    ----------
    all_trends : dict[str, list]
        ``{dwi_type: [trend_dict, ...]}`` as returned by JSON-parsing
        ``trends_json`` for each DWI type.

    Returns
    -------
    dict[str, dict[str, tuple[str, str]]]
        ``{norm_key: {dwi_type: (direction, description)}}``
        The ``norm_key`` is the canonical series name.  Each DWI type
        maps to its first matching direction and description.
    """
    merged: dict[str, dict[str, tuple[str, str]]] = {}
    # Also keep the best raw name per norm_key for display
    for dt, trends in all_trends.items():
        for t in trends:
            if not isinstance(t, dict):
                continue
            raw = t.get("series") or "overall"
            nk = _normalize_series_name(raw)
            direction = str(t.get("direction", ""))
            description = str(t.get("description", ""))
            merged.setdefault(nk, {})[dt] = (direction, description)
    return merged


# Keep a companion that returns display-friendly labels
def _best_display_name(
    all_trends: dict[str, list],
    norm_key: str,
) -> str:
    """Return the longest raw series name that normalises to *norm_key*."""
    best = norm_key
    for trends in all_trends.values():
        for t in trends:
            if not isinstance(t, dict):
                continue
            raw = t.get("series") or "overall"
            if _normalize_series_name(raw) == norm_key and len(raw) > len(best):
                best = raw
    return best


def _get_cohort_size(mat_data: dict | None) -> tuple[int, int, str]:
    """Extract (n_patients, n_timepoints, dwi_type) from MAT data.

    Returns (0, 0, "") if no cohort data is available.
    """
    if not mat_data:
        return 0, 0, ""
    for dt in DWI_TYPES:
        if dt in mat_data and "longitudinal" in mat_data[dt]:
            lon = mat_data[dt]["longitudinal"]
            n_pat = lon.get("num_patients", 0)
            n_tp = lon.get("num_timepoints", 0)
            if n_pat > 0:
                return n_pat, n_tp, dt  # type: ignore
    return 0, 0, ""


def _find_best_auc(log_data: dict | None, dwi_types_present: list[str]
                   ) -> tuple[dict | None, str]:
    """Find the best ROC/AUC result across all DWI types.

    Returns (best_roc_dict, dwi_type) or (None, "").
    """
    best_roc: dict | None = None
    best_dt = ""
    if not log_data:
        return None, ""
    for dt in dwi_types_present:
        dt_data = log_data.get(dt)
        if not dt_data:
            continue
        pred = dt_data.get("stats_predictive") or {}
        for roc in pred.get("roc_analyses", []):  # type: ignore
            auc = roc.get("auc", 0.0)  # type: ignore
            if best_roc is None or auc > best_roc.get("auc", 0.0):  # type: ignore
                best_roc = roc
                best_dt = dt
    return best_roc, best_dt


def _aggregate_dwi_statistics(
    log_data: dict | None,
    dwi_types_present: list[str],
    csv_data: dict | None = None,
) -> dict:
    """Aggregate key statistics across DWI types.

    Returns a dict with: total_sig, total_hrs, sig_hrs,
    best_auc_cards (list of (label, value) tuples).
    """
    total_sig = 0
    total_hrs = 0
    sig_hrs = 0
    auc_cards: list[tuple[str, float, str]] = []

    if log_data:
        for dt in dwi_types_present:
            dt_data = log_data.get(dt)
            if not dt_data:
                continue
            # Best AUC (track timepoint too)
            roc = dt_data.get("stats_predictive", {}).get("roc_analyses", [])
            best_roc_item = max(roc, key=lambda r: r.get("auc", 0), default=None) if roc else None
            best_auc = best_roc_item.get("auc", 0) if best_roc_item else 0
            best_tp = best_roc_item.get("timepoint", "") if best_roc_item else ""
            if best_auc > 0:
                auc_cards.append((dt, best_auc, best_tp))
            # GLME sig count
            sc = dt_data.get("stats_comparisons", {})
            total_sig += len([g for g in sc.get("glme_details", [])
                              if g["p"] < g["adj_alpha"]])
            # Hazard ratios
            hrs = dt_data.get("survival", {}).get("hazard_ratios", [])
            total_hrs += len(hrs)
            sig_hrs += len([hr for hr in hrs if hr.get("p", 1) < 0.05])

    return {
        "total_sig": total_sig,
        "total_hrs": total_hrs,
        "sig_hrs": sig_hrs,
        "auc_cards": auc_cards,
    }


def _compute_feature_overlap(
    log_data: dict | None,
    dwi_types_present: list[str],
) -> tuple[int, int]:
    """Compute cross-DWI feature overlap counts.

    Returns (total_shared, total_all) where total_shared is the number
    of features selected by >= 2 DWI types at any timepoint.
    """
    if not log_data or len(dwi_types_present) < 2:
        return 0, 0

    tp_features: dict[str, dict[str, list[str]]] = {}
    for dt in dwi_types_present:
        dt_data = log_data.get(dt)  # type: ignore
        if not dt_data:
            continue
        fs_list = dt_data.get("stats_predictive", {}).get("feature_selections", [])
        for fs in fs_list:
            tp = fs.get("timepoint", "?")
            tp_features.setdefault(tp, {})[dt] = fs.get("features", [])

    total_shared: int = 0
    total_all: int = 0
    for tp, dts_dict in tp_features.items():
        if len(dts_dict) < 2:
            continue
        all_feats: set[str] = set()
        for feats in dts_dict.values():
            all_feats.update(feats)
        total_all += len(all_feats)
        for feat in all_feats:
            if sum(1 for feats in dts_dict.values() if feat in feats) >= 2:
                total_shared = int(total_shared + 1)  # type: ignore

    return total_shared, total_all


def _aggregate_sanity_checks(
    log_data: dict | None,
    dwi_types_present: list[str],
) -> dict:
    """Aggregate sanity check flags across DWI types.

    Returns dict with: all_converged_count, total_conv_flags,
    total_dim_issues, sanity_types_checked.
    """
    all_converged_count: int = 0
    total_conv_flags: int = 0
    total_dim_issues: int = 0
    sanity_types_checked: int = 0

    if not log_data:
        return {
            "all_converged_count": 0,
            "total_conv_flags": 0,
            "total_dim_issues": 0,
            "sanity_types_checked": 0,
        }

    for dt in dwi_types_present:
        dt_data = log_data.get(dt)
        if not dt_data:
            continue
        san = dt_data.get("sanity_checks", {})
        if not san:
            continue
        sanity_types_checked = int(sanity_types_checked + 1)  # type: ignore
        if san.get("all_converged"):
            all_converged_count = int(all_converged_count + 1)  # type: ignore
        total_conv_flags = int(total_conv_flags + san.get("total_convergence", 0))  # type: ignore
        total_dim_issues = int(total_dim_issues + san.get("dim_mismatches", 0) + san.get("nan_dose_warnings", 0))  # type: ignore

    return {
        "all_converged_count": all_converged_count,
        "total_conv_flags": total_conv_flags,
        "total_dim_issues": total_dim_issues,
        "sanity_types_checked": sanity_types_checked,
    }


def _extract_significant_metrics(log_data: dict | None) -> dict:
    """Extract statistically significant metrics from all DWI types.

    Returns a dict with: sig_glme, sig_glme_details, sig_fdr_timepoints,
    sig_hr, best_roc, best_roc_dt, all_feature_selections, dosimetry, dosimetry_dt.
    """
    sig_glme: list[dict] = []
    sig_glme_details: list[dict] = []
    sig_fdr_timepoints: list[dict] = []
    sig_hr: list[dict] = []
    best_roc: dict | None = None
    best_roc_dt = ""
    all_feature_selections: list[dict] = []

    if log_data:
        for dt in DWI_TYPES:
            dt_data = log_data.get(dt)
            if not dt_data:
                continue
            # GLME interaction p-values
            stats_cmp = dt_data.get("stats_comparisons") or {}
            for p_val in stats_cmp.get("glme_interactions", []):  # type: ignore
                if p_val < 0.05:
                    sig_glme.append({"p": p_val, "dwi_type": dt})
            for detail in stats_cmp.get("glme_details", []):  # type: ignore
                if detail.get("p", 1.0) < detail.get("adj_alpha", 0.05):
                    sig_glme_details.append({**detail, "dwi_type": dt})
            for fdr_tp in stats_cmp.get("fdr_timepoints", []):  # type: ignore
                if fdr_tp.get("n_significant", 0) > 0:
                    sig_fdr_timepoints.append({**fdr_tp, "dwi_type": dt})
            # Survival hazard ratios
            surv = dt_data.get("survival") or {}
            for hr_entry in surv.get("hazard_ratios", []):  # type: ignore
                if hr_entry.get("p", 1.0) < 0.05:
                    sig_hr.append({**hr_entry, "dwi_type": dt})
            # Predictive model ROC/AUC
            pred = dt_data.get("stats_predictive") or {}
            for roc in pred.get("roc_analyses", []):  # type: ignore
                auc = roc.get("auc", 0.0)  # type: ignore
                if best_roc is None or auc > best_roc.get("auc", 0.0):  # type: ignore
                    best_roc = roc
                    best_roc_dt = dt
            for fs in pred.get("feature_selections", []):  # type: ignore
                if fs.get("features"):  # type: ignore
                    all_feature_selections.append({**fs, "dwi_type": dt})  # type: ignore

    return {
        "sig_glme": sig_glme,
        "sig_glme_details": sig_glme_details,
        "sig_fdr_timepoints": sig_fdr_timepoints,
        "sig_hr": sig_hr,
        "best_roc": best_roc,
        "best_roc_dt": best_roc_dt,
        "all_feature_selections": all_feature_selections,
    }


def _extract_dosimetry(mat_data: dict | None) -> tuple[dict, str]:
    """Extract first available dosimetry data from MAT data.

    Returns (dosimetry_dict, dwi_type) or ({}, "").
    """
    if not mat_data:
        return {}, ""
    for dt in DWI_TYPES:
        dosi = (mat_data.get(dt) or {}).get("dosimetry")
        if dosi:
            return dosi, dt
    return {}, ""


def _scalar_gy(val):
    """Return numeric Gy scalar from dict {mean,...} or plain float; None on failure."""
    if val is None:
        return None
    if isinstance(val, dict):
        m = val.get("mean")
        return m if (isinstance(m, (int, float)) and m == m) else None
    return val if (isinstance(val, (int, float)) and val == val) else None


def _compute_cross_dwi_trend_agreement(
    groups: dict | None,
    dwi_types_present: list[str],
) -> tuple[int, int, float]:
    """Compute cross-DWI trend agreement from Longitudinal_Mean_Metrics.

    Returns (n_agree, n_total_series, pct_agreement).
    """
    if not groups or "Longitudinal_Mean_Metrics" not in groups:
        return 0, 0, 0.0
    if len(dwi_types_present) < 2:
        return 0, 0, 0.0

    lmm = groups["Longitudinal_Mean_Metrics"]  # type: ignore
    all_trends: dict[str, list] = {}

    for dt, r in lmm.items():
        if dt == "Root":
            continue
        trends = _safe_json_load(r.get("trends_json", "[]"))
        if trends:
            all_trends[dt] = trends

    norm_map = _build_normalised_series_map(all_trends)

    n_agree: int = 0
    n_total: int = 0
    for norm_key, dt_entries in norm_map.items():
        dirs = {dt: d.lower() for dt, (d, _) in dt_entries.items() if d}
        if len(dirs) < 2:
            continue
        n_total += 1
        if len(set(dirs.values())) == 1:
            n_agree = int(n_agree + 1)  # type: ignore

    pct = (100 * float(n_agree) / float(n_total)) if n_total > 0 else 0.0  # type: ignore
    return n_agree, n_total, pct


def _compute_all_groups_trend_agreement(
    groups: dict | None,
    dwi_types_present: list[str],
) -> tuple[int, int, float]:
    """Compute cross-DWI trend agreement across ALL graph groups.

    Unlike :func:`_compute_cross_dwi_trend_agreement` which only checks
    Longitudinal_Mean_Metrics, this iterates over every graph group and
    counts agreement across all series with >= 2 DWI types.

    Returns (n_agree, n_total_series, pct_agreement).
    """
    if not groups or len(dwi_types_present) < 2:
        return 0, 0, 0.0

    n_agree: int = 0
    n_total_series: int = 0

    for base_name, dwi_dict in groups.items():  # type: ignore
        real = [t for t in dwi_dict if t != "Root"]
        if len(real) < 2:
            continue
        all_trends_dict: dict[str, list] = {}
        for dt_key in DWI_TYPES:
            if dt_key in dwi_dict:  # type: ignore
                trends = _safe_json_load(dwi_dict[dt_key].get("trends_json", "[]"))  # type: ignore
                if trends:
                    all_trends_dict[dt_key] = trends  # type: ignore
        if len(all_trends_dict) < 2:
            continue
        norm_map = _build_normalised_series_map(all_trends_dict)
        for norm_key, dt_entries in norm_map.items():
            dirs = {dt: d for dt, (d, _) in dt_entries.items()}
            if len(dirs) >= 2:
                n_total_series = int(n_total_series + 1)  # type: ignore
                if len(set(dirs.values())) == 1:
                    n_agree = int(n_agree + 1)  # type: ignore

    pct = (100 * float(n_agree) / float(n_total_series)) if n_total_series > 0 else 0.0  # type: ignore
    return n_agree, n_total_series, pct


def _extract_longitudinal_trend_consensus(
    groups: dict | None,
) -> tuple[str, str, list[str], list[str]]:
    """Extract D and f trend consensus from Longitudinal_Mean_Metrics.

    Parses trend directions from each DWI type's vision output and
    determines the consensus direction via keyword voting.

    Returns (d_consensus, f_consensus, d_trends_list, f_trends_list)
    where consensus values are 'increasing', 'decreasing', 'stable',
    or 'unknown'.
    """
    d_trends: list[str] = []
    f_trends: list[str] = []

    if not groups or "Longitudinal_Mean_Metrics" not in groups:
        return "unknown", "unknown", d_trends, f_trends

    for dt, r in groups["Longitudinal_Mean_Metrics"].items():  # type: ignore
        if dt == "Root":
            continue
        trends = _safe_json_load(r.get("trends_json", "[]"))
        for t in trends:
            if not isinstance(t, dict):
                continue
            series = t.get("series", "")
            direction = t.get("direction", "").lower()
            if series == "Mean D":
                d_trends.append(direction)
            elif series == "Mean f":
                f_trends.append(direction)

    # Use the same consensus logic as report_formatters._get_consensus
    from report.report_formatters import _get_consensus  # type: ignore
    return _get_consensus(d_trends), _get_consensus(f_trends), d_trends, f_trends
