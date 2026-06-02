#!/usr/bin/env python3
"""Format and data-extraction helpers for the interactive HTML report.

Split out from :mod:`report.generate_interactive_report`; the names are
re-exported there for backward compatibility.
"""

from __future__ import annotations

import html as _html


def _esc(text: str) -> str:
    """HTML-escape a string."""
    return _html.escape(str(text))


def _dwi_badge(dwi_type: str) -> str:
    """Return a colour-coded DWI badge span."""
    cls = {
        "Standard": "badge-standard", "dnCNN": "badge-dncnn",
        "IVIMnet": "badge-ivimnet",
    }.get(dwi_type, "badge-standard")
    return f'<span class="badge {cls}">{_esc(dwi_type)}</span>'


def _sig_class(p: float) -> str:
    """Return CSS class for a p-value significance level."""
    if p < 0.001:
        return "sig-3"
    if p < 0.01:
        return "sig-2"
    if p < 0.05:
        return "sig-1"
    return ""


def _trend_tag(direction: str) -> str:
    """Return a directional arrow badge for a trend."""
    d = direction.lower()
    if "increas" in d or "up" in d or "higher" in d or "rising" in d:
        cls, arrow = "trend-incr", "\u2191\u00a0"
    elif "decreas" in d or "down" in d or "lower" in d or "falling" in d or "drop" in d:
        cls, arrow = "trend-decr", "\u2193\u00a0"
    elif "flat" in d or "stable" in d or "constant" in d:
        cls, arrow = "trend-flat", "\u2192\u00a0"
    else:
        cls, arrow = "trend-nm", ""
    return f'<span class="trend-tag {cls}">{arrow}{_esc(direction)}</span>'


# ── Data extraction helpers ──────────────────────────────────────────────────

def _extract_patients(log_data: dict | None, mat_data: dict) -> list[dict]:
    """Build a list of patient dicts with per-timepoint metrics.

    Each patient dict has ``id``, ``outcome``, and ``timepoints`` (a list of
    dicts with ``dwi_type``, ``label``, and metric values).
    """
    patients: dict[str, dict] = {}

    # From log data: baseline metrics
    if log_data:
        for dwi_type, ddata in log_data.items():
            if not isinstance(ddata, dict):
                continue
            bl = ddata.get("baseline", {})
            for metric_name, metric_val in bl.items():
                if isinstance(metric_val, dict) and "per_patient" in metric_val:
                    for pid, val in metric_val["per_patient"].items():
                        if pid not in patients:
                            patients[pid] = {"id": pid, "outcome": "", "timepoints": []}
                        patients[pid]["timepoints"].append({
                            "dwi_type": dwi_type,
                            "label": "Baseline",
                            "metric": metric_name,
                            "value": val,
                        })

    # From MAT data: summary metrics
    for dwi_type, md in mat_data.items():
        if not isinstance(md, dict):
            continue
        summary = md.get("summary_metrics", {})
        if isinstance(summary, dict):
            for pid, metrics in summary.items():
                if pid not in patients:
                    patients[pid] = {"id": pid, "outcome": "", "timepoints": []}
                if isinstance(metrics, dict):
                    tp = {"dwi_type": dwi_type, "label": "Summary"}
                    tp.update({k: v for k, v in metrics.items()
                               if isinstance(v, (int, float)) and v == v})  # skip NaN
                    patients[pid]["timepoints"].append(tp)

        # Per-timepoint data
        longitudinal = md.get("longitudinal", {})
        if isinstance(longitudinal, dict):
            for pid, tp_list in longitudinal.items():
                if pid not in patients:
                    patients[pid] = {"id": pid, "outcome": "", "timepoints": []}
                if isinstance(tp_list, list):
                    for tp in tp_list:
                        if isinstance(tp, dict):
                            entry = {"dwi_type": dwi_type, "label": str(tp.get("timepoint", "?"))}
                            for k in ("adc_mean", "adc_median", "d_mean", "f_mean",
                                      "dstar_mean", "adc_volume", "d_volume"):
                                if k in tp and tp[k] is not None:
                                    try:
                                        v = float(tp[k])
                                        if v == v:  # skip NaN
                                            entry[k] = v
                                    except (ValueError, TypeError):
                                        pass
                            patients[pid]["timepoints"].append(entry)

    return sorted(patients.values(), key=lambda p: p["id"])


def _extract_core_methods(mat_data: dict) -> list[str]:
    """Extract the list of core method names from MAT data."""
    methods: set[str] = set()
    for md in mat_data.values():
        if not isinstance(md, dict):
            continue
        core = md.get("core_comparison", {})
        if isinstance(core, dict):
            for key in core:
                methods.add(str(key))
    return sorted(methods)
