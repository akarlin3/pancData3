#!/usr/bin/env python3
"""Generate an interactive HTML analysis report with client-side filtering.

Produces a standalone HTML file with embedded JavaScript that allows the
reader to:

- Filter by patient, DWI type, or core method via a sidebar
- Toggle between visualisations using tabbed panels
- Drill into individual scan results with expandable patient cards
- Compare metrics interactively with Chart.js bar/line charts
- Sort any table column by clicking its header
- Search across all content with a live search box

The report embeds all data as a JSON blob in a ``<script>`` tag and uses
Chart.js (loaded from CDN) for charting.  No server is required.

Usage:
    python generate_interactive_report.py [saved_files_path]
"""

from __future__ import annotations

import html as _html
import json
import sys
from datetime import datetime
from pathlib import Path

from tqdm import tqdm  # type: ignore

# Ensure analysis/ root is on sys.path so sibling packages are importable.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import (  # type: ignore
    DWI_TYPES,
    group_by_graph_name,
    load_graph_csv,
    resolve_folder,
    setup_utf8_stdout,
)
from parsers.parse_csv_results import parse_all_csvs  # type: ignore
from parsers.parse_log_metrics import parse_all_logs  # type: ignore
from report.interactive_constants import INTERACTIVE_CSS, INTERACTIVE_JS  # type: ignore

setup_utf8_stdout()


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


# ── Section builders ─────────────────────────────────────────────────────────

def _build_overview_section(
    log_data: dict | None, dwi_types: list[str], rows: list[dict],
    csv_data: dict | None, mat_data: dict, patients: list[dict],
) -> list[str]:
    """Build the overview section with stat cards."""
    h: list[str] = []
    h.append('<h2 id="overview">Overview</h2>')
    h.append('<div class="stat-grid">')

    # Patient count
    h.append(f'<div class="stat-card"><div class="label">Patients</div>'
             f'<div class="value">{len(patients)}</div>'
             f'<div class="sub">in dataset</div></div>')

    # DWI types
    h.append(f'<div class="stat-card"><div class="label">DWI Types</div>'
             f'<div class="value">{len(dwi_types)}</div>'
             f'<div class="sub">{", ".join(dwi_types)}</div></div>')

    # Graphs analysed
    if rows:
        h.append(f'<div class="stat-card"><div class="label">Graphs</div>'
                 f'<div class="value">{len(rows)}</div>'
                 f'<div class="sub">analysed</div></div>')

    # Best AUC from log data
    if log_data:
        best_auc = 0.0
        best_auc_dwi = ""
        for dwi, dd in log_data.items():
            if not isinstance(dd, dict):
                continue
            roc = dd.get("roc_auc", [])
            if isinstance(roc, list):
                for entry in roc:
                    if isinstance(entry, dict):
                        try:
                            auc_val = float(entry.get("auc", 0))
                            if auc_val > best_auc:
                                best_auc = auc_val
                                best_auc_dwi = dwi
                        except (ValueError, TypeError):
                            pass
        if best_auc > 0:
            h.append(f'<div class="stat-card"><div class="label">Best AUC</div>'
                     f'<div class="value">{best_auc:.3f}</div>'
                     f'<div class="sub">{_esc(best_auc_dwi)}</div></div>')

    # FDR significant count
    if csv_data and csv_data.get("fdr_global"):
        total_fdr = sum(len(v) for v in csv_data["fdr_global"].values() if isinstance(v, list))
        if total_fdr > 0:
            h.append(f'<div class="stat-card"><div class="label">FDR Significant</div>'
                     f'<div class="value">{total_fdr}</div>'
                     f'<div class="sub">surviving BH correction</div></div>')

    h.append('</div>')
    return h


def _build_patient_explorer(patients: list[dict], dwi_types: list[str]) -> list[str]:
    """Build the patient explorer section with expandable cards."""
    h: list[str] = []
    h.append('<h2 id="patients">Patient Explorer</h2>')
    h.append('<p class="meta">Click a patient card to expand and view per-timepoint metrics. '
             'Use sidebar filters to narrow the view.</p>')
    h.append('<div class="search-bar">'
             '<input type="text" id="search-input" placeholder="Search patients, metrics, values...">'
             '</div>')
    h.append('<div id="filter-count" class="filter-count"></div>')

    if not patients:
        h.append('<p class="meta">No patient data available.</p>')
        return h

    for p in patients:
        pid = _esc(p["id"])
        # Determine which DWI types this patient has data for
        pt_dwis = sorted(set(tp.get("dwi_type", "") for tp in p.get("timepoints", [])))
        dwi_badges = " ".join(_dwi_badge(d) for d in pt_dwis if d)

        h.append(f'<div class="patient-card" data-patient-id="{pid}">')
        h.append(f'<div class="patient-card-header">')
        h.append(f'<h3>{pid} {dwi_badges}</h3>')
        h.append('<span class="toggle-icon">\u25B6</span>')
        h.append('</div>')
        h.append('<div class="patient-card-body">')

        # Group timepoints by DWI type
        by_dwi: dict[str, list[dict]] = {}
        for tp in p.get("timepoints", []):
            dt = tp.get("dwi_type", "Unknown")
            by_dwi.setdefault(dt, []).append(tp)

        for dwi_type in dwi_types:
            tps = by_dwi.get(dwi_type, [])
            if not tps:
                continue
            h.append(f'<h4>{_dwi_badge(dwi_type)} — {len(tps)} timepoints</h4>')
            h.append('<div class="metric-grid">')
            for tp in tps:
                label = _esc(tp.get("label", "?"))
                for key in ("adc_mean", "adc_median", "d_mean", "f_mean",
                            "dstar_mean", "adc_volume", "d_volume", "metric", "value"):
                    val = tp.get(key)
                    if val is not None:
                        if isinstance(val, float):
                            display = f"{val:.6f}" if val < 0.01 else f"{val:.4f}"
                        else:
                            display = _esc(str(val))
                        h.append(f'<div class="metric-item">'
                                 f'<div class="m-label">{_esc(key)} ({label})</div>'
                                 f'<div class="m-value">{display}</div></div>')
            h.append('</div>')

        h.append('</div>')  # patient-card-body
        h.append('</div>')  # patient-card

    return h


def _build_visualizations_section(dwi_types: list[str]) -> list[str]:
    """Build the interactive visualizations section with tabs and charts."""
    h: list[str] = []
    h.append('<h2 id="visualizations">Interactive Visualizations</h2>')

    # Tab bar
    h.append('<div class="tab-bar" data-tab-group="viz">')
    h.append('<button class="tab-btn active" data-tab="longitudinal">Longitudinal</button>')
    h.append('<button class="tab-btn" data-tab="distribution">Distribution</button>')
    h.append('<button class="tab-btn" data-tab="comparison">DWI Comparison</button>')
    h.append('</div>')

    # Longitudinal tab
    h.append('<div class="tab-panel active" data-tab-group="viz" data-tab="longitudinal">')
    h.append('<p class="meta">Mean ADC values over time, grouped by DWI type. '
             'Toggle DWI types in the sidebar to compare.</p>')
    h.append('<div class="chart-container"><canvas id="chart-longitudinal"></canvas></div>')
    h.append('</div>')

    # Distribution tab
    h.append('<div class="tab-panel" data-tab-group="viz" data-tab="distribution">')
    h.append('<p class="meta">Distribution of mean ADC values across the filtered cohort.</p>')
    h.append('<div class="chart-container"><canvas id="chart-distribution"></canvas></div>')
    h.append('</div>')

    # Comparison tab
    h.append('<div class="tab-panel" data-tab-group="viz" data-tab="comparison">')
    h.append('<p class="meta">Side-by-side comparison of summary metrics across DWI types.</p>')
    h.append('<div class="chart-container"><canvas id="chart-dwi-compare"></canvas></div>')
    h.append('</div>')

    return h


def _build_significance_table(
    csv_data: dict | None, log_data: dict | None, dwi_types: list[str],
) -> list[str]:
    """Build a filterable significance table from CSV and log data."""
    h: list[str] = []
    h.append('<h2 id="significance">Statistical Significance</h2>')

    # FDR global table
    if csv_data and csv_data.get("fdr_global"):
        h.append('<h3>FDR Global Correction Results</h3>')
        h.append('<p class="meta">Metrics surviving Benjamini-Hochberg FDR correction. '
                 'Click column headers to sort.</p>')

        for dwi_type in DWI_TYPES:
            fdr = csv_data["fdr_global"].get(dwi_type, [])
            if not fdr:
                continue
            h.append(f'<h4>{_dwi_badge(dwi_type)} \u2014 {len(fdr)} tests</h4>')
            headers = list(fdr[0].keys())[:8]
            h.append('<table class="sortable"><thead><tr>')
            for hdr in headers:
                h.append(f'<th>{_esc(hdr)}</th>')
            h.append('</tr></thead><tbody>')
            for row in fdr:
                dwi = row.get("DWI_Type", dwi_type)
                h.append(f'<tr data-dwi="{_esc(dwi)}">')
                for hdr in headers:
                    val = str(row.get(hdr, ""))
                    try:
                        fval = float(val)
                        if "p" in hdr.lower() and 0 < fval < 0.05:
                            cls = _sig_class(fval)
                            h.append(f'<td class="{cls}">{_esc(val[:40])}</td>')
                            continue
                    except (ValueError, TypeError):
                        pass
                    h.append(f'<td>{_esc(val[:40])}</td>')
                h.append('</tr>')
            h.append('</tbody></table>')

    # Hazard ratios from log data
    if log_data:
        any_hr = False
        for dwi in dwi_types:
            dd = log_data.get(dwi, {})
            if not isinstance(dd, dict):
                continue
            hrs = dd.get("hazard_ratios", [])
            if isinstance(hrs, list) and hrs:
                if not any_hr:
                    h.append('<h3>Hazard Ratios</h3>')
                    any_hr = True
                h.append(f'<h4>{_dwi_badge(dwi)}</h4>')
                h.append('<table class="sortable"><thead><tr>'
                         '<th>Feature</th><th>HR</th><th>CI Low</th><th>CI High</th>'
                         '<th>p-value</th></tr></thead><tbody>')
                for hr in hrs:
                    if not isinstance(hr, dict):
                        continue
                    p = hr.get("p", 1.0)
                    cls = _sig_class(p) if isinstance(p, (int, float)) else ""
                    h.append(f'<tr data-dwi="{_esc(dwi)}">'
                             f'<td>{_esc(hr.get("feature", ""))}</td>'
                             f'<td>{hr.get("hr", "")}</td>'
                             f'<td>{hr.get("ci_lo", "")}</td>'
                             f'<td>{hr.get("ci_hi", "")}</td>'
                             f'<td class="{cls}">{p}</td></tr>')
                h.append('</tbody></table>')

    if not csv_data and not log_data:
        h.append('<p class="meta">No significance data available.</p>')

    return h


def _build_graph_explorer(rows: list[dict], groups: dict) -> list[str]:
    """Build the graph analysis explorer with filtering."""
    h: list[str] = []
    h.append('<h2 id="graphs">Graph Analysis</h2>')

    if not rows:
        h.append('<p class="meta">No graph analysis data available.</p>')
        return h

    h.append(f'<p class="meta">{len(rows)} graphs analysed across '
             f'{len(groups)} unique graph types.</p>')

    # Tab bar for graph view modes
    h.append('<div class="tab-bar" data-tab-group="graphs">')
    h.append('<button class="tab-btn active" data-tab="table">Table View</button>')
    h.append('<button class="tab-btn" data-tab="cross-dwi">Cross-DWI Comparison</button>')
    h.append('</div>')

    # Table view
    h.append('<div class="tab-panel active" data-tab-group="graphs" data-tab="table">')
    h.append('<table class="sortable"><thead><tr>'
             '<th>Graph</th><th>DWI Type</th><th>Axes</th>'
             '<th>Trends</th><th>Summary</th>'
             '</tr></thead><tbody>')
    for row in rows:
        dwi = row.get("dwi_type", "")
        graph_name = row.get("graph_name", row.get("file", ""))
        axes = row.get("axes_text", row.get("x_axis", ""))
        trends = row.get("trends_text", "")
        summary = row.get("full_summary", row.get("summary", ""))[:120]
        h.append(f'<tr data-dwi="{_esc(dwi)}">')
        h.append(f'<td><code>{_esc(graph_name)}</code></td>')
        h.append(f'<td>{_dwi_badge(dwi) if dwi else ""}</td>')
        h.append(f'<td class="meta">{_esc(axes)}</td>')
        h.append(f'<td>{_esc(trends)}</td>')
        h.append(f'<td class="meta">{_esc(summary)}</td>')
        h.append('</tr>')
    h.append('</tbody></table>')
    h.append('</div>')

    # Cross-DWI comparison view
    h.append('<div class="tab-panel" data-tab-group="graphs" data-tab="cross-dwi">')
    if groups:
        for gname, dwi_map in sorted(groups.items()):
            if len(dwi_map) < 2:
                continue  # Only show graphs present in multiple DWI types
            h.append(f'<h4>{_esc(gname)}</h4>')
            h.append('<div class="compare-panel">')
            for dwi_type, grow in sorted(dwi_map.items()):
                trends = grow.get("trends_text", "")
                summary = grow.get("full_summary", grow.get("summary", ""))[:200]
                h.append(f'<div class="compare-col">')
                h.append(f'<h4>{_dwi_badge(dwi_type)}</h4>')
                if trends:
                    h.append(f'<p><strong>Trends:</strong> {_esc(trends)}</p>')
                if summary:
                    h.append(f'<p class="meta">{_esc(summary)}</p>')
                h.append('</div>')
            h.append('</div>')
    else:
        h.append('<p class="meta">No cross-DWI comparison data available.</p>')
    h.append('</div>')

    return h


def _build_core_comparison(mat_data: dict) -> list[str]:
    """Build the core method comparison section from MAT data."""
    h: list[str] = []
    h.append('<h2 id="core-methods">Core Method Comparison</h2>')

    any_data = False
    for dwi_type, md in mat_data.items():
        if not isinstance(md, dict):
            continue
        core = md.get("core_comparison", {})
        if not isinstance(core, dict) or not core:
            continue
        any_data = True
        h.append(f'<h3>{_dwi_badge(dwi_type)}</h3>')

        # Dice matrix
        dice = core.get("dice_matrix")
        methods = core.get("method_names", [])
        if isinstance(dice, list) and isinstance(methods, list) and methods:
            h.append('<h4>Dice Similarity</h4>')
            h.append('<table class="sortable"><thead><tr><th></th>')
            for m in methods:
                h.append(f'<th>{_esc(m)}</th>')
            h.append('</tr></thead><tbody>')
            for i, row_vals in enumerate(dice):
                if not isinstance(row_vals, list):
                    continue
                h.append(f'<tr data-dwi="{_esc(dwi_type)}"><td><strong>{_esc(methods[i] if i < len(methods) else "?")}</strong></td>')
                for val in row_vals:
                    try:
                        fv = float(val)
                        h.append(f'<td>{fv:.3f}</td>')
                    except (ValueError, TypeError):
                        h.append(f'<td>{_esc(str(val))}</td>')
                h.append('</tr>')
            h.append('</tbody></table>')

        # Volume comparison
        volumes = core.get("volumes")
        if isinstance(volumes, dict) and volumes:
            h.append('<h4>Core Volumes</h4>')
            h.append('<table class="sortable"><thead><tr><th>Method</th><th>Volume</th></tr></thead><tbody>')
            for method, vol in volumes.items():
                h.append(f'<tr data-dwi="{_esc(dwi_type)}"><td>{_esc(method)}</td>'
                         f'<td>{vol}</td></tr>')
            h.append('</tbody></table>')

    if not any_data:
        h.append('<p class="meta">No core method comparison data available.</p>')

    return h


def _build_dosimetry_section(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build the dosimetry metrics section."""
    h: list[str] = []
    h.append('<h2 id="dosimetry">Dosimetry</h2>')

    any_data = False
    for dwi_type in dwi_types:
        md = mat_data.get(dwi_type, {})
        if not isinstance(md, dict):
            continue
        dosi = md.get("dosimetry", {})
        if not isinstance(dosi, dict) or not dosi:
            continue
        any_data = True
        h.append(f'<h3>{_dwi_badge(dwi_type)}</h3>')
        h.append('<div class="stat-grid">')
        for key, val in dosi.items():
            if isinstance(val, (int, float)):
                h.append(f'<div class="stat-card"><div class="label">{_esc(key)}</div>'
                         f'<div class="value">{val:.2f}</div></div>')
            elif isinstance(val, dict):
                for sub_key, sub_val in val.items():
                    if isinstance(sub_val, (int, float)):
                        h.append(f'<div class="stat-card"><div class="label">{_esc(key)}: {_esc(sub_key)}</div>'
                                 f'<div class="value">{sub_val:.2f}</div></div>')
        h.append('</div>')

    if not any_data:
        h.append('<p class="meta">No dosimetry data available.</p>')
    return h


# ── Main report assembly ─────────────────────────────────────────────────────

def generate_interactive_report(folder: Path) -> str:
    """Build the full interactive HTML report from all available data sources.

    Parameters
    ----------
    folder : Path
        Path to the ``saved_files_*`` output folder.

    Returns
    -------
    str
        Complete HTML document as a string.
    """
    # Extract timestamp from folder name.
    raw_ts = folder.name.replace("saved_files_", "")
    try:
        _dt = datetime.strptime(raw_ts, "%Y%m%d_%H%M%S")
        timestamp = _dt.strftime(f"%B {_dt.day}, %Y at %H:%M:%S")
    except ValueError:
        timestamp = raw_ts
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # ── Load data sources ──
    report_bar = tqdm(total=10, desc="Building interactive report", unit="step",
                      bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}] {postfix}")

    report_bar.set_postfix_str("loading CSVs", refresh=True)
    try:
        csv_data = parse_all_csvs(folder)
    except Exception:
        csv_data = None
    report_bar.update(1)

    report_bar.set_postfix_str("loading logs", refresh=True)
    try:
        log_data = parse_all_logs(folder)
    except Exception:
        log_data = None
    report_bar.update(1)

    report_bar.set_postfix_str("loading MAT data", refresh=True)
    mat_data: dict = {}
    for dt in DWI_TYPES:
        mat_path = folder / f"parsed_mat_metrics_{dt}.json"
        if mat_path.exists():
            try:
                with open(mat_path, "r", encoding="utf-8") as f:
                    mat_data[dt] = json.load(f)
            except Exception:
                pass
    report_bar.update(1)

    report_bar.set_postfix_str("loading graphs", refresh=True)
    rows = load_graph_csv(folder)
    groups = group_by_graph_name(rows) if rows else {}
    report_bar.update(1)

    # Detect DWI types present
    dwi_types_present: list[str] = []
    if log_data:
        for d in DWI_TYPES:
            if d in log_data:
                dwi_types_present.append(d)
    if not dwi_types_present:
        dwi_types_present = [d for d in DWI_TYPES if (folder / d).is_dir()]

    # Extract structured data
    report_bar.set_postfix_str("extracting patients", refresh=True)
    patients = _extract_patients(log_data, mat_data)
    core_methods = _extract_core_methods(mat_data)
    report_bar.update(1)

    # Prepare JSON data blob for JavaScript
    report_data = {
        "patients": patients,
        "dwi_types": dwi_types_present,
        "core_methods": core_methods,
        "graph_rows": [
            {k: v for k, v in r.items()
             if isinstance(v, (str, int, float, bool, type(None)))}
            for r in (rows or [])
        ],
    }

    # ── Assemble HTML ──
    h: list[str] = []
    h.append("<!DOCTYPE html>")
    h.append('<html lang="en">')
    h.append("<head>")
    h.append('<meta charset="utf-8">')
    h.append('<meta name="viewport" content="width=device-width, initial-scale=1">')
    h.append(f"<title>Interactive Report \u2014 {_esc(timestamp)}</title>")
    h.append(f"<style>{INTERACTIVE_CSS}</style>")
    # Chart.js from CDN
    h.append('<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.7/dist/chart.umd.min.js"></script>')
    # Inject data as JSON
    h.append("<script>")
    h.append("window.REPORT_DATA = ")
    h.append(json.dumps(report_data, default=str, ensure_ascii=False))
    h.append(";")
    h.append("</script>")
    h.append(f"<script>{INTERACTIVE_JS}</script>")
    h.append("</head>")
    h.append("<body>")

    # ── Sidebar ──
    h.append('<aside class="sidebar">')
    h.append("<h2>Filters</h2>")

    # DWI type filter
    h.append("<h3>DWI Type</h3>")
    h.append('<div id="dwi-filter-checks"></div>')

    # Patient filter
    h.append("<h3>Patient</h3>")
    h.append('<select id="patient-select" multiple size="5">')
    h.append('<option value="">All patients</option>')
    h.append("</select>")

    # Core method filter
    h.append('<div id="core-method-section">')
    h.append("<h3>Core Method</h3>")
    h.append('<div id="core-method-checks"></div>')
    h.append("</div>")

    # Navigation links
    h.append('<div class="nav-links">')
    h.append("<h3>Sections</h3>")
    nav_items = [
        ("overview", "Overview"),
        ("visualizations", "Visualizations"),
        ("patients", "Patient Explorer"),
        ("significance", "Significance"),
        ("graphs", "Graph Analysis"),
        ("core-methods", "Core Methods"),
        ("dosimetry", "Dosimetry"),
    ]
    for anchor, label in nav_items:
        h.append(f'<a href="#{anchor}">{label}</a>')
    h.append("</div>")

    h.append("</aside>")

    # ── Main content ──
    h.append('<main class="main-content">')
    h.append(f"<h1>Interactive Analysis Report \u2014 {_esc(timestamp)}</h1>")
    h.append(f'<div class="meta">Generated: {now} | DWI types: '
             f'{", ".join(dwi_types_present) or "None detected"}</div>')

    # Overview section
    report_bar.set_postfix_str("overview", refresh=True)
    h.extend(_build_overview_section(log_data, dwi_types_present, rows or [], csv_data, mat_data, patients))
    report_bar.update(1)

    # Visualizations section
    report_bar.set_postfix_str("visualizations", refresh=True)
    h.extend(_build_visualizations_section(dwi_types_present))
    report_bar.update(1)

    # Patient explorer
    report_bar.set_postfix_str("patient explorer", refresh=True)
    h.extend(_build_patient_explorer(patients, dwi_types_present))
    report_bar.update(1)

    # Significance tables
    report_bar.set_postfix_str("significance", refresh=True)
    h.extend(_build_significance_table(csv_data, log_data, dwi_types_present))

    # Graph analysis
    h.extend(_build_graph_explorer(rows or [], groups))

    # Core method comparison
    h.extend(_build_core_comparison(mat_data))

    # Dosimetry
    h.extend(_build_dosimetry_section(mat_data, dwi_types_present))
    report_bar.update(1)

    # Footer
    h.append(f"<footer>Interactive report generated by pancData3 on {now}</footer>")
    h.append("</main>")
    h.append("</body>")
    h.append("</html>")

    report_bar.set_postfix_str("complete", refresh=True)
    report_bar.close()

    return "\n".join(h)


def main():
    """CLI entry point: generate the interactive HTML report."""
    folder = resolve_folder(sys.argv)
    html_report = generate_interactive_report(folder)
    out_path = folder / "interactive_report.html"
    out_path.write_text(html_report, encoding="utf-8")
    print(f"Interactive report written to: {out_path}")
    print(f"  Length: {len(html_report)} characters")


if __name__ == "__main__":
    main()
