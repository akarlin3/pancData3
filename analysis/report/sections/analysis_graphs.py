"""Report sections: graph analysis (overview, issues, stats by type)."""

from __future__ import annotations

import json

from shared import (  # type: ignore
    extract_correlations,
    extract_pvalues,
    parse_dwi_info,
    safe_text,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _stat_card,
)

from shared import DWI_TYPES  # type: ignore


def _section_graph_overview(rows) -> list[str]:
    """Build the Graph Analysis Overview section.

    Displays two side-by-side summary tables:
    - Count of graphs by type (line, scatter, box, etc.)
    - Count of graphs by DWI type (Standard, dnCNN, IVIMnet, Root)

    Parameters
    ----------
    rows : list[dict]
        Vision CSV rows (may be empty, in which case the section is skipped).

    Returns
    -------
    list[str]
        HTML chunks for the graph overview section.
    """
    # ── 2. Graph Analysis Overview ──
    h: list[str] = []
    if rows:
        h.append(_h2("Graph Analysis Overview", "graph-overview"))

        type_counts: dict[str, int] = {}
        for r in rows:
            gt = r.get("graph_type", "unknown")
            type_counts[gt] = type_counts.get(gt, 0) + 1

        dwi_counts: dict[str, int] = {}
        for r in rows:
            dwi_type, _ = parse_dwi_info(r["file_path"])
            dwi_counts[dwi_type] = dwi_counts.get(dwi_type, 0) + 1

        h.append('<div style="display:grid;grid-template-columns:1fr 1fr;gap:1rem">')

        h.append("<div>")
        h.append("<h3>By Graph Type</h3>")
        h.append("<table><thead><tr><th>Type</th><th>Count</th></tr></thead><tbody>")
        for gt, cnt in sorted(type_counts.items(), key=lambda x: -x[1]):
            h.append(f"<tr><td>{_esc(gt)}</td><td>{cnt}</td></tr>")
        h.append("</tbody></table>")
        h.append("</div>")

        h.append("<div>")
        h.append("<h3>By DWI Type</h3>")
        h.append("<table><thead><tr><th>DWI Type</th><th>Graphs</th></tr></thead><tbody>")
        for dt in DWI_TYPES + ["Root"]:
            if dt in dwi_counts:
                h.append(f"<tr><td>{_dwi_badge(dt)}</td><td>{dwi_counts[dt]}</td></tr>")  # type: ignore
        h.append("</tbody></table>")
        h.append("</div>")

        h.append("</div>")

        # Per-graph-type statistical signal density
        h.append("<h3>Statistical Signal by Graph Type</h3>")
        h.append('<p class="meta">Shows which graph types contain the most '
                 'extractable statistical information (p-values, correlations).</p>')
        signal_data: list[tuple[str, int, int, int]] = []
        for gt, cnt in sorted(type_counts.items(), key=lambda x: -x[1]):
            gt_pvals = 0
            gt_corrs = 0
            for r in rows:
                if r.get("graph_type", "unknown") != gt:
                    continue
                all_text = safe_text(r, "summary", "trends_json", "inflection_points_json")
                gt_pvals += sum(1 for p, _ in extract_pvalues(all_text) if p < 0.05)
                gt_corrs += sum(1 for rv, _ in extract_correlations(all_text) if abs(rv) >= 0.3)
            if gt_pvals > 0 or gt_corrs > 0:
                signal_data.append((gt, cnt, gt_pvals, gt_corrs))

        if signal_data:
            h.append("<table><thead><tr>"
                     "<th>Graph Type</th><th>Count</th>"
                     "<th>Sig. p-values</th><th>Strong Correlations</th>"
                     "<th>Signal Density</th>"
                     "</tr></thead><tbody>")
            for gt, cnt, n_p, n_c in sorted(signal_data, key=lambda x: -(x[2] + x[3])):
                density = (n_p + n_c) / cnt if cnt > 0 else 0
                density_cls = "agree" if density >= 1.0 else ("" if density >= 0.5 else "")
                d_attr = f' class="{density_cls}"' if density_cls else ""
                h.append(
                    f"<tr><td><strong>{_esc(gt)}</strong></td>"
                    f"<td>{cnt}</td>"
                    f"<td>{n_p}</td>"
                    f"<td>{n_c}</td>"
                    f"<td{d_attr}>{density:.2f}/graph</td></tr>"
                )
            h.append("</tbody></table>")

    return h


def _section_graph_issues(rows) -> list[str]:
    """Build the Graph Issues section.

    Lists graphs that have quality issues detected by the vision model,
    grouped by DWI type.  Each issue is shown with the graph name and a
    description of the problem.

    Parameters
    ----------
    rows : list[dict]
        Vision CSV rows (may be empty, in which case the section is skipped).

    Returns
    -------
    list[str]
        HTML chunks for the graph issues section.
    """
    h: list[str] = []
    if not rows:
        return h

    # Collect rows that have issues or are error/unknown type.
    issue_rows = []
    for r in rows:
        issues_str = r.get("issues_json", "[]") or "[]"
        try:
            issues_list = json.loads(issues_str)
        except (json.JSONDecodeError, TypeError):
            issues_list = []
        graph_type = r.get("graph_type", "")
        # Also flag error/unknown graph types as issues.
        type_issues = []
        if graph_type == "error":
            type_issues.append("Graph analysis failed (API or processing error)")
        elif graph_type == "unknown":
            type_issues.append("Graph type could not be determined")
        all_issues = type_issues + (issues_list if isinstance(issues_list, list) else [])
        if all_issues:
            issue_rows.append((r, all_issues))

    if not issue_rows:
        return h

    h.append(_h2("Graph Issues", "graph-issues"))
    h.append(f"<p>{len(issue_rows)} of {len(rows)} graphs have detected quality issues.</p>")

    # Severity classification (change 7)
    _HIGH_KEYWORDS = {"cutoff", "overlap", "unreadable"}

    def _severity(r: dict, issues: list) -> str:
        gt = r.get("graph_type", "")
        if gt in ("error", "unknown"):
            return "Critical"
        combined = " ".join(str(iss).lower() for iss in issues)
        if any(kw in combined for kw in _HIGH_KEYWORDS):
            return "High"
        return "Low"

    n_critical = sum(1 for r, iss in issue_rows if _severity(r, iss) == "Critical")
    n_high = sum(1 for r, iss in issue_rows if _severity(r, iss) == "High")
    n_low = sum(1 for r, iss in issue_rows if _severity(r, iss) == "Low")

    # Severity summary stat cards
    sev_cards = []
    if n_critical:
        sev_cards.append(_stat_card("Critical", str(n_critical), "must re-run vision analysis"))
    if n_high:
        sev_cards.append(_stat_card("High", str(n_high), "readability/overlap issues"))
    if n_low:
        sev_cards.append(_stat_card("Low", str(n_low), "minor quality issues"))
    if sev_cards:
        h.append('<div class="stat-grid">')
        h.extend(sev_cards)
        h.append("</div>")

    # Legacy summary stat cards (kept for continuity)
    n_error = sum(1 for r, _ in issue_rows if r.get("graph_type") == "error")
    n_unknown = sum(1 for r, _ in issue_rows if r.get("graph_type") == "unknown")
    n_quality = len(issue_rows) - n_error - n_unknown
    cards = []
    if n_error:
        cards.append(_stat_card("Analysis Errors", str(n_error), "failed to analyse"))
    if n_unknown:
        cards.append(_stat_card("Unknown Type", str(n_unknown), "type not determined"))
    if n_quality:
        cards.append(_stat_card("Quality Issues", str(n_quality), "visual problems detected"))
    if cards:
        h.append('<div class="stat-grid">')
        h.extend(cards)
        h.append("</div>")

    # Render issues grouped by severity
    for severity_label in ("Critical", "High", "Low"):
        group = [(r, iss) for r, iss in issue_rows if _severity(r, iss) == severity_label]
        if not group:
            continue
        h.append(f"<h3>Severity: {_esc(severity_label)}</h3>")
        h.append("<table><thead><tr>"
                 "<th>#</th><th>DWI</th><th>Graph</th><th>Type</th><th>Issues</th>"
                 "</tr></thead><tbody>")
        for i, (r, issues) in enumerate(group, 1):
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            graph_type = r.get("graph_type", "")
            issues_html = "<ul>" + "".join(
                f"<li>{_esc(iss)}</li>" for iss in issues
            ) + "</ul>"
            h.append(
                f"<tr>"
                f"<td>{i}</td>"
                f"<td>{_dwi_badge(dwi_type)}</td>"
                f"<td>{_esc(base_name)}</td>"
                f"<td>{_esc(graph_type)}</td>"
                f"<td>{issues_html}</td>"
                f"</tr>"
            )
        h.append("</tbody></table>")
    return h


def _section_stats_by_graph_type(rows) -> list[str]:
    """Build the Statistics by Graph Type section.

    Aggregates significant p-values, correlations, and trend directions
    by graph type (box, line, scatter, heatmap, etc.), providing a
    bird's-eye view of which graph types carry the most statistical signal.

    Parameters
    ----------
    rows : list[dict]
        Vision CSV rows (may be empty).

    Returns
    -------
    list[str]
        HTML chunks for the statistics by graph type section.
    """
    h: list[str] = []
    if not rows:
        return h

    h.append(_h2("Statistics by Graph Type", "stats-by-type"))

    # Aggregate per graph type
    type_stats: dict[str, dict] = {}
    for r in rows:
        gt = r.get("graph_type", "unknown")
        if gt not in type_stats:
            type_stats[gt] = {
                "count": 0, "sig_p": 0, "nonsig_p": 0,
                "correlations": 0, "trends": {"increasing": 0, "decreasing": 0, "stable": 0, "other": 0},
            }
        ts = type_stats[gt]
        ts["count"] += 1

        # P-values
        all_text = safe_text(r, "summary", "trends_json", "inflection_points_json")
        for pval, _ in extract_pvalues(all_text):
            if pval < 0.05:
                ts["sig_p"] += 1
            else:
                ts["nonsig_p"] += 1

        # Correlations
        for rval, _ in extract_correlations(all_text):
            if abs(rval) >= 0.3:
                ts["correlations"] += 1

        # Trends
        trends_str = r.get("trends_json", "[]") or "[]"
        try:
            trends = json.loads(str(trends_str))
        except (json.JSONDecodeError, TypeError):
            trends = []
        if isinstance(trends, list):
            for t in trends:
                if not isinstance(t, dict):
                    continue
                d = str(t.get("direction", "")).lower()
                if "increas" in d or "up" in d or "higher" in d or "rising" in d:
                    ts["trends"]["increasing"] += 1
                elif "decreas" in d or "down" in d or "lower" in d or "falling" in d or "drop" in d:
                    ts["trends"]["decreasing"] += 1
                elif "flat" in d or "stable" in d or "constant" in d:
                    ts["trends"]["stable"] += 1
                else:
                    ts["trends"]["other"] += 1

    # Summary table
    h.append('<table class="table-compact"><thead><tr>'
             "<th>Graph Type</th><th>Count</th><th>Sig. p</th>"
             "<th>Non-sig. p</th><th>Corr.</th>"
             "<th>\u2191</th><th>\u2193</th><th>\u2192</th><th>Other</th>"
             "</tr></thead><tbody>")
    for gt in sorted(type_stats.keys(), key=lambda k: -type_stats[k]["count"]):
        ts = type_stats[gt]
        tr = ts["trends"]
        sig_cls = ' class="sig-1"' if ts["sig_p"] > 0 else ""
        h.append(f"<tr><td><strong>{_esc(gt)}</strong></td>"
                 f"<td>{ts['count']}</td>"
                 f"<td{sig_cls}>{ts['sig_p']}</td>"
                 f"<td>{ts['nonsig_p']}</td>"
                 f"<td>{ts['correlations']}</td>"
                 f"<td>{tr['increasing']}</td>"
                 f"<td>{tr['decreasing']}</td>"
                 f"<td>{tr['stable']}</td>"
                 f"<td>{tr['other']}</td></tr>")
    h.append("</tbody></table>")

    # Highlight graph types with the most signal
    most_sig = max(type_stats.items(), key=lambda x: x[1]["sig_p"], default=None)
    if most_sig and most_sig[1]["sig_p"] > 0:
        h.append(f'<div class="info-box"><strong>{_esc(most_sig[0])}</strong> graphs carry the most '
                 f'significant findings ({most_sig[1]["sig_p"]} p &lt; 0.05).</div>')

    return h
