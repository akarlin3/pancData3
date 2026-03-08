#!/usr/bin/env python3
"""Generate a comprehensive HTML analysis report.

Combines vision-based graph analysis (if available), direct CSV parsing, and
log file parsing into a single structured HTML report.

Usage:
    python generate_report.py [saved_files_path]
"""

from __future__ import annotations

import html
import json
import sys
from datetime import datetime
from pathlib import Path

import markdown

from shared import (
    DWI_TYPES,
    extract_correlations,
    extract_pvalues,
    group_by_graph_name,
    load_graph_csv,
    parse_dwi_info,
    resolve_folder,
    setup_utf8_stdout,
)
from parse_csv_results import parse_all_csvs
from parse_log_metrics import parse_all_logs

setup_utf8_stdout()

HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
  :root {{
    --bg: #ffffff;
    --fg: #1a1a2e;
    --accent: #0f3460;
    --accent-light: #e8eef6;
    --border: #d0d7de;
    --table-stripe: #f6f8fa;
    --sig-strong: #d32f2f;
    --sig-moderate: #e65100;
    --sig-mild: #f9a825;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                 "Helvetica Neue", Arial, sans-serif;
    line-height: 1.6;
    color: var(--fg);
    background: var(--bg);
    max-width: 1100px;
    margin: 0 auto;
    padding: 2rem 1.5rem;
  }}
  h1 {{
    font-size: 1.8rem;
    color: var(--accent);
    border-bottom: 3px solid var(--accent);
    padding-bottom: 0.5rem;
    margin-bottom: 1rem;
  }}
  h2 {{
    font-size: 1.4rem;
    color: var(--accent);
    border-bottom: 1px solid var(--border);
    padding-bottom: 0.3rem;
    margin-top: 2rem;
    margin-bottom: 0.8rem;
  }}
  h3 {{
    font-size: 1.15rem;
    margin-top: 1.5rem;
    margin-bottom: 0.5rem;
  }}
  p {{ margin-bottom: 0.6rem; }}
  code {{
    background: var(--accent-light);
    padding: 0.15em 0.4em;
    border-radius: 3px;
    font-size: 0.9em;
  }}
  table {{
    border-collapse: collapse;
    width: 100%;
    margin: 0.8rem 0 1.2rem;
    font-size: 0.92rem;
  }}
  th, td {{
    border: 1px solid var(--border);
    padding: 0.45rem 0.7rem;
    text-align: left;
  }}
  th {{
    background: var(--accent);
    color: #fff;
    font-weight: 600;
  }}
  tr:nth-child(even) {{ background: var(--table-stripe); }}
  tr:hover {{ background: var(--accent-light); }}
  ul, ol {{ margin: 0.5rem 0 0.5rem 1.5rem; }}
  li {{ margin-bottom: 0.25rem; }}
  strong {{ color: var(--accent); }}
  hr {{
    border: none;
    border-top: 1px solid var(--border);
    margin: 2rem 0;
  }}
  em {{ color: #555; }}
  .footer {{
    margin-top: 2rem;
    padding-top: 1rem;
    border-top: 1px solid var(--border);
    font-size: 0.85rem;
    color: #666;
  }}
</style>
</head>
<body>
{body}
</body>
</html>
"""


def _esc(text: str) -> str:
    """HTML-escape a string."""
    return html.escape(str(text))


def _sig_tag(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def _sig_class(p: float) -> str:
    """Return a CSS class name for the significance level."""
    if p < 0.001:
        return "sig-3"
    if p < 0.01:
        return "sig-2"
    if p < 0.05:
        return "sig-1"
    return ""


CSS = """\
:root {
    --bg: #ffffff; --fg: #1a1a2e; --muted: #64748b;
    --accent: #2563eb; --accent-light: #dbeafe;
    --border: #e2e8f0; --row-alt: #f8fafc;
    --green: #16a34a; --red: #dc2626; --amber: #d97706;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    color: var(--fg); background: var(--bg);
    line-height: 1.6; max-width: 1100px; margin: 0 auto; padding: 2rem 1.5rem;
}
h1 { font-size: 1.8rem; border-bottom: 3px solid var(--accent); padding-bottom: 0.5rem; margin-bottom: 1.5rem; }
h2 { font-size: 1.35rem; color: var(--accent); margin: 2rem 0 0.75rem; border-bottom: 1px solid var(--border); padding-bottom: 0.3rem; }
h3 { font-size: 1.1rem; margin: 1.25rem 0 0.5rem; }
p, ul { margin-bottom: 0.75rem; }
ul { padding-left: 1.5rem; }
.meta { color: var(--muted); font-size: 0.9rem; margin-bottom: 1.5rem; }
.meta span { display: inline-block; margin-right: 1.5rem; }
code { background: var(--row-alt); border: 1px solid var(--border); border-radius: 3px; padding: 0.1em 0.35em; font-size: 0.88em; }
table { width: 100%; border-collapse: collapse; margin: 0.75rem 0 1.25rem; font-size: 0.92rem; }
th { background: var(--accent); color: #fff; text-align: left; padding: 0.5rem 0.75rem; font-weight: 600; }
td { padding: 0.4rem 0.75rem; border-bottom: 1px solid var(--border); }
tr:nth-child(even) td { background: var(--row-alt); }
.sig-1 { color: var(--amber); font-weight: 600; }
.sig-2 { color: var(--red); font-weight: 600; }
.sig-3 { color: var(--red); font-weight: 700; }
.agree { color: var(--green); font-weight: 600; }
.differ { color: var(--red); font-weight: 700; }
.badge { display: inline-block; padding: 0.15em 0.5em; border-radius: 4px; font-size: 0.82rem; font-weight: 600; }
.badge-standard { background: #dbeafe; color: #1e40af; }
.badge-dncnn { background: #dcfce7; color: #166534; }
.badge-ivimnet { background: #fef3c7; color: #92400e; }
.badge-root { background: #f1f5f9; color: #475569; }
.summary-box { background: var(--accent-light); border-left: 4px solid var(--accent); padding: 1rem 1.25rem; border-radius: 0 6px 6px 0; margin: 1rem 0; }
footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid var(--border); color: var(--muted); font-size: 0.85rem; }
"""


def _dwi_badge(dwi_type: str) -> str:
    cls = {
        "Standard": "badge-standard", "dnCNN": "badge-dncnn",
        "IVIMnet": "badge-ivimnet",
    }.get(dwi_type, "badge-root")
    return f'<span class="badge {cls}">{_esc(dwi_type)}</span>'


def generate_report(folder: Path) -> str:
    """Build the full HTML report string."""
    timestamp = folder.name.replace("saved_files_", "")
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    try:
        csv_data = parse_all_csvs(folder)
    except Exception:
        csv_data = None

    try:
        log_data = parse_all_logs(folder)
    except Exception:
        log_data = None

    # Load MAT metrics (generated by parse_mat_metrics.py)
    mat_data = {}
    for dt in DWI_TYPES:
        mat_path = folder / f"parsed_mat_metrics_{dt}.json"
        if mat_path.exists():
            try:
                with open(mat_path, "r", encoding="utf-8") as f:
                    mat_data[dt] = json.load(f)
            except Exception as e:
                print(f"Error loading {mat_path.name}: {e}")

    # Load vision CSV if available
    rows = load_graph_csv(folder)
    groups = group_by_graph_name(rows) if rows else {}

    # Detect which DWI types are present
    dwi_types_present = []
    if log_data is not None:
        for d in DWI_TYPES:
            if d in log_data:
                dwi_types_present.append(d)
    if not dwi_types_present:
        dwi_types_present = [d for d in DWI_TYPES if (folder / d).is_dir()]

    h = []  # html chunks
    h.append("<!DOCTYPE html>")
    h.append('<html lang="en">')
    h.append("<head>")
    h.append('<meta charset="utf-8">')
    h.append('<meta name="viewport" content="width=device-width, initial-scale=1">')
    h.append(f"<title>Analysis Report \u2014 {_esc(timestamp)}</title>")
    h.append(f"<style>{CSS}</style>")
    h.append("</head>")
    h.append("<body>")

    # ── Header ──
    h.append(f"<h1>Analysis Report \u2014 {_esc(timestamp)}</h1>")
    h.append('<div class="meta">')
    h.append(f"<span><strong>Generated:</strong> {now}</span>")
    h.append(f"<span><strong>DWI types:</strong> {', '.join(dwi_types_present) or 'None detected'}</span>")
    if rows:
        h.append(f"<span><strong>Graphs analysed:</strong> {len(rows)}</span>")
    h.append("</div>")

    # ── 1. Executive Summary ──
    h.append("<h2>Executive Summary</h2>")
    h.append(f'<div class="summary-box">')
    h.append(f"<p>Pipeline run <code>{_esc(timestamp)}</code> processed "
             f"<strong>{len(dwi_types_present)}</strong> DWI type(s): "
             f"{', '.join(_dwi_badge(d) for d in dwi_types_present)}.</p>")
    if log_data:
        auc_items = []
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            roc = log_data[dwi_type].get("stats_predictive", {}).get("roc_analyses", [])
            best_auc = max((r.get("auc", 0) for r in roc), default=0)
            if best_auc > 0:
                auc_items.append(f"{_dwi_badge(dwi_type)} best AUC: <strong>{best_auc:.3f}</strong>")
        if auc_items:
            h.append("<ul>")
            for item in auc_items:
                h.append(f"<li>{item}</li>")
            h.append("</ul>")
    h.append("</div>")

    # ── 2. Graph Analysis Overview ──
    if rows:
        h.append("<h2>Graph Analysis Overview</h2>")

        # Count by type
        type_counts: dict[str, int] = {}
        for r in rows:
            gt = r.get("graph_type", "unknown")
            type_counts[gt] = type_counts.get(gt, 0) + 1

        h.append("<table><thead><tr><th>Graph Type</th><th>Count</th></tr></thead><tbody>")
        for gt, cnt in sorted(type_counts.items(), key=lambda x: -x[1]):
            h.append(f"<tr><td>{_esc(gt)}</td><td>{cnt}</td></tr>")
        h.append("</tbody></table>")

        # Count by DWI type
        dwi_counts: dict[str, int] = {}
        for r in rows:
            dwi_type, _ = parse_dwi_info(r["file_path"])
            dwi_counts[dwi_type] = dwi_counts.get(dwi_type, 0) + 1

        h.append("<table><thead><tr><th>DWI Type</th><th>Graphs</th></tr></thead><tbody>")
        for dt in DWI_TYPES + ["Root"]:
            if dt in dwi_counts:
                h.append(f"<tr><td>{_dwi_badge(dt)}</td><td>{dwi_counts[dt]}</td></tr>")
        h.append("</tbody></table>")

    # ── 3. Statistical Significance ──
    h.append("<h2>Statistical Significance</h2>")

    # From vision CSV
    sig_findings = []
    if rows:
        for r in rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            all_text = r["summary"] + " " + r["trends_json"] + " " + r["inflection_points_json"]
            for pval, context in extract_pvalues(all_text):
                if pval < 0.05:
                    sig_findings.append((pval, dwi_type, base_name, context))

    if sig_findings:
        sig_findings.sort()
        h.append("<h3>Vision-Extracted Significant Findings (p &lt; 0.05)</h3>")
        h.append("<table><thead><tr><th>p-value</th><th>Sig</th><th>DWI</th>"
                 "<th>Graph</th><th>Context</th></tr></thead><tbody>")
        for p, dwi, graph, ctx in sig_findings[:30]:
            ctx_clean = _esc(ctx.replace("\n", " ")[:100])
            cls = _sig_class(p)
            cls_attr = f' class="{cls}"' if cls else ""
            h.append(f"<tr><td{cls_attr}>{p:.4f}</td><td{cls_attr}>{_esc(_sig_tag(p))}</td>"
                     f"<td>{_dwi_badge(dwi)}</td><td>{_esc(graph)}</td>"
                     f"<td>{ctx_clean}</td></tr>")
        h.append("</tbody></table>")
        if len(sig_findings) > 30:
            h.append(f"<p><em>\u2026 and {len(sig_findings) - 30} more significant findings.</em></p>")

    # From direct CSV parsing
    if csv_data and csv_data.get("significant_metrics"):
        h.append("<h3>Pipeline CSV Significant Metrics</h3>")
        for dwi_type in DWI_TYPES:
            if dwi_type not in csv_data["significant_metrics"]:
                continue
            csv_rows = csv_data["significant_metrics"][dwi_type]
            h.append(f"<p>{_dwi_badge(dwi_type)} \u2014 {len(csv_rows)} significant metric(s)</p>")
            if csv_rows:
                headers = list(csv_rows[0].keys())[:6]
                h.append("<table><thead><tr>")
                for hdr in headers:
                    h.append(f"<th>{_esc(hdr)}</th>")
                h.append("</tr></thead><tbody>")
                for cr in csv_rows[:20]:
                    h.append("<tr>")
                    for hdr in headers:
                        h.append(f"<td>{_esc(str(cr.get(hdr, ''))[:30])}</td>")
                    h.append("</tr>")
                h.append("</tbody></table>")

    if not sig_findings and not (csv_data and csv_data.get("significant_metrics")):
        h.append("<p>No significant findings extracted.</p>")

    # ── 4. Cross-DWI Comparison ──
    if groups:
        h.append("<h2>Cross-DWI Comparison</h2>")

        priority_graphs = [
            "Dose_vs_Diffusion", "Longitudinal_Mean_Metrics",
            "Longitudinal_Mean_Metrics_ByOutcome", "Feature_BoxPlots",
            "Feature_Histograms", "core_method_dice_heatmap",
        ]

        for base_name in priority_graphs:
            if base_name not in groups:
                continue
            dwi_dict = groups[base_name]
            real = [t for t in dwi_dict if t != "Root"]
            if len(real) < 2:
                continue

            h.append(f"<h3>{_esc(base_name)}</h3>")

            all_trends: dict[str, list] = {}
            for dt in DWI_TYPES:
                if dt in dwi_dict:
                    try:
                        all_trends[dt] = json.loads(str(dwi_dict[dt].get("trends_json", "[]")))
                    except Exception:
                        pass

            if all_trends:
                all_series: set[str] = set()
                for trends in all_trends.values():
                    for t in trends:
                        all_series.add(t.get("series") or "overall")

                h.append("<table><thead><tr><th>Series</th>")
                for dt in DWI_TYPES:
                    h.append(f"<th>{_esc(dt)}</th>")
                h.append("<th>Agreement</th></tr></thead><tbody>")

                for series in sorted(all_series):
                    directions: dict[str, str] = {}
                    for dt in DWI_TYPES:
                        if dt not in all_trends:
                            continue
                        for t in all_trends[dt]:
                            if isinstance(t, dict):
                                if (t.get("series") or "overall") == series:
                                    directions[dt] = str(t.get("direction", ""))

                    if len(directions) >= 2:
                        vals = list(directions.values())
                        if len(set(vals)) == 1:
                            agree_html = '<span class="agree">AGREE</span>'
                        else:
                            agree_html = '<span class="differ">DIFFER</span>'
                        h.append(f"<tr><td>{_esc(series)}</td>")
                        for dt in DWI_TYPES:
                            h.append(f"<td>{_esc(directions.get(dt, '-'))}</td>")
                        h.append(f"<td>{agree_html}</td></tr>")

                h.append("</tbody></table>")

    # Cross-reference from CSV
    if csv_data and csv_data.get("cross_reference"):
        cross_ref_list = csv_data["cross_reference"]
        inconsistent = []
        if isinstance(cross_ref_list, list):
            for c in cross_ref_list:
                if isinstance(c, dict) and not c.get("consistent", True):
                    inconsistent.append(c)
        if inconsistent:
            h.append("<h3>Cross-DWI Significance Inconsistencies</h3>")
            h.append("<table><thead><tr><th>Metric</th><th>Timepoint</th>"
                     "<th>Significant In</th><th>Not Significant In</th></tr></thead><tbody>")
            for c in inconsistent:
                sig_in = c.get("significant_in", [])
                nsig_in = c.get("not_significant_in", [])
                
                sig_str = ", ".join(sig_in) if isinstance(sig_in, list) else "-"
                nsig_str = ", ".join(nsig_in) if isinstance(nsig_in, list) else "-"
                
                h.append(f"<tr><td>{_esc(str(c.get('metric', '')))}</td><td>{_esc(str(c.get('timepoint', '')))}</td>"
                         f"<td>{_esc(sig_str)}</td>"
                         f"<td>{_esc(nsig_str)}</td></tr>")
            h.append("</tbody></table>")

    # ── 5. Correlations ──
    if rows:
        h.append("<h2>Notable Correlations</h2>")
        corr_findings = []
        for r in rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            all_text = r["summary"] + " " + r["trends_json"]
            for rval, context in extract_correlations(all_text):
                if abs(rval) >= 0.3:
                    corr_findings.append((abs(rval), rval, dwi_type, base_name, context))

        if corr_findings:
            corr_findings.sort(reverse=True)
            h.append("<table><thead><tr><th>|r|</th><th>r</th><th>Strength</th>"
                     "<th>DWI</th><th>Graph</th></tr></thead><tbody>")
            for _, rval, dwi, graph, ctx in corr_findings[:20]:
                strength = "Strong" if abs(rval) >= 0.5 else "Moderate"
                h.append(f"<tr><td>{abs(rval):.2f}</td><td>{rval:+.2f}</td>"
                         f"<td>{strength}</td><td>{_dwi_badge(dwi)}</td>"
                         f"<td>{_esc(graph)}</td></tr>")
            h.append("</tbody></table>")
        else:
            h.append("<p>No notable correlations (|r| &ge; 0.3) found.</p>")

    # ── 6. Treatment Response ──
    if groups:
        has_longitudinal = any("Longitudinal" in bn for bn in groups)
        if has_longitudinal:
            h.append("<h2>Treatment Response \u2014 Longitudinal Trends</h2>")
            for base_name in sorted(groups.keys()):
                if "Longitudinal" not in base_name:
                    continue
                dwi_dict = groups[base_name]
                real = [t for t in dwi_dict if t != "Root"]
                if not real:
                    continue

                h.append(f"<h3>{_esc(base_name)}</h3>")
                h.append("<ul>")
                for dt in DWI_TYPES:
                    if dt not in dwi_dict:
                        continue
                    r = dwi_dict[dt]
                    summary = r["summary"]
                    if len(summary) > 200:
                        summary = summary[:200] + "..."
                    h.append(f"<li>{_dwi_badge(dt)}: {_esc(summary)}")

                    ips_str = r.get("inflection_points_json", "[]") if isinstance(r, dict) else "[]"
                    try:
                        ips = json.loads(str(ips_str))
                    except Exception:
                        ips = []
                        
                    if isinstance(ips, list) and ips:
                        h.append("<ul>")
                        for ip in ips:
                            if isinstance(ip, dict):
                                x = ip.get("approximate_x", "?")
                                h.append(f"<li>Inflection at x={_esc(str(x))}: {_esc(str(ip.get('description', '')))}</li>")
                        h.append("</ul>")
                    h.append("</li>")
                h.append("</ul>")

    # ── 7. Predictive Performance ──
    if log_data:
        h.append("<h2>Predictive Performance</h2>")

        has_roc = False
        roc_rows_html = []
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            roc = log_data[dwi_type].get("stats_predictive", {}).get("roc_analyses", [])
            for r in roc:
                has_roc = True
                auc_val = r.get("auc")
                sens_val = r.get("sensitivity")
                spec_val = r.get("specificity")
                auc = f"{auc_val:.3f}" if isinstance(auc_val, (int, float)) else "-"
                sens = f"{sens_val:.1f}%" if isinstance(sens_val, (int, float)) else "-"
                spec = f"{spec_val:.1f}%" if isinstance(spec_val, (int, float)) else "-"
                roc_rows_html.append(
                    f"<tr><td>{_dwi_badge(str(dwi_type))}</td><td>{_esc(str(r.get('timepoint', '-')))}</td>"
                    f"<td>{auc}</td><td>{sens}</td><td>{spec}</td></tr>"
                )

        if has_roc:
            h.append("<table><thead><tr><th>DWI</th><th>Timepoint</th><th>AUC</th>"
                     "<th>Sensitivity</th><th>Specificity</th></tr></thead><tbody>")
            h.extend(roc_rows_html)
            h.append("</tbody></table>")

        # Feature selections
        has_fs = False
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            fs = log_data[dwi_type].get("stats_predictive", {}).get("feature_selections", [])
            if fs:
                if not has_fs:
                    h.append("<h3>Selected Features (Elastic Net)</h3>")
                    has_fs = True
                h.append(f"<p>{_dwi_badge(dwi_type)}:</p><ul>")
                for sel in fs:
                    h.append(f"<li>{_esc(sel['timepoint'])} "
                             f"(\u03bb={sel['lambda']:.4f}): "
                             f"{_esc(', '.join(sel['features']))}</li>")
                h.append("</ul>")

        # Survival / Cox PH
        has_surv = False
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            sv = log_data[dwi_type].get("survival", {})
            hrs = sv.get("hazard_ratios", [])
            if hrs:
                if not has_surv:
                    h.append("<h3>Cox Proportional Hazards</h3>")
                    has_surv = True
                h.append(f"<p>{_dwi_badge(dwi_type)}:</p>")
                h.append("<table><thead><tr><th>Covariate</th><th>HR</th>"
                         "<th>95% CI</th><th>p</th></tr></thead><tbody>")
                for hr in hrs:
                    ci = f"[{hr.get('ci_lo', 0):.3f}, {hr.get('ci_hi', 0):.3f}]"
                    p_val = hr.get("p", 1.0)
                    sig = _sig_tag(p_val)
                    cls = _sig_class(p_val)
                    cls_attr = f' class="{cls}"' if cls else ""
                    h.append(f"<tr><td>{_esc(str(hr.get('covariate', '')))}</td>"
                             f"<td>{hr.get('hr', 0):.3f}</td><td>{ci}</td>"
                             f"<td{cls_attr}>{p_val:.4f} {_esc(sig)}</td></tr>")
                h.append("</tbody></table>")
                lrt = sv.get("global_lrt")
                if lrt:
                    h.append(f"<p>Global LRT: \u03c7\u00b2({lrt['df']}) = {lrt['chi2']:.2f}, "
                             f"p = {lrt['p']:.4f}</p>")

        if not has_roc and not has_fs and not has_surv:
            h.append("<p>No predictive performance data found in logs.</p>")

    # ── 7.5. Core Method & Dosimetry (MAT files) ──
    if mat_data:
        h.append("<h2>Supplemental Data (MAT Files)</h2>")
        
        # Dosimetry
        has_dos = any("dosimetry" in d for d in mat_data.values())
        if has_dos:
            h.append("<h3>Dosimetry (Target Coverage)</h3>")
            h.append("<table><thead><tr><th>DWI</th><th>D95 ADC (Gy)</th><th>V50 ADC (%)</th>"
                     "<th>D95 D (Gy)</th><th>V50 D (%)</th></tr></thead><tbody>")
            for dt in DWI_TYPES:
                if dt not in mat_data or "dosimetry" not in mat_data[dt]:
                    continue
                dos = mat_data[dt]["dosimetry"]
                if not dos:
                    continue
                h.append(
                    f"<tr><td>{_dwi_badge(dt)}</td>"
                    f"<td>{dos.get('d95_adc_mean', 'N/A'):.2f}</td>"
                    f"<td>{dos.get('v50_adc_mean', 'N/A') * 100:.1f}%</td>"
                    f"<td>{dos.get('d95_d_mean', 'N/A'):.2f}</td>"
                    f"<td>{dos.get('v50_d_mean', 'N/A') * 100:.1f}%</td></tr>"
                )
            h.append("</tbody></table>")

        # Core Methods
        has_core = any("core_method" in d for d in mat_data.values())
        if has_core:
            h.append("<h3>Core Method Comparison (Mean Dice)</h3>")
            for dt in DWI_TYPES:
                if dt not in mat_data or "core_method" not in mat_data[dt]:
                    continue
                core = mat_data[dt]["core_method"]
                if not core or not core.get("methods"):
                    continue
                h.append(f"<h4>{_dwi_badge(dt)}</h4>")
                methods = core["methods"]
                matrix = core["mean_dice_matrix"]
                
                # To keep it readable, we just show the top few methods against each other
                limit = min(6, len(methods))
                h.append("<table><thead><tr><th>Method</th>")
                for m in methods[:limit]:
                    h.append(f"<th>{_esc(m)}</th>")
                h.append("</tr></thead><tbody>")
                
                for i in range(limit):
                    h.append(f"<tr><td><strong>{_esc(methods[i])}</strong></td>")
                    for j in range(limit):
                        val = matrix[i][j]
                        if isinstance(val, (int, float)):
                            h.append(f"<td>{val:.2f}</td>")
                        else:
                            h.append("<td>-</td>")
                    h.append("</tr>")
                h.append("</tbody></table>")
                if len(methods) > limit:
                    h.append(f"<p><em>Showing {limit} of {len(methods)} methods. Full matrix available in .mat files.</em></p>")

    # ── 8. Appendix ──
    if rows:
        h.append("<h2>Appendix: All Graphs</h2>")
        h.append("<table><thead><tr><th>#</th><th>DWI</th><th>Type</th>"
                 "<th>Graph</th><th>Summary</th></tr></thead><tbody>")
        for i, r in enumerate(rows, 1):
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            summary = _esc(r["summary"][:80].replace("\n", " "))
            h.append(f"<tr><td>{i}</td><td>{_dwi_badge(dwi_type)}</td>"
                     f"<td>{_esc(r['graph_type'])}</td><td>{_esc(base_name)}</td>"
                     f"<td>{summary}\u2026</td></tr>")
        h.append("</tbody></table>")

    # Footer
    h.append(f"<footer>Report generated by pancData3 analysis suite on {now}</footer>")
    h.append("</body>")
    h.append("</html>")

    return "\n".join(h)


def markdown_to_html(md_text: str, title: str) -> str:
    """Convert Markdown text to a styled standalone HTML document."""
    body = markdown.markdown(
        md_text,
        extensions=["tables", "fenced_code", "toc"],
    )
    return HTML_TEMPLATE.format(title=title, body=body)


def main():
    folder = resolve_folder(sys.argv)
    md_report = generate_report(folder)

    timestamp = folder.name.replace("saved_files_", "")
    title = f"Analysis Report — {timestamp}"
    html_report = markdown_to_html(md_report, title)

    out_path = folder / "analysis_report.html"
    out_path.write_text(html_report, encoding="utf-8")
    print(f"Report written to: {out_path}")
    print(f"  Length: {len(html_report)} characters")


if __name__ == "__main__":
    main()
