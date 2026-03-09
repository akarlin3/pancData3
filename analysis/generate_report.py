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
import re
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


def _section(title: str, level: int = 2) -> str:
    """Return a Markdown-style section heading string (utility / test helper)."""
    return f"\n{'#' * level} {title}\n"


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
h4 { font-size: 0.98rem; margin: 1rem 0 0.4rem; color: var(--muted); }
p, ul { margin-bottom: 0.75rem; }
ul { padding-left: 1.5rem; }
.meta { color: var(--muted); font-size: 0.9rem; margin-bottom: 1.5rem; }
.meta span { display: inline-block; margin-right: 1.5rem; }
code { background: var(--row-alt); border: 1px solid var(--border); border-radius: 3px; padding: 0.1em 0.35em; font-size: 0.88em; }
table { width: 100%; border-collapse: collapse; margin: 0.75rem 0 1.25rem; font-size: 0.92rem; }
th { background: var(--accent); color: #fff; text-align: left; padding: 0.5rem 0.75rem; font-weight: 600; }
td { padding: 0.4rem 0.75rem; border-bottom: 1px solid var(--border); vertical-align: top; }
tr:nth-child(even) td { background: var(--row-alt); }
tr:hover td { background: var(--accent-light); }
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
.warn-box { background: #fef3c7; border-left: 4px solid var(--amber); padding: 0.75rem 1rem; border-radius: 0 6px 6px 0; margin: 0.75rem 0; font-size: 0.9rem; }
.info-box { background: #f0fdf4; border-left: 4px solid var(--green); padding: 0.75rem 1rem; border-radius: 0 6px 6px 0; margin: 0.75rem 0; font-size: 0.9rem; }
footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid var(--border); color: var(--muted); font-size: 0.85rem; }
nav.toc {
    position: sticky; top: 0; background: var(--bg);
    border-bottom: 2px solid var(--border); padding: 0.5rem 0 0.4rem;
    margin-bottom: 1.5rem; z-index: 100; overflow-x: auto;
    white-space: nowrap; font-size: 0.83rem;
}
nav.toc a { color: var(--accent); text-decoration: none; margin-right: 1.1rem; }
nav.toc a:hover { text-decoration: underline; }
details { margin: 0.3rem 0; }
details > summary {
    cursor: pointer; color: var(--accent); font-size: 0.88rem;
    padding: 0.2rem 0.1rem; list-style: none; user-select: none;
}
details > summary::before { content: "\25B6\00A0"; font-size: 0.7em; }
details[open] > summary::before { content: "\25BC\00A0"; }
details[open] > summary { margin-bottom: 0.4rem; }
.axis-info { font-size: 0.82rem; color: var(--muted); line-height: 1.4; }
.trend-tag {
    display: inline-block; padding: 0.1em 0.45em; border-radius: 3px;
    font-size: 0.8rem; margin: 0.1em 0.1em; border: 1px solid var(--border);
}
.trend-incr { background: #dcfce7; border-color: #86efac; color: #166534; }
.trend-decr { background: #fee2e2; border-color: #fca5a5; color: #991b1b; }
.trend-flat { background: var(--row-alt); }
.trend-nm   { background: #f5f3ff; border-color: #c4b5fd; color: #5b21b6; }
.full-summary { font-size: 0.88rem; color: #374151; margin-top: 0.25rem; }
.stat-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(220px, 1fr)); gap: 0.75rem; margin: 0.75rem 0 1rem; }
.stat-card { background: var(--row-alt); border: 1px solid var(--border); border-radius: 6px; padding: 0.75rem 1rem; }
.stat-card .label { font-size: 0.78rem; color: var(--muted); text-transform: uppercase; letter-spacing: 0.04em; }
.stat-card .value { font-size: 1.35rem; font-weight: 700; color: var(--accent); }
.stat-card .sub { font-size: 0.8rem; color: var(--muted); }
"""


def _dwi_badge(dwi_type: str) -> str:
    cls = {
        "Standard": "badge-standard", "dnCNN": "badge-dncnn",
        "IVIMnet": "badge-ivimnet",
    }.get(dwi_type, "badge-root")
    return f'<span class="badge {cls}">{_esc(dwi_type)}</span>'


def _trend_tag(direction: str) -> str:
    d = direction.lower()
    if "increas" in d or "up" in d or "higher" in d or "rising" in d:
        cls = "trend-incr"
        arrow = "\u2191\u00a0"
    elif "decreas" in d or "down" in d or "lower" in d or "falling" in d or "drop" in d:
        cls = "trend-decr"
        arrow = "\u2193\u00a0"
    elif "flat" in d or "stable" in d or "constant" in d:
        cls = "trend-flat"
        arrow = "\u2192\u00a0"
    else:
        cls = "trend-nm"
        arrow = ""
    return f'<span class="trend-tag {cls}">{arrow}{_esc(direction)}</span>'


# ── Navigation sections ────────────────────────────────────────────────────────

NAV_SECTIONS = [
    ("exec-summary", "Executive Summary"),
    ("data-quality", "Data Quality"),
    ("hypothesis", "Hypothesis"),
    ("graph-overview", "Graphs"),
    ("significance", "Statistics"),
    ("cross-dwi", "Cross-DWI"),
    ("fdr-global", "FDR Global"),
    ("correlations", "Correlations"),
    ("treatment", "Treatment"),
    ("predictive", "Predictive"),
    ("supplemental", "Supplemental"),
    ("appendix", "All Graphs"),
]


def _nav_bar() -> str:
    links = "".join(
        f'<a href="#{anchor}">{_esc(label)}</a>'
        for anchor, label in NAV_SECTIONS
    )
    return f'<nav class="toc">{links}</nav>'


def _h2(text: str, anchor: str) -> str:
    return f'<h2 id="{anchor}">{_esc(text)}</h2>'


def _stat_card(label: str, value: str, sub: str = "") -> str:
    sub_html = f'<div class="sub">{_esc(sub)}</div>' if sub else ""
    return (
        f'<div class="stat-card">'
        f'<div class="label">{_esc(label)}</div>'
        f'<div class="value">{_esc(value)}</div>'
        f'{sub_html}</div>'
    )


def _get_consensus(trend_list: list[str]) -> str:
    if not trend_list: return "unknown"
    increasers = sum(1 for x in trend_list if "increas" in x or "higher" in x or "up" in x)
    decreasers = sum(1 for x in trend_list if "decreas" in x or "lower" in x or "down" in x)
    if increasers > decreasers: return "increasing"
    if decreasers > increasers: return "decreasing"
    return "stable"


def _section_executive_summary(log_data, dwi_types_present, rows, csv_data, timestamp) -> list[str]:
    # ── 1. Executive Summary ──
    h: list[str] = []
    h.append(_h2("Executive Summary", "exec-summary"))
    h.append('<div class="summary-box">')
    h.append(f"<p>Pipeline run <code>{_esc(timestamp)}</code> processed "
             f"<strong>{len(dwi_types_present)}</strong> DWI type(s): "
             f"{', '.join(_dwi_badge(d) for d in dwi_types_present)}.</p>")

    # Stat cards row
    cards = []
    if rows:
        cards.append(_stat_card("Graphs Analysed", str(len(rows))))

    total_sig = 0
    if log_data:
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            roc = log_data[dwi_type].get("stats_predictive", {}).get("roc_analyses", [])
            best_auc = max((r.get("auc", 0) for r in roc), default=0)
            if best_auc > 0:
                cards.append(_stat_card(f"Best AUC ({dwi_type})", f"{best_auc:.3f}"))
            sc = log_data[dwi_type].get("stats_comparisons", {})
            total_sig += len([g for g in sc.get("glme_details", []) if g["p"] < g["adj_alpha"]])

    if total_sig > 0:
        cards.append(_stat_card("GLME Sig. Interactions", str(total_sig), "across all DWI types"))

    if csv_data and csv_data.get("significant_metrics"):
        n_csv_sig = sum(len(v) for v in csv_data["significant_metrics"].values())
        cards.append(_stat_card("CSV Sig. Metrics", str(n_csv_sig), "from pipeline CSVs"))

    if cards:
        h.append('<div class="stat-grid">')
        h.extend(cards)
        h.append("</div>")

    h.append("</div>")
    return h


def _section_hypothesis(groups) -> list[str]:
    # ── 1.5. Data-Driven Hypothesis ──
    h: list[str] = []
    h.append(_h2("Data-Driven Hypothesis", "hypothesis"))
    h.append('<div class="summary-box">')

    # Analyze trends for dynamic hypothesis
    d_trend_consensus = "unknown"
    f_trend_consensus = "unknown"

    vascular_specificity = ""
    cellular_specificity = ""
    vascular_inflections: list[str] = []
    cellular_inflections: list[str] = []

    if groups and "Longitudinal_Mean_Metrics" in groups:
        d_trends = []
        f_trends = []

        for dt, r in groups["Longitudinal_Mean_Metrics"].items():
            if dt == "Root": continue
            try:
                # 1. Parse Trends
                trends = json.loads(str(r.get("trends_json", "[]")))
                for t in trends:
                    if not isinstance(t, dict): continue
                    series = t.get("series", "")
                    direction = t.get("direction", "").lower()
                    if series == "Mean D":
                        d_trends.append(direction)
                    elif series == "Mean f":
                        f_trends.append(direction)

                # 2. Parse Inflection Points
                ips = json.loads(str(r.get("inflection_points_json", "[]")))
                for ip in ips:
                    if not isinstance(ip, dict): continue
                    x_val = float(ip.get("approximate_x", 0))
                    y_val = ip.get("approximate_y")
                    desc = ip.get("description", "").lower()

                    fx_label = f"Fx{int(x_val)}" if x_val > 0 else "baseline"

                    magnitude = ""
                    if y_val is not None:
                        try:
                            mag_val = abs(float(y_val))
                            if mag_val > 1:
                                magnitude = f" of ~{int(mag_val)}%"
                        except ValueError:
                            pass

                    if not magnitude:
                        pct_match = re.search(r'(\d+)%', desc)
                        if pct_match:
                            magnitude = f" of ~{pct_match.group(1)}%"

                    if "d*" in desc or "f" in desc or "vascular" in desc or "perfusion" in desc:
                        if "drop" in desc or "decrease" in desc or "decline" in desc:
                            vascular_inflections.append(f"a significant vascular drop{magnitude} observed around {fx_label}")
                        elif "increase" in desc or "rise" in desc or "peak" in desc:
                            vascular_inflections.append(f"a vascular increase{magnitude} observed around {fx_label}")

                    if "adc" in desc or " d " in desc or "diffusion" in desc:
                        if "increase" in desc or "rise" in desc or "peak" in desc:
                            cellular_inflections.append(f"a significant diffusion increase{magnitude} observed around {fx_label}")
                        elif "drop" in desc or "decrease" in desc or "decline" in desc:
                            cellular_inflections.append(f"a diffusion decrease{magnitude} observed around {fx_label}")

            except Exception:
                pass

        d_trend_consensus = _get_consensus(d_trends)
        f_trend_consensus = _get_consensus(f_trends)

        vascular_specificity = ""
        if vascular_inflections:
            unique_v = list(dict.fromkeys(vascular_inflections))
            vascular_specificity = f" Specifically, the data highlights {', and '.join(unique_v)}."

        cellular_specificity = ""
        if cellular_inflections:
            unique_c = list(dict.fromkeys(cellular_inflections))
            cellular_specificity = f" Specifically, the data highlights {', and '.join(unique_c)}."

    h.append("<p>Based on quantitative metrics and longitudinal trends extracted from the data, "
             "the following radiological-pathological hypothesis is proposed:</p>")
    h.append("<ul>")

    if d_trend_consensus == "increasing":
        h.append(f"<li><strong>Cellular Response (D, ADC):</strong> The data shows an increase in true diffusion (<em>D</em>) over time.{cellular_specificity} "
                 "This suggests that effective radiation therapy is inducing cellular necrosis and apoptosis, "
                 "leading to a breakdown of cell membranes. This expands the extracellular space, explaining the increased water mobility.</li>")
    elif d_trend_consensus == "decreasing":
        h.append(f"<li><strong>Cellular Response (D, ADC):</strong> The data shows a decrease in true diffusion (<em>D</em>) over time.{cellular_specificity} "
                 "This suggests limited cell kill or potential cellular swelling (cytotoxic edema), indicating a highly cellular, densely packed tumor resistant to therapy.</li>")
    else:
        h.append(f"<li><strong>Cellular Response (D, ADC):</strong> The data shows relatively stable or variable true diffusion (<em>D</em>) over time.{cellular_specificity} "
                 "This suggests a steady state between cellular destruction and tumor proliferation, or a timeline where major necrotic changes are not yet dominant.</li>")

    if f_trend_consensus == "decreasing":
        h.append(f"<li><strong>Vascular Response (f, D*):</strong> The data reveals drops in microcapillary perfusion fraction (<em>f</em>) and/or pseudo-diffusion (<em>D*</em>).{vascular_specificity} "
                 "These decreases indicate that radiation causes early endothelial damage and vascular regression, effectively cutting off the tumor's blood supply.</li>")
    elif f_trend_consensus == "increasing":
        h.append(f"<li><strong>Vascular Response (f, D*):</strong> The data shows an increase in microcapillary perfusion fraction (<em>f</em>).{vascular_specificity} "
                 "This suggests reactive angiogenesis, hyperemic inflammatory response, or a robust vascular supply aiding tumor survival and radiation resistance.</li>")
    else:
        h.append(f"<li><strong>Vascular Response (f, D*):</strong> The microcapillary perfusion fraction (<em>f</em>) remains relatively stable,{vascular_specificity} "
                 "suggesting that the tumor's vascular network has not been significantly altered or compromised by the treatment doses applied so far.</li>")

    outcome_specificity = ""
    if vascular_inflections or cellular_inflections:
        combined = list(dict.fromkeys(cellular_inflections + vascular_inflections))
        outcome_specificity = f" (supported by {', '.join(combined)})"

    if d_trend_consensus == "increasing" and f_trend_consensus == "decreasing":
        h.append(f"<li><strong>Outcome Trajectory:</strong> The combination of early decreases in perfusion coupled with "
                 f"subsequent increases in diffusion strongly supports the hypothesis that the tumor is responding effectively to treatment{outcome_specificity}, "
                 "correlating with long-term Local Control.</li>")
    else:
        h.append(f"<li><strong>Outcome Trajectory:</strong> The observed variation in diffusion and perfusion responses suggests "
                 f"a heterogeneous or limited overall treatment effect{outcome_specificity}. Tumors not exhibiting strong, concomitant increases in diffusion and decreases in perfusion "
                 "may remain at higher risk for Local Failure or harbor therapy-resistant sub-volumes.</li>")

    h.append("</ul>")
    h.append("</div>")
    return h


def _section_graph_overview(rows) -> list[str]:
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
                h.append(f"<tr><td>{_dwi_badge(dt)}</td><td>{dwi_counts[dt]}</td></tr>")
        h.append("</tbody></table>")
        h.append("</div>")

        h.append("</div>")
    return h


def _section_statistical_significance(rows, csv_data, log_data, dwi_types_present) -> list[str]:
    # ── 3. Statistical Significance ──
    h: list[str] = []
    h.append(_h2("Statistical Significance", "significance"))

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

    # GLME interaction details from logs
    if log_data:
        has_glme = any(
            log_data[d].get("stats_comparisons", {}).get("glme_details")
            for d in dwi_types_present if d in log_data
        )
        if has_glme:
            h.append("<h3>GLME Mixed-Effects Interaction Tests</h3>")
            h.append("<p><em>Per-metric interaction p-values from generalised linear mixed-effects models. "
                     "Adjusted alpha accounts for BH-FDR multiple comparison correction.</em></p>")

            for dwi_type in dwi_types_present:
                if dwi_type not in log_data:
                    continue
                sc = log_data[dwi_type].get("stats_comparisons", {})
                glme_details = sc.get("glme_details", [])
                glme_excl = sc.get("glme_excluded")
                if not glme_details:
                    continue

                h.append(f"<h4>{_dwi_badge(dwi_type)}</h4>")

                if glme_excl:
                    h.append(f'<div class="warn-box">'
                              f'\u26a0\ufe0f Competing-risk exclusion: '
                              f'<strong>{glme_excl["n_excluded"]}/{glme_excl["n_total"]}</strong> '
                              f'({glme_excl["pct"]:.1f}%) patients excluded from GLME.</div>')

                sig_details = [g for g in glme_details if g["p"] < g["adj_alpha"]]
                nsig_details = [g for g in glme_details if g["p"] >= g["adj_alpha"]]

                if sig_details:
                    h.append(f'<div class="info-box">'
                              f'\u2705 {len(sig_details)} of {len(glme_details)} metrics show significant '
                              f'interaction (p &lt; adjusted \u03b1).</div>')

                h.append("<table><thead><tr>"
                         "<th>Metric</th><th>p-value</th><th>Adj. \u03b1</th><th>Sig.</th>"
                         "</tr></thead><tbody>")
                for g in sorted(glme_details, key=lambda x: x["p"]):
                    is_sig = g["p"] < g["adj_alpha"]
                    cls = _sig_class(g["p"]) if is_sig else ""
                    cls_attr = f' class="{cls}"' if cls else ""
                    sig_mark = _sig_tag(g["p"]) if is_sig else ""
                    h.append(
                        f"<tr><td><code>{_esc(g['metric'])}</code></td>"
                        f"<td{cls_attr}>{g['p']:.4f}</td>"
                        f"<td>{g['adj_alpha']:.4f}</td>"
                        f"<td{cls_attr}>{_esc(sig_mark) or '\u2014'}</td></tr>"
                    )
                h.append("</tbody></table>")

        # FDR-significant timepoints
        has_fdr_tp = any(
            log_data[d].get("stats_comparisons", {}).get("fdr_timepoints")
            for d in dwi_types_present if d in log_data
        )
        if has_fdr_tp:
            h.append("<h3>FDR-Significant Metrics by Timepoint</h3>")
            h.append("<table><thead><tr><th>DWI</th><th>Timepoint</th>"
                     "<th>N Significant</th></tr></thead><tbody>")
            for dwi_type in dwi_types_present:
                if dwi_type not in log_data:
                    continue
                fdr_tps = log_data[dwi_type].get("stats_comparisons", {}).get("fdr_timepoints", [])
                for ftp in fdr_tps:
                    h.append(
                        f"<tr><td>{_dwi_badge(dwi_type)}</td>"
                        f"<td><code>{_esc(ftp['timepoint'])}</code></td>"
                        f"<td><strong>{ftp['n_significant']}</strong></td></tr>"
                    )
            h.append("</tbody></table>")

    if not sig_findings and not (csv_data and csv_data.get("significant_metrics")) and not log_data:
        h.append("<p>No significant findings extracted.</p>")
    return h


def _section_cross_dwi_comparison(groups, csv_data) -> list[str]:
    # ── 4. Cross-DWI Comparison ──
    h: list[str] = []
    if groups:
        h.append(_h2("Cross-DWI Comparison", "cross-dwi"))

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
                            d_str = directions.get(dt, "-")
                            cell = _trend_tag(d_str) if d_str != "-" else "-"
                            h.append(f"<td>{cell}</td>")
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
                h.append(f"<tr><td>{_esc(str(c.get('metric', '')))}</td>"
                         f"<td>{_esc(str(c.get('timepoint', '')))}</td>"
                         f"<td>{_esc(sig_str)}</td>"
                         f"<td><span class='differ'>{_esc(nsig_str)}</span></td></tr>")
            h.append("</tbody></table>")
    return h


def _section_correlations(rows) -> list[str]:
    # ── 5. Correlations ──
    h: list[str] = []
    if rows:
        h.append(_h2("Notable Correlations", "correlations"))
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
                     "<th>DWI</th><th>Graph</th><th>Context</th></tr></thead><tbody>")
            for _, rval, dwi, graph, ctx in corr_findings[:20]:
                strength = "Strong" if abs(rval) >= 0.5 else "Moderate"
                ctx_short = _esc(ctx.replace("\n", " ")[:80])
                h.append(f"<tr><td><strong>{abs(rval):.2f}</strong></td><td>{rval:+.2f}</td>"
                         f"<td>{strength}</td><td>{_dwi_badge(dwi)}</td>"
                         f"<td>{_esc(graph)}</td><td><em>{ctx_short}</em></td></tr>")
            h.append("</tbody></table>")
        else:
            h.append("<p>No notable correlations (|r| &ge; 0.3) found.</p>")
    return h


def _section_treatment_response(groups) -> list[str]:
    # ── 6. Treatment Response ──
    h: list[str] = []
    if groups:
        has_longitudinal = any("Longitudinal" in bn for bn in groups)
        if has_longitudinal:
            h.append(_h2("Treatment Response \u2014 Longitudinal Trends", "treatment"))
            for base_name in sorted(groups.keys()):
                if "Longitudinal" not in base_name:
                    continue
                dwi_dict = groups[base_name]
                real = [t for t in dwi_dict if t != "Root"]
                if not real:
                    continue

                h.append(f"<h3>{_esc(base_name)}</h3>")
                for dt in DWI_TYPES:
                    if dt not in dwi_dict:
                        continue
                    r = dwi_dict[dt]
                    summary = r.get("summary", "")

                    # Axis info
                    x_lbl = r.get("x_axis_label", "")
                    x_unit = r.get("x_axis_units", "")
                    y_lbl = r.get("y_axis_label", "")
                    y_unit = r.get("y_axis_units", "")
                    axis_parts = []
                    if x_lbl:
                        axis_parts.append(f"X: {x_lbl}" + (f" ({x_unit})" if x_unit else ""))
                    if y_lbl:
                        axis_parts.append(f"Y: {y_lbl}" + (f" ({y_unit})" if y_unit else ""))

                    h.append(f"<p>{_dwi_badge(dt)}")
                    if axis_parts:
                        h.append(f' <span class="axis-info">[{" | ".join(axis_parts)}]</span>')
                    h.append("</p>")

                    # Trends
                    trends_str = r.get("trends_json", "[]") if isinstance(r, dict) else "[]"
                    try:
                        trends = json.loads(str(trends_str))
                    except Exception:
                        trends = []
                    if isinstance(trends, list) and trends:
                        for t in trends:
                            if isinstance(t, dict):
                                series = t.get("series", "")
                                direction = t.get("direction", "")
                                desc = t.get("description", "")
                                h.append(f"<p style='margin-left:1rem'>")
                                if series:
                                    h.append(f"<code>{_esc(series)}</code> ")
                                h.append(_trend_tag(direction))
                                if desc:
                                    h.append(f" <em>{_esc(desc)}</em>")
                                h.append("</p>")

                    # Inflection points
                    ips_str = r.get("inflection_points_json", "[]") if isinstance(r, dict) else "[]"
                    try:
                        ips = json.loads(str(ips_str))
                    except Exception:
                        ips = []
                    if isinstance(ips, list) and ips:
                        h.append("<details><summary>Inflection points</summary><ul>")
                        for ip in ips:
                            if isinstance(ip, dict):
                                x = ip.get("approximate_x", "?")
                                y = ip.get("approximate_y", "")
                                desc_ip = ip.get("description", "")
                                coord = f"x={_esc(str(x))}"
                                if y is not None and y != "":
                                    coord += f", y={_esc(str(y))}"
                                h.append(f"<li>({coord}): {_esc(str(desc_ip))}</li>")
                        h.append("</ul></details>")

                    # Full summary in collapsible
                    if summary:
                        if len(summary) > 200:
                            h.append(f"<details><summary>Full summary</summary>"
                                     f'<p class="full-summary">{_esc(summary)}</p></details>')
                        else:
                            h.append(f'<p class="full-summary" style="margin-left:1rem">{_esc(summary)}</p>')
    return h


def _section_predictive_performance(log_data, dwi_types_present) -> list[str]:
    # ── 7. Predictive Performance ──
    h: list[str] = []
    if log_data:
        h.append(_h2("Predictive Performance", "predictive"))

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
                youden = r.get("youden_cutoff")
                auc = f"{auc_val:.3f}" if isinstance(auc_val, (int, float)) else "-"
                sens = f"{sens_val:.1f}%" if isinstance(sens_val, (int, float)) else "-"
                spec = f"{spec_val:.1f}%" if isinstance(spec_val, (int, float)) else "-"
                youden_str = f"{youden:.3f}" if isinstance(youden, (int, float)) else "-"
                roc_rows_html.append(
                    f"<tr><td>{_dwi_badge(str(dwi_type))}</td>"
                    f"<td>{_esc(str(r.get('timepoint', '-')))}</td>"
                    f"<td><strong>{auc}</strong></td>"
                    f"<td>{sens}</td><td>{spec}</td>"
                    f"<td>{youden_str}</td></tr>"
                )

        if has_roc:
            h.append("<h3>ROC / AUC Performance</h3>")
            h.append("<table><thead><tr><th>DWI</th><th>Timepoint</th><th>AUC</th>"
                     "<th>Sensitivity</th><th>Specificity</th><th>Youden Cutoff</th>"
                     "</tr></thead><tbody>")
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
                    features_html = ", ".join(
                        f"<code>{_esc(f)}</code>" for f in sel["features"]
                    )
                    h.append(f"<li><strong>{_esc(sel['timepoint'])}</strong> "
                             f"(\u03bb={sel['lambda']:.4f}): {features_html}</li>")
                h.append("</ul>")

        # Survival / Cox PH
        has_surv = False
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            sv = log_data[dwi_type].get("survival", {})
            hrs = sv.get("hazard_ratios", [])
            ipcw = sv.get("ipcw")
            if hrs:
                if not has_surv:
                    h.append("<h3>Cox Proportional Hazards</h3>")
                    has_surv = True
                h.append(f"<p>{_dwi_badge(dwi_type)}:")
                if ipcw:
                    h.append(f' <span class="meta">IPCW weights: '
                             f'[{ipcw["min_weight"]:.3f}, {ipcw["max_weight"]:.3f}]</span>')
                h.append("</p>")
                h.append("<table><thead><tr><th>Covariate</th><th>HR</th>"
                         "<th>95% CI</th><th>p</th></tr></thead><tbody>")
                for hr in sorted(hrs, key=lambda x: x.get("p", 1)):
                    ci = f"[{hr.get('ci_lo', 0):.3f}, {hr.get('ci_hi', 0):.3f}]"
                    p_val = hr.get("p", 1.0)
                    sig = _sig_tag(p_val)
                    cls = _sig_class(p_val)
                    cls_attr = f' class="{cls}"' if cls else ""
                    h.append(f"<tr><td><code>{_esc(str(hr.get('covariate', '')))}</code></td>"
                             f"<td>{hr.get('hr', 0):.3f}</td><td>{ci}</td>"
                             f"<td{cls_attr}>{p_val:.4f} {_esc(sig)}</td></tr>")
                h.append("</tbody></table>")
                lrt = sv.get("global_lrt")
                if lrt:
                    cls_lrt = _sig_class(lrt["p"])
                    h.append(f"<p>Global LRT: \u03c7\u00b2({lrt['df']}) = {lrt['chi2']:.2f}, "
                             f"p = <span{' class=' + repr(cls_lrt) if cls_lrt else ''}>{lrt['p']:.4f}</span></p>")

        if not has_roc and not has_fs and not has_surv:
            h.append("<p>No predictive performance data found in logs.</p>")
    return h


def _section_mat_data(mat_data) -> list[str]:
    # ── 7.5. Core Method & Dosimetry ──
    h: list[str] = []
    if mat_data:
        h.append(_h2("Supplemental Data (MAT Files)", "supplemental"))

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
                            # Colour-code: high Dice = green-ish
                            if i != j:
                                style = f" style=\"color: {'var(--green)' if val >= 0.7 else ('var(--amber)' if val >= 0.5 else 'var(--red)')}; font-weight:600\""
                            else:
                                style = ""
                            h.append(f"<td{style}>{val:.2f}</td>")
                        else:
                            h.append("<td>-</td>")
                    h.append("</tr>")
                h.append("</tbody></table>")
                if len(methods) > limit:
                    h.append(f"<p><em>Showing {limit} of {len(methods)} methods. "
                             f"Full matrix available in .mat files.</em></p>")
    return h


def _section_appendix(rows) -> list[str]:
    # ── 8. Appendix ──
    h: list[str] = []
    if rows:
        h.append(_h2("Appendix: All Graphs", "appendix"))
        h.append(f"<p>{len(rows)} graphs analysed. Expand each row for full details.</p>")
        h.append("<table><thead><tr>"
                 "<th>#</th><th>DWI</th><th>Type</th><th>Graph</th>"
                 "<th>Title</th><th>Axes</th><th>Trends</th><th>Summary</th>"
                 "</tr></thead><tbody>")
        for i, r in enumerate(rows, 1):
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            graph_title = r.get("graph_title", "") or ""
            graph_type = r.get("graph_type", "")

            # Axes summary
            x_lbl = r.get("x_axis_label", "")
            x_unit = r.get("x_axis_units", "")
            x_min = r.get("x_axis_range_min", "")
            x_max = r.get("x_axis_range_max", "")
            y_lbl = r.get("y_axis_label", "")
            y_unit = r.get("y_axis_units", "")
            y_min = r.get("y_axis_range_min", "")
            y_max = r.get("y_axis_range_max", "")
            c_lbl = r.get("color_axis_label", "")
            c_unit = r.get("color_axis_units", "")

            axis_lines = []
            if x_lbl:
                xr = f" [{x_min}\u2013{x_max}]" if x_min and x_max else ""
                xu = f" ({x_unit})" if x_unit else ""
                axis_lines.append(f"X: {x_lbl}{xu}{xr}")
            if y_lbl:
                yr = f" [{y_min}\u2013{y_max}]" if y_min and y_max else ""
                yu = f" ({y_unit})" if y_unit else ""
                axis_lines.append(f"Y: {y_lbl}{yu}{yr}")
            if c_lbl:
                cu = f" ({c_unit})" if c_unit else ""
                axis_lines.append(f"C: {c_lbl}{cu}")
            axes_cell = '<br>'.join(_esc(a) for a in axis_lines) if axis_lines else "\u2014"

            # Trends cell
            trends_str = r.get("trends_json", "[]") or "[]"
            try:
                trends_list = json.loads(trends_str)
            except Exception:
                trends_list = []
            trends_cell = ""
            if isinstance(trends_list, list) and trends_list:
                tags = []
                for t in trends_list[:4]:  # show max 4 trend tags
                    if isinstance(t, dict):
                        direction = t.get("direction", "")
                        series = t.get("series", "")
                        label = f"{series}: {direction}" if series else direction
                        tags.append(_trend_tag(label))
                trends_cell = "".join(tags)
                if len(trends_list) > 4:
                    trends_cell += f" <em>+{len(trends_list) - 4} more</em>"

            # Summary: short preview + collapsible full
            summary = r.get("summary", "") or ""
            summary_preview = _esc(summary[:100].replace("\n", " "))
            if len(summary) > 100:
                summary_cell = (
                    f"{summary_preview}\u2026"
                    f"<details><summary>more</summary>"
                    f'<p class="full-summary">{_esc(summary)}</p>'
                    f"</details>"
                )
            else:
                summary_cell = summary_preview

            h.append(
                f"<tr>"
                f"<td>{i}</td>"
                f"<td>{_dwi_badge(dwi_type)}</td>"
                f"<td>{_esc(graph_type)}</td>"
                f"<td>{_esc(base_name)}</td>"
                f'<td class="axis-info">{_esc(graph_title) if graph_title else "\u2014"}</td>'
                f'<td class="axis-info">{axes_cell}</td>'
                f"<td>{trends_cell}</td>"
                f"<td>{summary_cell}</td>"
                f"</tr>"
            )
        h.append("</tbody></table>")
    return h


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

    # ── Sticky nav ──
    h.append(_nav_bar())

    # ── Header ──
    h.append(f"<h1>Analysis Report \u2014 {_esc(timestamp)}</h1>")
    h.append('<div class="meta">')
    h.append(f"<span><strong>Generated:</strong> {now}</span>")
    h.append(f"<span><strong>DWI types:</strong> {', '.join(dwi_types_present) or 'None detected'}</span>")
    if rows:
        h.append(f"<span><strong>Graphs analysed:</strong> {len(rows)}</span>")
    h.append("</div>")

    h.extend(_section_executive_summary(log_data, dwi_types_present, rows, csv_data, timestamp))

    # ── 2. Data Quality ──
    h.append(_h2("Data Quality", "data-quality"))
    has_quality_data = False

    if log_data:
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            bl = log_data[dwi_type].get("baseline", {})
            outlier_flags = bl.get("outlier_flags", [])
            total_out = bl.get("total_outliers")
            baseline_exc = bl.get("baseline_exclusion")

            if not outlier_flags and not total_out and not baseline_exc:
                continue

            has_quality_data = True
            h.append(f"<h3>{_dwi_badge(dwi_type)} Baseline Quality</h3>")

            # Summary stats row
            qcards = []
            if total_out:
                pct = total_out.get("pct", 0)
                qcards.append(_stat_card(
                    "Outliers Removed",
                    f"{total_out['n_removed']}/{total_out['n_total']}",
                    f"{pct:.1f}% of samples"
                ))
            if baseline_exc:
                qcards.append(_stat_card(
                    "Baseline Excluded",
                    f"{baseline_exc['n_excluded']}/{baseline_exc['n_total']}",
                    "missing baseline data"
                ))
                lf_inc = baseline_exc.get("lf_rate_included")
                lf_exc = baseline_exc.get("lf_rate_excluded")
                if lf_inc is not None and lf_exc is not None:
                    qcards.append(_stat_card(
                        "LF Rate",
                        f"{lf_inc:.1f}% / {lf_exc:.1f}%",
                        "included / excluded"
                    ))
            if qcards:
                h.append('<div class="stat-grid">')
                h.extend(qcards)
                h.append("</div>")

            # Per-metric outlier flags table
            if outlier_flags:
                h.append("<h4>Per-Metric Outlier Flags</h4>")
                h.append("<table><thead><tr>"
                         "<th>Metric</th><th>Flagged</th>"
                         "<th>Local Failure</th><th>Local Control</th><th>Competing Risk</th>"
                         "</tr></thead><tbody>")
                for o in outlier_flags:
                    h.append(
                        f"<tr><td><code>{_esc(o['metric'])}</code></td>"
                        f"<td><strong>{o['n_flagged']}</strong></td>"
                        f"<td>{o['n_lf']}</td><td>{o['n_lc']}</td><td>{o['n_cr']}</td></tr>"
                    )
                h.append("</tbody></table>")

    if not has_quality_data:
        h.append('<p class="meta">No baseline quality data found in logs.</p>')

    h.extend(_section_hypothesis(groups))
    h.extend(_section_graph_overview(rows))
    h.extend(_section_statistical_significance(rows, csv_data, log_data, dwi_types_present))
    h.extend(_section_cross_dwi_comparison(groups, csv_data))

    # ── 7. FDR Global ──
    if csv_data and csv_data.get("fdr_global"):
        fdr_global = csv_data["fdr_global"]
        if fdr_global:
            h.append(_h2("FDR Global Correction Results", "fdr-global"))
            h.append("<p><em>Metrics surviving Benjamini-Hochberg FDR correction across the full "
                     "feature set, as exported by the pipeline.</em></p>")
            for dwi_type in DWI_TYPES:
                if dwi_type not in fdr_global:
                    continue
                fdr_rows = fdr_global[dwi_type]
                if not fdr_rows:
                    continue
                h.append(f"<h3>{_dwi_badge(dwi_type)} \u2014 {len(fdr_rows)} tests</h3>")
                headers = list(fdr_rows[0].keys())[:7]
                h.append("<table><thead><tr>")
                for hdr in headers:
                    h.append(f"<th>{_esc(hdr)}</th>")
                h.append("</tr></thead><tbody>")
                for fr in fdr_rows[:30]:
                    h.append("<tr>")
                    for hdr in headers:
                        val = str(fr.get(hdr, ""))
                        # Highlight p-value columns
                        try:
                            fval = float(val)
                            if "p" in hdr.lower() and 0 < fval < 0.05:
                                cls = _sig_class(fval)
                                h.append(f'<td class="{cls}">{_esc(val[:30])}</td>')
                                continue
                        except (ValueError, TypeError):
                            pass
                        h.append(f"<td>{_esc(val[:40])}</td>")
                    h.append("</tr>")
                h.append("</tbody></table>")
                if len(fdr_rows) > 30:
                    h.append(f"<p><em>\u2026 {len(fdr_rows) - 30} more rows. See FDR_Sig_Global.csv.</em></p>")
    else:
        h.append(_h2("FDR Global Correction Results", "fdr-global"))
        h.append('<p class="meta">No FDR_Sig_Global.csv files found.</p>')

    h.extend(_section_correlations(rows))
    h.extend(_section_treatment_response(groups))
    h.extend(_section_predictive_performance(log_data, dwi_types_present))
    h.extend(_section_mat_data(mat_data))
    h.extend(_section_appendix(rows))

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
    html_report = generate_report(folder)

    out_path = folder / "analysis_report.html"
    out_path.write_text(html_report, encoding="utf-8")
    print(f"Report written to: {out_path}")
    print(f"  Length: {len(html_report)} characters")


if __name__ == "__main__":
    main()
