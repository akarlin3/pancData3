#!/usr/bin/env python3
"""Generate a comprehensive HTML analysis report.

Combines vision-based graph analysis (if available), direct CSV parsing, and
log file parsing into a single structured HTML report.

Usage:
    python generate_report.py [saved_files_path]
"""

from __future__ import annotations

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


def _sig_tag(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def _section(title: str, level: int = 2) -> str:
    return f"\n{'#' * level} {title}\n"


def generate_report(folder: Path) -> str:
    """Build the full Markdown report string."""
    lines: list[str] = []
    timestamp = folder.name.replace("saved_files_", "")

    # Try to import sibling modules for direct parsing
    try:
        from parse_csv_results import parse_all_csvs
        csv_data = parse_all_csvs(folder)
    except Exception:
        csv_data = None

    try:
        from parse_log_metrics import parse_all_logs
        log_data = parse_all_logs(folder)
    except Exception:
        log_data = None

    # Load vision CSV if available
    rows = load_graph_csv(folder)
    groups = group_by_graph_name(rows) if rows else {}

    # Detect which DWI types are present
    dwi_types_present = [d for d in DWI_TYPES if (folder / d).is_dir()]

    # ── Header ──
    lines.append(f"# Analysis Report — {timestamp}")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"**Source folder:** `{folder}`")
    lines.append(f"**DWI types present:** {', '.join(dwi_types_present) or 'None detected'}")
    if rows:
        lines.append(f"**Total graphs analysed:** {len(rows)}")
    lines.append("")

    # ── 1. Executive Summary ──
    lines.append(_section("Executive Summary"))
    lines.append(f"This report summarises the analysis outputs from the pancData3 DWI pipeline run `{timestamp}`.")
    if dwi_types_present:
        lines.append(f"The pipeline processed **{len(dwi_types_present)}** DWI type(s): {', '.join(dwi_types_present)}.")

    # Quick stats from log data
    if log_data:
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            roc = log_data[dwi_type].get("stats_predictive", {}).get("roc_analyses", [])
            best_auc = max((r.get("auc", 0) for r in roc), default=0)
            if best_auc > 0:
                lines.append(f"- **{dwi_type}** best AUC: {best_auc:.3f}")

    lines.append("")

    # ── 2. Graph Analysis Overview ──
    if rows:
        lines.append(_section("Graph Analysis Overview"))

        # Count by type
        type_counts: dict[str, int] = {}
        for r in rows:
            gt = r.get("graph_type", "unknown")
            type_counts[gt] = type_counts.get(gt, 0) + 1

        lines.append("| Graph Type | Count |")
        lines.append("|---|---|")
        for gt, cnt in sorted(type_counts.items(), key=lambda x: -x[1]):
            lines.append(f"| {gt} | {cnt} |")
        lines.append("")

        # Count by DWI type
        dwi_counts: dict[str, int] = {}
        for r in rows:
            dwi_type, _ = parse_dwi_info(r["file_path"])
            dwi_counts[dwi_type] = dwi_counts.get(dwi_type, 0) + 1

        lines.append("| DWI Type | Graphs |")
        lines.append("|---|---|")
        for dt in DWI_TYPES + ["Root"]:
            if dt in dwi_counts:
                lines.append(f"| {dt} | {dwi_counts[dt]} |")
        lines.append("")

    # ── 3. Statistical Significance ──
    lines.append(_section("Statistical Significance"))

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
        lines.append("### Vision-Extracted Significant Findings (p < 0.05)\n")
        lines.append("| p-value | Sig | DWI | Graph | Context |")
        lines.append("|---|---|---|---|---|")
        for p, dwi, graph, ctx in sig_findings[:30]:  # cap at 30 rows
            ctx_clean = ctx.replace("|", "\\|").replace("\n", " ")[:100]
            lines.append(f"| {p:.4f} | {_sig_tag(p)} | {dwi} | {graph} | {ctx_clean} |")
        if len(sig_findings) > 30:
            lines.append(f"\n*... and {len(sig_findings) - 30} more significant findings.*\n")
        lines.append("")

    # From direct CSV parsing
    if csv_data and csv_data.get("significant_metrics"):
        lines.append("### Pipeline CSV Significant Metrics\n")
        for dwi_type in DWI_TYPES:
            if dwi_type not in csv_data["significant_metrics"]:
                continue
            csv_rows = csv_data["significant_metrics"][dwi_type]
            lines.append(f"**{dwi_type}** — {len(csv_rows)} significant metric(s)\n")
            if csv_rows:
                headers = list(csv_rows[0].keys())[:6]
                lines.append("| " + " | ".join(headers) + " |")
                lines.append("| " + " | ".join("---" for _ in headers) + " |")
                for cr in csv_rows[:20]:
                    vals = [str(cr.get(h, "")).replace("|", "\\|")[:30] for h in headers]
                    lines.append("| " + " | ".join(vals) + " |")
            lines.append("")

    if not sig_findings and not (csv_data and csv_data.get("significant_metrics")):
        lines.append("No significant findings extracted.\n")

    # ── 4. Cross-DWI Comparison ──
    if groups:
        lines.append(_section("Cross-DWI Comparison"))

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

            lines.append(f"\n### {base_name}\n")

            # Trend comparison
            all_trends: dict[str, list] = {}
            for dt in DWI_TYPES:
                if dt in dwi_dict:
                    all_trends[dt] = json.loads(dwi_dict[dt]["trends_json"])

            if all_trends:
                all_series: set[str] = set()
                for trends in all_trends.values():
                    for t in trends:
                        all_series.add(t.get("series") or "overall")

                lines.append("| Series | " + " | ".join(DWI_TYPES) + " | Agreement |")
                lines.append("| --- | " + " | ".join("---" for _ in DWI_TYPES) + " | --- |")

                for series in sorted(all_series):
                    directions: dict[str, str] = {}
                    for dt in DWI_TYPES:
                        if dt not in all_trends:
                            continue
                        for t in all_trends[dt]:
                            if (t.get("series") or "overall") == series:
                                directions[dt] = t["direction"]

                    if len(directions) >= 2:
                        vals = list(directions.values())
                        agree = "AGREE" if len(set(vals)) == 1 else "**DIFFER**"
                        row_vals = [directions.get(dt, "-") for dt in DWI_TYPES]
                        lines.append(f"| {series} | " + " | ".join(row_vals) + f" | {agree} |")

                lines.append("")

    # Cross-reference from CSV
    if csv_data and csv_data.get("cross_reference"):
        inconsistent = [c for c in csv_data["cross_reference"] if not c["consistent"]]
        if inconsistent:
            lines.append("### Cross-DWI Significance Inconsistencies\n")
            lines.append("| Metric | Timepoint | Significant In | Not Significant In |")
            lines.append("|---|---|---|---|")
            for c in inconsistent:
                lines.append(
                    f"| {c['metric']} | {c['timepoint']} "
                    f"| {', '.join(c['significant_in'])} "
                    f"| {', '.join(c['not_significant_in'])} |"
                )
            lines.append("")

    # ── 5. Correlations ──
    if rows:
        lines.append(_section("Notable Correlations"))
        corr_findings = []
        for r in rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            all_text = r["summary"] + " " + r["trends_json"]
            for rval, context in extract_correlations(all_text):
                if abs(rval) >= 0.3:
                    corr_findings.append((abs(rval), rval, dwi_type, base_name, context))

        if corr_findings:
            corr_findings.sort(reverse=True)
            lines.append("| |r| | r | Strength | DWI | Graph |")
            lines.append("|---|---|---|---|---|")
            for _, rval, dwi, graph, ctx in corr_findings[:20]:
                strength = "Strong" if abs(rval) >= 0.5 else "Moderate"
                lines.append(f"| {abs(rval):.2f} | {rval:+.2f} | {strength} | {dwi} | {graph} |")
            lines.append("")
        else:
            lines.append("No notable correlations (|r| >= 0.3) found.\n")

    # ── 6. Treatment Response ──
    if groups:
        lines.append(_section("Treatment Response — Longitudinal Trends"))
        for base_name in sorted(groups.keys()):
            if "Longitudinal" not in base_name:
                continue
            dwi_dict = groups[base_name]
            real = [t for t in dwi_dict if t != "Root"]
            if not real:
                continue

            lines.append(f"\n### {base_name}\n")
            for dt in DWI_TYPES:
                if dt not in dwi_dict:
                    continue
                r = dwi_dict[dt]
                summary = r["summary"]
                if len(summary) > 200:
                    summary = summary[:200] + "..."
                lines.append(f"- **{dt}**: {summary}")

                # Inflection points
                ips = json.loads(r["inflection_points_json"])
                for ip in ips:
                    x = ip.get("approximate_x", "?")
                    lines.append(f"  - Inflection at x={x}: {ip['description']}")

            lines.append("")

    # ── 7. Predictive Performance ──
    if log_data:
        lines.append(_section("Predictive Performance"))

        has_roc = False
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            roc = log_data[dwi_type].get("stats_predictive", {}).get("roc_analyses", [])
            if roc:
                if not has_roc:
                    lines.append("| DWI | Timepoint | AUC | Sensitivity | Specificity |")
                    lines.append("|---|---|---|---|---|")
                    has_roc = True
                for r in roc:
                    auc = f"{r['auc']:.3f}" if "auc" in r else "-"
                    sens = f"{r['sensitivity']:.1f}%" if "sensitivity" in r else "-"
                    spec = f"{r['specificity']:.1f}%" if "specificity" in r else "-"
                    lines.append(f"| {dwi_type} | {r['timepoint']} | {auc} | {sens} | {spec} |")

        if has_roc:
            lines.append("")

        # Feature selections
        has_fs = False
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            fs = log_data[dwi_type].get("stats_predictive", {}).get("feature_selections", [])
            if fs:
                if not has_fs:
                    lines.append("\n### Selected Features (Elastic Net)\n")
                    has_fs = True
                lines.append(f"**{dwi_type}:**\n")
                for sel in fs:
                    lines.append(f"- {sel['timepoint']} (lambda={sel['lambda']:.4f}): {', '.join(sel['features'])}")
                lines.append("")

        # Survival / Cox PH
        has_surv = False
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            sv = log_data[dwi_type].get("survival", {})
            hrs = sv.get("hazard_ratios", [])
            if hrs:
                if not has_surv:
                    lines.append("\n### Cox Proportional Hazards\n")
                    has_surv = True
                lines.append(f"**{dwi_type}:**\n")
                lines.append("| Covariate | HR | 95% CI | p |")
                lines.append("|---|---|---|---|")
                for hr in hrs:
                    ci = f"[{hr['ci_lo']:.3f}, {hr['ci_hi']:.3f}]"
                    sig = _sig_tag(hr["p"])
                    lines.append(f"| {hr['covariate']} | {hr['hr']:.3f} | {ci} | {hr['p']:.4f} {sig} |")
                lrt = sv.get("global_lrt")
                if lrt:
                    lines.append(f"\nGlobal LRT: chi2({lrt['df']}) = {lrt['chi2']:.2f}, p = {lrt['p']:.4f}\n")

        if not has_roc and not has_fs and not has_surv:
            lines.append("No predictive performance data found in logs.\n")

    # ── 8. Appendix ──
    if rows:
        lines.append(_section("Appendix: All Graphs"))
        lines.append("| # | DWI | Type | Graph | Summary |")
        lines.append("|---|---|---|---|---|")
        for i, r in enumerate(rows, 1):
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            summary = r["summary"][:80].replace("|", "\\|").replace("\n", " ")
            lines.append(f"| {i} | {dwi_type} | {r['graph_type']} | {base_name} | {summary}... |")
        lines.append("")

    # Footer
    lines.append("---")
    lines.append(f"*Report generated by pancData3 analysis suite on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*")

    return "\n".join(lines)


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
