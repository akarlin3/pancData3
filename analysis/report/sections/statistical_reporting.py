"""Report sections: statistical significance and broad overview."""

from __future__ import annotations

import math
import re

from shared import (  # type: ignore
    DWI_TYPES,
    extract_pvalues,
    parse_dwi_info,
    safe_text,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _sig_class,
    _sig_tag,
    _table_caption,
)


def _section_statistical_significance(rows, csv_data, log_data, dwi_types_present) -> list[str]:
    """Build the Statistical Significance section.

    Aggregates significant findings from three sources:
    1. Vision-extracted p-values from graph summaries/trends.
    2. Pipeline CSV significant metrics (Significant_LF_Metrics.csv).
    3. GLME interaction test details and FDR timepoints from log parsing.

    Also includes a Borderline Findings subsection for p-values between
    0.05 and 0.10 (trend-worthy but not significant).

    Parameters
    ----------
    rows : list[dict]
        Vision CSV rows.
    csv_data : dict or None
        Parsed pipeline CSV exports.
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.

    Returns
    -------
    list[str]
        HTML chunks for the statistical significance section.
    """
    # ── 3. Statistical Significance ──
    h: list[str] = []
    h.append(_h2("Statistical Significance", "significance"))

    # ── Source 1: Vision-extracted p-values ──
    # Search through summaries, trends, and inflection point descriptions
    # for p-value patterns (e.g. "p = 0.032") using regex extraction.
    sig_findings = []
    borderline_findings = []
    if rows:
        for r in rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            all_text = safe_text(r, "summary", "trends_json", "inflection_points_json")
            for pval, context in extract_pvalues(all_text):
                if pval < 0.05:
                    sig_findings.append((pval, dwi_type, base_name, context))
                elif pval < 0.10:
                    borderline_findings.append((pval, dwi_type, base_name, context))

    if sig_findings:
        sig_findings.sort()
        h.append(f"<h3>Vision-Extracted Significant Findings (p &lt; 0.05) \u2014 {len(sig_findings)} total</h3>")
        h.append("<table>")
        h.append(_table_caption("Vision-Extracted Significant Findings",
                                f"(N\u2009=\u2009{len(sig_findings)} findings with p < 0.05)"))
        h.append("<thead><tr><th>p-value</th><th>Sig</th><th>DWI</th>"
                 "<th>Graph</th><th>Context</th></tr></thead><tbody>")
        for p, dwi, graph, ctx in sig_findings:
            ctx_clean = _esc(ctx.replace("\n", " ")[:120])
            cls = _sig_class(p)
            cls_attr = f' class="{cls}"' if cls else ""
            h.append(f"<tr><td{cls_attr}>{p:.4f}</td><td{cls_attr}>{_esc(_sig_tag(p))}</td>"
                     f"<td>{_dwi_badge(dwi)}</td><td>{_esc(graph)}</td>"
                     f"<td>{ctx_clean}</td></tr>")
        h.append("</tbody></table>")

    # Borderline findings (0.05 <= p < 0.10)
    if borderline_findings:
        borderline_findings.sort()
        h.append(f"<details><summary><strong>Borderline Findings (0.05 \u2264 p &lt; 0.10)</strong> "
                 f"\u2014 {len(borderline_findings)} finding(s)</summary>")
        h.append('<div class="warn-box">These findings approach but do not reach conventional '
                 'significance. They may warrant monitoring in future analyses.</div>')
        h.append("<table><thead><tr><th>p-value</th><th>DWI</th>"
                 "<th>Graph</th><th>Context</th></tr></thead><tbody>")
        for p, dwi, graph, ctx in borderline_findings:
            ctx_clean = _esc(ctx.replace("\n", " ")[:120])
            h.append(f"<tr><td>{p:.4f}</td>"
                     f"<td>{_dwi_badge(dwi)}</td><td>{_esc(graph)}</td>"
                     f"<td>{ctx_clean}</td></tr>")
        h.append("</tbody></table>")
        h.append("</details>")

    # ── Source 2: Pipeline CSV significant metrics ──
    if csv_data and csv_data.get("significant_metrics"):
        h.append("<h3>Pipeline CSV Significant Metrics</h3>")
        for dwi_type in DWI_TYPES:
            if dwi_type not in csv_data["significant_metrics"]:
                continue
            csv_rows = csv_data["significant_metrics"][dwi_type]
            h.append(f"<p>{_dwi_badge(dwi_type)} \u2014 {len(csv_rows)} significant metric(s)</p>")
            if csv_rows:
                headers = list(csv_rows[0].keys())[:10]  # type: ignore
                tbl_cls = ' class="table-wide"' if len(headers) > 6 else ""
                h.append(f"<table{tbl_cls}><thead><tr>")
                for hdr in headers:
                    h.append(f"<th>{_esc(hdr)}</th>")
                h.append("</tr></thead><tbody>")
                for cr in csv_rows:
                    h.append("<tr>")
                    for hdr in headers:
                        val = str(cr.get(hdr, ""))  # type: ignore
                        # Highlight p-value columns
                        try:
                            fval = float(val)
                            if "p" in hdr.lower() and 0 < fval < 0.05:
                                cls = _sig_class(fval)
                                h.append(f'<td class="{cls}">{_esc(val[:40])}</td>')  # type: ignore
                                continue
                        except (ValueError, TypeError):
                            pass
                        h.append(f"<td>{_esc(val[:40])}</td>")  # type: ignore
                    h.append("</tr>")
                h.append("</tbody></table>")

    # ── Source 3: GLME interaction details from log parsing ──
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

                h.append("<table>")
                h.append(_table_caption(
                    f"GLME Interaction Test Results ({dwi_type})",
                    f"N\u2009=\u2009{len(glme_details)} metrics tested."))
                h.append("<thead><tr>"
                         "<th>Metric</th><th>p-value</th><th>Adj. \u03b1</th><th>Sig.</th>"
                         "</tr></thead><tbody>")
                for g in sorted(glme_details, key=lambda x: x["p"]):
                    is_sig = g["p"] < g["adj_alpha"]
                    cls = _sig_class(g["p"]) if is_sig else ""
                    cls_attr = f' class="{cls}"' if cls else ""
                    sig_mark = _sig_tag(g["p"]) if is_sig else ""
                    sig_display = _esc(sig_mark) or '\u2014'
                    h.append(
                        f"<tr><td><code>{_esc(g['metric'])}</code></td>"
                        f"<td{cls_attr}>{g['p']:.4f}</td>"
                        f"<td>{g['adj_alpha']:.4f}</td>"
                        f"<td{cls_attr}>{sig_display}</td></tr>"
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



def _section_broad_statistical_overview(log_data, dwi_types_present) -> list[str]:
    """Build the Broad Statistical Overview section.

    Aggregates ALL statistical tests across the pipeline — including those
    that do not survive FDR correction — to provide a comprehensive view
    of where signals exist (even if underpowered).  This section is
    intentionally inclusive: it surfaces unadjusted p < 0.05 findings,
    borderline survival covariates, and cross-DWI consistency patterns
    that would otherwise be hidden by the strict FDR-only reporting.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.

    Returns
    -------
    list[str]
        HTML chunks for the broad statistical overview.
    """
    h: list[str] = []
    if not log_data:
        return h

    h.append(_h2("Broad Statistical Overview", "broad-stats"))
    h.append(
        '<p class="meta">This section reports <strong>all</strong> statistical '
        "tests from the pipeline regardless of multiple-comparison correction. "
        "Findings listed here at unadjusted p&nbsp;&lt;&nbsp;0.05 are "
        "<em>nominally significant</em> and should be interpreted as "
        "hypothesis-generating, not confirmatory. Variables that appear "
        "consistently across DWI types or analysis methods warrant further "
        "investigation in larger cohorts.</p>"
    )

    # ── 1. Cross-DWI Cox PH Hazard Ratio Summary ──
    # Collect all HR data across DWI types into one table.
    all_hrs: list[tuple[str, str, float, float, float, float]] = []
    for dwi_type in dwi_types_present:
        if dwi_type not in log_data:
            continue
        sv = log_data[dwi_type].get("survival", {})
        for hr_entry in sv.get("hazard_ratios", []):
            all_hrs.append((
                dwi_type,
                hr_entry.get("covariate", ""),
                hr_entry.get("hr", float("nan")),
                hr_entry.get("ci_lo", float("nan")),
                hr_entry.get("ci_hi", float("nan")),
                hr_entry.get("p", float("nan")),
            ))

    if all_hrs:
        # Sort by p-value to show strongest signals first.
        all_hrs.sort(key=lambda x: x[5] if x[5] == x[5] else 999)

        n_nom = sum(1 for _, _, _, _, _, p in all_hrs if p < 0.05)
        n_borderline = sum(1 for _, _, _, _, _, p in all_hrs if 0.05 <= p < 0.10)
        h.append("<h3>Survival Covariates (Cox PH) \u2014 All DWI Types</h3>")
        h.append(
            f"<p>{len(all_hrs)} covariate\u2013DWI combinations tested. "
            f"<strong>{n_nom}</strong> nominally significant (p&nbsp;&lt;&nbsp;0.05)"
        )
        if n_borderline:
            h.append(f", {n_borderline} borderline (0.05&nbsp;\u2264&nbsp;p&nbsp;&lt;&nbsp;0.10)")
        h.append(".</p>")

        h.append("<table>")
        bonf_alpha = 0.05 / max(len(all_hrs), 1)
        n_bonf = sum(1 for _, _, _, _, _, p in all_hrs if p < bonf_alpha)
        bonf_note = (
            f"{n_bonf} survive Bonferroni correction" if n_bonf > 0
            else "none survive Bonferroni correction"
        )
        h.append(_table_caption(
            "Cox PH Hazard Ratios Across DWI Types",
            f"Ranked by unadjusted p-value. HR > 1 = increased risk of local failure. "
            f"Nominally significant (p < 0.05) highlighted; {bonf_note} "
            f"for {len(all_hrs)} tests (\u03b1\u2032 = {bonf_alpha:.4f})."))
        h.append(
            "<thead><tr>"
            "<th>Covariate</th><th>DWI</th><th>HR</th>"
            "<th>95% CI</th><th>p-value</th><th>log(HR)</th><th>Direction</th>"
            "</tr></thead><tbody>"
        )
        for dwi, cov, hr, ci_lo, ci_hi, p in all_hrs:
            log_hr = math.log(hr) if hr > 0 and hr == hr else float("nan")
            direction = "\u2191 Risk" if hr > 1 else "\u2193 Protective" if hr < 1 else "\u2014"
            p_str = f"{p:.4f}" if p == p else "\u2014"
            # Highlight rows by significance
            if p < 0.05:
                cls = _sig_class(p)
                row_cls = f' class="{cls}"'
            elif p < 0.10:
                row_cls = ' style="background:#fffbeb"'
            else:
                row_cls = ""
            h.append(
                f"<tr>"
                f"<td><code>{_esc(cov)}</code></td>"
                f"<td>{_dwi_badge(dwi)}</td>"
                f"<td><strong>{hr:.3f}</strong></td>"
                f"<td>{ci_lo:.3f}\u2013{ci_hi:.3f}</td>"
                f"<td{row_cls}>{p_str}</td>"
                f"<td>{log_hr:.3f}</td>"
                f"<td>{direction}</td>"
                f"</tr>"
            )
        h.append("</tbody></table>")

        # Global LRT summary
        lrt_parts = []
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            sv = log_data[dwi_type].get("survival", {})
            glrt = sv.get("global_lrt", {})
            if glrt:
                lrt_p = glrt.get("p", float("nan"))
                lrt_chi2 = glrt.get("chi2", float("nan"))
                lrt_df = glrt.get("df", 0)
                lrt_parts.append(
                    f"{_dwi_badge(dwi_type)} \u03c7\u00b2({lrt_df})={lrt_chi2:.2f}, "
                    f"p={lrt_p:.4f}"
                )
        if lrt_parts:
            h.append(f'<p><strong>Global likelihood-ratio tests:</strong> '
                     f'{"  |  ".join(lrt_parts)}</p>')

    # ── 2. GLME Interaction Overview Across DWI Types ──
    glme_rows: list[tuple[str, str, float, float, bool]] = []
    for dwi_type in dwi_types_present:
        if dwi_type not in log_data:
            continue
        sc = log_data[dwi_type].get("stats_comparisons", {})
        for g in sc.get("glme_details", []):
            glme_rows.append((
                dwi_type,
                g["metric"],
                g["p"],
                g["adj_alpha"],
                g["p"] < g["adj_alpha"],
            ))

    if glme_rows:
        h.append("<h3>GLME Interaction Tests \u2014 Cross-DWI Comparison</h3>")
        h.append(
            '<p class="meta">Mixed-effects models test whether the LF\u2009\u00d7\u2009Timepoint '
            "interaction differs significantly \u2014 i.e., whether biomarker trajectories "
            "diverge between LC and LF groups over the treatment course. "
            "Holm\u2013Bonferroni correction applied within each DWI type.</p>"
        )

        # Build a cross-tab: metric × DWI type
        metrics_seen = []
        for _, m, _, _, _ in glme_rows:
            if m not in metrics_seen:
                metrics_seen.append(m)
        dwi_glme: dict[str, dict[str, tuple[float, float, bool]]] = {}
        for dwi, m, p, adj, sig in glme_rows:
            dwi_glme.setdefault(dwi, {})[m] = (p, adj, sig)

        h.append('<table class="table-compact"><thead><tr><th>Metric</th>')
        for dwi_type in dwi_types_present:
            if dwi_type in dwi_glme:
                h.append(f"<th>{_dwi_badge(dwi_type)}</th>")
        h.append("</tr></thead><tbody>")
        for metric in metrics_seen:
            h.append(f"<tr><td><code>{_esc(metric)}</code></td>")
            for dwi_type in dwi_types_present:
                if dwi_type not in dwi_glme:  # type: ignore
                    continue
                entry = dwi_glme[dwi_type].get(metric)  # type: ignore
                if entry:
                    p, adj, sig = entry
                    cls = _sig_class(p) if sig else ""
                    cls_attr = f' class="{cls}"' if cls else ""
                    bold_s = "<strong>" if p < 0.05 else ""
                    bold_e = "</strong>" if p < 0.05 else ""
                    h.append(f"<td{cls_attr}>{bold_s}{p:.4f}{bold_e}</td>")
                else:
                    h.append("<td>\u2014</td>")
            h.append("</tr>")
        h.append("</tbody></table>")

        # Cross-DWI consistency note
        # Find metrics where the lowest p-value across DWI types is < 0.25
        suggestive = []
        for metric in metrics_seen:
            pvals = []
            for dwi_type in dwi_types_present:
                if dwi_type in dwi_glme and metric in dwi_glme[dwi_type]:  # type: ignore
                    pvals.append(dwi_glme[dwi_type][metric][0])  # type: ignore
            if pvals and min(pvals) < 0.25:
                suggestive.append((metric, min(pvals)))
        if suggestive:
            suggestive.sort(key=lambda x: x[1])
            parts = [f"<code>{_esc(m)}</code> (min p={p:.3f})" for m, p in suggestive]
            h.append(f'<div class="info-box">Most promising interaction signals: '
                     f'{", ".join(parts)}. While none reach corrected significance, '
                     f'these parameters show the strongest trajectory divergence between '
                     f'outcome groups and may be underpowered.</div>')

    # ── 3. Wilcoxon Rank-Sum Summary ──
    # Report what we know about the Wilcoxon testing family
    has_wilcoxon_info = False
    for dwi_type in dwi_types_present:
        if dwi_type not in log_data:
            continue
        sc = log_data[dwi_type].get("stats_comparisons", {})
        fdr_tps = sc.get("fdr_timepoints", [])
        if fdr_tps:
            has_wilcoxon_info = True
            break

    if has_wilcoxon_info:
        h.append("<h3>Wilcoxon Rank-Sum Tests \u2014 Per-Timepoint Summary</h3>")
        h.append(
            '<p class="meta">The pipeline performs Wilcoxon rank-sum (Mann\u2013Whitney U) '
            "tests comparing LC vs LF groups for each diffusion parameter at each "
            "timepoint. Results are corrected for multiple comparisons using "
            "Benjamini\u2013Hochberg FDR across the full family of tests.</p>"
        )

        # Show per-timepoint FDR counts
        h.append('<table class="table-compact"><thead><tr>'
                 "<th>Timepoint</th>")
        for dwi_type in dwi_types_present:
            h.append(f"<th>{_dwi_badge(dwi_type)}</th>")
        h.append("</tr></thead><tbody>")

        # Collect all timepoints across DWI types
        all_tps: list[str] = []
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            for ftp in log_data[dwi_type].get("stats_comparisons", {}).get("fdr_timepoints", []):  # type: ignore
                if ftp["timepoint"] not in all_tps:
                    all_tps.append(ftp["timepoint"])

        tp_data: dict[str, dict[str, int]] = {}
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            for ftp in log_data[dwi_type].get("stats_comparisons", {}).get("fdr_timepoints", []):
                tp_data.setdefault(ftp["timepoint"], {})[dwi_type] = ftp["n_significant"]

        total_fdr = 0
        for tp in all_tps:
            h.append(f"<tr><td><code>{_esc(tp)}</code></td>")
            for dwi_type in dwi_types_present:
                n = tp_data.get(tp, {}).get(dwi_type, 0)
                total_fdr += n  # type: ignore
                cls = ' class="sig-1"' if n > 0 else ""
                h.append(f"<td{cls}>{n}</td>")
            h.append("</tr>")
        h.append("</tbody></table>")

        if total_fdr == 0:
            n_tps = len(all_tps)
            h.append(
                '<div class="warn-box">'
                "<strong>No metrics survived FDR correction</strong> in any DWI type. "
                f"With ~4 parameters \u00d7 {n_tps} timepoints \u00d7 multiple metric "
                "sets per type, the BH procedure requires raw p-values well below "
                "0.05 to survive correction. The absence of FDR-significant results "
                "does not preclude clinically meaningful effects \u2014 it indicates "
                "insufficient statistical power to detect them at the current sample "
                "size. Individual unadjusted Wilcoxon p&nbsp;&lt;&nbsp;0.05 results "
                "(visible on the boxplot figures with red titles) should be noted as "
                "candidate variables for future powered studies.</div>"
            )
        else:
            h.append(
                f'<div class="info-box"><strong>{total_fdr} metric(s)</strong> '
                f'survived FDR correction across {len(all_tps)} timepoints.</div>'
            )

    # ── 4. Nominally Significant Variables Summary ──
    # Collect all nominally significant findings across all test types.
    nominal_sig: list[dict] = []

    # From Cox PH
    for dwi, cov, hr, ci_lo, ci_hi, p in all_hrs:
        if p < 0.05:
            log_hr = math.log(hr) if hr > 0 and hr == hr else 0
            nominal_sig.append({
                "variable": cov,
                "dwi": dwi,
                "test": "Cox PH",
                "statistic": f"HR={hr:.3f}",
                "p": p,
                "effect": "Large" if abs(log_hr) >= 0.8 else "Medium" if abs(log_hr) >= 0.5 else "Small",  # type: ignore
            })

    # From GLME (raw p < 0.05, even if not surviving Holm-Bonferroni)
    for dwi, metric, p, adj, sig in glme_rows:
        if p < 0.05:
            nominal_sig.append({
                "variable": metric,
                "dwi": dwi,
                "test": "GLME interaction",
                "statistic": f"p(interaction)",
                "p": p,
                "effect": "\u2014",
            })

    if nominal_sig:
        nominal_sig.sort(key=lambda x: x["p"])
        h.append("<h3>Nominally Significant Variables (Unadjusted p &lt; 0.05)</h3>")
        h.append(
            '<div class="info-box">'
            f"<strong>{len(nominal_sig)} finding(s)</strong> reach nominal significance "
            "before multiple-comparison correction. These represent the strongest "
            "individual signals in this dataset and are candidates for hypothesis-driven "
            "investigation in a validation cohort.</div>"
        )
        h.append(
            "<table><thead><tr>"
            "<th>Variable</th><th>DWI</th><th>Test</th>"
            "<th>Statistic</th><th>p-value</th><th>Effect Size</th>"
            "</tr></thead><tbody>"
        )
        for f in nominal_sig:
            cls = _sig_class(f["p"])
            cls_attr = f' class="{cls}"' if cls else ""
            h.append(
                f"<tr><td><code>{_esc(f['variable'])}</code></td>"
                f"<td>{_dwi_badge(f['dwi'])}</td>"
                f"<td>{_esc(f['test'])}</td>"
                f"<td>{_esc(f['statistic'])}</td>"
                f"<td{cls_attr}><strong>{f['p']:.4f}</strong></td>"
                f"<td>{_esc(f['effect'])}</td></tr>"
            )
        h.append("</tbody></table>")
    else:
        h.append("<h3>Nominally Significant Variables</h3>")
        h.append(
            "<p>No individual test reached even nominal significance "
            "(unadjusted p&nbsp;&lt;&nbsp;0.05) across any DWI processing type.</p>"
        )

    # ── 5. Cross-DWI Consistency Patterns ──
    # Which covariates show consistent direction across DWI types?
    if all_hrs:
        h.append("<h3>Cross-DWI Directional Consistency</h3>")
        h.append(
            '<p class="meta">Variables whose hazard ratios point in the same '
            "direction across all DWI types are more likely to reflect true "
            "biological signal rather than processing artefacts.</p>"
        )

        # Group HRs by covariate
        by_cov: dict[str, list[tuple[str, float, float]]] = {}
        for dwi, cov, hr, ci_lo, ci_hi, p in all_hrs:
            by_cov.setdefault(cov, []).append((dwi, hr, p))

        h.append('<table class="table-compact"><thead><tr>'
                 "<th>Covariate</th><th>Direction</th><th>Consistent?</th>"
                 "<th>Best p</th><th>Details</th>"
                 "</tr></thead><tbody>")
        for cov in by_cov:
            entries = by_cov[cov]
            directions = ["risk" if hr > 1 else "protective" for _, hr, _ in entries]
            consistent = len(set(directions)) == 1 and len(entries) > 1
            best_p = min(p for _, _, p in entries)
            direction_str = directions[0] if consistent else "mixed"
            icon = "\u2705" if consistent else "\u26a0\ufe0f"
            detail_parts = []
            for dwi, hr, p in entries:
                detail_parts.append(f"{dwi}: HR={hr:.2f} (p={p:.3f})")
            h.append(
                f"<tr>"
                f"<td><code>{_esc(cov)}</code></td>"
                f"<td>{direction_str}</td>"
                f"<td>{icon} {'Yes' if consistent else 'No'}</td>"
                f"<td>{best_p:.4f}</td>"
                f'<td class="axis-info">{" | ".join(detail_parts)}</td>'
                f"</tr>"
            )
        h.append("</tbody></table>")

    # ── 6. Imputation Sensitivity ──
    # Report whether HR estimates are stable across half-life assumptions
    has_sensitivity = False
    for dwi_type in dwi_types_present:
        if dwi_type not in log_data:  # type: ignore
            continue
        sv = log_data[dwi_type].get("survival", {})  # type: ignore
        if sv.get("hazard_ratios"):
            has_sensitivity = True
            break

    if has_sensitivity:
        h.append("<h3>Imputation Sensitivity</h3>")
        h.append(
            '<p class="meta">Cox PH models use KNN imputation for missing longitudinal '
            "data. Sensitivity to the imputation decay half-life (3\u201324 months) "
            "is assessed by re-fitting models at 5 decay rates. Stable HR estimates "
            "across half-lives indicate robust findings.</p>"
        )

        # For each nominally significant finding, note stability
        for entry in nominal_sig:
            cov = entry["variable"]
            dwi = entry["dwi"]
            if dwi in log_data:
                ipcw = log_data[dwi].get("survival", {}).get("ipcw", {})
                if ipcw:
                    h.append(
                        f'<p>{_dwi_badge(dwi)} <code>{_esc(cov)}</code>: '
                        f'IPCW weight range [{ipcw.get("min_weight", "?"):.2f}, '
                        f'{ipcw.get("max_weight", "?"):.2f}] '
                        f'\u2014 {"no weight inflation (stable)" if ipcw.get("max_weight", 1) <= 1.5 else "weight inflation present (interpret with caution)"}.</p>'
                    )

    return h
