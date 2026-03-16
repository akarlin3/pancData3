"""Report sections: effect size analysis and multiple comparisons."""

from __future__ import annotations

from shared import (  # type: ignore
    DWI_TYPES,
)
from report.report_formatters import (  # type: ignore
    _cite,
    _dwi_badge,
    _effect_size_class,
    _effect_size_label,
    _esc,
    _forest_plot_cell,
    _h2,
    _sig_class,
    _sig_tag,
    _stat_card,
    _table_caption,
)


def _section_effect_sizes(log_data, dwi_types_present, csv_data) -> list[str]:
    """Build the Effect Size and Sample Size Reporting section.

    Reports effect sizes alongside p-values to distinguish statistical
    from clinical significance, following best practices for biomedical
    publication.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    csv_data : dict or None
        Parsed pipeline CSV exports.

    Returns
    -------
    list[str]
        HTML chunks for the effect size section.
    """
    h: list[str] = []
    h.append(_h2("Effect Size Analysis", "effect-sizes"))
    h.append(
        f'<p class="meta">Effect sizes provide a measure of practical significance '
        f'independent of sample size{_cite("cohen_d")}. Hazard ratios are reported '
        f'with 95% confidence '
        f'intervals; HR > 1 indicates increased risk of the endpoint. The log(HR) is '
        f'used as an approximate effect size metric (|log(HR)| \u2265 0.5 \u2248 '
        f'medium effect).</p>'
    )

    has_data = False

    if log_data:
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            sv = log_data[dwi_type].get("survival", {})
            hrs = sv.get("hazard_ratios", [])
            if not hrs:
                continue

            has_data = True
            h.append(f"<h3>{_dwi_badge(dwi_type)} Hazard Ratio Effect Sizes</h3>")
            h.append("<table>")
            h.append(_table_caption(
                f"Hazard Ratio Effect Sizes ({dwi_type})",
                f"Effect size classification per Cohen (1988): "
                f"small (|log HR| < 0.5), medium (0.5\u20130.8), large (\u2265 0.8)."))
            h.append("<thead><tr>"
                     "<th>Covariate</th><th>HR</th><th>95% CI</th>"
                     "<th>p-value</th><th>log(HR)</th><th>Effect</th>"
                     "<th>Forest Plot</th>"
                     "</tr></thead><tbody>")

            import math
            for hr in sorted(hrs, key=lambda x: x.get("p", 1)):
                hr_val = hr.get("hr", 1.0)
                ci_lo = hr.get("ci_lo", hr_val)
                ci_hi = hr.get("ci_hi", hr_val)
                p_val = hr.get("p", 1.0)
                log_hr = math.log(hr_val) if hr_val > 0 else 0
                eff_cls = _effect_size_class(log_hr)
                eff_lbl = _effect_size_label(log_hr)
                sig = _sig_tag(p_val)
                p_cls = _sig_class(p_val)
                p_attr = f' class="{p_cls}"' if p_cls else ""

                forest = _forest_plot_cell(hr_val, ci_lo, ci_hi, p_val)

                h.append(
                    f"<tr><td><code>{_esc(str(hr.get('covariate', '')))}</code></td>"
                    f"<td><strong>{hr_val:.3f}</strong></td>"
                    f"<td class=\"ci-text\">[{ci_lo:.3f}, {ci_hi:.3f}]</td>"
                    f"<td{p_attr}>{p_val:.4f} {_esc(sig)}</td>"
                    f"<td>{log_hr:+.3f}</td>"
                    f'<td class="{eff_cls}">{_esc(eff_lbl)}</td>'
                    f"<td>{forest}</td></tr>"
                )
            h.append("</tbody></table>")

            # CI width commentary
            narrow_ci = [hr for hr in hrs if hr.get("ci_hi", 1) - hr.get("ci_lo", 1) < 0.5]
            wide_ci = [hr for hr in hrs if hr.get("ci_hi", 1) - hr.get("ci_lo", 1) >= 1.5]
            if wide_ci:
                h.append(f'<div class="warn-box">{len(wide_ci)} covariate(s) have wide '
                         f'confidence intervals (CI width \u2265 1.5), suggesting imprecise '
                         f'effect estimates likely due to limited sample size. '
                         f'For CIs wider than \u00b10.5 log(HR), consider: '
                         f'(1) increasing cohort size (target n \u2265 50 for primary endpoint), '
                         f'(2) reducing the number of covariates tested, or '
                         f'(3) pre-specifying primary endpoints before analysis.</div>')
            if narrow_ci:
                h.append(f'<div class="info-box">{len(narrow_ci)} covariate(s) have narrow '
                         f'CIs (width < 0.5), indicating precise effect estimation.</div>')

            # ── Competing-risk interpretation note ──
            h.append(
                '<div class="info-box"><strong>Competing-risk interpretation:</strong> '
                "Hazard ratios shown are cause-specific (competing-risk adjusted via IPCW "
                "weighting). Standard Cox HRs without competing-risk adjustment would "
                "typically overestimate risk in the presence of competing events (e.g., "
                "death without local failure).</div>"
            )

        # ── ROC/AUC effect size interpretation ──
        auc_data = []
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            roc = log_data[dwi_type].get("stats_predictive", {}).get("roc_analyses", [])
            for r_item in roc:
                auc_val = r_item.get("auc", 0)
                if auc_val > 0:
                    auc_data.append((dwi_type, r_item.get("timepoint", "?"), auc_val))

        if auc_data:
            has_data = True
            h.append("<h3>Discriminative Performance Interpretation</h3>")
            h.append(
                f'<p class="meta">AUC interpretation{_cite("hosmer_lemeshow")}: '
                f'0.5 = no discrimination, '
                f'0.6\u20130.7 = poor, 0.7\u20130.8 = acceptable, '
                f'0.8\u20130.9 = excellent, > 0.9 = outstanding.</p>'
            )
            h.append("<table>")
            h.append(_table_caption(
                "Discriminative Performance by Timepoint",
                "Classification per Hosmer & Lemeshow (2000)."))
            h.append("<thead><tr><th>DWI</th><th>Timepoint</th><th>AUC</th>"
                     "<th>Discrimination</th></tr></thead><tbody>")
            for dwi, tp, auc_val in sorted(auc_data, key=lambda x: -x[2]):
                if auc_val >= 0.9:
                    disc = "Outstanding"
                    cls = "effect-lg"
                elif auc_val >= 0.8:
                    disc = "Excellent"
                    cls = "effect-lg"
                elif auc_val >= 0.7:
                    disc = "Acceptable"
                    cls = "effect-md"
                elif auc_val >= 0.6:
                    disc = "Poor"
                    cls = "effect-sm"
                else:
                    disc = "No discrimination"
                    cls = "effect-sm"
                h.append(
                    f"<tr><td>{_dwi_badge(dwi)}</td>"
                    f"<td>{_esc(str(tp))}</td>"
                    f"<td><strong>{auc_val:.3f}</strong></td>"
                    f'<td class="{cls}">{_esc(disc)}</td></tr>'
                )
            h.append("</tbody></table>")

            # ── AUC benchmark ──
            best_auc_val = max((a for _, _, a in auc_data), default=0.0)
            bench_note = ""
            if best_auc_val >= 0.80:
                bench_note = (
                    f" The observed peak AUC of {best_auc_val:.3f} <strong>exceeds</strong> "
                    "most published DWI benchmarks and warrants prospective validation."
                )
            elif best_auc_val >= 0.70:
                bench_note = (
                    f" The observed peak AUC of {best_auc_val:.3f} is <strong>within</strong> "
                    "the published DWI benchmark range and is generally considered clinically useful."
                )
            elif best_auc_val > 0:
                bench_note = (
                    f" The observed peak AUC of {best_auc_val:.3f} is <strong>below</strong> "
                    "the typical benchmark range; further optimisation or larger cohorts may be needed."
                )
            h.append(
                '<div class="info-box"><strong>Literature benchmark for DWI biomarker AUC in '
                "pancreatic RT response:</strong> published studies report AUC 0.68\u20130.82 "
                "for ADC-based response prediction. AUC \u2265 0.70 is generally considered "
                "clinically useful. AUC \u2265 0.80 exceeds most published DWI benchmarks and "
                f"warrants prospective validation.{bench_note}</div>"
            )

    if not has_data:
        h.append('<p class="meta">No effect size data available (requires survival or predictive log data).</p>')

    return h



def _section_multiple_comparisons(log_data, dwi_types_present, csv_data) -> list[str]:
    """Build the Multiple Comparisons Correction Summary section.

    Provides a comprehensive overview of raw vs adjusted p-values
    across all metrics and DWI types, showing the impact of FDR
    correction.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    csv_data : dict or None
        Parsed pipeline CSV exports.

    Returns
    -------
    list[str]
        HTML chunks for the multiple comparisons section.
    """
    h: list[str] = []
    h.append(_h2("Multiple Comparison Correction Summary", "mult-comp"))

    # ── FDR methodology disclosure ──
    h.append(
        '<div class="info-box"><strong>FDR Methodology:</strong> '
        "Multiple comparisons corrected using the Benjamini\u2013Hochberg (BH) procedure. "
        "For <em>m</em> tests, the adjusted significance threshold for the <em>i</em>-th "
        "ranked p-value is \u03b1(<em>i</em>) = (<em>i</em>/<em>m</em>) \u00d7 0.05. "
        "This controls the false discovery rate (FDR) \u2014 the expected proportion of "
        "false positives among rejected hypotheses \u2014 at 5%.</div>"
    )

    has_data = False

    if log_data:
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            sc = log_data[dwi_type].get("stats_comparisons", {})
            glme_details = sc.get("glme_details", [])
            if not glme_details:
                continue

            has_data = True
            n_total = len(glme_details)
            n_raw_sig = len([g for g in glme_details if g["p"] < 0.05])
            n_fdr_sig = len([g for g in glme_details if g["p"] < g["adj_alpha"]])

            h.append(f"<h3>{_dwi_badge(dwi_type)} \u2014 BH-FDR Correction Impact</h3>")
            h.append('<div class="stat-grid">')
            h.append(_stat_card("Total Tests", str(n_total)))
            h.append(_stat_card("Raw Significant", f"{n_raw_sig}/{n_total}",
                                f"at \u03b1 = 0.05"))
            h.append(_stat_card("FDR Significant", f"{n_fdr_sig}/{n_total}",
                                f"at adjusted \u03b1"))
            lost = n_raw_sig - n_fdr_sig
            if lost > 0:
                h.append(_stat_card("Lost to Correction", str(lost),
                                    "would be false discoveries"))
            h.append("</div>")

            # ── Expected vs observed rejections ──
            expected_fp = round(0.05 * n_total)
            excess = n_raw_sig - expected_fp
            excess_str = (f"+{excess}" if excess > 0 else str(excess))
            h.append(
                f'<div class="info-box">Under the null hypothesis (no true effects), '
                f"approximately 5% of {n_total} test(s) are expected to be falsely "
                f"significant by chance. "
                f"<strong>Expected false discoveries:</strong> {expected_fp}. "
                f"<strong>Observed pre-FDR significant:</strong> {n_raw_sig}. "
                f"<strong>Excess over expected:</strong> {excess_str} "
                f"({'positive excess suggests genuine signal' if excess > 0 else 'no excess above chance level'}).</div>"
            )

            # ── Alternative correction methods ──
            bonferroni_thresh = 0.05 / n_total if n_total > 0 else 0.05
            n_bonferroni_sig = len([g for g in glme_details if g["p"] < bonferroni_thresh])
            h.append(
                f"<details><summary><strong>Alternative correction methods</strong></summary>"
                f'<div class="methods-box"><ul>'
                f"<li><strong>Bonferroni:</strong> \u03b1&nbsp;=&nbsp;0.05/{n_total}"
                f"&nbsp;=&nbsp;{bonferroni_thresh:.4f} (most conservative; controls "
                f"family-wise error rate). Surviving metrics: "
                f"<strong>{n_bonferroni_sig}/{n_total}</strong>.</li>"
                f"<li><strong>BH-FDR (used here):</strong> Controls the expected proportion "
                f"of false positives among rejected hypotheses. Surviving metrics: "
                f"<strong>{n_fdr_sig}/{n_total}</strong>.</li>"
                f"</ul></div></details>"
            )

            # Detailed table: raw p-value vs adjusted threshold
            h.append("<table>")
            h.append(_table_caption(
                f"BH-FDR Correction Detail ({dwi_type})",
                f"Raw vs adjusted significance for {n_total} metrics."))
            h.append("<thead><tr>"
                     "<th>Metric</th><th>Raw p-value</th>"
                     "<th>BH Adj. \u03b1</th><th>Raw Sig?</th>"
                     "<th>FDR Sig?</th><th>Status</th>"
                     "</tr></thead><tbody>")
            for g in sorted(glme_details, key=lambda x: x["p"]):
                raw_sig = g["p"] < 0.05
                fdr_sig = g["p"] < g["adj_alpha"]
                raw_mark = "\u2713" if raw_sig else "\u2717"
                fdr_mark = "\u2713" if fdr_sig else "\u2717"

                if fdr_sig:
                    status = "Confirmed"
                    status_cls = "agree"
                elif raw_sig and not fdr_sig:
                    status = "Rejected by FDR"
                    status_cls = "differ"
                else:
                    status = "Not significant"
                    status_cls = ""

                p_cls = _sig_class(g["p"]) if fdr_sig else ""
                p_attr = f' class="{p_cls}"' if p_cls else ""
                s_attr = f' class="{status_cls}"' if status_cls else ""
                h.append(
                    f"<tr><td><code>{_esc(g['metric'])}</code></td>"
                    f"<td{p_attr}>{g['p']:.4f}</td>"
                    f"<td>{g['adj_alpha']:.4f}</td>"
                    f"<td>{raw_mark}</td><td>{fdr_mark}</td>"
                    f"<td{s_attr}><strong>{_esc(status)}</strong></td></tr>"
                )
            h.append("</tbody></table>")

    # Also show FDR global from CSV
    if csv_data and csv_data.get("fdr_global"):
        fdr_global = csv_data["fdr_global"]
        n_total_fdr = sum(len(v) for v in fdr_global.values())
        if n_total_fdr > 0:
            has_data = True
            h.append("<h3>Global FDR-Surviving Metrics</h3>")
            h.append(f"<p><strong>{n_total_fdr}</strong> metric(s) survived global "
                     f"Benjamini\u2013Hochberg correction across the full feature "
                     f"set, distributed as follows:</p>")
            h.append("<ul>")
            for dt in DWI_TYPES:
                if dt in fdr_global and fdr_global[dt]:
                    h.append(f"<li>{_dwi_badge(dt)}: {len(fdr_global[dt])} metric(s)</li>")
            h.append("</ul>")

    if not has_data:
        h.append('<p class="meta">No multiple comparison data available.</p>')

    return h
