"""Report section: predictive performance (ROC/AUC, elastic net, Cox PH)."""

from __future__ import annotations

from report.report_formatters import (  # type: ignore
    _cite,
    _dwi_badge,
    _esc,
    _h2,
    _sig_class,
    _sig_tag,
    _table_caption,
)


def _section_predictive_performance(log_data, dwi_types_present) -> list[str]:
    """Build the Predictive Performance section.

    Displays three sub-sections from parsed log data:
    1. **ROC/AUC Performance** -- per-timepoint AUC, sensitivity,
       specificity, and Youden optimal cutoff.
    2. **Selected Features (Elastic Net)** -- features retained by
       elastic-net regularisation at each timepoint.
    3. **Cox Proportional Hazards** -- hazard ratios, 95% CIs, p-values,
       IPCW weights, and global likelihood-ratio test.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.

    Returns
    -------
    list[str]
        HTML chunks for the predictive performance section.
    """
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
                # AUC discrimination rating (Hosmer–Lemeshow scale)
                disc = "-"
                if isinstance(auc_val, (int, float)):
                    if auc_val >= 0.9:
                        disc = "Outstanding"
                    elif auc_val >= 0.8:
                        disc = "Excellent"
                    elif auc_val >= 0.7:
                        disc = "Acceptable"
                    elif auc_val >= 0.6:
                        disc = "Poor"
                    else:
                        disc = "No discrimination"
                roc_rows_html.append(
                    f"<tr><td>{_dwi_badge(str(dwi_type))}</td>"
                    f"<td>{_esc(str(r.get('timepoint', '-')))}</td>"
                    f"<td><strong>{auc}</strong></td>"
                    f"<td>{disc}</td>"
                    f"<td>{sens}</td><td>{spec}</td>"
                    f"<td>{youden_str}</td></tr>"
                )

        if has_roc:
            h.append("<h3>ROC / AUC Performance</h3>")
            h.append("<table>")
            h.append(_table_caption(
                "Receiver Operating Characteristic Performance",
                "AUC from nested LOOCV with patient-stratified folds."))
            h.append("<thead><tr><th>DWI</th><th>Timepoint</th><th>AUC</th>"
                     "<th>Discrimination</th>"
                     "<th>Sensitivity</th><th>Specificity</th><th>Youden Cutoff</th>"
                     "</tr></thead><tbody>")
            h.extend(roc_rows_html)
            h.append("</tbody></table>")

        # Feature selections
        has_fs = False
        all_lambdas: dict[str, list] = {}
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            fs = log_data[dwi_type].get("stats_predictive", {}).get("feature_selections", [])
            if fs:
                if not has_fs:
                    h.append("<h3>Selected Features (Elastic Net)</h3>")
                    has_fs = True
                h.append(f"<p>{_dwi_badge(dwi_type)}:</p>")
                h.append("<table><thead><tr><th>Timepoint</th><th>\u03bb</th>"
                         "<th># Features</th><th>Selected Features</th></tr></thead><tbody>")
                for sel in fs:
                    features_html = ", ".join(
                        f"<code>{_esc(f)}</code>" for f in sel["features"]
                    )
                    h.append(f"<tr><td><strong>{_esc(sel['timepoint'])}</strong></td>"
                             f"<td>{sel['lambda']:.4f}</td>"
                             f"<td>{len(sel['features'])}</td>"
                             f"<td>{features_html}</td></tr>")
                    all_lambdas.setdefault(dwi_type, []).append(
                        (sel["timepoint"], sel["lambda"], len(sel["features"]))
                    )
                h.append("</tbody></table>")

        # Firth penalised-likelihood refits
        has_firth = False
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            firth = log_data[dwi_type].get("stats_predictive", {}).get("firth_refits", [])
            if firth:
                if not has_firth:
                    h.append("<h3>Firth Penalised-Likelihood Refits</h3>")
                    h.append(f'<p class="meta">Firth bias reduction{_cite("firth")} was applied '
                             "to address separation or convergence issues in the logistic model.</p>")
                    has_firth = True
                firth_parts = []
                for fr in firth:
                    tp = _esc(fr["timepoint"])
                    nf = fr["n_features"]
                    firth_parts.append(f"<code>{tp}</code> ({nf} features)")
                h.append(f"<p>{_dwi_badge(dwi_type)}: Firth refit successful at "
                         f"{', '.join(firth_parts)}</p>")

        # Lambda trend analysis
        if all_lambdas:
            h.append("<h4>Regularisation Trend</h4>")
            h.append('<p class="meta">Higher \u03bb indicates stronger regularisation '
                     '(fewer features retained). Trend across timepoints reflects '
                     'evolving discriminative power.</p>')
            for dwi_type, entries in all_lambdas.items():
                if len(entries) >= 2:
                    first_lam = entries[0][1]
                    last_lam = entries[-1][1]
                    if last_lam > first_lam * 1.2:
                        trend_desc = "increasing (tighter regularisation at later timepoints)"
                    elif last_lam < first_lam * 0.8:
                        trend_desc = "decreasing (more features viable at later timepoints)"
                    else:
                        trend_desc = "stable across timepoints"
                    h.append(f"<p>{_dwi_badge(dwi_type)}: \u03bb trend is <strong>{trend_desc}</strong> "
                             f"({entries[0][0]}: {first_lam:.4f} \u2192 {entries[-1][0]}: {last_lam:.4f})</p>")

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
                h.append("<table>")
                h.append(_table_caption(
                    f"Cause-Specific Cox Proportional Hazards ({dwi_type})",
                    "HR > 1 indicates increased hazard for local failure."))
                h.append("<thead><tr><th>Covariate</th><th>HR</th>"
                         "<th>95% CI</th><th>Direction</th><th>p</th></tr></thead><tbody>")
                for hr in sorted(hrs, key=lambda x: x.get("p", 1)):
                    ci = f"[{hr.get('ci_lo', 0):.3f}, {hr.get('ci_hi', 0):.3f}]"
                    p_val = hr.get("p", 1.0)
                    sig = _sig_tag(p_val)
                    cls = _sig_class(p_val)
                    cls_attr = f' class="{cls}"' if cls else ""
                    hr_val = hr.get('hr', 1.0)
                    if hr_val > 1:
                        dir_html = '<span class="differ">Risk \u2191</span>'
                    elif hr_val < 1:
                        dir_html = '<span class="agree">Protective \u2193</span>'
                    else:
                        dir_html = 'Neutral'
                    h.append(f"<tr><td><code>{_esc(str(hr.get('covariate', '')))}</code></td>"
                             f"<td>{hr_val:.3f}</td><td>{ci}</td>"
                             f"<td>{dir_html}</td>"
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
