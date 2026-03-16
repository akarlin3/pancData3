"""Report sections: manuscript-ready outputs and predictive performance."""

from __future__ import annotations

import math

from report.report_formatters import (  # type: ignore
    _cite,
    _copy_button,
    _dwi_badge,
    _esc,
    _h2,
    _manuscript_sentence,
    _sig_class,
    _sig_tag,
    _table_caption,
)
from report.sections._helpers import (  # type: ignore
    _compute_all_groups_trend_agreement,
    _extract_dosimetry,
    _extract_longitudinal_trend_consensus,
    _get_cohort_size,
    _scalar_gy,
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



def _section_manuscript_ready_findings(
    log_data, dwi_types_present, csv_data, mat_data, groups
) -> list[str]:
    """Build a 'Key Findings for Manuscript' section with copyable sentences.

    Each sentence is dynamically generated from the actual analysis data
    and formatted for direct insertion into a Results section.  Copy
    buttons allow one-click extraction of individual sentences.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    csv_data : dict or None
        Parsed pipeline CSV exports.
    mat_data : dict
        Parsed MAT file metrics.
    groups : dict
        Grouped vision graph data.

    Returns
    -------
    list[str]
        HTML chunks for the key manuscript findings section.
    """
    h: list[str] = []
    sentences: list[str] = []

    # ── Cohort description sentence ──
    n_patients, n_timepoints, _ = _get_cohort_size(mat_data)
    if n_patients > 0:
        sentences.append(
            f"A total of {n_patients} patients with pancreatic cancer "
            f"were included, with DWI acquired at {n_timepoints} timepoints "
            f"during the course of radiotherapy."
        )

    # ── Baseline exclusion sentence ──
    if log_data:
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            bl = log_data[dt].get("baseline", {})
            exc = bl.get("baseline_exclusion")
            if exc and exc.get("n_excluded", 0) > 0:
                sentences.append(
                    f"Of these, {exc['n_excluded']} of {exc['n_total']} "
                    f"patients were excluded due to missing baseline data, "
                    f"leaving {exc['n_total'] - exc['n_excluded']} patients "
                    f"for baseline analysis."
                )
                break

    # ── GLME significance sentence ──
    glme_excluded = None
    if log_data:
        total_glme_tested = 0
        total_glme_sig = 0
        sig_metric_names: list[str] = []
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            sc = log_data[dt].get("stats_comparisons", {})
            details = sc.get("glme_details", [])
            total_glme_tested += len(details)
            for g in details:
                if g["p"] < g["adj_alpha"]:
                    total_glme_sig += 1  # type: ignore
                    if g["metric"] not in sig_metric_names:
                        sig_metric_names.append(g["metric"])
            if sc.get("glme_excluded") and glme_excluded is None:
                glme_excluded = sc["glme_excluded"]
        if total_glme_tested > 0:
            exc_str = ""
            if glme_excluded:
                n_analysed = glme_excluded["n_total"] - glme_excluded["n_excluded"]  # type: ignore
                exc_str = (
                    f"Among {n_analysed} evaluable patients "
                    f"(excluding {glme_excluded['n_excluded']} with "  # type: ignore
                    f"competing-risk events), "
                )
            if total_glme_sig > 0:
                metric_str = ", ".join(sig_metric_names[:3])  # type: ignore
                sentences.append(
                    f"{exc_str}GLME interaction testing identified "
                    f"{total_glme_sig} "
                    f"of {total_glme_tested} metrics with significant "
                    f"time-by-outcome interaction effects after BH-FDR "
                    f"correction ({metric_str})."
                )
            else:
                sentences.append(
                    f"{exc_str}none of the {total_glme_tested} metrics tested "
                    f"showed significant time-by-outcome interaction "
                    f"effects after BH-FDR correction (all p > adjusted alpha)."
                )

    # ── FDR global sentence ──
    if csv_data and csv_data.get("fdr_global"):
        n_fdr = sum(len(v) for v in csv_data["fdr_global"].values())
        if n_fdr > 0:
            sentences.append(
                f"A total of {n_fdr} metric(s) survived global "
                f"Benjamini-Hochberg FDR correction across the full "
                f"feature set."
            )

    # ── Best AUC sentence ──
    if log_data:
        best_auc = 0.0
        best_type = ""
        best_tp = ""
        best_sens = None
        best_spec = None
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            roc = log_data[dt].get("stats_predictive", {}).get("roc_analyses", [])
            for r_item in roc:
                a = r_item.get("auc", 0)
                if a > best_auc:
                    best_auc = a
                    best_type = dt
                    best_tp = r_item.get("timepoint", "")
                    best_sens = r_item.get("sensitivity")
                    best_spec = r_item.get("specificity")
        if best_auc > 0:  # type: ignore
            parts = [f"AUC = {best_auc:.3f}"]
            if best_sens is not None:
                parts.append(f"sensitivity = {best_sens:.1f}%")
            if best_spec is not None:
                parts.append(f"specificity = {best_spec:.1f}%")
            sentences.append(
                f"The elastic-net logistic regression model achieved "
                f"peak discriminative performance at {best_tp} "
                f"using {best_type} DWI ({', '.join(parts)})."
            )

    # ── Cox PH hazard ratio sentences ──
    if log_data:
        import math as _math
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            hrs = log_data[dt].get("survival", {}).get("hazard_ratios", [])
            sig_hrs = [hr for hr in hrs if hr.get("p", 1) < 0.05]
            for hr_item in sorted(sig_hrs, key=lambda x: x["p"])[:2]:  # type: ignore
                hr_val = hr_item["hr"]
                ci_lo = hr_item.get("ci_lo")
                ci_hi = hr_item.get("ci_hi")
                p_val = hr_item["p"]
                cov = hr_item["covariate"]
                direction = "increased" if hr_val > 1 else "decreased"
                ci_str = ""
                if ci_lo is not None and ci_hi is not None:
                    ci_str = f"; 95% CI [{ci_lo:.2f}, {ci_hi:.2f}]"
                # Effect size classification based on log(HR)
                log_hr = _math.log(hr_val) if hr_val > 0 else 0
                abs_log_hr = abs(log_hr)  # type: ignore
                if abs_log_hr >= 0.8:
                    eff_label = "large effect"
                elif abs_log_hr >= 0.5:
                    eff_label = "medium effect"
                else:
                    eff_label = "small effect"
                sentences.append(
                    f"In cause-specific Cox regression, {cov} was "
                    f"significantly associated with {direction} hazard "
                    f"of local failure (HR = {hr_val:.3f}{ci_str}; "
                    f"p = {p_val:.4f}; {eff_label})."
                )

    # ── Dosimetry sentence ──
    dosi, _ = _extract_dosimetry(mat_data)
    if dosi:
        d95 = _scalar_gy(dosi.get("d95_adc_mean"))
        v50 = _scalar_gy(dosi.get("v50_adc_mean"))
        if d95 is not None:
            v50_pct = (v50 * 100 if v50 is not None and v50 <= 1.0
                       else v50) if v50 is not None else None
            parts = [f"D95 = {d95:.1f} Gy"]
            if v50_pct is not None:
                parts.append(f"V50 = {v50_pct:.0f}%")
            sentences.append(
                f"Dosimetric analysis of ADC-defined resistant "
                f"sub-volumes showed {', '.join(parts)}."
            )

    # ── Longitudinal trend sentence ──
    d_cons, f_cons, _, _ = _extract_longitudinal_trend_consensus(groups)
    if d_cons != "unknown" or f_cons != "unknown":
        d_word = {"increasing": "increased", "decreasing": "decreased"}.get(
            d_cons, "remained stable"
        )
        f_word = {"increasing": "increased", "decreasing": "decreased"}.get(
            f_cons, "remained stable"
        )
        sentences.append(
            f"Longitudinally, true diffusion (D) {d_word} while "
            f"perfusion fraction (f) {f_word} over the course of "
            f"treatment, consistent with "
            + (
                "therapy-induced cellular necrosis and vascular regression."
                if d_cons == "increasing" and f_cons == "decreasing"
                else "a heterogeneous treatment response pattern."
            )
        )

    # ── Cross-DWI agreement sentence ──
    n_agree, n_total_series, pct = _compute_all_groups_trend_agreement(
        groups, dwi_types_present
    )
    if n_total_series > 0:
        sentences.append(
            f"Cross-DWI-type trend agreement was {pct:.0f}% "
            f"({n_agree}/{n_total_series} data series), "
            + (
                "supporting robustness of findings across "
                "acquisition strategies."
                if pct >= 70
                else "suggesting some acquisition-dependent variability "
                "in observed trends."
            )
        )

    if not sentences:
        return h

    h.append(_h2("Key Findings for Manuscript", "manuscript-findings"))
    h.append(
        '<p class="meta">Publication-ready sentences generated from the '
        "analysis data. Click <strong>Copy</strong> to copy individual "
        "sentences for direct insertion into your manuscript's Results "
        "section.</p>"
    )

    for sent in sentences:
        h.append(_manuscript_sentence(sent))

    # Copy-all button for the entire block
    all_text = " ".join(sentences)
    h.append(
        f'<div style="margin-top:0.75rem" id="all-findings" '
        f'data-copy="{_esc(all_text)}">'
        f'{_copy_button("all-findings")} '
        f'<span class="meta">Copy all findings as a paragraph</span>'
        f"</div>"
    )

    return h



def _section_results_draft(
    log_data, dwi_types_present, csv_data, mat_data, groups
) -> list[str]:
    """Build a draft Results section with publication-ready paragraphs.

    Generates a structured Results section draft that follows the standard
    organisation of a biomedical research paper: Cohort, Baseline, Group
    Comparisons, Predictive Modelling, Survival Analysis, Cross-DWI.
    Each paragraph is individually copyable.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    csv_data : dict or None
        Parsed pipeline CSV exports.
    mat_data : dict
        Parsed MAT file metrics.
    groups : dict
        Grouped vision graph data.

    Returns
    -------
    list[str]
        HTML chunks for the draft results section.
    """
    h: list[str] = []
    paragraphs: list[tuple[str, str]] = []  # (subsection_title, paragraph_text)

    # ── 1. Cohort paragraph ──
    n_patients, n_timepoints, _ = _get_cohort_size(mat_data)

    baseline_exc = None
    if log_data:
        for dt in dwi_types_present:
            if dt not in log_data:  # type: ignore
                continue
            bl = log_data[dt].get("baseline", {})  # type: ignore
            if bl.get("baseline_exclusion"):
                baseline_exc = bl["baseline_exclusion"]
                break

    if n_patients > 0:
        parts = [
            f"A total of {n_patients} patients with histologically confirmed "
            f"pancreatic cancer who underwent concurrent chemoradiotherapy were "
            f"included in this analysis. DWI was acquired at {n_timepoints} "
            f"timepoints during the course of treatment."
        ]
        if baseline_exc and baseline_exc.get("n_excluded", 0) > 0:
            n_inc = baseline_exc["n_total"] - baseline_exc["n_excluded"]
            parts.append(
                f" {baseline_exc['n_excluded']} patients were excluded due to "
                f"missing baseline DWI data, leaving {n_inc} patients for "
                f"baseline analysis."
            )
        paragraphs.append(("Patient Cohort", "".join(parts)))

    # ── 2. Group comparisons paragraph ──
    if log_data:
        total_tested = 0
        total_sig = 0
        sig_metrics: list[str] = []
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            sc = log_data[dt].get("stats_comparisons", {})
            details = sc.get("glme_details", [])
            total_tested += len(details)
            for g in details:
                if g["p"] < g["adj_alpha"]:
                    total_sig += 1  # type: ignore
                    if g["metric"] not in sig_metrics:
                        sig_metrics.append(g["metric"])

        if total_tested > 0:
            glme_exc = None
            for dt in dwi_types_present:
                if dt in log_data:
                    glme_exc = log_data[dt].get("stats_comparisons", {}).get("glme_excluded")
                    if glme_exc:
                        break

            parts = []
            if glme_exc:
                parts.append(
                    f"After excluding {glme_exc['n_excluded']} of "  # type: ignore
                    f"{glme_exc['n_total']} patients ({glme_exc['pct']:.1f}%) "  # type: ignore
                    f"with competing-risk events, "
                )
            else:
                parts.append("")

            if total_sig > 0:
                metric_str = ", ".join(sig_metrics[:5])  # type: ignore
                parts.append(
                    f"GLME interaction testing identified {total_sig} of "
                    f"{total_tested} metrics with significant time\u00d7outcome "
                    f"interaction effects after Benjamini\u2013Hochberg FDR "
                    f"correction ({metric_str}). "
                )
            else:
                parts.append(
                    f"None of the {total_tested} metrics tested demonstrated "
                    f"significant time\u00d7outcome interaction effects after "
                    f"Benjamini\u2013Hochberg FDR correction. "
                )

            # Add FDR global if available
            if csv_data and csv_data.get("fdr_global"):
                n_fdr = sum(len(v) for v in csv_data["fdr_global"].values())
                if n_fdr > 0:
                    parts.append(
                        f"Globally, {n_fdr} metric(s) survived FDR correction "
                        f"across the full feature set."
                    )

            paragraphs.append(("Group Comparisons", "".join(parts)))

    # ── 3. Predictive modelling paragraph ──
    if log_data:
        all_roc: list[tuple[str, dict]] = []
        all_features: list[tuple[str, dict]] = []
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            sp = log_data[dt].get("stats_predictive", {})
            for r_item in sp.get("roc_analyses", []):
                all_roc.append((dt, r_item))
            for fs in sp.get("feature_selections", []):
                all_features.append((dt, fs))

        if all_roc:
            # Find best AUC
            best_dt, best_roc = max(all_roc, key=lambda x: x[1].get("auc", 0))
            best_auc = best_roc.get("auc", 0)
            best_tp = best_roc.get("timepoint", "")
            best_sens = best_roc.get("sensitivity")
            best_spec = best_roc.get("specificity")

            parts = [
                f"Elastic-net regularised logistic regression with LOOCV "
                f"achieved peak discriminative performance at {best_tp} "
                f"using {best_dt} DWI (AUC = {best_auc:.3f}"
            ]
            if best_sens is not None and best_spec is not None:
                parts.append(
                    f"; sensitivity = {best_sens:.1f}%, "
                    f"specificity = {best_spec:.1f}%"
                )
            parts.append("). ")

            # Feature counts across timepoints
            if all_features:
                unique_feats: set[str] = set()
                for _, fs in all_features:
                    unique_feats.update(fs.get("features", []))
                parts.append(
                    f"Across all timepoints and DWI types, {len(unique_feats)} "
                    f"unique features were selected by the elastic-net model. "
                )

            # Cross-DWI AUC comparison
            if len(dwi_types_present) >= 2:
                best_per_type: dict[str, float] = {}
                for dt, r_item in all_roc:
                    a = r_item.get("auc", 0)
                    if dt not in best_per_type or a > best_per_type[dt]:
                        best_per_type[dt] = a
                if len(best_per_type) >= 2:
                    auc_strs = [f"{dt}: {a:.3f}" for dt, a in
                                sorted(best_per_type.items(), key=lambda x: -x[1])]
                    parts.append(
                        f"Peak AUC by DWI type: {'; '.join(auc_strs)}."
                    )

            paragraphs.append(("Predictive Modelling", "".join(parts)))

    # ── 4. Survival analysis paragraph ──
    if log_data:
        all_hrs: list[tuple[str, dict]] = []
        global_lrt = None
        ipcw_used = False
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            sv = log_data[dt].get("survival", {})
            for hr_item in sv.get("hazard_ratios", []):
                all_hrs.append((dt, hr_item))
            if sv.get("global_lrt") and global_lrt is None:
                global_lrt = sv["global_lrt"]
            if sv.get("ipcw"):
                ipcw_used = True

        if all_hrs:
            sig_hrs = [(dt, hr) for dt, hr in all_hrs if hr.get("p", 1) < 0.05]
            parts = ["Cause-specific Cox proportional hazards modelling"]
            if ipcw_used:
                parts.append(" with IPCW adjustment for competing risks")
            if global_lrt:
                parts.append(
                    f" yielded a global model fit of "
                    f"\u03c7\u00b2({global_lrt['df']}) = {global_lrt['chi2']:.2f}, "  # type: ignore
                    f"p = {global_lrt['p']:.4f}. "  # type: ignore
                )
            else:
                parts.append(". ")

            if sig_hrs:
                for dt, hr_item in sorted(sig_hrs, key=lambda x: x[1]["p"])[:3]:  # type: ignore
                    hr_val = hr_item["hr"]
                    ci_lo = hr_item.get("ci_lo")
                    ci_hi = hr_item.get("ci_hi")
                    direction = "increased" if hr_val > 1 else "decreased"
                    ci_str = ""
                    if ci_lo is not None and ci_hi is not None:
                        ci_str = f"; 95% CI [{ci_lo:.2f}, {ci_hi:.2f}]"
                    parts.append(
                        f"{hr_item['covariate']} was associated with "
                        f"{direction} hazard of local failure "
                        f"(HR = {hr_val:.3f}{ci_str}; p = {hr_item['p']:.4f}). "
                    )
            else:
                parts.append(
                    f"None of the {len(all_hrs)} covariates tested achieved "
                    f"statistical significance at the 0.05 level. "
                )

            paragraphs.append(("Survival Analysis", "".join(parts)))

    # ── 5. Cross-DWI comparison paragraph ──
    n_agree, n_total_series, pct = _compute_all_groups_trend_agreement(
        groups, dwi_types_present
    )
    if n_total_series > 0:
        robust = "high" if pct >= 80 else "moderate" if pct >= 60 else "limited"
        paragraphs.append(("Cross-DWI Comparison",
            f"Trend agreement across the {len(dwi_types_present)} DWI "
            f"processing strategies (Standard, DnCNN, IVIMnet) was "
            f"{pct:.0f}% ({n_agree}/{n_total_series} data series), "
            f"indicating {robust} robustness of the observed biomarker "
            f"trajectories to the choice of acquisition and "
            f"post-processing method."
        ))

    # ── 6. Dosimetry paragraph ──
    dosi, _ = _extract_dosimetry(mat_data)
    if dosi:
        d95 = _scalar_gy(dosi.get("d95_adc_mean"))
        v50 = _scalar_gy(dosi.get("v50_adc_mean"))
        if d95 is not None:
            v50_pct = (v50 * 100 if v50 is not None and v50 <= 1.0
                       else v50) if v50 is not None else None
            coverage = "adequate" if d95 >= 45.0 else "suboptimal"
            parts = [
                f"Dosimetric analysis of ADC-defined resistant "
                f"sub-volumes demonstrated {coverage} target coverage "
                f"(D95 = {d95:.1f} Gy"
            ]
            if v50_pct is not None:
                parts.append(f"; V50 = {v50_pct:.0f}%")
            parts.append(
                "), suggesting that diffusion-defined resistant "
                "regions may receive insufficient dose with current "
                "treatment plans."
                if d95 < 45.0
                else "), indicating that current treatment plans "
                "provide adequate coverage of diffusion-defined "
                "resistant sub-volumes."
            )
            paragraphs.append(("Dosimetric Analysis", "".join(parts)))

    if not paragraphs:
        return h

    h.append(_h2("Draft Results Section", "results-draft"))
    h.append(
        '<p class="meta">Auto-generated Results section draft following '
        "standard biomedical publication structure. Each subsection is "
        "individually copyable for direct insertion into your manuscript. "
        "Review and adapt the language before submission.</p>"
    )

    # Copy-all button for the entire results section
    all_text = "\n\n".join(
        f"{title}\n\n{text}" for title, text in paragraphs
    )
    h.append(
        f'<div style="margin-bottom:1rem" id="results-draft-all" '
        f'data-copy="{_esc(all_text)}">'
        f'{_copy_button("results-draft-all")} '
        f'<span class="meta">Copy entire Results draft</span>'
        f"</div>"
    )

    for title, text in paragraphs:
        h.append(f"<h3>{_esc(title)}</h3>")
        h.append(_manuscript_sentence(text))

    return h
