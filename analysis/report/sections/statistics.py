"""Report sections: statistics group."""

from __future__ import annotations

import re

from shared import (  # type: ignore
    DWI_TYPES,
)
from report_formatters import (  # type: ignore
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



def _section_model_diagnostics(log_data, dwi_types_present, mat_data) -> list[str]:
    """Build the Model Diagnostics section.

    Reports on model assumptions, convergence, and potential issues
    that should be disclosed in a publication.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    mat_data : dict
        Parsed MAT file metrics.

    Returns
    -------
    list[str]
        HTML chunks for the model diagnostics section.
    """
    h: list[str] = []
    h.append(_h2("Model Diagnostics", "model-diag"))

    diagnostics_found = False

    if log_data:
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue

            issues: list[str] = []

            # IPCW weight range check
            sv = log_data[dwi_type].get("survival", {})
            ipcw = sv.get("ipcw")
            if ipcw:
                max_w = ipcw.get("max_weight", 1.0)
                min_w = ipcw.get("min_weight", 1.0)
                if max_w > 5.0:
                    issues.append(
                        f"IPCW maximum weight = {max_w:.2f} (>5.0). Extreme weights "
                        f"may indicate near-violation of the positivity assumption "
                        f"or sparse censoring strata. Consider weight truncation."
                    )
                elif max_w > 2.0:
                    issues.append(
                        f"IPCW weight range [{min_w:.2f}, {max_w:.2f}] is moderately "
                        f"dispersed. Monitor for undue influence of high-weight observations."
                    )
                else:
                    issues.append(
                        f"IPCW weight range [{min_w:.2f}, {max_w:.2f}] is well-behaved, "
                        f"suggesting adequate censoring model specification."
                    )

            # Competing-risk exclusions
            sc = log_data[dwi_type].get("stats_comparisons", {})
            glme_excl = sc.get("glme_excluded")
            if glme_excl:
                pct = glme_excl.get("pct", 0)
                if pct > 20:
                    issues.append(
                        f"High competing-risk exclusion rate: {pct:.1f}% of patients. "
                        f"This substantially reduces effective sample size and may "
                        f"introduce selection bias."
                    )
                elif pct > 10:
                    issues.append(
                        f"Moderate competing-risk exclusion: {pct:.1f}% of patients. "
                        f"Results should be interpreted in the context of the remaining cohort."
                    )

            # Baseline exclusions
            bl = log_data[dwi_type].get("baseline", {})
            baseline_exc = bl.get("baseline_exclusion")
            if baseline_exc:
                n_exc = baseline_exc.get("n_excluded", 0)
                n_tot = baseline_exc.get("n_total", 1)
                exc_pct = 100 * n_exc / n_tot if n_tot > 0 else 0
                if exc_pct > 15:
                    issues.append(
                        f"High baseline missingness: {n_exc}/{n_tot} ({exc_pct:.1f}%) "
                        f"patients excluded for missing baseline data. Assess whether "
                        f"missingness is informative (MCAR vs MAR vs MNAR)."
                    )
                lf_inc = baseline_exc.get("lf_rate_included")
                lf_exc = baseline_exc.get("lf_rate_excluded")
                if lf_inc is not None and lf_exc is not None:
                    diff = abs(lf_inc - lf_exc)
                    if diff > 15:
                        issues.append(
                            f"LF rate differs substantially between included ({lf_inc:.1f}%) "
                            f"and excluded ({lf_exc:.1f}%) patients (\u0394 = {diff:.1f}pp). "
                            f"This suggests potentially informative missingness (MAR/MNAR)."
                        )

            # Outlier removal impact
            total_out = bl.get("total_outliers")
            if total_out:
                out_pct = total_out.get("pct", 0)
                if out_pct > 10:
                    issues.append(
                        f"Outlier removal rate: {out_pct:.1f}%. High removal rates "
                        f"may bias results; consider sensitivity analysis with and "
                        f"without outlier removal."
                    )

            # Lambda stability for elastic net
            fs = log_data[dwi_type].get("stats_predictive", {}).get("feature_selections", [])
            if len(fs) >= 2:
                lambdas = [s["lambda"] for s in fs]
                lam_range = max(lambdas) / min(lambdas) if min(lambdas) > 0 else float("inf")
                if lam_range > 10:
                    issues.append(
                        f"Elastic net \u03bb varies by >{lam_range:.0f}\u00d7 across timepoints "
                        f"({min(lambdas):.4f}\u2013{max(lambdas):.4f}), suggesting "
                        f"substantial variation in signal-to-noise ratio over time."
                    )

            if issues:
                diagnostics_found = True
                h.append(f"<h3>{_dwi_badge(dwi_type)}</h3>")
                for issue in issues:
                    is_warning = any(w in issue.lower() for w in ["high", "extreme", "violation", "substantially", "bias"])
                    box_cls = "warn-box" if is_warning else "diag-box"
                    h.append(f'<div class="{box_cls}">{_esc(issue)}</div>')

    # General assumptions note
    h.append("<h3>Assumptions and Caveats</h3>")
    h.append('<div class="methods-box">')
    h.append("<ul>")
    h.append(
        "<li><strong>Proportional hazards:</strong> Cox models assume that the hazard "
        "ratio is constant over time. Violations (e.g., time-varying effects) are "
        "not explicitly tested here. Schoenfeld residual tests should be performed "
        "for definitive assessment.</li>"
    )
    h.append(
        "<li><strong>\u26a0 Proportional hazards (PH) assumption not formally tested:</strong> "
        "PH assumption was not formally tested in these outputs. Schoenfeld residuals or "
        "log\u2013log plots are recommended to verify PH before interpreting Cox model "
        "hazard ratios as time-constant effects.</li>"
    )
    h.append(
        "<li><strong>Non-parametric tests:</strong> Wilcoxon rank-sum tests make no "
        "distributional assumptions but are less powerful than parametric alternatives "
        "when normality holds. With small sample sizes, this trade-off favours robustness.</li>"
    )
    h.append(
        "<li><strong>Multiple testing:</strong> BH-FDR controls the expected false "
        "discovery proportion but does not control the family-wise error rate. For "
        "confirmatory analyses, Bonferroni or Holm correction may be more appropriate.</li>"
    )
    h.append(
        "<li><strong>LOOCV bias-variance:</strong> Leave-one-out cross-validation "
        "provides nearly unbiased performance estimates but with high variance. "
        "AUC confidence intervals from LOOCV should be interpreted cautiously.</li>"
    )
    h.append(
        "<li><strong>Feature selection stability:</strong> Elastic net feature "
        "selection may vary across folds. Features consistently selected across "
        "multiple timepoints and DWI types carry higher confidence.</li>"
    )
    h.append(
        "<li><strong>Missing data imputation:</strong> Missing data were handled by "
        "k-nearest-neighbor (KNN) imputation with strict temporal-leakage bounds "
        "(training patients only; future timepoints excluded from imputation reference set). "
        "Formal imputation quality assessment (MICE convergence, predictive mean matching) "
        "is not available from current outputs. Visual inspection of imputed distributions "
        "is recommended.</li>"
    )
    h.append(
        "<li><strong>Elastic net \u03bb selection:</strong> Elastic net \u03bb was selected "
        "by 5-fold cross-validation minimizing binomial deviance. The 1-standard-error (1SE) "
        "rule was not applied; minimum-error \u03bb was used, which may slightly overfit "
        "compared to the 1SE alternative.</li>"
    )
    h.append("</ul>")
    h.append("</div>")
    diagnostics_found = True

    return h



def _section_sensitivity_analysis(log_data, dwi_types_present, mat_data) -> list[str]:
    """Build a Sensitivity Analysis section for publication robustness.

    Evaluates the robustness of key findings by examining:
    - Outlier removal impact on sample size
    - Events-per-variable (EPV) ratio for predictive models
    - Competing-risk exclusion impact on GLME analysis
    - Concordance of significant findings across DWI types

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    mat_data : dict
        Parsed MAT file metrics.

    Returns
    -------
    list[str]
        HTML chunks for the sensitivity analysis section.
    """
    h: list[str] = []
    items: list[tuple[str, str, str]] = []  # (title, detail, box_class)

    if not log_data:
        return h

    # 1. Events-per-variable (EPV) check for predictive models
    # EPV < 10 is a classic overfitting risk indicator
    n_patients = 0
    if mat_data:
        for dt in DWI_TYPES:
            if dt in mat_data and "longitudinal" in mat_data[dt]:  # type: ignore
                n = mat_data[dt]["longitudinal"].get("num_patients", 0)  # type: ignore
                if n > n_patients:
                    n_patients = n

    for dt in dwi_types_present:
        if dt not in log_data:
            continue
        fs_list = log_data[dt].get("stats_predictive", {}).get("feature_selections", [])
        bl = log_data[dt].get("baseline", {})
        baseline_exc = bl.get("baseline_exclusion")
        n_eff = n_patients
        if baseline_exc:
            n_eff = baseline_exc.get("n_total", n_patients) - baseline_exc.get("n_excluded", 0)

        for fs in fs_list:
            n_feat = len(fs.get("features", []))
            if n_feat > 0 and n_eff > 0:
                # Assume ~40% event rate for LF (typical in this cohort)
                n_events_approx = max(int(n_eff * 0.4), 1)  # type: ignore
                epv = n_events_approx / n_feat
                if epv < 5:
                    items.append((
                        f"Low EPV ({dt}, {fs.get('timepoint', '?')})",
                        f"Events-per-variable ratio \u2248 {epv:.1f} "
                        f"({n_events_approx} estimated events / {n_feat} features). "
                        f"EPV < 10 increases overfitting risk; EPV < 5 suggests "
                        f"model estimates may be unreliable. LOOCV partially "
                        f"mitigates this, but external validation is essential.",
                        "warn-box"
                    ))
                elif epv < 10:
                    items.append((
                        f"Marginal EPV ({dt}, {fs.get('timepoint', '?')})",
                        f"Events-per-variable ratio \u2248 {epv:.1f}. "
                        f"This is below the recommended minimum of 10 for "
                        f"logistic regression. Coefficient estimates may be "
                        f"unstable across resampling.",
                        "diag-box"
                    ))

    # 2. Cross-DWI concordance of significant GLME metrics
    sig_metrics_by_dwi: dict[str, set[str]] = {}
    for dt in dwi_types_present:
        if dt not in log_data:
            continue
        sc = log_data[dt].get("stats_comparisons", {})
        sig_set = set()
        for g in sc.get("glme_details", []):
            if g["p"] < g["adj_alpha"]:
                sig_set.add(g["metric"])
        if sig_set:
            sig_metrics_by_dwi[dt] = sig_set

    if len(sig_metrics_by_dwi) >= 2:
        all_sig = set()
        for s in sig_metrics_by_dwi.values():
            all_sig.update(s)
        shared_sig = set.intersection(*sig_metrics_by_dwi.values()) if sig_metrics_by_dwi else set()  # type: ignore
        if all_sig:
            n_shared = len(shared_sig)
            n_total = len(all_sig)
            if n_shared == 0:
                items.append((
                    "No cross-DWI GLME concordance",
                    f"Of {n_total} metric(s) reaching FDR-adjusted significance, "
                    f"none are significant across all DWI types. "
                    f"This suggests findings may be processing-specific rather "
                    f"than reflecting robust biological signal.",
                    "warn-box"
                ))
            elif n_shared < n_total:
                items.append((
                    "Partial cross-DWI GLME concordance",
                    f"{n_shared} of {n_total} significant metric(s) are concordant "
                    f"across all DWI types ({', '.join(sorted(shared_sig))}). "
                    f"Concordant metrics have higher confidence for clinical validity.",
                    "info-box"
                ))
            else:
                items.append((
                    "Full cross-DWI GLME concordance",
                    f"All {n_total} significant metric(s) are concordant across "
                    f"DWI types, supporting robust biological signal.",
                    "info-box"
                ))

    # 3. Hazard ratio stability check: flag HRs with very wide CIs
    for dt in dwi_types_present:
        if dt not in log_data:
            continue
        hrs = log_data[dt].get("survival", {}).get("hazard_ratios", [])
        unstable = [hr for hr in hrs
                    if hr.get("ci_hi", 1) / max(hr.get("ci_lo", 1), 0.001) > 10]
        if unstable:
            names = ", ".join(hr.get("covariate", "?") for hr in unstable[:3])  # type: ignore
            items.append((
                f"Unstable HR estimates ({dt})",
                f"{len(unstable)} covariate(s) have CI ratio > 10 "
                f"({names}), indicating very imprecise effect estimates "
                f"likely due to small sample size or near-separation.",
                "warn-box"
            ))

    if not items:
        return h

    h.append(_h2("Sensitivity Analysis", "sensitivity"))
    h.append(
        '<p class="meta">Robustness checks to evaluate the reliability of key '
        'findings. Issues flagged here should be acknowledged as limitations '
        'in publication and addressed in future validation studies.</p>'
    )
    for title, detail, box_cls in items:
        h.append(f'<div class="{box_cls}"><strong>{_esc(title)}:</strong> {detail}</div>')

    return h



def _section_power_analysis(log_data, dwi_types_present, mat_data) -> list[str]:
    """Build the Statistical Power Commentary section.

    Provides context on statistical power based on the observed cohort
    size and effect sizes, helping readers interpret non-significant
    results appropriately.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    mat_data : dict
        Parsed MAT file metrics.

    Returns
    -------
    list[str]
        HTML chunks for the power analysis section.
    """
    import math

    h: list[str] = []

    # Determine cohort size
    n_patients = 0
    if mat_data:
        for dt in DWI_TYPES:
            if dt in mat_data and "longitudinal" in mat_data[dt]:  # type: ignore
                n = mat_data[dt]["longitudinal"].get("num_patients", 0)  # type: ignore
                if n > n_patients:
                    n_patients = n

    if n_patients == 0:
        return h

    h.append(_h2("Statistical Power Context", "power"))
    h.append(
        '<p class="meta">Post-hoc power commentary based on observed cohort '
        'size and effect sizes. These are approximate guidelines, not formal '
        'power calculations (which require pre-specified effect sizes and '
        'significance levels).</p>'
    )

    # Compute minimum detectable effect sizes for common tests
    # Using rule-of-thumb: for two-sided test at alpha=0.05, 80% power,
    # Wilcoxon requires ~n per group; Cohen's d detectable ≈ 2.8/sqrt(n/2)
    n_per_group = n_patients // 2  # approximate split
    if n_per_group > 0:
        min_cohen_d = 2.8 / math.sqrt(n_per_group)
        # Minimum detectable HR (approximate): HR = exp(d * pi/sqrt(3))
        # where d is the standardised effect size on the log-hazard scale
        min_hr = math.exp(min_cohen_d * 0.5)

        h.append('<div class="methods-box">')
        h.append(f"<h3>Approximate Detectable Effect Sizes (n = {n_patients})</h3>")
        h.append(
            f"<p>With approximately {n_per_group} patients per outcome group "
            f"and conventional settings (\u03b1 = 0.05, 80% power):</p>"
        )
        h.append("<ul>")
        h.append(
            f"<li><strong>Group comparison (Wilcoxon):</strong> Minimum "
            f"detectable standardised effect size d \u2248 {min_cohen_d:.2f} "
            f"({'large' if min_cohen_d >= 0.8 else 'medium' if min_cohen_d >= 0.5 else 'small'} "
            f"effect). "
            f"{'The study is underpowered for detecting small to medium effects.' if min_cohen_d > 0.5 else 'The study has adequate power for medium effects.'}"
            f"</li>"
        )
        h.append(
            f"<li><strong>Survival analysis (Cox PH):</strong> Minimum "
            f"detectable hazard ratio HR \u2248 {min_hr:.2f} (or 1/{min_hr:.2f} = "
            f"{1/min_hr:.2f} for protective effects). Hazard ratios closer to 1.0 "
            f"than this threshold cannot be reliably detected.</li>"
        )

        n_tests_approx: int = 0
        if log_data:
            for dt in dwi_types_present:
                if dt in log_data:
                    sc = log_data[dt].get("stats_comparisons", {})
                    n_tests_approx = int(n_tests_approx + len(sc.get("glme_details", [])))  # type: ignore
        if n_tests_approx > 0:
            # After BH-FDR, effective alpha is approximately alpha * k / m
            # where k is rank; for the median test, effective alpha ≈ 0.025
            adj_min_d = 2.8 / math.sqrt(n_per_group) * 1.15  # ~15% penalty
            h.append(
                f"<li><strong>After FDR correction ({n_tests_approx} tests):</strong> "
                f"The effective significance threshold is reduced, requiring "
                f"slightly larger effects (d \u2248 {adj_min_d:.2f}) for detection. "
                f"Non-significant results after FDR correction should not be "
                f"interpreted as evidence of no effect.</li>"
            )
        h.append("</ul>")

        # Classify observed effect sizes
        if log_data:
            detectable_hrs = []
            undetectable_hrs = []
            for dt in dwi_types_present:
                if dt not in log_data:
                    continue
                hrs = log_data[dt].get("survival", {}).get("hazard_ratios", [])
                for hr_item in hrs:
                    hr_val = hr_item.get("hr", 1.0)
                    log_hr = abs(math.log(hr_val)) if hr_val > 0 else 0
                    if log_hr >= math.log(min_hr):
                        detectable_hrs.append(hr_item)
                    else:
                        undetectable_hrs.append(hr_item)

            if detectable_hrs or undetectable_hrs:
                h.append("<h4>Observed vs Detectable Effect Sizes</h4>")
                if detectable_hrs:
                    h.append(
                        f'<div class="info-box">'
                        f'{len(detectable_hrs)} covariate(s) have effect sizes '
                        f'within the study\'s detection range (|log HR| \u2265 '
                        f'{math.log(min_hr):.2f}).</div>'
                    )
                if undetectable_hrs:
                    h.append(
                        f'<div class="warn-box">'
                        f'{len(undetectable_hrs)} covariate(s) have effect sizes '
                        f'below the study\'s detection threshold. Non-significant '
                        f'p-values for these covariates may reflect insufficient '
                        f'power rather than absence of effect.</div>'
                    )

        h.append("</div>")

    return h


