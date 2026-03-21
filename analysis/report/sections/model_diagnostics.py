"""Report sections: model diagnostics and sensitivity analysis."""

from __future__ import annotations

from shared import (  # type: ignore
    DWI_TYPES,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


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

    # ── Model Calibration ──
    if log_data:
        cal_found = False
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            sp = log_data[dwi_type].get("stats_predictive", {})
            cal = sp.get("calibration")
            if not cal:
                continue

            if not cal_found:
                h.append("<h3>Model Calibration</h3>")
                h.append(
                    '<p class="meta">Calibration assesses whether predicted probabilities '
                    'match observed event rates. A well-calibrated model has Brier score '
                    '&lt; 0.25 and Hosmer\u2013Lemeshow p &gt; 0.05.</p>'
                )
                cal_found = True
            diagnostics_found = True

            brier = cal.get("brier_score", float("nan"))
            hl_p = cal.get("hl_p", float("nan"))
            cal_slope = cal.get("calibration_slope", float("nan"))

            brier_cls = "agree" if isinstance(brier, (int, float)) and brier < 0.25 else "differ"
            hl_cls = "agree" if isinstance(hl_p, (int, float)) and hl_p > 0.05 else "differ"

            h.append(f"<h4>{_dwi_badge(dwi_type)}</h4>")
            h.append("<table><thead><tr>"
                     "<th>Metric</th><th>Value</th><th>Interpretation</th>"
                     "</tr></thead><tbody>")

            brier_interp = "Useful discrimination" if brier < 0.25 else "Poor discrimination"
            hl_interp = "Adequate calibration" if hl_p > 0.05 else "Poor calibration"

            h.append(
                f'<tr><td>Brier Score</td><td class="{brier_cls}">'
                f'<strong>{brier:.4f}</strong></td>'
                f'<td>{brier_interp}</td></tr>'
            )
            h.append(
                f'<tr><td>Hosmer\u2013Lemeshow p</td><td class="{hl_cls}">'
                f'<strong>{hl_p:.4f}</strong></td>'
                f'<td>{hl_interp}</td></tr>'
            )
            if isinstance(cal_slope, (int, float)) and not (cal_slope != cal_slope):
                slope_interp = "Perfect" if abs(cal_slope - 1.0) < 0.1 else (
                    "Overfitting" if cal_slope < 0.8 else "Acceptable")
                h.append(
                    f'<tr><td>Calibration Slope</td>'
                    f'<td><strong>{cal_slope:.3f}</strong></td>'
                    f'<td>{slope_interp}</td></tr>'
                )
            h.append("</tbody></table>")

    # ── Decision Curve Analysis ──
    if log_data:
        dca_found = False
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            sp = log_data[dwi_type].get("stats_predictive", {})
            dca = sp.get("decision_curve")
            if not dca:
                continue

            if not dca_found:
                h.append("<h3>Decision Curve Analysis</h3>")
                h.append(
                    '<p class="meta">Decision curve analysis evaluates clinical '
                    'utility by computing net benefit across a range of threshold '
                    'probabilities. A model provides clinical value at thresholds '
                    'where its net benefit exceeds both the treat-all and treat-none '
                    'default strategies.</p>'
                )
                dca_found = True
            diagnostics_found = True

            useful_lo = dca.get("useful_range_lo", float("nan"))
            useful_hi = dca.get("useful_range_hi", float("nan"))
            h.append(f"<h4>{_dwi_badge(dwi_type)}</h4>")
            if isinstance(useful_lo, (int, float)) and isinstance(useful_hi, (int, float)):
                h.append(
                    f'<p>The model provides positive net benefit over default strategies '
                    f'at threshold probabilities between <strong>{useful_lo:.2f}</strong> '
                    f'and <strong>{useful_hi:.2f}</strong>, indicating clinical utility '
                    f'for patients with this range of prior risk estimates.</p>'
                )
            else:
                h.append(
                    '<p>Decision curve analysis data available but threshold range '
                    'could not be extracted from log outputs.</p>'
                )

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



def build_dca_section(saved_files: dict) -> str:
    """Build a standalone Decision Curve Analysis section.

    Renders a net benefit vs threshold probability table and a brief
    clinical interpretation paragraph from parsed pipeline output.

    Parameters
    ----------
    saved_files : dict
        Parsed pipeline output (keyed by DWI type at top level, with
        ``stats_predictive.decision_curve`` nested within).

    Returns
    -------
    str
        HTML string for the DCA section.  Empty string if no DCA data
        is present.
    """
    if not saved_files:
        return ""

    h: list[str] = []
    dca_found = False

    for dwi_type, dwi_data in saved_files.items():
        if not isinstance(dwi_data, dict):
            continue
        sp = dwi_data.get("stats_predictive", {})
        if not isinstance(sp, dict):
            continue
        dca = sp.get("decision_curve")
        if not dca or not isinstance(dca, dict):
            continue

        if not dca_found:
            h.append(_h2("Decision Curve Analysis", "dca"))
            h.append(
                '<p class="meta">Decision curve analysis evaluates the clinical '
                'utility of a predictive model by computing net benefit across a '
                'range of threshold probabilities. A model adds value at thresholds '
                'where its net benefit exceeds both the treat-all and treat-none '
                'default strategies.</p>'
            )
            dca_found = True

        h.append(f"<h3>{_dwi_badge(dwi_type)}</h3>")

        # Net benefit table
        thresholds = dca.get("thresholds", [])
        net_benefits = dca.get("net_benefits", [])
        treat_all = dca.get("treat_all", [])

        if thresholds and net_benefits:
            h.append(
                "<table><thead><tr>"
                "<th>Threshold Probability</th>"
                "<th>Model Net Benefit</th>"
                "<th>Treat All</th>"
                "<th>Advantage</th>"
                "</tr></thead><tbody>"
            )
            for i, thr in enumerate(thresholds):
                nb = net_benefits[i] if i < len(net_benefits) else None
                ta = treat_all[i] if i < len(treat_all) else None
                nb_str = f"{nb:.4f}" if isinstance(nb, (int, float)) else "N/A"
                ta_str = f"{ta:.4f}" if isinstance(ta, (int, float)) else "N/A"
                if isinstance(nb, (int, float)) and isinstance(ta, (int, float)):
                    adv = nb - ta
                    adv_cls = "agree" if adv > 0 else "differ"
                    adv_str = f'<span class="{adv_cls}">{adv:+.4f}</span>'
                else:
                    adv_str = "N/A"
                h.append(
                    f"<tr><td>{_esc(str(thr))}</td>"
                    f"<td><strong>{nb_str}</strong></td>"
                    f"<td>{ta_str}</td>"
                    f"<td>{adv_str}</td></tr>"
                )
            h.append("</tbody></table>")

        # Clinical interpretation
        useful_lo = dca.get("useful_range_lo")
        useful_hi = dca.get("useful_range_hi")
        if (isinstance(useful_lo, (int, float))
                and isinstance(useful_hi, (int, float))):
            h.append(
                f'<p>The model provides positive net benefit over default '
                f'strategies at threshold probabilities between '
                f'<strong>{useful_lo:.2f}</strong> and '
                f'<strong>{useful_hi:.2f}</strong>. This range represents the '
                f'clinical decision space where using the model to guide '
                f'treatment decisions is expected to produce better outcomes '
                f'than either treating all patients or treating none.</p>'
            )
        else:
            h.append(
                '<p class="meta">Decision curve data available but the useful '
                'threshold range could not be determined.</p>'
            )

    if not dca_found:
        return ""
    return "\n".join(h)


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
