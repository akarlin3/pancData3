"""Report sections: model diagnostics."""

from __future__ import annotations

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
