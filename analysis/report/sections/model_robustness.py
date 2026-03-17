"""Report section: model robustness (imputation sensitivity & time-varying Cox).

Renders v2.1-dev pipeline outputs:
- Imputation method AUC comparison table
- Time-varying Cox HR summary for PH-violating covariates
"""

from __future__ import annotations

from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _stat_card,
    _table_caption,
)


def _section_model_robustness(log_data, dwi_types_present, mat_data) -> list[str]:
    """Build the Model Robustness section.

    Displays imputation sensitivity AUC comparison and time-varying Cox
    HR summaries extracted from v2.1-dev MATLAB pipeline logs.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics (keyed by DWI type).
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    mat_data : dict
        Parsed MAT file metrics (keyed by DWI type).

    Returns
    -------
    list[str]
        HTML chunks for the model robustness section.
    """
    h: list[str] = []
    has_content = False

    # ── Imputation Sensitivity ──
    imp_html = _build_imputation_sensitivity(log_data, dwi_types_present)
    if imp_html:
        has_content = True

    # ── Time-Varying Cox ──
    tv_html = _build_time_varying_cox(log_data, dwi_types_present)
    if tv_html:
        has_content = True

    if not has_content:
        return h

    h.append(_h2("Model Robustness", "model-robustness"))
    h.append(
        '<p class="meta">Robustness analyses from v2.1-dev pipeline outputs: '
        "imputation sensitivity compares prediction stability across imputation "
        "strategies; time-varying Cox models address proportional hazards "
        "violations.</p>"
    )

    h.extend(imp_html)
    h.extend(tv_html)

    return h


def _build_imputation_sensitivity(log_data, dwi_types_present) -> list[str]:
    """Build imputation sensitivity comparison sub-section."""
    h: list[str] = []
    if not log_data:
        return h

    found = False
    for dwi_type in dwi_types_present:
        if dwi_type not in log_data:
            continue
        sp = log_data[dwi_type].get("stats_predictive", {})
        imp = sp.get("imputation_sensitivity", [])
        if not imp:
            continue

        if not found:
            h.append("<h3>Imputation Sensitivity</h3>")
            h.append(
                '<p class="meta">AUC comparison across four imputation strategies '
                "(KNN, LOCF, Mean, Linear Interpolation). Consistent AUC values "
                "indicate that predictive performance is robust to imputation "
                "choice.</p>"
            )
            found = True

        h.append(f"<h4>{_dwi_badge(dwi_type)}</h4>")

        # Summary stat cards
        aucs = [e["auc"] for e in imp if isinstance(e.get("auc"), (int, float)) and e["auc"] == e["auc"]]
        if aucs:
            best = max(aucs)
            worst = min(aucs)
            spread = best - worst
            best_method = next(e["method"] for e in imp if e.get("auc") == best)

            cards = [
                _stat_card("Best AUC", f"{best:.3f}", best_method),
                _stat_card("AUC Spread", f"{spread:.3f}",
                           "stable" if spread < 0.05 else "variable"),
            ]
            h.append('<div class="stat-grid">')
            h.extend(cards)
            h.append("</div>")

        # Comparison table
        h.append(_table_caption(
            "Imputation Sensitivity",
            f"AUC and imputed value counts for {_esc(dwi_type)}."
        ))
        h.append(
            "<table><thead><tr>"
            "<th>Method</th><th>AUC</th><th>N Imputed</th>"
            "</tr></thead><tbody>"
        )
        for entry in imp:
            auc_val = entry.get("auc")
            if isinstance(auc_val, float) and auc_val == auc_val:
                auc_str = f"{auc_val:.3f}"
            else:
                auc_str = "N/A"
            h.append(
                f"<tr><td>{_esc(entry['method'])}</td>"
                f"<td><strong>{auc_str}</strong></td>"
                f"<td>{entry.get('n_imputed', 0)}</td></tr>"
            )
        h.append("</tbody></table>")

        # Interpretation
        if aucs and spread < 0.05:
            h.append(
                '<div class="info-box">Predictive performance is robust across '
                "imputation methods (AUC spread &lt; 0.05).</div>"
            )
        elif aucs and spread >= 0.10:
            h.append(
                '<div class="warn-box">Substantial AUC variation across imputation '
                f"methods (spread = {spread:.3f}). Results may be sensitive to "
                "imputation choice. Consider reporting results for multiple "
                "strategies.</div>"
            )

    return h


def _build_time_varying_cox(log_data, dwi_types_present) -> list[str]:
    """Build time-varying Cox model sub-section."""
    h: list[str] = []
    if not log_data:
        return h

    found = False
    for dwi_type in dwi_types_present:
        if dwi_type not in log_data:
            continue
        sv = log_data[dwi_type].get("survival", {})
        tv = sv.get("time_varying_cox")
        if not tv:
            continue

        if not found:
            h.append("<h3>Time-Varying Cox Models</h3>")
            h.append(
                '<p class="meta">When Schoenfeld residual tests detect '
                "proportional hazards (PH) violations, extended Cox models with "
                "covariate \u00d7 log(time) interactions are fit to capture "
                "time-varying effects. A significant interaction term indicates "
                "that the hazard ratio changes over time.</p>"
            )
            found = True

        violated = tv.get("violated_covariates", [])
        h.append(f"<h4>{_dwi_badge(dwi_type)} \u2014 PH violations: "
                 f"{', '.join(_esc(v) for v in violated) or 'none'}</h4>")

        if tv.get("stratified_by"):
            h.append(
                f'<p>Stratified Cox model applied (stratified by '
                f'<code>{_esc(tv["stratified_by"])}</code>).</p>'
            )

        models = tv.get("interaction_models", [])
        if models:
            h.append(_table_caption(
                "Time-Varying Cox Interactions",
                f"Extended Cox model results for {_esc(dwi_type)}."
            ))
            h.append(
                "<table><thead><tr>"
                "<th>Covariate</th>"
                "<th>Base Coef</th><th>Base p</th>"
                "<th>Interaction Coef</th><th>Interaction p</th>"
                "<th>Interpretation</th>"
                "</tr></thead><tbody>"
            )
            for im in models:
                base_c = im.get("base_coef")
                base_p = im.get("base_p")
                int_c = im["interaction_coef"]
                int_p = im["interaction_p"]

                base_c_str = f"{base_c:.4f}" if isinstance(base_c, (int, float)) else "N/A"
                base_p_str = f"{base_p:.4f}" if isinstance(base_p, (int, float)) else "N/A"

                # Interpret the direction of the time-varying effect
                if int_p < 0.05:
                    if int_c > 0:
                        interp = "HR increases over time"
                    else:
                        interp = "HR decreases over time"
                    interp_cls = "differ"
                else:
                    interp = "No significant time variation"
                    interp_cls = "agree"

                h.append(
                    f"<tr><td><code>{_esc(im['covariate'])}</code></td>"
                    f"<td>{base_c_str}</td>"
                    f"<td>{base_p_str}</td>"
                    f"<td><strong>{int_c:.4f}</strong></td>"
                    f"<td><strong>{int_p:.4f}</strong></td>"
                    f'<td><span class="{interp_cls}">{interp}</span></td></tr>'
                )
            h.append("</tbody></table>")
        elif violated:
            h.append(
                '<p class="meta">PH violations detected but extended Cox '
                "interaction models were not successfully fit.</p>"
            )

    return h
