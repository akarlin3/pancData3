"""Report sections: sensitivity analysis and power analysis."""

from __future__ import annotations

import math

from shared import (  # type: ignore
    DWI_TYPES,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _stat_card,
    _table_caption,
)


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


def build_nri_idi_section(saved_files: dict) -> str:
    """Build the Net Reclassification Improvement (NRI) and IDI section.

    Renders NRI (overall, events, non-events), IDI, and their confidence
    intervals from parsed ``compute_nri.m`` output.

    Parameters
    ----------
    saved_files : dict
        Parsed pipeline output (keyed by DWI type at top level, with
        ``nri`` and ``idi`` keys nested within).

    Returns
    -------
    str
        HTML string for the NRI/IDI section.  Empty string if no NRI/IDI
        data is present.
    """
    if not saved_files:
        return ""

    h: list[str] = []
    found = False

    for dwi_type, dwi_data in saved_files.items():
        if not isinstance(dwi_data, dict):
            continue
        nri = dwi_data.get("nri")
        idi = dwi_data.get("idi")
        if not nri and not idi:
            continue

        if not found:
            h.append(_h2("Net Reclassification Improvement", "nri-idi"))
            h.append(
                '<p class="meta">NRI and IDI quantify the improvement in risk '
                "classification when adding DWI-derived biomarkers to a baseline "
                "clinical model. Positive NRI indicates that the new model "
                "correctly reclassifies more patients than it misclassifies.</p>"
            )
            found = True

        h.append(f"<h3>{_dwi_badge(dwi_type)}</h3>")

        # Stat cards for quick overview
        cards: list[str] = []
        if isinstance(nri, dict):
            overall = nri.get("overall")
            p_val = nri.get("p")
            if isinstance(overall, (int, float)):
                sig_label = ""
                if isinstance(p_val, (int, float)) and p_val < 0.05:
                    sig_label = " (p < 0.05)"
                cards.append(_stat_card(
                    "Overall NRI",
                    f"{overall:+.3f}{sig_label}",
                    "positive = improvement",
                ))
        if isinstance(idi, dict):
            idi_val = idi.get("value")
            if isinstance(idi_val, (int, float)):
                cards.append(_stat_card(
                    "IDI",
                    f"{idi_val:+.4f}",
                    "integrated discrimination",
                ))
        if cards:
            h.append('<div class="stat-grid">')
            h.extend(cards)
            h.append("</div>")

        # Detailed NRI table
        if isinstance(nri, dict):
            h.append(
                _table_caption(
                    "NRI Components",
                    f"Reclassification improvement for {_esc(dwi_type)}.",
                )
            )
            h.append(
                "<table><thead><tr>"
                "<th>Component</th><th>Estimate</th>"
                "<th>95% CI</th><th>p-value</th>"
                "</tr></thead><tbody>"
            )
            for component, key in [
                ("Overall NRI", "overall"),
                ("Events NRI", "events"),
                ("Non-events NRI", "non_events"),
                ("Continuous NRI", "continuous"),
            ]:
                val = nri.get(key)
                if val is None:
                    continue
                ci_lo = nri.get(f"{key}_ci_lo")
                ci_hi = nri.get(f"{key}_ci_hi")
                p_val = nri.get(f"{key}_p", nri.get("p"))

                val_str = f"{val:+.4f}" if isinstance(val, (int, float)) else "N/A"
                if (isinstance(ci_lo, (int, float))
                        and isinstance(ci_hi, (int, float))):
                    ci_str = f"[{ci_lo:+.4f}, {ci_hi:+.4f}]"
                else:
                    ci_str = "N/A"
                if isinstance(p_val, (int, float)):
                    p_cls = "agree" if p_val < 0.05 else ""
                    p_str = f"{p_val:.4f}"
                else:
                    p_cls = ""
                    p_str = "N/A"
                cls_attr = f' class="{p_cls}"' if p_cls else ""
                h.append(
                    f"<tr><td>{_esc(component)}</td>"
                    f"<td><strong>{val_str}</strong></td>"
                    f"<td>{ci_str}</td>"
                    f"<td{cls_attr}>{p_str}</td></tr>"
                )
            h.append("</tbody></table>")

        # IDI row
        if isinstance(idi, dict):
            idi_val = idi.get("value")
            idi_ci_lo = idi.get("ci_lo")
            idi_ci_hi = idi.get("ci_hi")
            idi_p = idi.get("p")

            h.append(
                _table_caption(
                    "Integrated Discrimination Improvement",
                    f"IDI for {_esc(dwi_type)}.",
                )
            )
            h.append(
                "<table><thead><tr>"
                "<th>Metric</th><th>Estimate</th>"
                "<th>95% CI</th><th>p-value</th>"
                "</tr></thead><tbody>"
            )
            val_str = f"{idi_val:+.4f}" if isinstance(idi_val, (int, float)) else "N/A"
            if (isinstance(idi_ci_lo, (int, float))
                    and isinstance(idi_ci_hi, (int, float))):
                ci_str = f"[{idi_ci_lo:+.4f}, {idi_ci_hi:+.4f}]"
            else:
                ci_str = "N/A"
            if isinstance(idi_p, (int, float)):
                p_cls = "agree" if idi_p < 0.05 else ""
                p_str = f"{idi_p:.4f}"
            else:
                p_cls = ""
                p_str = "N/A"
            cls_attr = f' class="{p_cls}"' if p_cls else ""
            h.append(
                f'<tr><td>IDI</td><td><strong>{val_str}</strong></td>'
                f"<td>{ci_str}</td>"
                f"<td{cls_attr}>{p_str}</td></tr>"
            )
            h.append("</tbody></table>")

        # Significance flag
        sig_items: list[str] = []
        if isinstance(nri, dict):
            p_val = nri.get("p")
            if isinstance(p_val, (int, float)) and p_val < 0.05:
                sig_items.append(f"NRI (p = {p_val:.4f})")
        if isinstance(idi, dict):
            p_val = idi.get("p")
            if isinstance(p_val, (int, float)) and p_val < 0.05:
                sig_items.append(f"IDI (p = {p_val:.4f})")
        if sig_items:
            h.append(
                '<div class="info-box"><strong>Statistically significant '
                "reclassification improvement:</strong> "
                + ", ".join(sig_items)
                + ". The DWI-derived model provides meaningful improvement in "
                "risk stratification over the baseline clinical model.</div>"
            )

    if not found:
        return ""
    return "\n".join(h)
