"""Report section: draft Results section for manuscript preparation."""

from __future__ import annotations

from report.report_formatters import (  # type: ignore
    _copy_button,
    _dwi_badge,
    _esc,
    _h2,
    _manuscript_sentence,
)
from report.sections._helpers import (  # type: ignore
    _compute_all_groups_trend_agreement,
    _extract_dosimetry,
    _get_cohort_size,
    _scalar_gy,
)


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
