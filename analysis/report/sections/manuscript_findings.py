"""Report section: manuscript-ready key findings with copyable sentences."""

from __future__ import annotations

import math

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
    _extract_longitudinal_trend_consensus,
    _get_cohort_size,
    _scalar_gy,
)


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
