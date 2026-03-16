"""Report section: executive summary."""

from __future__ import annotations

from shared import (  # type: ignore
    extract_correlations,
    safe_text,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _stat_card,
)
from report.sections._helpers import (  # type: ignore
    _aggregate_dwi_statistics,
    _aggregate_sanity_checks,
    _compute_feature_overlap,
    _find_best_auc,
    _get_cohort_size,
)


def _section_executive_summary(log_data, dwi_types_present, rows, csv_data, timestamp, mat_data=None) -> list[str]:
    """Build the Executive Summary section.

    Displays a summary box with DWI type badges, stat cards for graph
    count, best AUC per DWI type, total significant GLME interactions,
    CSV-derived significant metric count, hazard ratio summary,
    FDR-surviving metric count, correlation count, and cohort info.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics (from :func:`parse_log_metrics.parse_all_logs`).
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    rows : list[dict]
        Vision CSV rows (may be empty).
    csv_data : dict or None
        Parsed pipeline CSV exports.
    timestamp : str
        Run timestamp string for display.
    mat_data : dict or None
        Parsed MAT file metrics (optional).

    Returns
    -------
    list[str]
        HTML chunks for the executive summary section.
    """
    # ── 1. Executive Summary ──
    h: list[str] = []
    h.append(_h2("Executive Summary", "exec-summary"))
    h.append('<div class="abstract-box">')
    h.append(f"<h4>Study Overview</h4>")
    h.append(f"<p>Pipeline run <code>{_esc(timestamp)}</code> processed "
             f"<strong>{len(dwi_types_present)}</strong> DWI type(s): "
             f"{', '.join(_dwi_badge(d) for d in dwi_types_present)}.</p>")

    # Stat cards row
    # Build stat cards for the summary grid.
    cards = []
    if rows:
        cards.append(_stat_card("Graphs Analysed", str(len(rows))))

    stats = _aggregate_dwi_statistics(log_data, dwi_types_present)
    total_sig = stats["total_sig"]
    total_hrs = stats["total_hrs"]
    sig_hrs = stats["sig_hrs"]
    for dt_label, auc_val, auc_tp in stats["auc_cards"]:
        if auc_val >= 0.80:
            disc_label = "excellent discrimination"
        elif auc_val >= 0.70:
            disc_label = "acceptable discrimination"
        else:
            disc_label = "limited discrimination"
        tp_str = f"{dt_label} / {auc_tp} \u2014 {disc_label}" if auc_tp else f"{dt_label} \u2014 {disc_label}"
        cards.append(_stat_card(f"Best AUC ({dt_label})", f"{auc_val:.3f}", tp_str))

    if total_sig > 0:
        cards.append(_stat_card("GLME Sig. Interactions", str(total_sig),
                                "significant interactions (direction not parsed from logs; see Effect Sizes section)"))

    if csv_data and csv_data.get("significant_metrics"):
        n_csv_sig = sum(len(v) for v in csv_data["significant_metrics"].values())
        cards.append(_stat_card("CSV Sig. Metrics", str(n_csv_sig), "from pipeline CSVs"))

    if csv_data and csv_data.get("fdr_global"):
        n_fdr = sum(len(v) for v in csv_data["fdr_global"].values())
        if n_fdr > 0:
            # Build per-DWI breakdown for subtitle
            per_dwi_parts = []
            for dt in dwi_types_present:
                dt_fdr = csv_data["fdr_global"].get(dt, [])
                if dt_fdr:
                    per_dwi_parts.append(f"{dt}: {len(dt_fdr)}")
            fdr_subtitle = " | ".join(per_dwi_parts) if per_dwi_parts else "metrics after BH correction"
            cards.append(_stat_card("FDR-Surviving", str(n_fdr), fdr_subtitle))

    if total_hrs > 0:
        cards.append(_stat_card("Hazard Ratios", f"{sig_hrs}/{total_hrs}", "significant covariates"))

    # Count notable correlations from vision data
    if rows:
        n_corr = 0
        for r in rows:
            all_text = safe_text(r, "summary", "trends_json")
            for rval, _ in extract_correlations(all_text):
                if abs(rval) >= 0.3:
                    n_corr += 1  # type: ignore
        if n_corr > 0:
            cards.append(_stat_card("Notable Correlations", str(n_corr), "|r| \u2265 0.3"))

    # Cross-DWI feature overlap
    total_shared, total_all = _compute_feature_overlap(log_data, dwi_types_present)
    if total_all > 0:
        cards.append(_stat_card("Feature Overlap",
                                f"{total_shared}/{total_all}",
                                "shared across DWI types"))

    # Cohort info from MAT data
    n_pat, n_tp, _ = _get_cohort_size(mat_data)
    if n_pat > 0:
        cards.append(_stat_card("Cohort", f"{n_pat} patients", f"{n_tp} timepoints"))

    # Sanity check summary (convergence + alignment)
    san_agg = _aggregate_sanity_checks(log_data, dwi_types_present)
    if san_agg["sanity_types_checked"] > 0:
        if (san_agg["all_converged_count"] == san_agg["sanity_types_checked"]
                and san_agg["total_dim_issues"] == 0):
            cards.append(_stat_card("Data Quality", "All Passed",
                                    "convergence + alignment"))
        else:
            issue_parts = []
            if san_agg["total_conv_flags"] > 0:
                issue_parts.append(f"{san_agg['total_conv_flags']} conv.")
            if san_agg["total_dim_issues"] > 0:
                issue_parts.append(f"{san_agg['total_dim_issues']} align.")
            cards.append(_stat_card("Data Quality",
                                    ", ".join(issue_parts) if issue_parts else "Passed",
                                    "flags across DWI types"))

    if cards:
        h.append('<div class="stat-grid">')
        h.extend(cards)
        h.append("</div>")

    # ── Structured abstract subsections ──
    h.append("<h4>Objective</h4>")
    h.append("<p>To evaluate the predictive value of diffusion-weighted MRI (DWI) "
             "biomarkers\u2014including ADC, IVIM-derived D, f, and D*\u2014for "
             "treatment response assessment in pancreatic cancer patients receiving "
             "radiotherapy, using Standard, DnCNN-denoised, and IVIMnet-processed "
             "DWI acquisition strategies.</p>")

    h.append("<h4>Methods</h4>")
    methods_parts = []
    if n_pat > 0:
        methods_parts.append(
            f"Longitudinal DWI data from {n_pat} patients across "
            f"{n_tp} timepoints were analysed."
        )
    methods_parts.append(
        "Statistical analysis included Wilcoxon rank-sum tests with Benjamini\u2013Hochberg "
        "FDR correction, generalised linear mixed-effects (GLME) interaction models, "
        "elastic-net regularised logistic regression with LOOCV, and cause-specific "
        "Cox proportional hazards models with IPCW for competing risks."
    )
    h.append(f"<p>{' '.join(methods_parts)}</p>")

    h.append("<h4>Key Results</h4>")
    result_bullets = []
    if total_sig > 0:
        result_bullets.append(f"{total_sig} GLME interaction(s) achieved significance after FDR correction")
    if csv_data and csv_data.get("fdr_global"):
        n_fdr = sum(len(v) for v in csv_data["fdr_global"].values())
        if n_fdr > 0:
            result_bullets.append(f"{n_fdr} metric(s) survived Benjamini\u2013Hochberg FDR correction globally")
    if sig_hrs > 0:
        result_bullets.append(f"{sig_hrs} of {total_hrs} Cox PH covariate(s) were statistically significant (p < 0.05)")
    # Best AUC across all types
    best_roc_overall, best_auc_type = _find_best_auc(log_data, dwi_types_present)
    best_overall_auc = best_roc_overall.get("auc", 0.0) if best_roc_overall else 0.0
    if best_overall_auc > 0:
        result_bullets.append(f"Peak discriminative performance: AUC = {best_overall_auc:.3f} ({best_auc_type})")
    # Cross-DWI trend agreement key result (for abstract)
    if log_data and len(dwi_types_present) >= 2:
        # Recompute from features overlap data already gathered above
        if total_all > 0 and total_shared > 0:
            result_bullets.append(
                f"{total_shared} of {total_all} elastic-net-selected features "
                f"({100 * total_shared / total_all:.0f}%) are shared across DWI types"
            )

    # Sanity check key result
    if san_agg["sanity_types_checked"] > 0:
        if (san_agg["all_converged_count"] == san_agg["sanity_types_checked"]
                and san_agg["total_conv_flags"] == 0):
            result_bullets.append("All voxel-level model fits converged successfully across DWI types")
        elif san_agg["total_conv_flags"] > 0:
            result_bullets.append(f"{san_agg['total_conv_flags']} convergence flag(s) raised across DWI types (see Data Completeness)")
    if result_bullets:
        h.append("<ul>")
        for rb in result_bullets:
            h.append(f"<li>{_esc(rb)}</li>")
        h.append("</ul>")
    else:
        h.append("<p>See detailed sections below for full results.</p>")

    h.append("<h4>Conclusions</h4>")
    if best_overall_auc >= 0.75:
        h.append("<p>DWI-derived biomarkers demonstrate promising discriminative ability for "
                 "treatment response prediction. Cross-DWI-type analysis enables assessment "
                 "of robustness across acquisition and post-processing strategies. "
                 "Longitudinal trends and inflection-point analysis provide mechanistic "
                 "insight into the temporal dynamics of tumour response.</p>")
    else:
        h.append("<p>DWI-derived biomarkers show preliminary potential for treatment response "
                 "assessment. Further investigation with larger cohorts is warranted to "
                 "confirm clinical utility across acquisition strategies.</p>")

    h.append("</div>")
    return h
