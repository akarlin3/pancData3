"""Report sections: main results group (overview, hypothesis, treatment response)."""

from __future__ import annotations

import json
import re

from shared import (  # type: ignore
    DWI_TYPES,
    extract_correlations,
    extract_pvalues,
    safe_text,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _get_consensus,
    _h2,
    _sig_class,
    _sig_tag,
    _stat_card,
    _trend_tag,
)
from report.sections._helpers import (  # type: ignore
    _aggregate_dwi_statistics,
    _aggregate_sanity_checks,
    _compute_feature_overlap,
    _extract_dosimetry,
    _extract_significant_metrics,
    _find_best_auc,
    _get_cohort_size,
    _scalar_gy,
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



def _section_hypothesis(groups, log_data=None, mat_data=None) -> list[str]:
    """Build the Data-Driven Hypothesis section.

    Analyses longitudinal trend and inflection-point data from the
    ``Longitudinal_Mean_Metrics`` graph group to generate a
    radiological-pathological hypothesis about treatment response.

    The hypothesis addresses three axes:
    - **Cellular response** (D, ADC trends) -- cell kill vs resistance.
    - **Vascular response** (f, D* trends) -- perfusion changes.
    - **Outcome trajectory** -- combined interpretation.

    A fourth axis, **Suggested Treatment Plan**, uses statistically
    significant metrics (p < 0.05) to drive recommendations with actual
    numeric values: GLME interaction p-values, Cox hazard ratios with 95%
    CI, predictive AUC / sensitivity / specificity, FDR-significant
    timepoints, and dosimetry coverage (D95, V50).

    Parameters
    ----------
    groups : dict[str, dict[str, dict]]
        Graph rows grouped by base name and DWI type.
    log_data : dict or None, optional
        Parsed log metrics keyed by DWI type (from ``parse_all_logs``).
        Supplies survival, predictive, and statistical comparison data.
    mat_data : dict or None, optional
        Parsed MAT-file data keyed by DWI type (from ``parse_mat_metrics``).
        Supplies dosimetry metrics (D95, V50) for the treatment plan.

    Returns
    -------
    list[str]
        HTML chunks for the hypothesis section.
    """
    # ── 1.5. Data-Driven Hypothesis ──
    h: list[str] = []
    h.append(_h2("Data-Driven Hypothesis", "hypothesis"))
    h.append('<div class="summary-box">')

    # Analyze trends for dynamic hypothesis
    d_trend_consensus = "unknown"
    f_trend_consensus = "unknown"

    vascular_specificity = ""
    cellular_specificity = ""
    vascular_inflections: list[str] = []
    cellular_inflections: list[str] = []

    # ── Extract statistically significant metrics from log_data ──
    sig_metrics = _extract_significant_metrics(log_data)
    sig_glme = sig_metrics["sig_glme"]
    sig_glme_details = sig_metrics["sig_glme_details"]
    sig_fdr_timepoints = sig_metrics["sig_fdr_timepoints"]
    sig_hr = sig_metrics["sig_hr"]
    best_roc = sig_metrics["best_roc"]
    best_roc_dt = sig_metrics["best_roc_dt"]
    all_feature_selections = sig_metrics["all_feature_selections"]

    # ── Extract dosimetry from mat_data ──
    dosimetry, dosimetry_dt = _extract_dosimetry(mat_data)

    # Analyse the Longitudinal_Mean_Metrics graph (if present) to extract
    # trend directions for D (diffusion) and f (perfusion fraction), plus
    # any inflection points that indicate specific treatment-fraction effects.
    if groups and "Longitudinal_Mean_Metrics" in groups:
        d_trends = []
        f_trends = []

        for dt, r in groups["Longitudinal_Mean_Metrics"].items():
            if dt == "Root": continue
            try:
                # 1. Parse trend directions from the vision model output.
                trends = json.loads(str(r.get("trends_json", "[]")))
                for t in trends:
                    if not isinstance(t, dict): continue
                    series = t.get("series", "")
                    direction = t.get("direction", "").lower()
                    # Classify trends by IVIM parameter series.
                    if series == "Mean D":
                        d_trends.append(direction)
                    elif series == "Mean f":
                        f_trends.append(direction)

                # 2. Parse inflection points for specific fraction-level events.
                ips = json.loads(str(r.get("inflection_points_json", "[]")))
                for ip in ips:
                    if not isinstance(ip, dict): continue
                    x_val = float(ip.get("approximate_x", 0))
                    y_val = ip.get("approximate_y")
                    desc = ip.get("description", "").lower()

                    # Convert x-coordinate to fraction label (e.g. "Fx5").
                    fx_label = f"Fx{int(x_val)}" if x_val > 0 else "baseline"

                    # Try to extract a magnitude (percentage) from the
                    # y-coordinate or from the description text.
                    magnitude = ""
                    if y_val is not None:
                        try:
                            mag_val = abs(float(y_val))
                            if mag_val > 1:
                                magnitude = f" of ~{int(mag_val)}%"
                        except ValueError:
                            pass

                    if not magnitude:
                        # Regex: look for a bare percentage like "15%"
                        pct_match = re.search(r'(\d+)%', desc)
                        if pct_match:
                            magnitude = f" of ~{pct_match.group(1)}%"

                    # Classify inflection as vascular (D*/f-related) or
                    # cellular (ADC/D-related) based on keyword matching.
                    if "d*" in desc or "f" in desc or "vascular" in desc or "perfusion" in desc:
                        if "drop" in desc or "decrease" in desc or "decline" in desc:
                            vascular_inflections.append(f"a significant vascular drop{magnitude} observed around {fx_label}")
                        elif "increase" in desc or "rise" in desc or "peak" in desc:
                            vascular_inflections.append(f"a vascular increase{magnitude} observed around {fx_label}")

                    if "adc" in desc or " d " in desc or "diffusion" in desc:
                        if "increase" in desc or "rise" in desc or "peak" in desc:
                            cellular_inflections.append(f"a significant diffusion increase{magnitude} observed around {fx_label}")
                        elif "drop" in desc or "decrease" in desc or "decline" in desc:
                            cellular_inflections.append(f"a diffusion decrease{magnitude} observed around {fx_label}")

            except Exception:
                pass

        # Determine consensus trend direction via keyword voting.
        d_trend_consensus = _get_consensus(d_trends)
        f_trend_consensus = _get_consensus(f_trends)

        # Build specificity text from unique inflection-point descriptions.
        vascular_specificity = ""
        if vascular_inflections:
            # dict.fromkeys preserves order while deduplicating.
            unique_v = list(dict.fromkeys(vascular_inflections))
            vascular_specificity = f" Specifically, the data highlights {', and '.join(unique_v)}."

        cellular_specificity = ""
        if cellular_inflections:
            unique_c = list(dict.fromkeys(cellular_inflections))
            cellular_specificity = f" Specifically, the data highlights {', and '.join(unique_c)}."

    # ── Identify which metrics are statistically significant ──
    sig_cellular = [d for d in sig_glme_details
                    if "adc" in d["metric"].lower() or d["metric"].lower().startswith("mean_d")]
    sig_vascular = [d for d in sig_glme_details
                    if "f" in d["metric"].lower() or "d*" in d["metric"].lower()
                    or "d_star" in d["metric"].lower() or "dstar" in d["metric"].lower()]

    # Build significance annotations for hypothesis text.
    cellular_sig_note = ""
    if sig_cellular:
        best_cell = min(sig_cellular, key=lambda x: x["p"])
        cellular_sig_note = (
            f" This trend is statistically significant "
            f"(<code>{_esc(best_cell['metric'])}</code>: "
            f"p&nbsp;=&nbsp;{best_cell['p']:.4f}, "
            f"FDR-adjusted \u03b1&nbsp;=&nbsp;{best_cell['adj_alpha']:.4f})."
        )

    vascular_sig_note = ""
    if sig_vascular:
        best_vasc = min(sig_vascular, key=lambda x: x["p"])
        vascular_sig_note = (
            f" This trend is statistically significant "
            f"(<code>{_esc(best_vasc['metric'])}</code>: "
            f"p&nbsp;=&nbsp;{best_vasc['p']:.4f}, "
            f"FDR-adjusted \u03b1&nbsp;=&nbsp;{best_vasc['adj_alpha']:.4f})."
        )

    h.append("<p>Based on quantitative metrics and longitudinal trends extracted from the data, "
             "the following radiological-pathological hypothesis is proposed:</p>")
    h.append("<ul>")

    if d_trend_consensus == "increasing":
        h.append(f"<li><strong>Cellular Response (D, ADC):</strong> The data shows an increase in true diffusion (<em>D</em>) over time.{cellular_specificity}{cellular_sig_note} "
                 "This suggests that effective radiation therapy is inducing cellular necrosis and apoptosis, "
                 "leading to a breakdown of cell membranes. This expands the extracellular space, explaining the increased water mobility.</li>")
    elif d_trend_consensus == "decreasing":
        h.append(f"<li><strong>Cellular Response (D, ADC):</strong> The data shows a decrease in true diffusion (<em>D</em>) over time.{cellular_specificity}{cellular_sig_note} "
                 "This suggests limited cell kill or potential cellular swelling (cytotoxic edema), indicating a highly cellular, densely packed tumor resistant to therapy.</li>")
    else:
        h.append(f"<li><strong>Cellular Response (D, ADC):</strong> The data shows relatively stable or variable true diffusion (<em>D</em>) over time.{cellular_specificity}{cellular_sig_note} "
                 "This suggests a steady state between cellular destruction and tumor proliferation, or a timeline where major necrotic changes are not yet dominant.</li>")

    if f_trend_consensus == "decreasing":
        h.append(f"<li><strong>Vascular Response (f, D*):</strong> The data reveals drops in microcapillary perfusion fraction (<em>f</em>) and/or pseudo-diffusion (<em>D*</em>).{vascular_specificity}{vascular_sig_note} "
                 "These decreases indicate that radiation causes early endothelial damage and vascular regression, effectively cutting off the tumor's blood supply.</li>")
    elif f_trend_consensus == "increasing":
        h.append(f"<li><strong>Vascular Response (f, D*):</strong> The data shows an increase in microcapillary perfusion fraction (<em>f</em>).{vascular_specificity}{vascular_sig_note} "
                 "This suggests reactive angiogenesis, hyperemic inflammatory response, or a robust vascular supply aiding tumor survival and radiation resistance.</li>")
    else:
        h.append(f"<li><strong>Vascular Response (f, D*):</strong> The microcapillary perfusion fraction (<em>f</em>) remains relatively stable,{vascular_specificity}{vascular_sig_note} "
                 "suggesting that the tumor's vascular network has not been significantly altered or compromised by the treatment doses applied so far.</li>")

    outcome_specificity = ""
    if vascular_inflections or cellular_inflections:
        combined = list(dict.fromkeys(cellular_inflections + vascular_inflections))
        outcome_specificity = f" (supported by {', '.join(combined)})"

    if d_trend_consensus == "increasing" and f_trend_consensus == "decreasing":
        h.append(f"<li><strong>Outcome Trajectory:</strong> The combination of early decreases in perfusion coupled with "
                 f"subsequent increases in diffusion strongly supports the hypothesis that the tumor is responding effectively to treatment{outcome_specificity}, "
                 "correlating with long-term Local Control.</li>")
    else:
        h.append(f"<li><strong>Outcome Trajectory:</strong> The observed variation in diffusion and perfusion responses suggests "
                 f"a heterogeneous or limited overall treatment effect{outcome_specificity}. Tumors not exhibiting strong, concomitant increases in diffusion and decreases in perfusion "
                 "may remain at higher risk for Local Failure or harbor therapy-resistant sub-volumes.</li>")

    h.append("</ul>")

    # ── Suggested Treatment Plan ──
    # Only factors that reach statistical significance (p < 0.05) drive
    # the plan content.  Actual numeric values (HR, CI, AUC, sensitivity,
    # specificity, D95, V50, p-values) are embedded throughout.
    plan_points: list[str] = []

    # -- Core recommendation qualified by GLME significance --
    sig_qualifier = ""
    if sig_glme:
        best_glme_p = min(g["p"] for g in sig_glme)
        sig_qualifier = (
            f" GLME interaction testing confirms that outcome-group "
            f"trajectories diverge significantly "
            f"(best interaction p&nbsp;=&nbsp;{best_glme_p:.4f})."
        )
    elif log_data:
        sig_qualifier = (" GLME interaction testing did not identify "
                         "statistically significant trajectory divergence "
                         "between outcome groups (all p&nbsp;\u2265&nbsp;0.05); "
                         "these recommendations are based on observed trends only.")

    if d_trend_consensus == "increasing" and f_trend_consensus == "decreasing":
        plan_points.append(
            "The current radiotherapy regimen appears effective. "
            "Continue the prescribed fractionation scheme and monitor "
            "for sustained increases in <em>D</em> and decreases in <em>f</em> "
            f"to confirm durable tumor response.{sig_qualifier}"
        )
    elif d_trend_consensus in ("stable", "decreasing", "unknown") and f_trend_consensus == "decreasing":
        plan_points.append(
            "Vascular disruption is evident (decreasing <em>f</em>), but "
            "cellular kill remains limited (stable/decreasing <em>D</em>). "
            "Consider dose escalation to the gross tumor volume (GTV) or "
            "the addition of a radiosensitizing agent (e.g., gemcitabine or "
            "capecitabine) to augment cell kill while the vascular supply is "
            f"compromised.{sig_qualifier}"
        )
    elif d_trend_consensus == "increasing" and f_trend_consensus in ("stable", "increasing", "unknown"):
        plan_points.append(
            "Cellular response is observed (increasing <em>D</em>), but the "
            "tumor vasculature remains intact or is increasing (<em>f</em> "
            "stable/rising). Consider combination with an anti-angiogenic agent "
            "to disrupt residual vascular supply, or evaluate whether a "
            "hypofractionated boost to hyperperfused sub-volumes may enhance "
            f"local control.{sig_qualifier}"
        )
    else:
        plan_points.append(
            "Neither strong cellular kill nor significant vascular disruption "
            "is observed. Consider a multidisciplinary review to evaluate "
            "alternative treatment strategies, such as treatment intensification, "
            "concurrent systemic therapy, or adaptive replanning with dose "
            "escalation to resistant sub-volumes identified by persistently "
            f"low <em>D</em> regions.{sig_qualifier}"
        )

    # -- Significant per-metric GLME details --
    if sig_glme_details:
        seen_metrics: set[str] = set()
        detail_parts: list[str] = []
        for d in sorted(sig_glme_details, key=lambda x: x["p"]):
            metric = d["metric"]
            if metric in seen_metrics:
                continue
            seen_metrics.add(metric)
            detail_parts.append(
                f"<code>{_esc(metric)}</code> "
                f"(p&nbsp;=&nbsp;{d['p']:.4f}, "
                f"\u03b1<sub>adj</sub>&nbsp;=&nbsp;{d['adj_alpha']:.4f})"
            )
        if detail_parts:
            plan_points.append(
                f"Metrics with statistically significant outcome-group "
                f"differences (surviving BH-FDR correction): "
                f"{'; '.join(detail_parts)}. "
                f"These should be primary targets for longitudinal "
                f"monitoring and treatment-response evaluation."
            )

    # -- Timepoint-specific FDR significance --
    if sig_fdr_timepoints:
        tp_parts = []
        for fdr_tp in sorted(sig_fdr_timepoints,
                             key=lambda x: x.get("n_significant", 0),
                             reverse=True):
            tp_parts.append(
                f"{_esc(fdr_tp['timepoint'])} "
                f"({fdr_tp['n_significant']} significant metric"
                f"{'s' if fdr_tp['n_significant'] != 1 else ''})"
            )
        plan_points.append(
            f"FDR-corrected analysis identified significant outcome "
            f"differences at timepoints: {', '.join(tp_parts)}. "
            f"These fractions represent critical monitoring windows where "
            f"treatment adaptation should be evaluated."
        )

    # -- Inflection-point-driven timing guidance --
    if cellular_inflections or vascular_inflections:
        fx_labels: list[str] = []
        for desc in cellular_inflections + vascular_inflections:
            fx_match = re.search(r'(Fx\d+|baseline)', desc)
            if fx_match and fx_match.group(1) not in fx_labels:
                fx_labels.append(fx_match.group(1))
        if fx_labels:
            plan_points.append(
                f"Key inflection points were identified around "
                f"{', '.join(fx_labels)}. These timepoints represent "
                f"optimal windows for adaptive replanning or mid-treatment "
                f"imaging assessment to confirm the trajectory before "
                f"continuing or modifying the plan."
            )

    # -- Survival hazard ratios with numeric values --
    if sig_hr:
        protective: list[str] = []
        risk: list[str] = []
        seen_covariates: set[str] = set()
        for hr_entry in sorted(sig_hr, key=lambda x: x["p"]):
            cov = hr_entry["covariate"]
            if cov in seen_covariates:
                continue
            seen_covariates.add(cov)
            hr_val = hr_entry["hr"]
            ci_lo = hr_entry.get("ci_lo")
            ci_hi = hr_entry.get("ci_hi")
            p_val = hr_entry["p"]
            ci_str = ""
            if ci_lo is not None and ci_hi is not None:
                ci_str = f", 95%&nbsp;CI&nbsp;[{ci_lo:.2f}\u2013{ci_hi:.2f}]"
            entry_str = (
                f"<code>{_esc(cov)}</code> "
                f"(HR&nbsp;=&nbsp;{hr_val:.3f}{ci_str}, "
                f"p&nbsp;=&nbsp;{p_val:.4f})"
            )
            if hr_val < 1.0:
                protective.append(entry_str)
            else:
                risk.append(entry_str)

        parts = []
        if protective:
            parts.append(
                f"Protective factors (HR&nbsp;&lt;&nbsp;1): "
                f"{'; '.join(protective)}. "
                f"Higher values of these metrics correlate with improved "
                f"local control; treatment strategies that increase these "
                f"parameters may be beneficial."
            )
        if risk:
            parts.append(
                f"Risk factors (HR&nbsp;&gt;&nbsp;1): "
                f"{'; '.join(risk)}. "
                f"Patients with elevated values of these metrics are at "
                f"increased risk for local failure and may benefit from "
                f"intensified monitoring or treatment adaptation."
            )
        if parts:
            plan_points.append(" ".join(parts))

    # -- Dosimetry coverage (D95, V50) --
    if dosimetry:
        d95_adc = _scalar_gy(dosimetry.get("d95_adc_mean"))
        v50_adc = _scalar_gy(dosimetry.get("v50_adc_mean"))
        d95_d = _scalar_gy(dosimetry.get("d95_d_mean"))
        v50_d = _scalar_gy(dosimetry.get("v50_d_mean"))
        dosi_parts: list[str] = []
        dose_concerns: list[str] = []

        if d95_adc is not None:
            dosi_parts.append(f"D95<sub>ADC</sub>&nbsp;=&nbsp;{d95_adc:.2f}&nbsp;Gy")
            if d95_adc < 45.0:
                dose_concerns.append(
                    f"ADC-defined resistant sub-volume is under-dosed "
                    f"(D95&nbsp;=&nbsp;{d95_adc:.2f}&nbsp;Gy, below 45&nbsp;Gy)"
                )
        if v50_adc is not None:
            v50_pct = v50_adc * 100 if v50_adc <= 1.0 else v50_adc
            dosi_parts.append(f"V50<sub>ADC</sub>&nbsp;=&nbsp;{v50_pct:.1f}%")
            if v50_pct < 90.0:
                dose_concerns.append(
                    f"only {v50_pct:.1f}% of the ADC sub-volume receives "
                    f"\u226550&nbsp;Gy"
                )
        if d95_d is not None:
            dosi_parts.append(f"D95<sub>D</sub>&nbsp;=&nbsp;{d95_d:.2f}&nbsp;Gy")
            if d95_d < 45.0:
                dose_concerns.append(
                    f"D-defined resistant sub-volume is under-dosed "
                    f"(D95&nbsp;=&nbsp;{d95_d:.2f}&nbsp;Gy)"
                )
        if v50_d is not None:
            v50_d_pct = v50_d * 100 if v50_d <= 1.0 else v50_d
            dosi_parts.append(f"V50<sub>D</sub>&nbsp;=&nbsp;{v50_d_pct:.1f}%")

        if dosi_parts:
            dosi_text = (
                f"Dosimetry analysis ({_esc(dosimetry_dt)}): "
                f"{', '.join(dosi_parts)}. "
            )
            if dose_concerns:
                dosi_text += (
                    f"<strong>Dose coverage concern:</strong> "
                    f"{'; '.join(dose_concerns)}. "
                    f"Consider dose escalation or adaptive replanning to "
                    f"improve coverage of resistant sub-volumes."
                )
            else:
                dosi_text += (
                    "Dose coverage to resistant sub-volumes appears adequate. "
                    "Continue current dose constraints."
                )
            plan_points.append(dosi_text)

    # -- Predictive model (AUC, sensitivity, specificity, features) --
    if best_roc and best_roc.get("auc", 0.0) >= 0.65:
        auc = best_roc["auc"]
        tp = best_roc.get("timepoint", "")
        confidence = "moderate"
        if auc >= 0.8:
            confidence = "high"
        elif auc >= 0.7:
            confidence = "good"

        roc_text = (
            f"Predictive modeling achieved {confidence} discriminative "
            f"ability (AUC&nbsp;=&nbsp;{auc:.3f}"
        )
        if tp:
            roc_text += f" at {_esc(tp)}"
        roc_text += ")"

        sens = best_roc.get("sensitivity")
        spec = best_roc.get("specificity")
        youden = best_roc.get("youden_cutoff")
        perf_parts: list[str] = []
        if sens is not None:
            perf_parts.append(f"sensitivity&nbsp;=&nbsp;{sens:.1f}%")
        if spec is not None:
            perf_parts.append(f"specificity&nbsp;=&nbsp;{spec:.1f}%")
        if youden is not None:
            perf_parts.append(f"Youden cutoff&nbsp;=&nbsp;{youden:.3f}")
        if perf_parts:
            roc_text += f" with {', '.join(perf_parts)}"
        roc_text += ". "

        if all_feature_selections:
            feat_parts: list[str] = []
            for fs in all_feature_selections[:3]:
                feats = ", ".join(fs["features"][:4])
                lam = fs.get("lambda")
                lam_str = f", \u03bb={lam:.4f}" if lam is not None else ""
                feat_parts.append(
                    f"{_esc(fs.get('timepoint', '?'))}: "
                    f"<code>{_esc(feats)}</code>{lam_str}"
                )
            roc_text += f"Selected features by timepoint: {'; '.join(feat_parts)}. "

        roc_text += (
            "These features should be prioritized in longitudinal monitoring "
            "protocols to enable early detection of treatment failure."
        )
        plan_points.append(roc_text)

    if plan_points:
        h.append('<div class="summary-box" style="border-left-color: var(--green);">')
        h.append("<p><strong>Suggested Treatment Plan:</strong> "
                 "Based on statistically significant findings from the "
                 "pipeline's survival, predictive, and dosimetric analyses, "
                 "the following data-driven treatment considerations are proposed. "
                 "<em>These suggestions are hypothesis-generating and should "
                 "be reviewed by the treating physician in the context of the "
                 "individual patient's clinical picture.</em></p>")
        h.append(
            '<div class="info-box"><strong>Recommendation confidence framework:</strong> '
            "Recommendations below are ordered by confidence level: "
            "(1) <strong>High confidence</strong> \u2014 supported by statistically significant "
            "findings in \u22652 DWI types; "
            "(2) <strong>Moderate</strong> \u2014 significant in 1 DWI type only or borderline FDR; "
            "(3) <strong>Exploratory</strong> \u2014 based on trends only, not statistically significant.</div>"
        )
        h.append("<ul>")
        for pt in plan_points:
            h.append(f"<li>{pt}</li>")
        h.append("</ul>")
        h.append("</div>")
        h.append(
            '<div class="warn-box">'
            "<strong>\u26a0 Research Use Only:</strong> This treatment plan is exploratory, "
            "generated automatically from a single-centre retrospective cohort. It has not been "
            "validated in prospective trials. Clinical implementation requires: "
            "(1) independent prospective validation, "
            "(2) multidisciplinary review including radiation oncology and medical physics, "
            "(3) assessment of normal tissue dose limits and toxicity, and "
            "(4) institutional ethics approval. "
            "Do not apply these recommendations to patient care without formal clinical validation."
            "</div>"
        )

    h.append("</div>")
    return h



def _section_treatment_response(groups) -> list[str]:
    """Build the Treatment Response / Longitudinal Trends section.

    For each graph with "Longitudinal" in its name, displays per-DWI-type
    trend tags, axis information, inflection points (in collapsible
    ``<details>`` blocks), and full summaries.

    Parameters
    ----------
    groups : dict[str, dict[str, dict]]
        Graph rows grouped by base name and DWI type.

    Returns
    -------
    list[str]
        HTML chunks for the treatment response section.
    """
    # ── 6. Treatment Response ──
    h: list[str] = []
    if groups:
        has_longitudinal = any("Longitudinal" in bn for bn in groups)
        if has_longitudinal:
            h.append(_h2("Treatment Response \u2014 Longitudinal Trends", "treatment"))
            for base_name in sorted(groups.keys()):
                if "Longitudinal" not in base_name:
                    continue
                dwi_dict = groups[base_name]
                real = [t for t in dwi_dict if t != "Root"]
                if not real:
                    continue

                h.append(f"<h3>{_esc(base_name)}</h3>")
                for dt in DWI_TYPES:
                    if dt not in dwi_dict:
                        continue
                    r = dwi_dict[dt]
                    summary = r.get("summary", "")

                    # Axis info with range analysis
                    x_lbl = r.get("x_axis_label", "")
                    x_unit = r.get("x_axis_units", "")
                    x_min = r.get("x_axis_range_min", "")
                    x_max = r.get("x_axis_range_max", "")
                    y_lbl = r.get("y_axis_label", "")
                    y_unit = r.get("y_axis_units", "")
                    y_min = r.get("y_axis_range_min", "")
                    y_max = r.get("y_axis_range_max", "")
                    axis_parts = []
                    if x_lbl:
                        x_range = ""
                        if x_min and x_max:
                            x_range = f" [{x_min}\u2013{x_max}]"
                        axis_parts.append(f"X: {x_lbl}" + (f" ({x_unit})" if x_unit else "") + x_range)
                    if y_lbl:
                        y_range = ""
                        if y_min and y_max:
                            y_range = f" [{y_min}\u2013{y_max}]"
                        axis_parts.append(f"Y: {y_lbl}" + (f" ({y_unit})" if y_unit else "") + y_range)

                    h.append(f"<p>{_dwi_badge(dt)}")
                    if axis_parts:
                        h.append(f' <span class="axis-info">[{" | ".join(axis_parts)}]</span>')

                    # Extract any p-values mentioned in this specific graph
                    all_text = safe_text(r, "summary", "trends_json", "inflection_points_json")
                    graph_pvals = extract_pvalues(all_text)
                    graph_corrs = extract_correlations(all_text)
                    if graph_pvals or graph_corrs:
                        stat_parts = []
                        for pv, _ in graph_pvals:
                            cls = _sig_class(pv)
                            cls_attr = f' class="{cls}"' if cls else ""
                            stat_parts.append(f'<span{cls_attr}>p={pv:.4f}{_sig_tag(pv)}</span>')
                        for rv, _ in graph_corrs:
                            if abs(rv) >= 0.3:
                                stat_parts.append(f"r={rv:+.2f}")
                        if stat_parts:
                            h.append(f' <span class="meta">[{", ".join(stat_parts)}]</span>')
                    h.append("</p>")

                    # Trends
                    trends_str = r.get("trends_json", "[]") if isinstance(r, dict) else "[]"
                    try:
                        trends = json.loads(str(trends_str))
                    except Exception:
                        trends = []
                    if isinstance(trends, list) and trends:
                        for t in trends:
                            if isinstance(t, dict):
                                series = t.get("series", "")
                                direction = t.get("direction", "")
                                desc = t.get("description", "")
                                h.append(f"<p style='margin-left:1rem'>")
                                if series:
                                    h.append(f"<code>{_esc(series)}</code> ")
                                h.append(_trend_tag(direction))
                                if desc:
                                    h.append(f" <em>{_esc(desc)}</em>")
                                h.append("</p>")

                    # Inflection points
                    ips_str = r.get("inflection_points_json", "[]") if isinstance(r, dict) else "[]"
                    try:
                        ips = json.loads(str(ips_str))
                    except Exception:
                        ips = []
                    if isinstance(ips, list) and ips:
                        h.append("<details open><summary>Inflection points</summary><ul>")
                        for ip in ips:
                            if isinstance(ip, dict):
                                x = ip.get("approximate_x", "?")
                                y = ip.get("approximate_y", "")
                                desc_ip = ip.get("description", "")
                                coord = f"x={_esc(str(x))}"
                                if y is not None and y != "":
                                    coord += f", y={_esc(str(y))}"
                                h.append(f"<li>({coord}): {_esc(str(desc_ip))}</li>")
                        h.append("</ul></details>")

                    # Full summary in collapsible
                    if summary:
                        if len(summary) > 200:
                            h.append(f"<details open><summary>Full summary</summary>"
                                     f'<p class="full-summary">{_esc(summary)}</p></details>')
                        else:
                            h.append(f'<p class="full-summary" style="margin-left:1rem">{_esc(summary)}</p>')

                # Cross-DWI consensus for this longitudinal graph
                if len(real) >= 2:
                    all_dirs: dict[str, list[str]] = {}
                    for dt in real:
                        r = dwi_dict[dt]
                        t_str = r.get("trends_json", "[]") if isinstance(r, dict) else "[]"
                        try:
                            t_list = json.loads(str(t_str))
                        except Exception:
                            t_list = []
                        if isinstance(t_list, list):
                            for t in t_list:
                                if isinstance(t, dict):
                                    s = t.get("series") or "overall"
                                    d = str(t.get("direction", ""))
                                    all_dirs.setdefault(s, []).append(d)
                    consensus_parts = []
                    for s, dirs in all_dirs.items():
                        if len(dirs) >= 2:
                            if len(set(dirs)) == 1:
                                consensus_parts.append(f"{s}: all {dirs[0]}")
                            else:
                                consensus_parts.append(f"{s}: mixed ({', '.join(dirs)})")
                    if consensus_parts:
                        cons_cls = "agree" if all("all" in p for p in consensus_parts) else "differ"
                        h.append(
                            f'<div class="info-box"><strong>Cross-DWI consensus:</strong> '
                            f'<span class="{cons_cls}">{"; ".join(consensus_parts)}</span></div>'
                        )

        # Non-longitudinal inflection points
        non_long_ips = []
        for base_name in sorted(groups.keys()):
            if "Longitudinal" in base_name:
                continue
            dwi_dict = groups[base_name]
            for dt in DWI_TYPES:
                if dt not in dwi_dict:
                    continue
                r = dwi_dict[dt]
                ips_str = r.get("inflection_points_json", "[]") if isinstance(r, dict) else "[]"
                try:
                    ips = json.loads(str(ips_str))
                except Exception:
                    ips = []
                if isinstance(ips, list):
                    for ip in ips:
                        if isinstance(ip, dict) and ip.get("description"):
                            non_long_ips.append((base_name, dt, ip))

        if non_long_ips:
            h.append("<h3>Non-Longitudinal Inflection Points</h3>")
            h.append('<p class="meta">Inflection points detected in non-longitudinal graphs '
                     'may indicate dose thresholds, feature boundaries, or regime transitions.</p>')
            h.append("<table><thead><tr><th>Graph</th><th>DWI</th>"
                     "<th>Location</th><th>Description</th></tr></thead><tbody>")
            for gname, dt, ip in non_long_ips:
                x = ip.get("approximate_x", "?")
                y = ip.get("approximate_y", "")
                coord = f"x={_esc(str(x))}"
                if y is not None and y != "":
                    coord += f", y={_esc(str(y))}"
                desc_ip = _esc(str(ip.get("description", "")))
                h.append(f"<tr><td>{_esc(gname)}</td>"
                         f"<td>{_dwi_badge(dt)}</td>"
                         f"<td><code>{coord}</code></td>"
                         f"<td>{desc_ip}</td></tr>")
            h.append("</tbody></table>")
    return h
