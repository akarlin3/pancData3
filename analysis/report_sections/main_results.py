"""Report sections: main results group."""

from __future__ import annotations

import json
import re

from shared import (
    DWI_TYPES,
    extract_correlations,
    extract_pvalues,
    parse_dwi_info,
)
from report_formatters import (
    _cite,
    _copy_button,
    _dwi_badge,
    _esc,
    _get_consensus,
    _h2,
    _manuscript_sentence,
    _sig_class,
    _sig_tag,
    _stat_card,
    _table_caption,
    _trend_tag,
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

    total_sig = 0
    total_hrs = 0
    sig_hrs = 0
    if log_data:
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            # Find the best AUC across all timepoints for this DWI type.
            roc = log_data[dwi_type].get("stats_predictive", {}).get("roc_analyses", [])
            best_auc = max((r.get("auc", 0) for r in roc), default=0)
            if best_auc > 0:
                cards.append(_stat_card(f"Best AUC ({dwi_type})", f"{best_auc:.3f}"))
            # Count GLME metrics that pass their BH-adjusted significance threshold.
            sc = log_data[dwi_type].get("stats_comparisons", {})
            total_sig += len([g for g in sc.get("glme_details", []) if g["p"] < g["adj_alpha"]])
            # Count hazard ratios
            hrs = log_data[dwi_type].get("survival", {}).get("hazard_ratios", [])
            total_hrs += len(hrs)
            sig_hrs += len([hr for hr in hrs if hr.get("p", 1) < 0.05])

    if total_sig > 0:
        cards.append(_stat_card("GLME Sig. Interactions", str(total_sig), "across all DWI types"))

    if csv_data and csv_data.get("significant_metrics"):
        n_csv_sig = sum(len(v) for v in csv_data["significant_metrics"].values())
        cards.append(_stat_card("CSV Sig. Metrics", str(n_csv_sig), "from pipeline CSVs"))

    if csv_data and csv_data.get("fdr_global"):
        n_fdr = sum(len(v) for v in csv_data["fdr_global"].values())
        if n_fdr > 0:
            cards.append(_stat_card("FDR-Surviving", str(n_fdr), "metrics after BH correction"))

    if total_hrs > 0:
        cards.append(_stat_card("Hazard Ratios", f"{sig_hrs}/{total_hrs}", "significant covariates"))

    # Count notable correlations from vision data
    if rows:
        n_corr = 0
        for r in rows:
            all_text = r.get("summary", "") + " " + r.get("trends_json", "")
            for rval, _ in extract_correlations(all_text):
                if abs(rval) >= 0.3:
                    n_corr += 1
        if n_corr > 0:
            cards.append(_stat_card("Notable Correlations", str(n_corr), "|r| \u2265 0.3"))

    # Cross-DWI feature overlap
    if log_data and len(dwi_types_present) >= 2:
        tp_features: dict[str, dict[str, list[str]]] = {}
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            fs_list = log_data[dt].get("stats_predictive", {}).get("feature_selections", [])
            for fs in fs_list:
                tp = fs.get("timepoint", "?")
                tp_features.setdefault(tp, {})[dt] = fs.get("features", [])
        total_shared = 0
        total_all = 0
        for tp, dts_dict in tp_features.items():
            if len(dts_dict) < 2:
                continue
            all_feats: set[str] = set()
            for feats in dts_dict.values():
                all_feats.update(feats)
            total_all += len(all_feats)
            for feat in all_feats:
                if sum(1 for feats in dts_dict.values() if feat in feats) >= 2:
                    total_shared += 1
        if total_all > 0:
            cards.append(_stat_card("Feature Overlap",
                                    f"{total_shared}/{total_all}",
                                    "shared across DWI types"))

    # Cohort info from MAT data
    if mat_data:
        for dt in DWI_TYPES:
            if dt in mat_data and "longitudinal" in mat_data[dt]:
                lon = mat_data[dt]["longitudinal"]
                n_pat = lon.get("num_patients", 0)
                n_tp = lon.get("num_timepoints", 0)
                if n_pat > 0:
                    cards.append(_stat_card("Cohort", f"{n_pat} patients", f"{n_tp} timepoints"))
                    break  # only show once

    # Sanity check summary (convergence + alignment)
    if log_data:
        all_converged_count = 0
        total_conv_flags = 0
        total_dim_issues = 0
        sanity_types_checked = 0
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            san = log_data[dt].get("sanity_checks", {})
            if not san:
                continue
            sanity_types_checked += 1
            if san.get("all_converged"):
                all_converged_count += 1
            total_conv_flags += san.get("total_convergence", 0)
            total_dim_issues += san.get("dim_mismatches", 0) + san.get("nan_dose_warnings", 0)
        if sanity_types_checked > 0:
            if all_converged_count == sanity_types_checked and total_dim_issues == 0:
                cards.append(_stat_card("Data Quality", "All Passed",
                                        "convergence + alignment"))
            else:
                issue_parts = []
                if total_conv_flags > 0:
                    issue_parts.append(f"{total_conv_flags} conv.")
                if total_dim_issues > 0:
                    issue_parts.append(f"{total_dim_issues} align.")
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
    if mat_data:
        for dt in DWI_TYPES:
            if dt in mat_data and "longitudinal" in mat_data[dt]:
                lon = mat_data[dt]["longitudinal"]
                n_pat = lon.get("num_patients", 0)
                n_tp = lon.get("num_timepoints", 0)
                if n_pat > 0:
                    methods_parts.append(
                        f"Longitudinal DWI data from {n_pat} patients across "
                        f"{n_tp} timepoints were analysed."
                    )
                    break
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
    best_overall_auc = 0
    best_auc_type = ""
    if log_data:
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            roc = log_data[dwi_type].get("stats_predictive", {}).get("roc_analyses", [])
            for r_item in roc:
                a = r_item.get("auc", 0)
                if a > best_overall_auc:
                    best_overall_auc = a
                    best_auc_type = dwi_type
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
    if log_data:
        san_all_ok = True
        san_total_flags = 0
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            san = log_data[dt].get("sanity_checks", {})
            if san and not san.get("all_converged"):
                san_all_ok = False
            san_total_flags += san.get("total_convergence", 0) if san else 0
        if san_all_ok and san_total_flags == 0:
            result_bullets.append("All voxel-level model fits converged successfully across DWI types")
        elif san_total_flags > 0:
            result_bullets.append(f"{san_total_flags} convergence flag(s) raised across DWI types (see Data Completeness)")
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
    # These gate which factors are highlighted in the hypothesis and plan.
    sig_glme: list[dict] = []            # GLME interactions with p < 0.05
    sig_glme_details: list[dict] = []    # Per-metric rows with p < adj_alpha
    sig_fdr_timepoints: list[dict] = []  # Timepoints with FDR-sig metrics
    sig_hr: list[dict] = []              # Significant hazard ratios
    best_roc: dict | None = None         # Best ROC analysis result
    best_roc_dt = ""                     # DWI type of best ROC
    all_feature_selections: list[dict] = []  # Per-timepoint elastic net
    if log_data:
        for dt in DWI_TYPES:
            dt_data = log_data.get(dt)
            if not dt_data:
                continue
            # GLME interaction p-values
            stats_cmp = dt_data.get("stats_comparisons") or {}
            for p_val in stats_cmp.get("glme_interactions", []):
                if p_val < 0.05:
                    sig_glme.append({"p": p_val, "dwi_type": dt})
            for detail in stats_cmp.get("glme_details", []):
                if detail.get("p", 1.0) < detail.get("adj_alpha", 0.05):
                    sig_glme_details.append({**detail, "dwi_type": dt})
            for fdr_tp in stats_cmp.get("fdr_timepoints", []):
                if fdr_tp.get("n_significant", 0) > 0:
                    sig_fdr_timepoints.append({**fdr_tp, "dwi_type": dt})
            # Survival hazard ratios
            surv = dt_data.get("survival") or {}
            for hr_entry in surv.get("hazard_ratios", []):
                if hr_entry.get("p", 1.0) < 0.05:
                    sig_hr.append({**hr_entry, "dwi_type": dt})
            # Predictive model ROC/AUC
            pred = dt_data.get("stats_predictive") or {}
            for roc in pred.get("roc_analyses", []):
                auc = roc.get("auc", 0.0)
                if best_roc is None or auc > best_roc.get("auc", 0.0):
                    best_roc = roc
                    best_roc_dt = dt
            for fs in pred.get("feature_selections", []):
                if fs.get("features"):
                    all_feature_selections.append({**fs, "dwi_type": dt})

    # ── Extract dosimetry from mat_data ──
    dosimetry: dict = {}
    dosimetry_dt = ""
    if mat_data:
        for dt in DWI_TYPES:
            dosi = (mat_data.get(dt) or {}).get("dosimetry")
            if dosi and not dosimetry:
                dosimetry = dosi
                dosimetry_dt = dt

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
        d95_adc = dosimetry.get("d95_adc_mean")
        v50_adc = dosimetry.get("v50_adc_mean")
        d95_d = dosimetry.get("d95_d_mean")
        v50_d = dosimetry.get("v50_d_mean")
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
        h.append("<ul>")
        for pt in plan_points:
            h.append(f"<li>{pt}</li>")
        h.append("</ul>")
        h.append("</div>")

    h.append("</div>")
    return h



def _section_statistical_significance(rows, csv_data, log_data, dwi_types_present) -> list[str]:
    """Build the Statistical Significance section.

    Aggregates significant findings from three sources:
    1. Vision-extracted p-values from graph summaries/trends.
    2. Pipeline CSV significant metrics (Significant_LF_Metrics.csv).
    3. GLME interaction test details and FDR timepoints from log parsing.

    Also includes a Borderline Findings subsection for p-values between
    0.05 and 0.10 (trend-worthy but not significant).

    Parameters
    ----------
    rows : list[dict]
        Vision CSV rows.
    csv_data : dict or None
        Parsed pipeline CSV exports.
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.

    Returns
    -------
    list[str]
        HTML chunks for the statistical significance section.
    """
    # ── 3. Statistical Significance ──
    h: list[str] = []
    h.append(_h2("Statistical Significance", "significance"))

    # ── Source 1: Vision-extracted p-values ──
    # Search through summaries, trends, and inflection point descriptions
    # for p-value patterns (e.g. "p = 0.032") using regex extraction.
    sig_findings = []
    borderline_findings = []
    if rows:
        for r in rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            all_text = r["summary"] + " " + r["trends_json"] + " " + r["inflection_points_json"]
            for pval, context in extract_pvalues(all_text):
                if pval < 0.05:
                    sig_findings.append((pval, dwi_type, base_name, context))
                elif pval < 0.10:
                    borderline_findings.append((pval, dwi_type, base_name, context))

    if sig_findings:
        sig_findings.sort()
        h.append(f"<h3>Vision-Extracted Significant Findings (p &lt; 0.05) \u2014 {len(sig_findings)} total</h3>")
        h.append("<table>")
        h.append(_table_caption("Vision-Extracted Significant Findings",
                                f"(N\u2009=\u2009{len(sig_findings)} findings with p < 0.05)"))
        h.append("<thead><tr><th>p-value</th><th>Sig</th><th>DWI</th>"
                 "<th>Graph</th><th>Context</th></tr></thead><tbody>")
        for p, dwi, graph, ctx in sig_findings:
            ctx_clean = _esc(ctx.replace("\n", " ")[:120])
            cls = _sig_class(p)
            cls_attr = f' class="{cls}"' if cls else ""
            h.append(f"<tr><td{cls_attr}>{p:.4f}</td><td{cls_attr}>{_esc(_sig_tag(p))}</td>"
                     f"<td>{_dwi_badge(dwi)}</td><td>{_esc(graph)}</td>"
                     f"<td>{ctx_clean}</td></tr>")
        h.append("</tbody></table>")

    # Borderline findings (0.05 <= p < 0.10)
    if borderline_findings:
        borderline_findings.sort()
        h.append(f"<details><summary><strong>Borderline Findings (0.05 \u2264 p &lt; 0.10)</strong> "
                 f"\u2014 {len(borderline_findings)} finding(s)</summary>")
        h.append('<div class="warn-box">These findings approach but do not reach conventional '
                 'significance. They may warrant monitoring in future analyses.</div>')
        h.append("<table><thead><tr><th>p-value</th><th>DWI</th>"
                 "<th>Graph</th><th>Context</th></tr></thead><tbody>")
        for p, dwi, graph, ctx in borderline_findings:
            ctx_clean = _esc(ctx.replace("\n", " ")[:120])
            h.append(f"<tr><td>{p:.4f}</td>"
                     f"<td>{_dwi_badge(dwi)}</td><td>{_esc(graph)}</td>"
                     f"<td>{ctx_clean}</td></tr>")
        h.append("</tbody></table>")
        h.append("</details>")

    # ── Source 2: Pipeline CSV significant metrics ──
    if csv_data and csv_data.get("significant_metrics"):
        h.append("<h3>Pipeline CSV Significant Metrics</h3>")
        for dwi_type in DWI_TYPES:
            if dwi_type not in csv_data["significant_metrics"]:
                continue
            csv_rows = csv_data["significant_metrics"][dwi_type]
            h.append(f"<p>{_dwi_badge(dwi_type)} \u2014 {len(csv_rows)} significant metric(s)</p>")
            if csv_rows:
                headers = list(csv_rows[0].keys())[:10]
                h.append("<table><thead><tr>")
                for hdr in headers:
                    h.append(f"<th>{_esc(hdr)}</th>")
                h.append("</tr></thead><tbody>")
                for cr in csv_rows:
                    h.append("<tr>")
                    for hdr in headers:
                        val = str(cr.get(hdr, ""))
                        # Highlight p-value columns
                        try:
                            fval = float(val)
                            if "p" in hdr.lower() and 0 < fval < 0.05:
                                cls = _sig_class(fval)
                                h.append(f'<td class="{cls}">{_esc(val[:40])}</td>')
                                continue
                        except (ValueError, TypeError):
                            pass
                        h.append(f"<td>{_esc(val[:40])}</td>")
                    h.append("</tr>")
                h.append("</tbody></table>")

    # ── Source 3: GLME interaction details from log parsing ──
    if log_data:
        has_glme = any(
            log_data[d].get("stats_comparisons", {}).get("glme_details")
            for d in dwi_types_present if d in log_data
        )
        if has_glme:
            h.append("<h3>GLME Mixed-Effects Interaction Tests</h3>")
            h.append("<p><em>Per-metric interaction p-values from generalised linear mixed-effects models. "
                     "Adjusted alpha accounts for BH-FDR multiple comparison correction.</em></p>")

            for dwi_type in dwi_types_present:
                if dwi_type not in log_data:
                    continue
                sc = log_data[dwi_type].get("stats_comparisons", {})
                glme_details = sc.get("glme_details", [])
                glme_excl = sc.get("glme_excluded")
                if not glme_details:
                    continue

                h.append(f"<h4>{_dwi_badge(dwi_type)}</h4>")

                if glme_excl:
                    h.append(f'<div class="warn-box">'
                              f'\u26a0\ufe0f Competing-risk exclusion: '
                              f'<strong>{glme_excl["n_excluded"]}/{glme_excl["n_total"]}</strong> '
                              f'({glme_excl["pct"]:.1f}%) patients excluded from GLME.</div>')

                sig_details = [g for g in glme_details if g["p"] < g["adj_alpha"]]
                nsig_details = [g for g in glme_details if g["p"] >= g["adj_alpha"]]

                if sig_details:
                    h.append(f'<div class="info-box">'
                              f'\u2705 {len(sig_details)} of {len(glme_details)} metrics show significant '
                              f'interaction (p &lt; adjusted \u03b1).</div>')

                h.append("<table>")
                h.append(_table_caption(
                    f"GLME Interaction Test Results ({dwi_type})",
                    f"N\u2009=\u2009{len(glme_details)} metrics tested."))
                h.append("<thead><tr>"
                         "<th>Metric</th><th>p-value</th><th>Adj. \u03b1</th><th>Sig.</th>"
                         "</tr></thead><tbody>")
                for g in sorted(glme_details, key=lambda x: x["p"]):
                    is_sig = g["p"] < g["adj_alpha"]
                    cls = _sig_class(g["p"]) if is_sig else ""
                    cls_attr = f' class="{cls}"' if cls else ""
                    sig_mark = _sig_tag(g["p"]) if is_sig else ""
                    sig_display = _esc(sig_mark) or '\u2014'
                    h.append(
                        f"<tr><td><code>{_esc(g['metric'])}</code></td>"
                        f"<td{cls_attr}>{g['p']:.4f}</td>"
                        f"<td>{g['adj_alpha']:.4f}</td>"
                        f"<td{cls_attr}>{sig_display}</td></tr>"
                    )
                h.append("</tbody></table>")

        # FDR-significant timepoints
        has_fdr_tp = any(
            log_data[d].get("stats_comparisons", {}).get("fdr_timepoints")
            for d in dwi_types_present if d in log_data
        )
        if has_fdr_tp:
            h.append("<h3>FDR-Significant Metrics by Timepoint</h3>")
            h.append("<table><thead><tr><th>DWI</th><th>Timepoint</th>"
                     "<th>N Significant</th></tr></thead><tbody>")
            for dwi_type in dwi_types_present:
                if dwi_type not in log_data:
                    continue
                fdr_tps = log_data[dwi_type].get("stats_comparisons", {}).get("fdr_timepoints", [])
                for ftp in fdr_tps:
                    h.append(
                        f"<tr><td>{_dwi_badge(dwi_type)}</td>"
                        f"<td><code>{_esc(ftp['timepoint'])}</code></td>"
                        f"<td><strong>{ftp['n_significant']}</strong></td></tr>"
                    )
            h.append("</tbody></table>")

    if not sig_findings and not (csv_data and csv_data.get("significant_metrics")) and not log_data:
        h.append("<p>No significant findings extracted.</p>")
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
                    all_text = (r.get("summary", "") + " " +
                                r.get("trends_json", "") + " " +
                                r.get("inflection_points_json", ""))
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
                            h.append(f"<details><summary>Full summary</summary>"
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
    n_patients = 0
    n_timepoints = 0
    if mat_data:
        for dt in DWI_TYPES:
            if dt in mat_data and "longitudinal" in mat_data[dt]:
                lon = mat_data[dt]["longitudinal"]
                n = lon.get("num_patients", 0)
                tp = lon.get("num_timepoints", 0)
                if n > n_patients:
                    n_patients = n
                    n_timepoints = tp
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
                    total_glme_sig += 1
                    if g["metric"] not in sig_metric_names:
                        sig_metric_names.append(g["metric"])
            if sc.get("glme_excluded") and glme_excluded is None:
                glme_excluded = sc["glme_excluded"]
        if total_glme_tested > 0:
            exc_str = ""
            if glme_excluded:
                n_analysed = glme_excluded["n_total"] - glme_excluded["n_excluded"]
                exc_str = (
                    f"Among {n_analysed} evaluable patients "
                    f"(excluding {glme_excluded['n_excluded']} with "
                    f"competing-risk events), "
                )
            if total_glme_sig > 0:
                metric_str = ", ".join(sig_metric_names[:3])
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
        if best_auc > 0:
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
            for hr_item in sorted(sig_hrs, key=lambda x: x["p"])[:2]:
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
                abs_log_hr = abs(log_hr)
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
    if mat_data:
        for dt in DWI_TYPES:
            dosi = (mat_data.get(dt) or {}).get("dosimetry")
            if dosi:
                d95 = dosi.get("d95_adc_mean")
                v50 = dosi.get("v50_adc_mean")
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
                break

    # ── Longitudinal trend sentence ──
    if groups and "Longitudinal_Mean_Metrics" in groups:
        d_trends = []
        f_trends = []
        for dt, r in groups["Longitudinal_Mean_Metrics"].items():
            if dt == "Root":
                continue
            try:
                trends = json.loads(str(r.get("trends_json", "[]")))
                for t in trends:
                    if isinstance(t, dict):
                        series = t.get("series", "")
                        direction = t.get("direction", "").lower()
                        if series == "Mean D":
                            d_trends.append(direction)
                        elif series == "Mean f":
                            f_trends.append(direction)
            except Exception:
                pass
        d_cons = _get_consensus(d_trends)
        f_cons = _get_consensus(f_trends)
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
    if groups:
        n_agree = 0
        n_total_series = 0
        for base_name, dwi_dict in groups.items():
            real = [t for t in dwi_dict if t != "Root"]
            if len(real) < 2:
                continue
            all_trends_dict: dict[str, list] = {}
            for dt_key in DWI_TYPES:
                if dt_key in dwi_dict:
                    try:
                        all_trends_dict[dt_key] = json.loads(
                            str(dwi_dict[dt_key].get("trends_json", "[]"))
                        )
                    except Exception:
                        pass
            if len(all_trends_dict) >= 2:
                all_series: set[str] = set()
                for trends in all_trends_dict.values():
                    for t in trends:
                        if isinstance(t, dict):
                            all_series.add(t.get("series") or "overall")
                for series in all_series:
                    dirs: dict[str, str] = {}
                    for dt_key, trends in all_trends_dict.items():
                        for t in trends:
                            if isinstance(t, dict) and (t.get("series") or "overall") == series:
                                dirs[dt_key] = str(t.get("direction", ""))
                    if len(dirs) >= 2:
                        n_total_series += 1
                        if len(set(dirs.values())) == 1:
                            n_agree += 1
        if n_total_series > 0:
            pct = 100 * n_agree / n_total_series
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
    n_patients = 0
    n_timepoints = 0
    if mat_data:
        for dt in DWI_TYPES:
            if dt in mat_data and "longitudinal" in mat_data[dt]:
                lon = mat_data[dt]["longitudinal"]
                n = lon.get("num_patients", 0)
                tp = lon.get("num_timepoints", 0)
                if n > n_patients:
                    n_patients = n
                    n_timepoints = tp

    baseline_exc = None
    if log_data:
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            bl = log_data[dt].get("baseline", {})
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
                    total_sig += 1
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
                    f"After excluding {glme_exc['n_excluded']} of "
                    f"{glme_exc['n_total']} patients ({glme_exc['pct']:.1f}%) "
                    f"with competing-risk events, "
                )
            else:
                parts.append("")

            if total_sig > 0:
                metric_str = ", ".join(sig_metrics[:5])
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
                    f"\u03c7\u00b2({global_lrt['df']}) = {global_lrt['chi2']:.2f}, "
                    f"p = {global_lrt['p']:.4f}. "
                )
            else:
                parts.append(". ")

            if sig_hrs:
                for dt, hr_item in sorted(sig_hrs, key=lambda x: x[1]["p"])[:3]:
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
    if groups and len(dwi_types_present) >= 2:
        n_agree = 0
        n_total_series = 0
        for base_name, dwi_dict in groups.items():
            real = [t for t in dwi_dict if t != "Root"]
            if len(real) < 2:
                continue
            all_trends_dict: dict[str, list] = {}
            for dt_key in DWI_TYPES:
                if dt_key in dwi_dict:
                    try:
                        all_trends_dict[dt_key] = json.loads(
                            str(dwi_dict[dt_key].get("trends_json", "[]"))
                        )
                    except Exception:
                        pass
            if len(all_trends_dict) >= 2:
                all_series: set[str] = set()
                for trends in all_trends_dict.values():
                    for t in trends:
                        if isinstance(t, dict):
                            all_series.add(t.get("series") or "overall")
                for series in all_series:
                    dirs: dict[str, str] = {}
                    for dt_key, trends in all_trends_dict.items():
                        for t in trends:
                            if isinstance(t, dict) and (t.get("series") or "overall") == series:
                                dirs[dt_key] = str(t.get("direction", ""))
                    if len(dirs) >= 2:
                        n_total_series += 1
                        if len(set(dirs.values())) == 1:
                            n_agree += 1
        if n_total_series > 0:
            pct = 100 * n_agree / n_total_series
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
    if mat_data:
        for dt in DWI_TYPES:
            dosi = (mat_data.get(dt) or {}).get("dosimetry")
            if dosi:
                d95 = dosi.get("d95_adc_mean")
                v50 = dosi.get("v50_adc_mean")
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
                break

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


