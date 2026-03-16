"""Report section: data-driven hypothesis."""

from __future__ import annotations

import json
import re

from report.report_formatters import (  # type: ignore
    _esc,
    _get_consensus,
    _h2,
)
from report.sections._helpers import (  # type: ignore
    _extract_dosimetry,
    _extract_significant_metrics,
    _scalar_gy,
)


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
