"""Section builders for the HTML analysis report.

Each public function in this module corresponds to one major section of the
HTML report.  They accept pre-loaded data (log dicts, CSV dicts, vision rows,
grouped graph dicts) and return a ``list[str]`` of HTML chunks that are
concatenated by :func:`generate_report.generate_report`.

Sections are assembled in the order they appear in the report:

1. Executive Summary
2. Cohort Overview
3. Data-Driven Hypothesis
4. Graph Analysis Overview
5. Statistics by Graph Type
6. Statistical Significance (with Borderline Findings)
7. Cross-DWI Comparison (all graph groups)
8. Notable Correlations (grouped by strength)
9. Treatment Response (Longitudinal Trends + non-longitudinal inflection points)
10. Predictive Performance (ROC/AUC, features, Cox PH, lambda trends)
11. Supplemental Data (MAT files, full core method matrix)
12. Appendix: All Graphs
"""

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
    _dwi_badge,
    _effect_size_class,
    _effect_size_label,
    _esc,
    _forest_plot_cell,
    _get_consensus,
    _h2,
    _sig_class,
    _sig_tag,
    _stat_card,
    _table_caption,
    _trend_tag,
)


def _section_publication_header() -> list[str]:
    """Build a publication metadata block with author/institution placeholders.

    This section provides fill-in-the-blank placeholders for manuscript
    metadata that cannot be auto-generated (author names, affiliations,
    corresponding author, IRB approval number).

    Returns
    -------
    list[str]
        HTML chunks for the publication metadata block.
    """
    h: list[str] = []
    h.append('<div class="pub-meta">')
    h.append("<h4>Manuscript Information</h4>")
    h.append('<p><strong>Title:</strong> <span class="placeholder">'
             "Diffusion-Weighted MRI Biomarkers for Treatment Response "
             "Assessment in Pancreatic Cancer: A Multi-Strategy Analysis"
             "</span></p>")
    h.append('<p><strong>Authors:</strong> <span class="placeholder">'
             "[Author 1, Author 2, ... ]</span></p>")
    h.append('<p><strong>Affiliations:</strong> <span class="placeholder">'
             "Department of Medical Physics, Memorial Sloan Kettering Cancer Center, "
             "New York, NY</span></p>")
    h.append('<p><strong>Corresponding Author:</strong> <span class="placeholder">'
             "[Name, email]</span></p>")
    h.append('<p><strong>IRB Approval:</strong> <span class="placeholder">'
             "[Protocol number]</span></p>")
    h.append("</div>")
    return h


def _section_data_availability() -> list[str]:
    """Build the Data Availability Statement section.

    Returns
    -------
    list[str]
        HTML chunks for the data availability section.
    """
    h: list[str] = []
    h.append(_h2("Data Availability", "data-availability"))
    h.append('<div class="methods-box">')
    h.append(
        "<p>The clinical imaging data supporting this study are subject to "
        "institutional review board (IRB) restrictions and cannot be made "
        "publicly available due to patient privacy requirements under HIPAA. "
        "De-identified summary statistics and analysis code are available "
        "from the corresponding author upon reasonable request.</p>"
    )
    h.append(
        "<p>The analysis pipeline source code (pancData3) is available "
        "under the MIT License. Post-hoc analysis scripts and report "
        "generation tools are included in the repository.</p>"
    )
    h.append("</div>")
    return h


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
        from shared import extract_correlations
        n_corr = 0
        for r in rows:
            all_text = r.get("summary", "") + " " + r.get("trends_json", "")
            for rval, _ in extract_correlations(all_text):
                if abs(rval) >= 0.3:
                    n_corr += 1
        if n_corr > 0:
            cards.append(_stat_card("Notable Correlations", str(n_corr), "|r| \u2265 0.3"))

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


def _section_cohort_overview(mat_data, log_data, dwi_types_present) -> list[str]:
    """Build the Cohort Overview section.

    Displays patient counts, timepoint counts, and data quality
    summary across all DWI types from MAT file and log data.

    Parameters
    ----------
    mat_data : dict
        Parsed MAT file metrics.
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.

    Returns
    -------
    list[str]
        HTML chunks for the cohort overview section.
    """
    h: list[str] = []
    has_data = False

    # Check if we have any cohort data to show
    if mat_data:
        for dt in DWI_TYPES:
            if dt in mat_data and "longitudinal" in mat_data[dt]:
                has_data = True
                break

    if log_data:
        for dt in dwi_types_present:
            if dt in log_data:
                bl = log_data[dt].get("baseline", {})
                if bl.get("total_outliers") or bl.get("baseline_exclusion"):
                    has_data = True
                    break

    if not has_data:
        return h

    h.append(_h2("Cohort Overview", "cohort"))

    # Longitudinal data summary from MAT files
    if mat_data:
        cohort_rows = []
        for dt in DWI_TYPES:
            if dt not in mat_data or "longitudinal" not in mat_data[dt]:
                continue
            lon = mat_data[dt]["longitudinal"]
            n_pat = lon.get("num_patients", 0)
            n_tp = lon.get("num_timepoints", 0)
            if n_pat > 0 or n_tp > 0:
                cohort_rows.append((dt, n_pat, n_tp))

        if cohort_rows:
            h.append("<h3>Longitudinal Data Dimensions</h3>")
            h.append("<table><thead><tr><th>DWI Type</th><th>Patients</th>"
                     "<th>Timepoints</th></tr></thead><tbody>")
            for dt, n_pat, n_tp in cohort_rows:
                h.append(f"<tr><td>{_dwi_badge(dt)}</td>"
                         f"<td><strong>{n_pat}</strong></td>"
                         f"<td>{n_tp}</td></tr>")
            h.append("</tbody></table>")

    # Cross-DWI data quality summary
    if log_data:
        quality_rows = []
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            bl = log_data[dt].get("baseline", {})
            total_out = bl.get("total_outliers")
            baseline_exc = bl.get("baseline_exclusion")
            if total_out or baseline_exc:
                out_str = f"{total_out['n_removed']}/{total_out['n_total']} ({total_out['pct']:.1f}%)" if total_out else "\u2014"
                exc_str = f"{baseline_exc['n_excluded']}/{baseline_exc['n_total']}" if baseline_exc else "\u2014"
                lf_inc = baseline_exc.get("lf_rate_included") if baseline_exc else None
                lf_exc = baseline_exc.get("lf_rate_excluded") if baseline_exc else None
                lf_str = f"{lf_inc:.1f}% / {lf_exc:.1f}%" if lf_inc is not None and lf_exc is not None else "\u2014"
                quality_rows.append((dt, out_str, exc_str, lf_str))

        if quality_rows:
            h.append("<h3>Data Quality Summary</h3>")
            h.append("<table><thead><tr><th>DWI Type</th><th>Outliers Removed</th>"
                     "<th>Baseline Excluded</th><th>LF Rate (incl/excl)</th>"
                     "</tr></thead><tbody>")
            for dt, out_str, exc_str, lf_str in quality_rows:
                h.append(f"<tr><td>{_dwi_badge(dt)}</td>"
                         f"<td>{_esc(out_str)}</td>"
                         f"<td>{_esc(exc_str)}</td>"
                         f"<td>{_esc(lf_str)}</td></tr>")
            h.append("</tbody></table>")

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


def _section_graph_overview(rows) -> list[str]:
    """Build the Graph Analysis Overview section.

    Displays two side-by-side summary tables:
    - Count of graphs by type (line, scatter, box, etc.)
    - Count of graphs by DWI type (Standard, dnCNN, IVIMnet, Root)

    Parameters
    ----------
    rows : list[dict]
        Vision CSV rows (may be empty, in which case the section is skipped).

    Returns
    -------
    list[str]
        HTML chunks for the graph overview section.
    """
    # ── 2. Graph Analysis Overview ──
    h: list[str] = []
    if rows:
        h.append(_h2("Graph Analysis Overview", "graph-overview"))

        type_counts: dict[str, int] = {}
        for r in rows:
            gt = r.get("graph_type", "unknown")
            type_counts[gt] = type_counts.get(gt, 0) + 1

        dwi_counts: dict[str, int] = {}
        for r in rows:
            dwi_type, _ = parse_dwi_info(r["file_path"])
            dwi_counts[dwi_type] = dwi_counts.get(dwi_type, 0) + 1

        h.append('<div style="display:grid;grid-template-columns:1fr 1fr;gap:1rem">')

        h.append("<div>")
        h.append("<h3>By Graph Type</h3>")
        h.append("<table><thead><tr><th>Type</th><th>Count</th></tr></thead><tbody>")
        for gt, cnt in sorted(type_counts.items(), key=lambda x: -x[1]):
            h.append(f"<tr><td>{_esc(gt)}</td><td>{cnt}</td></tr>")
        h.append("</tbody></table>")
        h.append("</div>")

        h.append("<div>")
        h.append("<h3>By DWI Type</h3>")
        h.append("<table><thead><tr><th>DWI Type</th><th>Graphs</th></tr></thead><tbody>")
        for dt in DWI_TYPES + ["Root"]:
            if dt in dwi_counts:
                h.append(f"<tr><td>{_dwi_badge(dt)}</td><td>{dwi_counts[dt]}</td></tr>")
        h.append("</tbody></table>")
        h.append("</div>")

        h.append("</div>")
    return h


def _section_stats_by_graph_type(rows) -> list[str]:
    """Build the Statistics by Graph Type section.

    Aggregates significant p-values, correlations, and trend directions
    by graph type (box, line, scatter, heatmap, etc.), providing a
    bird's-eye view of which graph types carry the most statistical signal.

    Parameters
    ----------
    rows : list[dict]
        Vision CSV rows (may be empty).

    Returns
    -------
    list[str]
        HTML chunks for the statistics by graph type section.
    """
    h: list[str] = []
    if not rows:
        return h

    h.append(_h2("Statistics by Graph Type", "stats-by-type"))

    # Aggregate per graph type
    type_stats: dict[str, dict] = {}
    for r in rows:
        gt = r.get("graph_type", "unknown")
        if gt not in type_stats:
            type_stats[gt] = {
                "count": 0, "sig_p": 0, "nonsig_p": 0,
                "correlations": 0, "trends": {"increasing": 0, "decreasing": 0, "stable": 0, "other": 0},
            }
        ts = type_stats[gt]
        ts["count"] += 1

        # P-values
        all_text = r.get("summary", "") + " " + r.get("trends_json", "") + " " + r.get("inflection_points_json", "")
        for pval, _ in extract_pvalues(all_text):
            if pval < 0.05:
                ts["sig_p"] += 1
            else:
                ts["nonsig_p"] += 1

        # Correlations
        for rval, _ in extract_correlations(all_text):
            if abs(rval) >= 0.3:
                ts["correlations"] += 1

        # Trends
        trends_str = r.get("trends_json", "[]") or "[]"
        try:
            trends = json.loads(str(trends_str))
        except Exception:
            trends = []
        if isinstance(trends, list):
            for t in trends:
                if not isinstance(t, dict):
                    continue
                d = str(t.get("direction", "")).lower()
                if "increas" in d or "up" in d or "higher" in d or "rising" in d:
                    ts["trends"]["increasing"] += 1
                elif "decreas" in d or "down" in d or "lower" in d or "falling" in d or "drop" in d:
                    ts["trends"]["decreasing"] += 1
                elif "flat" in d or "stable" in d or "constant" in d:
                    ts["trends"]["stable"] += 1
                else:
                    ts["trends"]["other"] += 1

    # Summary table
    h.append("<table><thead><tr>"
             "<th>Graph Type</th><th>Count</th><th>Sig. p-values</th>"
             "<th>Non-sig. p</th><th>Correlations</th>"
             "<th>\u2191 Incr.</th><th>\u2193 Decr.</th><th>\u2192 Stable</th><th>Other</th>"
             "</tr></thead><tbody>")
    for gt in sorted(type_stats.keys(), key=lambda k: -type_stats[k]["count"]):
        ts = type_stats[gt]
        tr = ts["trends"]
        sig_cls = ' class="sig-1"' if ts["sig_p"] > 0 else ""
        h.append(f"<tr><td><strong>{_esc(gt)}</strong></td>"
                 f"<td>{ts['count']}</td>"
                 f"<td{sig_cls}>{ts['sig_p']}</td>"
                 f"<td>{ts['nonsig_p']}</td>"
                 f"<td>{ts['correlations']}</td>"
                 f"<td>{tr['increasing']}</td>"
                 f"<td>{tr['decreasing']}</td>"
                 f"<td>{tr['stable']}</td>"
                 f"<td>{tr['other']}</td></tr>")
    h.append("</tbody></table>")

    # Highlight graph types with the most signal
    most_sig = max(type_stats.items(), key=lambda x: x[1]["sig_p"], default=None)
    if most_sig and most_sig[1]["sig_p"] > 0:
        h.append(f'<div class="info-box"><strong>{_esc(most_sig[0])}</strong> graphs carry the most '
                 f'significant findings ({most_sig[1]["sig_p"]} p &lt; 0.05).</div>')

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


def _section_cross_dwi_comparison(groups, csv_data) -> list[str]:
    """Build the Cross-DWI Comparison section.

    For priority graphs that exist in multiple DWI types, displays a
    trend-direction comparison table showing whether Standard, dnCNN,
    and IVIMnet agree or differ on each data series' direction.

    Also includes a significance inconsistency table from the CSV
    cross-reference (metrics significant in some DWI types but not others).

    Parameters
    ----------
    groups : dict[str, dict[str, dict]]
        Graph rows grouped by base name and DWI type.
    csv_data : dict or None
        Parsed pipeline CSV exports (contains ``cross_reference`` key).

    Returns
    -------
    list[str]
        HTML chunks for the cross-DWI comparison section.
    """
    # ── 4. Cross-DWI Comparison ──
    h: list[str] = []
    if groups:
        h.append(_h2("Cross-DWI Comparison", "cross-dwi"))

        priority_graphs = [
            "Dose_vs_Diffusion", "Longitudinal_Mean_Metrics",
            "Longitudinal_Mean_Metrics_ByOutcome", "Feature_BoxPlots",
            "Feature_Histograms", "core_method_dice_heatmap",
        ]

        # Build the list of all graph groups with 2+ DWI types, priority first
        all_comparable = []
        seen = set()
        for base_name in priority_graphs:
            if base_name in groups:
                real = [t for t in groups[base_name] if t != "Root"]
                if len(real) >= 2:
                    all_comparable.append((base_name, True))
                    seen.add(base_name)
        for base_name in sorted(groups.keys()):
            if base_name not in seen:
                real = [t for t in groups[base_name] if t != "Root"]
                if len(real) >= 2:
                    all_comparable.append((base_name, False))

        # Summary counts
        n_agree = 0
        n_differ = 0

        for base_name, is_priority in all_comparable:
            dwi_dict = groups[base_name]

            priority_label = ' <span class="badge badge-standard">priority</span>' if is_priority else ""
            h.append(f"<h3>{_esc(base_name)}{priority_label}</h3>")

            all_trends: dict[str, list] = {}
            for dt in DWI_TYPES:
                if dt in dwi_dict:
                    try:
                        all_trends[dt] = json.loads(str(dwi_dict[dt].get("trends_json", "[]")))
                    except Exception:
                        pass

            if all_trends:
                all_series: set[str] = set()
                for trends in all_trends.values():
                    for t in trends:
                        all_series.add(t.get("series") or "overall")

                h.append("<table><thead><tr><th>Series</th>")
                for dt in DWI_TYPES:
                    h.append(f"<th>{_esc(dt)}</th>")
                h.append("<th>Agreement</th></tr></thead><tbody>")

                for series in sorted(all_series):
                    directions: dict[str, str] = {}
                    for dt in DWI_TYPES:
                        if dt not in all_trends:
                            continue
                        for t in all_trends[dt]:
                            if isinstance(t, dict):
                                if (t.get("series") or "overall") == series:
                                    directions[dt] = str(t.get("direction", ""))

                    if len(directions) >= 2:
                        vals = list(directions.values())
                        if len(set(vals)) == 1:
                            agree_html = '<span class="agree">AGREE</span>'
                            n_agree += 1
                        else:
                            agree_html = '<span class="differ">DIFFER</span>'
                            n_differ += 1
                        h.append(f"<tr><td>{_esc(series)}</td>")
                        for dt in DWI_TYPES:
                            d_str = directions.get(dt, "-")
                            cell = _trend_tag(d_str) if d_str != "-" else "-"
                            h.append(f"<td>{cell}</td>")
                        h.append(f"<td>{agree_html}</td></tr>")

                h.append("</tbody></table>")

        # Overall agreement summary
        if n_agree + n_differ > 0:
            pct_agree = 100 * n_agree / (n_agree + n_differ)
            cls = "agree" if pct_agree >= 70 else ("differ" if pct_agree < 50 else "")
            cls_attr = f' class="{cls}"' if cls else ""
            h.append(f'<div class="summary-box"><strong>Cross-DWI Agreement:</strong> '
                     f'<span{cls_attr}>{n_agree}/{n_agree + n_differ} series agree '
                     f'({pct_agree:.0f}%)</span>, '
                     f'{n_differ} differ across {len(all_comparable)} graph group(s).</div>')

    # Cross-reference from CSV
    if csv_data and csv_data.get("cross_reference"):
        cross_ref_list = csv_data["cross_reference"]
        inconsistent = []
        if isinstance(cross_ref_list, list):
            for c in cross_ref_list:
                if isinstance(c, dict) and not c.get("consistent", True):
                    inconsistent.append(c)
        if inconsistent:
            h.append("<h3>Cross-DWI Significance Inconsistencies</h3>")
            h.append("<table><thead><tr><th>Metric</th><th>Timepoint</th>"
                     "<th>Significant In</th><th>Not Significant In</th></tr></thead><tbody>")
            for c in inconsistent:
                sig_in = c.get("significant_in", [])
                nsig_in = c.get("not_significant_in", [])
                sig_str = ", ".join(sig_in) if isinstance(sig_in, list) else "-"
                nsig_str = ", ".join(nsig_in) if isinstance(nsig_in, list) else "-"
                h.append(f"<tr><td>{_esc(str(c.get('metric', '')))}</td>"
                         f"<td>{_esc(str(c.get('timepoint', '')))}</td>"
                         f"<td>{_esc(sig_str)}</td>"
                         f"<td><span class='differ'>{_esc(nsig_str)}</span></td></tr>")
            h.append("</tbody></table>")
    return h


def _section_correlations(rows) -> list[str]:
    """Build the Notable Correlations section.

    Extracts correlation coefficients (r, rs, r-squared) from vision
    graph summaries and trend descriptions, filtering for |r| >= 0.3.
    Results are sorted by absolute correlation strength (descending).

    Parameters
    ----------
    rows : list[dict]
        Vision CSV rows.

    Returns
    -------
    list[str]
        HTML chunks for the correlations section.
    """
    # ── 5. Correlations ──
    h: list[str] = []
    if rows:
        h.append(_h2("Notable Correlations", "correlations"))
        corr_findings = []
        for r in rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            all_text = r["summary"] + " " + r["trends_json"]
            for rval, context in extract_correlations(all_text):
                if abs(rval) >= 0.3:
                    corr_findings.append((abs(rval), rval, dwi_type, base_name, context))

        if corr_findings:
            corr_findings.sort(reverse=True)

            # Split into strong (|r| >= 0.5) and moderate (0.3 <= |r| < 0.5)
            strong = [(a, r, d, g, c) for a, r, d, g, c in corr_findings if a >= 0.5]
            moderate = [(a, r, d, g, c) for a, r, d, g, c in corr_findings if a < 0.5]

            h.append(f"<p>{len(corr_findings)} correlation(s) with |r| \u2265 0.3 "
                     f"({len(strong)} strong, {len(moderate)} moderate).</p>")

            if strong:
                h.append("<h3>Strong Correlations (|r| \u2265 0.5)</h3>")
                h.append("<table><thead><tr><th>|r|</th><th>r</th><th>Strength</th>"
                         "<th>DWI</th><th>Graph</th><th>Context</th></tr></thead><tbody>")
                for _, rval, dwi, graph, ctx in strong:
                    ctx_short = _esc(ctx.replace("\n", " ")[:120])
                    h.append(f'<tr><td class="sig-2"><strong>{abs(rval):.2f}</strong></td>'
                             f"<td>{rval:+.2f}</td>"
                             f"<td>Strong</td><td>{_dwi_badge(dwi)}</td>"
                             f"<td>{_esc(graph)}</td><td><em>{ctx_short}</em></td></tr>")
                h.append("</tbody></table>")

            if moderate:
                h.append("<h3>Moderate Correlations (0.3 \u2264 |r| &lt; 0.5)</h3>")
                h.append("<table><thead><tr><th>|r|</th><th>r</th><th>Strength</th>"
                         "<th>DWI</th><th>Graph</th><th>Context</th></tr></thead><tbody>")
                for _, rval, dwi, graph, ctx in moderate:
                    ctx_short = _esc(ctx.replace("\n", " ")[:120])
                    h.append(f"<tr><td><strong>{abs(rval):.2f}</strong></td><td>{rval:+.2f}</td>"
                             f"<td>Moderate</td><td>{_dwi_badge(dwi)}</td>"
                             f"<td>{_esc(graph)}</td><td><em>{ctx_short}</em></td></tr>")
                h.append("</tbody></table>")
        else:
            h.append("<p>No notable correlations (|r| &ge; 0.3) found.</p>")
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

                    # Axis info
                    x_lbl = r.get("x_axis_label", "")
                    x_unit = r.get("x_axis_units", "")
                    y_lbl = r.get("y_axis_label", "")
                    y_unit = r.get("y_axis_units", "")
                    axis_parts = []
                    if x_lbl:
                        axis_parts.append(f"X: {x_lbl}" + (f" ({x_unit})" if x_unit else ""))
                    if y_lbl:
                        axis_parts.append(f"Y: {y_lbl}" + (f" ({y_unit})" if y_unit else ""))

                    h.append(f"<p>{_dwi_badge(dt)}")
                    if axis_parts:
                        h.append(f' <span class="axis-info">[{" | ".join(axis_parts)}]</span>')
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
                roc_rows_html.append(
                    f"<tr><td>{_dwi_badge(str(dwi_type))}</td>"
                    f"<td>{_esc(str(r.get('timepoint', '-')))}</td>"
                    f"<td><strong>{auc}</strong></td>"
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
                         "<th>95% CI</th><th>p</th></tr></thead><tbody>")
                for hr in sorted(hrs, key=lambda x: x.get("p", 1)):
                    ci = f"[{hr.get('ci_lo', 0):.3f}, {hr.get('ci_hi', 0):.3f}]"
                    p_val = hr.get("p", 1.0)
                    sig = _sig_tag(p_val)
                    cls = _sig_class(p_val)
                    cls_attr = f' class="{cls}"' if cls else ""
                    h.append(f"<tr><td><code>{_esc(str(hr.get('covariate', '')))}</code></td>"
                             f"<td>{hr.get('hr', 0):.3f}</td><td>{ci}</td>"
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


def _section_mat_data(mat_data) -> list[str]:
    """Build the Supplemental Data (MAT Files) section.

    Displays two sub-sections from parsed MAT-file JSON:
    1. **Dosimetry** -- mean D95 and V50 values for ADC and D sub-volumes.
    2. **Core Method Comparison** -- truncated mean-Dice heatmap matrix
       showing agreement between tumor core delineation methods.

    Parameters
    ----------
    mat_data : dict
        Mapping of DWI type to parsed MAT metrics dict (from
        ``parsed_mat_metrics_{dwi}.json``).

    Returns
    -------
    list[str]
        HTML chunks for the supplemental data section.
    """
    # ── 7.5. Core Method & Dosimetry ──
    h: list[str] = []
    if mat_data:
        h.append(_h2("Supplemental Data (MAT Files)", "supplemental"))

        has_dos = any("dosimetry" in d for d in mat_data.values())
        if has_dos:
            h.append("<h3>Dosimetry (Target Coverage)</h3>")
            h.append("<table><thead><tr><th>DWI</th><th>D95 ADC (Gy)</th><th>V50 ADC (%)</th>"
                     "<th>D95 D (Gy)</th><th>V50 D (%)</th></tr></thead><tbody>")
            for dt in DWI_TYPES:
                if dt not in mat_data or "dosimetry" not in mat_data[dt]:
                    continue
                dos = mat_data[dt]["dosimetry"]
                if not dos:
                    continue
                h.append(
                    f"<tr><td>{_dwi_badge(dt)}</td>"
                    f"<td>{dos.get('d95_adc_mean', 'N/A'):.2f}</td>"
                    f"<td>{dos.get('v50_adc_mean', 'N/A') * 100:.1f}%</td>"
                    f"<td>{dos.get('d95_d_mean', 'N/A'):.2f}</td>"
                    f"<td>{dos.get('v50_d_mean', 'N/A') * 100:.1f}%</td></tr>"
                )
            h.append("</tbody></table>")

        has_core = any("core_method" in d for d in mat_data.values())
        if has_core:
            h.append("<h3>Core Method Comparison (Mean Dice)</h3>")
            for dt in DWI_TYPES:
                if dt not in mat_data or "core_method" not in mat_data[dt]:
                    continue
                core = mat_data[dt]["core_method"]
                if not core or not core.get("methods"):
                    continue
                h.append(f"<h4>{_dwi_badge(dt)} \u2014 {len(core['methods'])} methods</h4>")
                methods = core["methods"]
                matrix = core["mean_dice_matrix"]
                n = len(methods)
                h.append('<div style="overflow-x:auto">')
                h.append("<table><thead><tr><th>Method</th>")
                for m in methods:
                    h.append(f"<th>{_esc(m)}</th>")
                h.append("</tr></thead><tbody>")
                for i in range(n):
                    h.append(f"<tr><td><strong>{_esc(methods[i])}</strong></td>")
                    for j in range(n):
                        val = matrix[i][j] if i < len(matrix) and j < len(matrix[i]) else None
                        if isinstance(val, (int, float)):
                            # Colour-code: high Dice = green-ish
                            if i != j:
                                style = f" style=\"color: {'var(--green)' if val >= 0.7 else ('var(--amber)' if val >= 0.5 else 'var(--red)')}; font-weight:600\""
                            else:
                                style = ""
                            h.append(f"<td{style}>{val:.2f}</td>")
                        else:
                            h.append("<td>-</td>")
                    h.append("</tr>")
                h.append("</tbody></table>")
                h.append("</div>")

                # Summary statistics for core method agreement
                off_diag = []
                for i in range(n):
                    for j in range(i + 1, n):
                        if i < len(matrix) and j < len(matrix[i]):
                            val = matrix[i][j]
                            if isinstance(val, (int, float)) and val > 0:
                                off_diag.append(val)
                if off_diag:
                    avg_dice = sum(off_diag) / len(off_diag)
                    min_dice = min(off_diag)
                    max_dice = max(off_diag)
                    h.append(f'<div class="info-box">Mean pairwise Dice: <strong>{avg_dice:.3f}</strong> '
                             f'(range: {min_dice:.3f}\u2013{max_dice:.3f} across {len(off_diag)} pairs)</div>')
    return h


def _section_methods(dwi_types_present, mat_data, log_data) -> list[str]:
    """Build the Methods section describing statistical methodology.

    Provides a publication-ready description of all statistical and
    analytical methods used in the pipeline, including IVIM modelling,
    group comparisons, multiple comparison correction, survival analysis,
    and predictive modelling.

    Parameters
    ----------
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    mat_data : dict
        Parsed MAT file metrics.
    log_data : dict or None
        Parsed log metrics.

    Returns
    -------
    list[str]
        HTML chunks for the methods section.
    """
    h: list[str] = []
    h.append(_h2("Statistical Methods", "methods"))
    h.append('<div class="methods-box">')

    # ── DWI Acquisition & Processing ──
    h.append("<h3>DWI Acquisition and Processing</h3>")
    h.append(
        "<p>Diffusion-weighted images were acquired and processed using three "
        "complementary strategies: <strong>Standard</strong> (conventional DWI), "
        f"<strong>DnCNN</strong> (deep learning denoised using a convolutional neural "
        f"network){_cite('dncnn')}, and <strong>IVIMnet</strong> (deep learning IVIM "
        "parameter estimation). "
        "For each strategy, apparent diffusion coefficient (ADC) maps were computed via "
        "mono-exponential fitting, and intravoxel incoherent motion (IVIM) parameters "
        "\u2014 true diffusion coefficient (<em>D</em>), perfusion fraction (<em>f</em>), "
        "and pseudo-diffusion coefficient (<em>D*</em>) \u2014 were estimated using "
        f"segmented and Bayesian fitting approaches{_cite('ivim')}.</p>"
    )

    # ── Tumour Delineation ──
    h.append("<h3>Tumour Sub-volume Delineation</h3>")
    h.append(
        "<p>Tumour core sub-volumes were identified using configurable delineation "
        "methods (default: ADC thresholding). Eleven methods were compared pairwise "
        f"using Dice similarity coefficient{_cite('dice')} and Hausdorff distance, "
        "including threshold-based "
        "(ADC, D, D\u00b7f intersection), clustering-based (Otsu, GMM, k-means, spectral), "
        "region-based (region growing, active contours), percentile-based, and functional "
        "diffusion map (fDM) approaches.</p>"
    )

    # ── Group Comparisons ──
    h.append("<h3>Group Comparisons</h3>")
    h.append(
        "<p>Differences between treatment outcome groups (Local Failure vs Local Control) "
        f"were assessed using the <strong>Wilcoxon rank-sum test</strong> "
        f"(Mann\u2013Whitney U){_cite('wilcoxon')}, "
        "a non-parametric test appropriate for small sample sizes and non-normally distributed "
        "DWI-derived biomarkers. Tests were performed independently at each imaging timepoint.</p>"
    )

    # ── GLME ──
    h.append("<h3>Mixed-Effects Modelling</h3>")
    h.append(
        "<p><strong>Generalised linear mixed-effects models (GLME)</strong> were used to test "
        "for time\u00d7outcome interaction effects, with patient as a random intercept to account "
        "for repeated measures. This approach tests whether the trajectory of each DWI metric "
        "differs significantly between outcome groups over time, while accounting for "
        "within-patient correlation.</p>"
    )

    # ── Multiple Comparisons ──
    h.append("<h3>Multiple Comparison Correction</h3>")
    h.append(
        "<p>To control the false discovery rate across the large number of metrics tested, "
        f"the <strong>Benjamini\u2013Hochberg (BH) procedure</strong>{_cite('bh_fdr')} "
        "was applied. Each "
        "metric\u2019s p-value was compared to an individually adjusted significance "
        "threshold (\u03b1<sub>adj</sub> = 0.05 \u00d7 rank / total tests), rather than a "
        "fixed \u03b1 = 0.05. This controls the expected proportion of false discoveries "
        "among rejected hypotheses at 5%.</p>"
    )

    # ── Predictive Modelling ──
    h.append("<h3>Predictive Modelling</h3>")
    pred_text_parts = [
        f"<p><strong>Elastic-net regularised logistic regression</strong>{_cite('elastic_net')} "
        "(mixing parameter "
        "\u03b1 = 0.5) was used for binary outcome prediction at each timepoint. "
        "The optimal regularisation parameter (\u03bb) was selected via 5-fold "
        "cross-validation with patient-stratified folds to prevent data leakage.",
    ]
    # Check if LOOCV was used
    if log_data:
        for dt in dwi_types_present:
            if dt in log_data:
                roc = log_data[dt].get("stats_predictive", {}).get("roc_analyses", [])
                if roc:
                    pred_text_parts.append(
                        " Discriminative performance was evaluated using leave-one-out "
                        "cross-validation (LOOCV) to generate unbiased out-of-fold risk "
                        "scores, reported as area under the receiver operating "
                        "characteristic curve (AUC)."
                    )
                    break
    pred_text_parts.append(
        " Feature collinearity was addressed by pruning highly correlated "
        "features (|r| > 0.8) prior to model fitting, retaining the feature with "
        "higher univariate AUC.</p>"
    )
    h.append("".join(pred_text_parts))

    # ── Survival Analysis ──
    h.append("<h3>Survival Analysis</h3>")
    ipcw_used = False
    if log_data:
        for dt in dwi_types_present:
            if dt in log_data and log_data[dt].get("survival", {}).get("ipcw"):
                ipcw_used = True
                break
    surv_text = (
        f"<p><strong>Cause-specific Cox proportional hazards models</strong>"
        f"{_cite('cox_ph')} were "
        "used to estimate hazard ratios (HR) with 95% confidence intervals for "
        "DWI-derived covariates. To account for competing risks (non-tumour-related "
        "mortality), "
    )
    if ipcw_used:
        surv_text += (
            f"<strong>inverse probability of censoring weighting (IPCW)</strong>"
            f"{_cite('ipcw')} "
            "was applied to adjust for informative censoring bias. "
        )
    else:
        surv_text += "competing-risk patients were excluded from the analysis. "
    surv_text += (
        "Model significance was assessed using the global likelihood ratio test (LRT). "
        "Where separation or convergence issues arose, Firth\u2019s penalised likelihood "
        f"method was used as a bias-reduction technique{_cite('firth')}.</p>"
    )
    h.append(surv_text)

    # ── Software ──
    h.append("<h3>Software and Reproducibility</h3>")
    h.append(
        "<p>All pipeline computations were performed in MATLAB (R2021a+) with the "
        "Statistics and Machine Learning Toolbox and Image Processing Toolbox. "
        "Post-hoc analysis and report generation used Python 3.12+. "
        "Parallel processing was limited to 2 workers with deterministic checkpointing "
        "for reproducibility. All analyses used patient-stratified splits with explicit "
        "temporal leakage prevention in imputation, scaling, and cross-validation.</p>"
    )

    h.append("</div>")
    return h


def _section_effect_sizes(log_data, dwi_types_present, csv_data) -> list[str]:
    """Build the Effect Size and Sample Size Reporting section.

    Reports effect sizes alongside p-values to distinguish statistical
    from clinical significance, following best practices for biomedical
    publication.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    csv_data : dict or None
        Parsed pipeline CSV exports.

    Returns
    -------
    list[str]
        HTML chunks for the effect size section.
    """
    h: list[str] = []
    h.append(_h2("Effect Size Analysis", "effect-sizes"))
    h.append(
        f'<p class="meta">Effect sizes provide a measure of practical significance '
        f'independent of sample size{_cite("cohen_d")}. Hazard ratios are reported '
        f'with 95% confidence '
        f'intervals; HR > 1 indicates increased risk of the endpoint. The log(HR) is '
        f'used as an approximate effect size metric (|log(HR)| \u2265 0.5 \u2248 '
        f'medium effect).</p>'
    )

    has_data = False

    if log_data:
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            sv = log_data[dwi_type].get("survival", {})
            hrs = sv.get("hazard_ratios", [])
            if not hrs:
                continue

            has_data = True
            h.append(f"<h3>{_dwi_badge(dwi_type)} Hazard Ratio Effect Sizes</h3>")
            h.append("<table>")
            h.append(_table_caption(
                f"Hazard Ratio Effect Sizes ({dwi_type})",
                f"Effect size classification per Cohen (1988): "
                f"small (|log HR| < 0.5), medium (0.5\u20130.8), large (\u2265 0.8)."))
            h.append("<thead><tr>"
                     "<th>Covariate</th><th>HR</th><th>95% CI</th>"
                     "<th>p-value</th><th>log(HR)</th><th>Effect</th>"
                     "<th>Forest Plot</th>"
                     "</tr></thead><tbody>")

            import math
            for hr in sorted(hrs, key=lambda x: x.get("p", 1)):
                hr_val = hr.get("hr", 1.0)
                ci_lo = hr.get("ci_lo", hr_val)
                ci_hi = hr.get("ci_hi", hr_val)
                p_val = hr.get("p", 1.0)
                log_hr = math.log(hr_val) if hr_val > 0 else 0
                eff_cls = _effect_size_class(log_hr)
                eff_lbl = _effect_size_label(log_hr)
                sig = _sig_tag(p_val)
                p_cls = _sig_class(p_val)
                p_attr = f' class="{p_cls}"' if p_cls else ""

                forest = _forest_plot_cell(hr_val, ci_lo, ci_hi, p_val)

                h.append(
                    f"<tr><td><code>{_esc(str(hr.get('covariate', '')))}</code></td>"
                    f"<td><strong>{hr_val:.3f}</strong></td>"
                    f"<td class=\"ci-text\">[{ci_lo:.3f}, {ci_hi:.3f}]</td>"
                    f"<td{p_attr}>{p_val:.4f} {_esc(sig)}</td>"
                    f"<td>{log_hr:+.3f}</td>"
                    f'<td class="{eff_cls}">{_esc(eff_lbl)}</td>'
                    f"<td>{forest}</td></tr>"
                )
            h.append("</tbody></table>")

            # CI width commentary
            narrow_ci = [hr for hr in hrs if hr.get("ci_hi", 1) - hr.get("ci_lo", 1) < 0.5]
            wide_ci = [hr for hr in hrs if hr.get("ci_hi", 1) - hr.get("ci_lo", 1) >= 1.5]
            if wide_ci:
                h.append(f'<div class="warn-box">{len(wide_ci)} covariate(s) have wide '
                         f'confidence intervals (CI width \u2265 1.5), suggesting imprecise '
                         f'effect estimates likely due to limited sample size.</div>')
            if narrow_ci:
                h.append(f'<div class="info-box">{len(narrow_ci)} covariate(s) have narrow '
                         f'CIs (width < 0.5), indicating precise effect estimation.</div>')

        # ── ROC/AUC effect size interpretation ──
        auc_data = []
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            roc = log_data[dwi_type].get("stats_predictive", {}).get("roc_analyses", [])
            for r_item in roc:
                auc_val = r_item.get("auc", 0)
                if auc_val > 0:
                    auc_data.append((dwi_type, r_item.get("timepoint", "?"), auc_val))

        if auc_data:
            has_data = True
            h.append("<h3>Discriminative Performance Interpretation</h3>")
            h.append(
                f'<p class="meta">AUC interpretation{_cite("hosmer_lemeshow")}: '
                f'0.5 = no discrimination, '
                f'0.6\u20130.7 = poor, 0.7\u20130.8 = acceptable, '
                f'0.8\u20130.9 = excellent, > 0.9 = outstanding.</p>'
            )
            h.append("<table>")
            h.append(_table_caption(
                "Discriminative Performance by Timepoint",
                "Classification per Hosmer & Lemeshow (2000)."))
            h.append("<thead><tr><th>DWI</th><th>Timepoint</th><th>AUC</th>"
                     "<th>Discrimination</th></tr></thead><tbody>")
            for dwi, tp, auc_val in sorted(auc_data, key=lambda x: -x[2]):
                if auc_val >= 0.9:
                    disc = "Outstanding"
                    cls = "effect-lg"
                elif auc_val >= 0.8:
                    disc = "Excellent"
                    cls = "effect-lg"
                elif auc_val >= 0.7:
                    disc = "Acceptable"
                    cls = "effect-md"
                elif auc_val >= 0.6:
                    disc = "Poor"
                    cls = "effect-sm"
                else:
                    disc = "No discrimination"
                    cls = "effect-sm"
                h.append(
                    f"<tr><td>{_dwi_badge(dwi)}</td>"
                    f"<td>{_esc(str(tp))}</td>"
                    f"<td><strong>{auc_val:.3f}</strong></td>"
                    f'<td class="{cls}">{_esc(disc)}</td></tr>'
                )
            h.append("</tbody></table>")

    if not has_data:
        h.append('<p class="meta">No effect size data available (requires survival or predictive log data).</p>')

    return h


def _section_multiple_comparisons(log_data, dwi_types_present, csv_data) -> list[str]:
    """Build the Multiple Comparisons Correction Summary section.

    Provides a comprehensive overview of raw vs adjusted p-values
    across all metrics and DWI types, showing the impact of FDR
    correction.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    csv_data : dict or None
        Parsed pipeline CSV exports.

    Returns
    -------
    list[str]
        HTML chunks for the multiple comparisons section.
    """
    h: list[str] = []
    h.append(_h2("Multiple Comparison Correction Summary", "mult-comp"))

    has_data = False

    if log_data:
        for dwi_type in dwi_types_present:
            if dwi_type not in log_data:
                continue
            sc = log_data[dwi_type].get("stats_comparisons", {})
            glme_details = sc.get("glme_details", [])
            if not glme_details:
                continue

            has_data = True
            n_total = len(glme_details)
            n_raw_sig = len([g for g in glme_details if g["p"] < 0.05])
            n_fdr_sig = len([g for g in glme_details if g["p"] < g["adj_alpha"]])

            h.append(f"<h3>{_dwi_badge(dwi_type)} \u2014 BH-FDR Correction Impact</h3>")
            h.append('<div class="stat-grid">')
            h.append(_stat_card("Total Tests", str(n_total)))
            h.append(_stat_card("Raw Significant", f"{n_raw_sig}/{n_total}",
                                f"at \u03b1 = 0.05"))
            h.append(_stat_card("FDR Significant", f"{n_fdr_sig}/{n_total}",
                                f"at adjusted \u03b1"))
            lost = n_raw_sig - n_fdr_sig
            if lost > 0:
                h.append(_stat_card("Lost to Correction", str(lost),
                                    "would be false discoveries"))
            h.append("</div>")

            # Detailed table: raw p-value vs adjusted threshold
            h.append("<table>")
            h.append(_table_caption(
                f"BH-FDR Correction Detail ({dwi_type})",
                f"Raw vs adjusted significance for {n_total} metrics."))
            h.append("<thead><tr>"
                     "<th>Metric</th><th>Raw p-value</th>"
                     "<th>BH Adj. \u03b1</th><th>Raw Sig?</th>"
                     "<th>FDR Sig?</th><th>Status</th>"
                     "</tr></thead><tbody>")
            for g in sorted(glme_details, key=lambda x: x["p"]):
                raw_sig = g["p"] < 0.05
                fdr_sig = g["p"] < g["adj_alpha"]
                raw_mark = "\u2713" if raw_sig else "\u2717"
                fdr_mark = "\u2713" if fdr_sig else "\u2717"

                if fdr_sig:
                    status = "Confirmed"
                    status_cls = "agree"
                elif raw_sig and not fdr_sig:
                    status = "Rejected by FDR"
                    status_cls = "differ"
                else:
                    status = "Not significant"
                    status_cls = ""

                p_cls = _sig_class(g["p"]) if fdr_sig else ""
                p_attr = f' class="{p_cls}"' if p_cls else ""
                s_attr = f' class="{status_cls}"' if status_cls else ""
                h.append(
                    f"<tr><td><code>{_esc(g['metric'])}</code></td>"
                    f"<td{p_attr}>{g['p']:.4f}</td>"
                    f"<td>{g['adj_alpha']:.4f}</td>"
                    f"<td>{raw_mark}</td><td>{fdr_mark}</td>"
                    f"<td{s_attr}><strong>{_esc(status)}</strong></td></tr>"
                )
            h.append("</tbody></table>")

    # Also show FDR global from CSV
    if csv_data and csv_data.get("fdr_global"):
        fdr_global = csv_data["fdr_global"]
        n_total_fdr = sum(len(v) for v in fdr_global.values())
        if n_total_fdr > 0:
            has_data = True
            h.append("<h3>Global FDR-Surviving Metrics</h3>")
            h.append(f"<p><strong>{n_total_fdr}</strong> metric(s) survived global "
                     f"Benjamini\u2013Hochberg correction across the full feature "
                     f"set, distributed as follows:</p>")
            h.append("<ul>")
            for dt in DWI_TYPES:
                if dt in fdr_global and fdr_global[dt]:
                    h.append(f"<li>{_dwi_badge(dt)}: {len(fdr_global[dt])} metric(s)</li>")
            h.append("</ul>")

    if not has_data:
        h.append('<p class="meta">No multiple comparison data available.</p>')

    return h


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
    h.append("</ul>")
    h.append("</div>")
    diagnostics_found = True

    return h


def _section_limitations(log_data, dwi_types_present, mat_data) -> list[str]:
    """Build the Study Limitations section.

    Generates a contextual limitations discussion based on the actual
    data characteristics observed (sample size, missingness, etc.).

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
        HTML chunks for the limitations section.
    """
    h: list[str] = []
    h.append(_h2("Limitations", "limitations"))

    limitations = []

    # Sample size
    n_patients = 0
    if mat_data:
        for dt in DWI_TYPES:
            if dt in mat_data and "longitudinal" in mat_data[dt]:
                n = mat_data[dt]["longitudinal"].get("num_patients", 0)
                if n > n_patients:
                    n_patients = n

    if n_patients > 0 and n_patients < 50:
        limitations.append(
            f"<strong>Small sample size (n = {n_patients}):</strong> The limited "
            f"cohort size restricts statistical power, increases the risk of "
            f"overfitting in predictive models, and limits generalisability. "
            f"Effect size estimates may be imprecise, as reflected by wide "
            f"confidence intervals in Cox regression."
        )
    elif n_patients > 0:
        limitations.append(
            f"<strong>Moderate sample size (n = {n_patients}):</strong> While "
            f"adequate for exploratory analyses, the cohort size may limit "
            f"power for detecting small effect sizes and reduces the "
            f"reliability of subgroup analyses."
        )

    limitations.append(
        "<strong>Single-institution cohort:</strong> All data were acquired at a "
        "single institution, which may introduce scanner-specific bias, "
        "protocol-dependent effects, and limit external validity. Multi-centre "
        "validation is needed to confirm generalisability."
    )

    limitations.append(
        "<strong>Retrospective design:</strong> This analysis is retrospective "
        "in nature. Unmeasured confounders (performance status, genetic subtypes, "
        "concurrent systemic therapy variations) may influence outcomes and were "
        "not controlled for in the current analysis."
    )

    # Check for missing data issues
    if log_data:
        for dt in dwi_types_present:
            if dt in log_data:
                bl = log_data[dt].get("baseline", {})
                exc = bl.get("baseline_exclusion")
                if exc and exc.get("n_excluded", 0) > 0:
                    limitations.append(
                        f"<strong>Missing data:</strong> "
                        f"{exc['n_excluded']}/{exc['n_total']} patients were "
                        f"excluded due to missing baseline data. If missingness "
                        f"is non-random (e.g., sicker patients less likely to "
                        f"complete baseline imaging), results may be biased towards "
                        f"a healthier sub-population."
                    )
                    break

    limitations.append(
        "<strong>DWI-specific limitations:</strong> IVIM parameter estimation "
        "is sensitive to the choice of b-values, number of signal averages, "
        "and fitting algorithm. DnCNN denoising may alter the noise distribution "
        "in ways that affect downstream parameter estimation, and IVIMnet "
        "predictions depend on the training set composition."
    )

    limitations.append(
        "<strong>Tumour delineation:</strong> GTV contours were propagated using "
        "deformable image registration, which may introduce geometric errors, "
        "particularly in regions of large anatomical deformation (e.g., due to "
        "bowel gas motion or tumour shrinkage)."
    )

    limitations.append(
        "<strong>Multiple comparisons:</strong> Despite BH-FDR correction, "
        "the large number of metrics tested across multiple timepoints and "
        "DWI types increases the cumulative risk of spurious findings. "
        "Results should be interpreted as hypothesis-generating rather "
        "than confirmatory."
    )

    h.append('<ul class="limitation-list">')
    for lim in limitations:
        h.append(f"<li>{lim}</li>")
    h.append("</ul>")

    return h


def _section_conclusions(log_data, dwi_types_present, csv_data, mat_data, groups) -> list[str]:
    """Build the Conclusions section.

    Synthesises key findings from all data sources into a structured
    conclusions paragraph suitable for publication.

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
        HTML chunks for the conclusions section.
    """
    h: list[str] = []
    h.append(_h2("Conclusions", "conclusions"))
    h.append('<div class="conclusion-box">')

    findings = []

    # 1. Key significant biomarkers
    n_fdr = 0
    if csv_data and csv_data.get("fdr_global"):
        n_fdr = sum(len(v) for v in csv_data["fdr_global"].values())
    total_glme_sig = 0
    if log_data:
        for dt in dwi_types_present:
            if dt in log_data:
                sc = log_data[dt].get("stats_comparisons", {})
                total_glme_sig += len([g for g in sc.get("glme_details", []) if g["p"] < g["adj_alpha"]])

    if total_glme_sig > 0 or n_fdr > 0:
        parts = []
        if total_glme_sig > 0:
            parts.append(f"{total_glme_sig} metric(s) demonstrating significant "
                         f"time\u00d7outcome interaction effects in GLME models")
        if n_fdr > 0:
            parts.append(f"{n_fdr} metric(s) surviving global FDR correction")
        findings.append(
            f"DWI-derived biomarkers show statistically significant associations "
            f"with treatment outcome, with {' and '.join(parts)}."
        )

    # 2. Predictive performance
    best_auc = 0
    best_type = ""
    best_tp = ""
    if log_data:
        for dt in dwi_types_present:
            if dt in log_data:
                roc = log_data[dt].get("stats_predictive", {}).get("roc_analyses", [])
                for r_item in roc:
                    a = r_item.get("auc", 0)
                    if a > best_auc:
                        best_auc = a
                        best_type = dt
                        best_tp = r_item.get("timepoint", "")

    if best_auc > 0:
        if best_auc >= 0.8:
            disc = "excellent"
        elif best_auc >= 0.7:
            disc = "acceptable"
        else:
            disc = "limited"
        findings.append(
            f"Elastic-net regularised logistic regression achieved {disc} "
            f"discriminative performance (AUC = {best_auc:.3f}, "
            f"{best_type} at {best_tp}), supporting the potential of "
            f"DWI biomarkers for early treatment response prediction."
        )

    # 3. Cox PH
    sig_covs = []
    if log_data:
        for dt in dwi_types_present:
            if dt in log_data:
                hrs = log_data[dt].get("survival", {}).get("hazard_ratios", [])
                for hr_item in hrs:
                    if hr_item.get("p", 1) < 0.05:
                        sig_covs.append((dt, hr_item.get("covariate", "?"),
                                         hr_item.get("hr", 1), hr_item.get("p", 1)))
    if sig_covs:
        cov_list = ", ".join(f"{c[1]} (HR={c[2]:.2f}, p={c[3]:.3f})" for c in sig_covs[:3])
        findings.append(
            f"Cause-specific Cox regression identified significant prognostic "
            f"covariates: {cov_list}."
        )

    # 4. Cross-DWI agreement
    if groups:
        n_agree = 0
        n_total = 0
        for base_name, dwi_dict in groups.items():
            real = [t for t in dwi_dict if t != "Root"]
            if len(real) < 2:
                continue
            all_trends_dict: dict[str, list] = {}
            for dt in DWI_TYPES:
                if dt in dwi_dict:
                    try:
                        all_trends_dict[dt] = json.loads(str(dwi_dict[dt].get("trends_json", "[]")))
                    except Exception:
                        pass
            if len(all_trends_dict) >= 2:
                all_series: set[str] = set()
                for trends in all_trends_dict.values():
                    for t in trends:
                        if isinstance(t, dict):
                            all_series.add(t.get("series") or "overall")
                for series in all_series:
                    directions: dict[str, str] = {}
                    for dt_key, trends in all_trends_dict.items():
                        for t in trends:
                            if isinstance(t, dict) and (t.get("series") or "overall") == series:
                                directions[dt_key] = str(t.get("direction", ""))
                    if len(directions) >= 2:
                        n_total += 1
                        if len(set(directions.values())) == 1:
                            n_agree += 1
        if n_total > 0:
            pct = 100 * n_agree / n_total
            findings.append(
                f"Cross-DWI-type trend agreement is {pct:.0f}% ({n_agree}/{n_total} "
                f"series), {'supporting' if pct >= 70 else 'suggesting limited'} "
                f"robustness of findings across acquisition strategies."
            )

    # 5. Hypothesis direction
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
        if d_cons == "increasing" and f_cons == "decreasing":
            findings.append(
                "Longitudinal trends show the canonical response pattern "
                "(increasing diffusion, decreasing perfusion), consistent "
                "with therapy-induced cellular necrosis and vascular regression."
            )

    if findings:
        h.append("<ol>")
        for f in findings:
            h.append(f"<li>{_esc(f)}</li>")
        h.append("</ol>")
    else:
        h.append("<p>Detailed findings are presented in the sections above. "
                 "Overall, the analysis demonstrates the feasibility of "
                 "multi-parametric DWI analysis for treatment response assessment "
                 "in pancreatic cancer.</p>")

    # Clinical significance statement
    h.append(
        '<div class="summary-box" style="border-left-color: var(--green);">'
        "<p><strong>Clinical Significance:</strong> "
        "These findings suggest that longitudinal DWI biomarkers may enable "
        "non-invasive, early identification of patients at risk for local "
        "failure during radiotherapy. If validated prospectively, this could "
        "support adaptive treatment strategies\u2014such as dose escalation "
        "to resistant sub-volumes or early intensification of systemic "
        "therapy\u2014within the existing fractionation schedule. The "
        "cross-DWI-type analysis demonstrates that key findings are robust "
        "to the choice of post-processing strategy, increasing confidence "
        "in their clinical applicability.</p></div>"
    )

    h.append(
        "<p><strong>Future directions:</strong> Prospective validation in an "
        "independent multi-centre cohort is warranted. Investigation of "
        "radiomics and texture features, time-dependent covariates in Cox "
        "models, and deep learning\u2013based outcome prediction may further "
        "improve prognostic accuracy. Additionally, integration of "
        "circulating tumour DNA (ctDNA) and PET imaging biomarkers with "
        "DWI-derived metrics could enable multi-modal prediction models with "
        "improved sensitivity for early treatment failure detection.</p>"
    )

    h.append("</div>")
    return h


def _section_appendix(rows) -> list[str]:
    """Build the Appendix: All Graphs section.

    Lists every analysed graph in a detailed table with columns for
    DWI type, graph type, title, axis details, trend tags, and a
    collapsible full summary.

    Parameters
    ----------
    rows : list[dict]
        Vision CSV rows (may be empty, in which case the section is skipped).

    Returns
    -------
    list[str]
        HTML chunks for the appendix section.
    """
    # ── 8. Appendix ──
    h: list[str] = []
    if rows:
        h.append(_h2("Appendix: All Graphs", "appendix"))
        h.append(f"<p>{len(rows)} graphs analysed. Expand each row for full details.</p>")
        h.append("<table><thead><tr>"
                 "<th>#</th><th>DWI</th><th>Type</th><th>Graph</th>"
                 "<th>Title</th><th>Axes</th><th>Trends</th><th>Summary</th>"
                 "</tr></thead><tbody>")
        for i, r in enumerate(rows, 1):
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            graph_title = r.get("graph_title", "") or ""
            graph_type = r.get("graph_type", "")

            # Axes summary
            x_lbl = r.get("x_axis_label", "")
            x_unit = r.get("x_axis_units", "")
            x_min = r.get("x_axis_range_min", "")
            x_max = r.get("x_axis_range_max", "")
            y_lbl = r.get("y_axis_label", "")
            y_unit = r.get("y_axis_units", "")
            y_min = r.get("y_axis_range_min", "")
            y_max = r.get("y_axis_range_max", "")
            c_lbl = r.get("color_axis_label", "")
            c_unit = r.get("color_axis_units", "")

            axis_lines = []
            if x_lbl:
                xr = f" [{x_min}\u2013{x_max}]" if x_min and x_max else ""
                xu = f" ({x_unit})" if x_unit else ""
                axis_lines.append(f"X: {x_lbl}{xu}{xr}")
            if y_lbl:
                yr = f" [{y_min}\u2013{y_max}]" if y_min and y_max else ""
                yu = f" ({y_unit})" if y_unit else ""
                axis_lines.append(f"Y: {y_lbl}{yu}{yr}")
            if c_lbl:
                cu = f" ({c_unit})" if c_unit else ""
                axis_lines.append(f"C: {c_lbl}{cu}")
            axes_cell = '<br>'.join(_esc(a) for a in axis_lines) if axis_lines else "\u2014"

            # Trends cell
            trends_str = r.get("trends_json", "[]") or "[]"
            try:
                trends_list = json.loads(trends_str)
            except Exception:
                trends_list = []
            trends_cell = ""
            if isinstance(trends_list, list) and trends_list:
                tags = []
                for t in trends_list:
                    if isinstance(t, dict):
                        direction = t.get("direction", "")
                        series = t.get("series", "")
                        label = f"{series}: {direction}" if series else direction
                        tags.append(_trend_tag(label))
                trends_cell = "".join(tags)

            # Summary: short preview + collapsible full
            summary = r.get("summary", "") or ""
            summary_preview = _esc(summary[:100].replace("\n", " "))
            if len(summary) > 100:
                summary_cell = (
                    f"{summary_preview}\u2026"
                    f"<details><summary>more</summary>"
                    f'<p class="full-summary">{_esc(summary)}</p>'
                    f"</details>"
                )
            else:
                summary_cell = summary_preview

            title_display = _esc(graph_title) if graph_title else "\u2014"
            h.append(
                f"<tr>"
                f"<td>{i}</td>"
                f"<td>{_dwi_badge(dwi_type)}</td>"
                f"<td>{_esc(graph_type)}</td>"
                f"<td>{_esc(base_name)}</td>"
                f'<td class="axis-info">{title_display}</td>'
                f'<td class="axis-info">{axes_cell}</td>'
                f"<td>{trends_cell}</td>"
                f"<td>{summary_cell}</td>"
                f"</tr>"
            )
        h.append("</tbody></table>")
    return h
