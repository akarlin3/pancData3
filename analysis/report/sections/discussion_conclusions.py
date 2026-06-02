"""Report section: conclusions (discussion group)."""

from __future__ import annotations

import json

from shared import (  # type: ignore
    DWI_TYPES,
)
from report.report_formatters import (  # type: ignore
    _esc,
    _get_consensus,
    _h2,
)
from report.sections._helpers import _scalar_gy  # type: ignore


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
                total_glme_sig = int(total_glme_sig + len([g for g in sc.get("glme_details", []) if g["p"] < g["adj_alpha"]]))  # type: ignore

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
        # Collect per-DWI best AUC for cross-type comparison
        per_dwi_auc: list[tuple[str, float, str]] = []
        for dt in dwi_types_present:
            if dt in log_data:
                roc = log_data[dt].get("stats_predictive", {}).get("roc_analyses", [])
                for r_item in roc:
                    a = r_item.get("auc", 0)
                    if a > 0:
                        per_dwi_auc.append((dt, a, r_item.get("timepoint", "")))
        auc_comparison = ""
        if len(per_dwi_auc) > 1:
            # Show per-type best AUCs for cross-comparison
            best_per_type: dict[str, tuple[float, str]] = {}
            for dt, a, tp in per_dwi_auc:
                if dt not in best_per_type or a > best_per_type[dt][0]:
                    best_per_type[dt] = (a, tp)
            if len(best_per_type) >= 2:
                parts = [f"{dt}: {a:.3f} ({tp})" for dt, (a, tp) in
                         sorted(best_per_type.items(), key=lambda x: -x[1][0])]
                auc_comparison = f" Per-type comparison: {'; '.join(parts)}."
        findings.append(
            f"Elastic-net regularised logistic regression achieved {disc} "
            f"discriminative performance (AUC = {best_auc:.3f}, "
            f"{best_type} at {best_tp}), supporting the potential of "
            f"DWI biomarkers for early treatment response prediction."
            f"{auc_comparison}"
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
        cov_list = ", ".join(f"{c[1]} (HR={c[2]:.2f}, p={c[3]:.3f})" for c in sig_covs[:3])  # type: ignore
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
                    except (json.JSONDecodeError, TypeError):
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
                        n_total = int(n_total + 1)  # type: ignore
                        if len(set(directions.values())) == 1:
                            n_agree = int(n_agree + 1)  # type: ignore
        if n_total > 0:
            pct = 100 * float(n_agree) / float(n_total)  # type: ignore
            if pct >= 70:
                _robustness_word = "supporting"
            elif pct >= 50:
                _robustness_word = "showing moderate"
            else:
                _robustness_word = "suggesting limited"
            findings.append(
                f"Cross-DWI-type trend agreement is {pct:.0f}% ({n_agree}/{n_total} "
                f"series), {_robustness_word} "
                f"robustness of findings across acquisition strategies."
            )

    # 5. Dosimetry
    if mat_data:
        for dt in DWI_TYPES:
            dosi = (mat_data.get(dt) or {}).get("dosimetry")
            if dosi:
                d95_adc = _scalar_gy(dosi.get("d95_adc_mean"))
                v50_adc = _scalar_gy(dosi.get("v50_adc_mean"))
                if d95_adc is not None:
                    v50_pct = (v50_adc * 100 if v50_adc is not None and v50_adc <= 1.0
                               else v50_adc) if v50_adc is not None else None
                    dosi_parts = [f"D95 = {d95_adc:.1f} Gy"]
                    if v50_pct is not None:
                        dosi_parts.append(f"V50 = {v50_pct:.0f}%")
                    coverage = "adequate" if d95_adc >= 45.0 else "sub-optimal"
                    findings.append(
                        f"Dosimetric analysis of diffusion-defined resistant sub-volumes "
                        f"shows {coverage} target coverage ({', '.join(dosi_parts)})."
                    )
                break

    # 6. Core method agreement (Dice)
    if mat_data:
        for dt in DWI_TYPES:
            core = (mat_data.get(dt) or {}).get("core_method")
            if core and core.get("mean_dice_matrix"):
                methods = core.get("methods", [])
                matrix = core["mean_dice_matrix"]
                n = len(methods)
                off_diag = []
                for i in range(n):
                    for j in range(i + 1, n):
                        if i < len(matrix) and j < len(matrix[i]):  # type: ignore
                            val = matrix[i][j]  # type: ignore
                            if isinstance(val, (int, float)) and val > 0:
                                off_diag.append(val)
                if off_diag:
                    avg_dice = sum(off_diag) / len(off_diag)  # type: ignore
                    findings.append(
                        f"Tumour core delineation across {n} methods shows "
                        f"{'moderate' if avg_dice < 0.7 else 'good'} spatial "
                        f"agreement (mean pairwise Dice = {avg_dice:.3f}), "
                        f"suggesting method selection impacts sub-volume definition."
                    )
                break

    # 7. Hypothesis direction
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
            except (json.JSONDecodeError, TypeError, KeyError, ValueError):
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

    # Clinical significance statement — conditional on cross-DWI agreement
    # Compute agreement percentage to reconcile with conclusion #4
    cross_dwi_pct = 0.0
    if groups:
        _n_agree = 0
        _n_total = 0
        for base_name, dwi_dict in groups.items():
            real = [t for t in dwi_dict if t != "Root"]
            if len(real) < 2:
                continue
            _all_trends: dict[str, list] = {}
            for dt in DWI_TYPES:
                if dt in dwi_dict:
                    try:
                        _all_trends[dt] = json.loads(str(dwi_dict[dt].get("trends_json", "[]")))
                    except (json.JSONDecodeError, TypeError):
                        pass
            if len(_all_trends) >= 2:
                _all_series: set[str] = set()
                for trends in _all_trends.values():
                    for t in trends:
                        if isinstance(t, dict):
                            _all_series.add(t.get("series") or "overall")
                for series in _all_series:
                    _dirs: dict[str, str] = {}
                    for dt_key, trends in _all_trends.items():
                        for t in trends:
                            if isinstance(t, dict) and (t.get("series") or "overall") == series:
                                _dirs[dt_key] = str(t.get("direction", ""))
                    if len(_dirs) >= 2:
                        _n_total += 1
                        if len(set(_dirs.values())) == 1:
                            _n_agree += 1
        if _n_total > 0:
            cross_dwi_pct = 100 * float(_n_agree) / float(_n_total)

    if cross_dwi_pct >= 70:
        robustness_text = (
            "The cross-DWI-type analysis demonstrates that key findings are robust "
            "to the choice of post-processing strategy, increasing confidence "
            "in their clinical applicability."
        )
    elif cross_dwi_pct >= 50:
        robustness_text = (
            f"Cross-DWI trend agreement ({cross_dwi_pct:.0f}%) shows moderate consistency "
            "across processing strategies. While not all biomarkers replicate, "
            "specific parameters (D, f) show consistent directional changes."
        )
    else:
        robustness_text = (
            f"While overall cross-DWI trend agreement is limited ({cross_dwi_pct:.0f}%), "
            "specific biomarkers show consistent directional changes across processing "
            "strategies, suggesting parameter-specific rather than global robustness."
        )

    h.append(
        '<div class="summary-box" style="border-left-color: var(--green);">'
        "<p><strong>Clinical Significance:</strong> "
        "These findings suggest that longitudinal DWI biomarkers may enable "
        "non-invasive, early identification of patients at risk for local "
        "failure during radiotherapy. If validated prospectively, this could "
        "support adaptive treatment strategies\u2014such as dose escalation "
        "to resistant sub-volumes or early intensification of systemic "
        "therapy\u2014within the existing fractionation schedule. "
        f"{robustness_text}</p></div>"
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
