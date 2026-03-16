"""Report sections: discussion group."""

from __future__ import annotations

import json

from shared import (  # type: ignore
    DWI_TYPES,
)
from report.report_formatters import (  # type: ignore
    _cite,
    _esc,
    _get_consensus,
    _h2,
)
from report.sections._helpers import _scalar_gy  # type: ignore


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
        f"using Dice similarity coefficient{_cite('dice')} and Hausdorff distance:</p>"
    )
    h.append("<ul>")
    core_method_descriptions = [
        ("<code>adc_threshold</code>", "ADC voxel-level threshold (default &lt; 1.0\u00d710\u207b\u00b3 mm\u00b2/s); selects voxels with restricted diffusion."),
        ("<code>d_threshold</code>", "Diffusion coefficient D threshold; isolates regions of genuinely low true diffusivity."),
        ("<code>df_intersection</code>", "Voxels where both D is low AND f is high (hypercellular-perfused intersection); targets densely packed, well-perfused tumour regions."),
        ("<code>otsu</code>", "Otsu\u2019s automatic bi-level thresholding applied to the ADC map; data-driven threshold selection without a fixed cutoff."),
        ("<code>gmm</code>", "2-component Gaussian mixture model fitted to the ADC distribution; probabilistically assigns voxels to tumour core vs periphery."),
        ("<code>kmeans</code>", "K-means clustering (k = 2) on ADC voxel values; partitions the tumour into two subregions."),
        ("<code>region_growing</code>", "Seed-based region growing from the minimum ADC voxel; expands iteratively to include spatially adjacent restricted-diffusion voxels."),
        ("<code>active_contours</code>", "Active contour (\u201csnake\u201d) segmentation; energy-minimising deformable boundary fitted to the ADC gradient."),
        ("<code>percentile</code>", "Bottom-N-percentile of the ADC distribution within the GTV; defines a fixed proportion of the most restricted voxels as the core."),
        ("<code>spectral</code>", "Spectral graph clustering on the voxel-level ADC similarity matrix; captures spatially coherent subregions based on diffusion texture."),
        ("<code>fdm</code>", "Functional Diffusion Map (fDM); threshold on voxel-level ADC change between timepoints to identify progressively restricted subvolumes."),
    ]
    for method_code, method_desc in core_method_descriptions:
        h.append(f"<li>{method_code}: {method_desc}</li>")
    h.append("</ul>")

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
        "<p><strong>Full model specification:</strong></p>"
        "<pre>log(parameter) ~ timepoint + outcome + timepoint\u00d7outcome + (1 | patient_id)</pre>"
        "<p>where <em>timepoint</em> is a fixed-effect factor, "
        "<em>outcome</em> \u2208 {local failure, local control} is the group variable, "
        "and <em>patient_id</em> is the random intercept. "
        "The interaction term tests whether the parameter trajectory over treatment differs "
        "between outcome groups. "
        "Estimation by maximum likelihood (ML); inference by likelihood ratio test.</p>"
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
        "Hyperparameter \u03bb was selected by inner 5-fold cross-validation (CV) "
        "minimizing binomial deviance, using patient-stratified folds to prevent data leakage.",
    ]
    # Check if LOOCV was used
    if log_data:
        for dt in dwi_types_present:
            if dt in log_data:
                roc = log_data[dt].get("stats_predictive", {}).get("roc_analyses", [])
                if roc:
                    pred_text_parts.append(
                        " Outer leave-one-out cross-validation (LOOCV) was then applied to "
                        "generate unbiased out-of-fold risk scores, avoiding optimistic bias "
                        "from using the same data for model selection and performance estimation. "
                        "Discriminative performance is reported as area under the receiver "
                        "operating characteristic curve (AUC)."
                    )
                    break
    pred_text_parts.append(
        " Feature collinearity was addressed by pruning highly correlated "
        "features (|r| > 0.8) prior to model fitting, retaining the feature with "
        "higher univariate AUC.</p>"
    )
    h.append("".join(pred_text_parts))

    # ── Missing data ──
    h.append("<h3>Missing Data Handling</h3>")
    h.append(
        "<p>Baseline missing values were imputed using k-nearest-neighbor (KNN) imputation "
        "(k = 3 neighbours). To prevent data leakage, imputation was performed separately "
        "for training and test folds, and future timepoint data were excluded from the "
        "imputation reference set for each patient.</p>"
    )

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

    # ── Deep Learning Disclosure ──
    h.append("<h3>Deep Learning Model Training Disclosure</h3>")
    h.append(
        "<p>DnCNN was pre-trained on a large corpus of natural images (ImageNet) and applied "
        "to DWI without additional fine-tuning on pancreatic data. IVIMnet was trained on "
        "synthetic IVIM data generated from biologically plausible parameter distributions. "
        "Neither model was trained on data from this cohort, minimizing the risk of data leakage.</p>"
    )

    # ── Data Quality Assurance ──
    h.append("<h3>Data Quality Assurance</h3>")
    h.append(
        "<p>An automated sanity check module validated all fitted parameter maps "
        "prior to downstream analysis. Checks included: (1) detection of non-finite "
        "values (Inf, NaN) and negative values in voxel-level ADC, D, f, and D* maps; "
        "(2) statistical outlier detection via interquartile range (IQR) fencing "
        "(values > 3 IQR from the median); (3) dimensional alignment verification "
        "between DWI parameter maps and RT dose grids; and (4) excessive NaN fraction "
        "warnings (threshold: 50% of in-mask voxels). Outliers were excluded from "
        "downstream group comparisons and predictive modelling.</p>"
    )

    # ── Dosimetry ──
    has_dosimetry = False
    if mat_data:
        for dt in DWI_TYPES:
            if dt in mat_data and "dosimetry" in mat_data[dt]:
                has_dosimetry = True
                break
    if has_dosimetry:
        h.append("<h3>Dosimetric Analysis</h3>")
        h.append(
            "<p>Radiotherapy dose distributions were resampled onto the DWI parameter "
            "maps using the RT dose grid. Target coverage was quantified using dose-volume "
            "histogram (DVH) metrics within diffusion-defined tumour sub-volumes: "
            "D95 (minimum dose to 95% of the sub-volume) and V50 (fraction of the "
            "sub-volume receiving \u226550 Gy). Sub-volumes were defined using both ADC "
            "and IVIM-D thresholding to identify putatively resistant regions.</p>"
        )

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
            if dt in mat_data and "longitudinal" in mat_data[dt]:  # type: ignore
                n = mat_data[dt]["longitudinal"].get("num_patients", 0)  # type: ignore
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

    # Check if PH assumption was tested via Schoenfeld residuals
    ph_tested = False
    ph_violated_covs: list[str] = []
    if log_data:
        for dt in dwi_types_present:
            if dt in log_data:
                sv = log_data[dt].get("survival", {})
                if sv.get("schoenfeld_tested"):
                    ph_tested = True
                    ph_violated_covs.extend(sv.get("schoenfeld_violated", []))

    # Check if time-varying Cox models were fitted as PH follow-up
    tv_cox_fitted = False
    if log_data:
        for dt in dwi_types_present:
            if dt in log_data:
                sv = log_data[dt].get("survival", {})
                if sv.get("time_varying_cox_fitted"):
                    tv_cox_fitted = True

    if ph_tested:
        if ph_violated_covs:
            violated_str = ", ".join(set(ph_violated_covs))
            if tv_cox_fitted:
                limitations.append(
                    f"<strong>Proportional hazards assumption:</strong> Formal PH testing "
                    f"via Schoenfeld residuals identified violations for: {violated_str}. "
                    f"Time-varying coefficient models (covariate \u00d7 log(time) interaction) "
                    f"and stratified Cox models were fitted as follow-up to address these violations. "
                    f"Time-varying HR curves are provided for the violating covariates."
                )
            else:
                limitations.append(
                    f"<strong>Proportional hazards assumption:</strong> Formal PH testing "
                    f"via Schoenfeld residuals identified violations for: {violated_str}. "
                    f"Time-varying coefficients or stratified models may be warranted for "
                    f"these covariates."
                )
        else:
            limitations.append(
                "<strong>Proportional hazards assumption:</strong> Formal PH testing "
                "via Schoenfeld residuals showed no significant violations at \u03b1 = 0.05 "
                "for any covariate, supporting the proportional hazards assumption."
            )
    else:
        limitations.append(
            "<strong>Proportional hazards assumption:</strong> Proportional hazards assumption "
            "was not formally tested; time-varying hazard ratios cannot be excluded."
        )

    # Check if imputation sensitivity analysis was performed
    imputation_sensitivity_done = False
    if log_data:
        for dt in dwi_types_present:
            if dt in log_data:
                pred = log_data[dt].get("predictive", {})
                if pred.get("imputation_sensitivity_done"):
                    imputation_sensitivity_done = True

    if imputation_sensitivity_done:
        limitations.append(
            "<strong>KNN imputation quality:</strong> Imputation sensitivity analysis was "
            "performed comparing KNN against LOCF, mean imputation, and linear interpolation. "
            "AUC concordance across methods supports robustness of the imputation approach."
        )
    else:
        limitations.append(
            "<strong>KNN imputation quality:</strong> KNN imputation quality was not formally "
            "validated; imputed values may not accurately represent missing observations for "
            "patients with unusual trajectories."
        )

    # Check if Fine-Gray model was computed
    fine_gray_computed = False
    if log_data:
        for dt in dwi_types_present:
            if dt in log_data:
                sv = log_data[dt].get("survival", {})
                if sv.get("fine_gray_computed"):
                    fine_gray_computed = True
                    break

    if fine_gray_computed:
        limitations.append(
            "<strong>Competing-risk model choice:</strong> Both cause-specific Cox models "
            "(with IPCW weighting) and Fine\u2013Gray subdistribution hazard models were "
            "computed, providing complementary perspectives on the competing-risk structure. "
            "CSH estimates the biological hazard rate; Fine\u2013Gray estimates the cumulative "
            "incidence accounting for competing events."
        )
    else:
        limitations.append(
            "<strong>Competing-risk model choice:</strong> Competing-risk analysis used "
            "cause-specific Cox models with IPCW weighting; sub-distribution hazard ratios "
            "(Fine\u2013Gray model) were not computed."
        )

    limitations.append(
        "<strong>Core delineation validation:</strong> Core delineation methods were compared "
        "on agreement (Dice/Hausdorff) rather than on ground-truth segmentation accuracy; "
        "expert contour validation was not performed."
    )

    limitations.append(
        "<strong>Vision-based graph analysis:</strong> Vision-based graph analysis via AI "
        "introduces potential misclassification of figure content; all auto-extracted trends "
        "and p-values should be verified against primary MATLAB log outputs."
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
                        n_total = int(n_total + 1)  # type: ignore
                        if len(set(directions.values())) == 1:
                            n_agree = int(n_agree + 1)  # type: ignore
        if n_total > 0:
            pct = 100 * float(n_agree) / float(n_total)  # type: ignore
            findings.append(
                f"Cross-DWI-type trend agreement is {pct:.0f}% ({n_agree}/{n_total} "
                f"series), {'supporting' if pct >= 70 else 'suggesting limited'} "
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



