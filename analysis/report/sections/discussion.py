"""Report sections: discussion group.

The ``_section_limitations`` and ``_section_conclusions`` builders live in
sibling modules (``discussion_limitations``, ``discussion_conclusions``) and
are re-exported here for backward compatibility.
"""

from __future__ import annotations

from shared import (  # type: ignore
    DWI_TYPES,
)
from report.report_formatters import (  # type: ignore
    _cite,
    _h2,
)

# Re-export sibling discussion sections for backward compatibility.
from report.sections.discussion_limitations import _section_limitations  # noqa: F401  # type: ignore
from report.sections.discussion_conclusions import _section_conclusions  # noqa: F401  # type: ignore


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
