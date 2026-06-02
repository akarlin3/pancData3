"""Report section: study limitations (discussion group)."""

from __future__ import annotations

from shared import (  # type: ignore
    DWI_TYPES,
)
from report.report_formatters import (  # type: ignore
    _h2,
)


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
