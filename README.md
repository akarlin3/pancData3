# pancData

Analysis pipeline for pancreatic DWI (Diffusion-Weighted Imaging) data, including IVIM model fitting, ADC computation, and correlation with RT dose maps.

## Scripts

### `load_dwi_data_code/load_dwi_data_forAvery.m`
Main pipeline script. Loads DICOM DWI data, fits ADC and IVIM models, extracts dose metrics, computes longitudinal biomarkers, and performs statistical analyses (ANOVA, ROC, Cox regression).

### `load_dwi_data_code/sanity_checks.m`
Data validation script ("Understand the Data"). Run after loading the workspace produced by the main pipeline. Performs:
- **Verify Convergence** — checks all voxel-level ADC and IVIM fit outputs for `Inf`, `NaN`, and negative diffusion coefficients.
- **Identify Missingness & Outliers** — summarises `NaN` patterns across patients/fractions and flags biomarker values that are extreme outliers (>3 IQR from the cohort median).
- **Data Alignment** — confirms that RT Dose and DWI volumes are spatially consistent by comparing voxel counts and checking for excessive `NaN` fractions inside the GTV mask.

### `load_dwi_data_code/visualize_results.m`
Visualization script ("Visualizing It"). Run after loading the workspace produced by the main pipeline. Generates:
- **Parameter Maps overlaid on Anatomy** — ADC maps displayed on top of b=0 anatomical MRI slices with GTV contours for representative patients.
- **Distributions of Extracted Features** — histograms and box plots of baseline ADC, D, f, and D* grouped by clinical outcome (Local Control vs Local Failure).
- **Scatter Plots for Dose Correlation** — RT Dose (mean and D95) plotted against diffusion metrics with Spearman correlation and linear trend lines.

### `load_dwi_data_code/metrics.m`
Statistical analysis script ("Quantitative Metrics & Outcome Analysis"). Run after loading the workspace produced by the main pipeline. Performs:
- **Repeatability** — within-subject coefficient of variation (wCV) for ADC, D, f, and D* across repeat Fx1 scans.
- **Clinical Outcome Loading** — reads the master spreadsheet to extract local failure status and time-to-event data for each patient.
- **Longitudinal Visualization** — absolute and percent-change plots of diffusion biomarkers across fractions.
- **DVH Sub-volume Analysis** — dose–volume histogram metrics for DWI-defined sub-volumes.
- **Univariate Analysis** — ANOVA box plots comparing Local Control vs Local Failure at each fraction.
- **Multiple-Comparison Corrections** — Benjamini-Hochberg FDR (full metric sweep) and Holm-Bonferroni FWER control.
- **Feature Selection** — Elastic Net (α = 0.5) regularized logistic regression with 5-fold cross-validation and correlation pre-filtering (|r| > 0.8).
- **ROC Analysis** — Youden's J optimal cutoff and AUC computation.
- **Multivariable Logistic Regression** — with Firth penalized fallback for sparse data.
- **Survival Analysis** — Leave-One-Out cross-validated Kaplan-Meier curves stratified by imaging biomarker risk group.
- **Parametric Maps** — representative ADC map panels for responder vs non-responder patients.

### `load_dwi_data_code/test.m`
MATLAB test script for the statistical and analytical procedures in `metrics.m`. Run with:
```matlab
test
```
Covers:
- **Benjamini-Hochberg FDR** — q-value computation, monotonicity, capping at 1, and edge cases.
- **Holm-Bonferroni FWER** — threshold formula, step-down early-stop behavior, and boundary conditions.
- **Correlation Pre-filtering** — drop/keep logic for |r| > 0.8, including negative correlation and chain-dropping.
- **LOOCV Kaplan-Meier** — no data-leakage verification, training set size, and risk-assignment correctness.
- **Elastic Net** — confirms α = 0.5 is used (not pure LASSO).
- **All-Fractions Loop** — verifies the analysis iterates over all fractions (2 : nTp), not a hardcoded subset.
- **Targeted FDR Removal** — confirms post-hoc cherry-picked sections have been removed.
- **ADC Fitting via fit_adc_mono** — confirms `fit_adc_mono` is called instead of inline OLS.

### `load_dwi_data_code/dependencies/`
Utility functions used by the pipeline: NIfTI loaders, IVIM model fitting (segmented and Bayesian), ADC fitting, RT dose sampling, and DVH computation.
