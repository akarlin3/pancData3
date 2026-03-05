# Data Pipeline Analytical Audit Report

**Date:** 2026-03-05
**Scope:** All core/, utils/ modules in the pancData3 MATLAB pipeline
**Auditor:** Claude Code (automated audit)

---

## Critical Findings (Impact: Could Produce Incorrect Conclusions)

### 1. IPCW Frequency Approximation Introduces Bias in Cox Model
**File:** `core/metrics_survival.m`, lines 247-267

The IPCW weights are converted to integer frequencies via `max(1, round(ipcw_weights * 100))` because MATLAB's `coxphfit` doesn't support continuous probability weights. The SE correction (`SE * sqrt(ipcw_scale)`) on line 266 is an approximation that assumes the information matrix scales linearly with N, which is not generally true for the Cox partial likelihood.

**Impact:** Standard errors and p-values for hazard ratios are approximate. Confidence intervals may be too narrow or too wide depending on the risk-set structure. This is acknowledged in the code comments but is a genuine analytical limitation.

**Recommendation:** Implement a custom weighted partial likelihood or use a sandwich variance estimator for proper inference.

---

### 2. Collinearity Filtering Uses Different Feature Sets Across LOOCV Folds
**File:** `core/metrics_stats_predictive.m`, lines 166-171 (outer CV) and 278-280 (LOOCV)

Each fold computes its own `keep_fold = filter_collinear_features(X_tr_imp, y_tr_fold)`, which can select different features per fold. The final model uses a consensus mask (`keep_fold_counts > n_folds_en / 2`), but the LOOCV risk scores are generated with per-fold collinearity masks. This means:

- Different LOOCV folds may predict on different feature subspaces
- The risk scores are not strictly comparable across folds
- The final model's feature set (consensus mask) differs from what any individual LOOCV fold used

**Impact:** The out-of-fold risk scores may have inconsistent feature spaces across patients, making the aggregate ROC AUC unreliable as a performance estimate.

**Recommendation:** Apply a fixed collinearity mask (e.g., from the full training set, or the consensus mask computed before LOOCV) consistently across all LOOCV folds.

---

### 3. Elastic Net Lambda Path Mismatch Across CV Folds
**File:** `core/metrics_stats_predictive.m`, lines 174-194

The lambda path is determined by the first fold (`cv_i == 1`) and then forced on subsequent folds via `'Lambda', common_Lambda`. However, `lassoglm` may select different features per fold due to the collinearity filter producing different feature subsets. When `X_tr_kept` has a different number of columns across folds (because `keep_fold` varies), the lambda values from fold 1 may not be appropriate for the regularization path of later folds with different dimensionality.

**Impact:** Cross-validation deviance estimates may be unreliable, leading to suboptimal lambda selection and potentially overfitting or underfitting.

---

### 4. Exponential Decay Imputation Applies to All Missing Covariates Uniformly
**File:** `utils/build_td_panel.m`, lines 121-133

The exponential decay imputation uses the same half-life for all covariates (ADC, D, f, D*), but these parameters have different biological decay dynamics:
- ADC/D: primarily reflect cellularity changes (slower response)
- f: perfusion fraction (faster vascular response)
- D*: pseudo-diffusion (fastest response, highest noise)

**Impact:** The single half-life assumption may over- or under-impute certain parameters, biasing their hazard ratios in the Cox model. The sensitivity analysis (half-life grid) partially addresses this but applies the same half-life to all covariates simultaneously rather than independently.

---

### 5. Exponential Decay Imputation Decays Toward Zero
**File:** `utils/build_td_panel.m`, line 132

```matlab
cov_row(decay_mask) = orig_X(decay_mask) .* exp(-lambda_decay * dt_per_feat(decay_mask));
```

The decay model assumes all features decay toward zero. For DWI parameters like ADC, the pre-treatment value is a positive physiological quantity. After treatment, values may return toward baseline, not zero. For patients with missing late-timepoint scans and long follow-up, imputed diffusion parameters will be artificially low, potentially creating a spurious association between low parameter values and events (since patients with events are more likely to have missing late scans).

**Impact:** Systematic downward bias in imputed late-timepoint covariates, particularly for parameters with high physiological floor values (ADC, D).

**Recommendation:** Decay toward patient-specific baseline rather than zero: `baseline + (last_observed - baseline) * exp(-lambda * dt)`.

---

### 6. Zero-Variance Feature Handling Inconsistent Between Scaling Modes
**File:** `utils/scale_td_panel.m`, lines 123-124 vs lines 50-55

In `baseline` mode (line 55), zero-variance features are correctly zeroed out. In `per_week` mode (line 123-124), zero-variance features have `sd_col` clamped to 1 instead, producing `(x - mu) / 1`, which preserves arbitrary centered values. Similarly, single-observation weeks (line 127-128) get `sd_col = 1`.

```matlab
% per_week mode (line 123-124):
if sd_col == 0
    sd_col = 1;  % zero-variance: clamp to 1
end
% vs baseline mode (line 50-55):
if sig == 0 || isnan(sig)
    X_td_scaled(:, fi) = 0;  % zero-variance: zero out
end
```

**Impact:** Zero-variance features may appear to carry discriminative information in per_week mode (the default), artificially inflating their importance in downstream Cox regression.

**Recommendation:** Use the baseline mode's approach — zero out constant features — in both modes.

---

### 7. Landmark Day Selection May Be Wrong for Non-Standard Protocols
**File:** `core/metrics_survival.m`, lines 120-122

```matlab
n_on_tx = max(length(td_scan_days) - 1, 1);
landmark_idx = n_on_tx;
landmark_day = td_scan_days(landmark_idx);
```

This assumes the last element of `td_scan_days` is the post-treatment scan. If `td_scan_days` has a different structure (e.g., two post-treatment scans), the landmark would be set incorrectly, potentially including post-treatment intervals in the analysis set while claiming to exclude them.

**Impact:** Could introduce immortal time bias if scan_days doesn't follow the assumed structure.

---

### 8. DVH Metrics Silently Overwritten Across Repeat Scans
**File:** `core/load_dwi_data.m`, lines 432-437

```matlab
pat_dmean_gtvp(1,fi) = scan_result.dmean_gtvp;
pat_d95_gtvp(1,fi) = scan_result.d95_gtvp;
pat_v50gy_gtvp(1,fi) = scan_result.v50gy_gtvp;
```

These arrays are indexed by `(1, fi)` with no repeat dimension `rpi`. When multiple repeat scans exist for a fraction, only the **last** repeat's dose metrics survive. Contrast with `pat_adc_mean(1,fi,rpi)` on line 426 which correctly includes the repeat index.

**Impact:** For fractions with multiple repeat acquisitions, DVH metrics (mean dose, D95, V50Gy) are silently lost for all but the final repeat. Dose-parameter correlation analyses for multi-repeat fractions use incorrect dose values.

**Recommendation:** Either index by `rpi` (like the imaging metrics) or explicitly average across repeats.

---

### 9. Dose Resampling Called with Empty/Invalid Inputs
**File:** `core/process_single_scan.m`, lines 87-110

Two related issues:

(a) When dose DICOM directory exists but contains no `.dcm` files, `rtdosefile = ''` (line 102) but execution falls through to `sample_rtdose_on_image(b0list, '')` (line 110) with no guard. This will either crash or produce a corrupt dose NIfTI that persists on disk and is loaded in subsequent runs.

(b) `b0list = cell(1)` (line 87) creates `{[]}` (a cell containing one empty element), not an empty cell. If no b=0 DICOM images are found, `b0list` remains `{[]}` — an invalid reference image list passed to `sample_rtdose_on_image`.

**Impact:** Silent dose map corruption when dose directory is present but empty, or when DWI data lacks b=0 images. The corrupted NIfTI is cached on disk, persisting across pipeline reruns.

**Recommendation:** Add guard `if b0count > 0 && ~isempty(rtdosefile)` before line 110, and initialize `b0list = cell(0,1)`.

---

### 10. ROC Analysis Includes Competing-Risk Patients as Negatives
**File:** `core/metrics_stats_predictive.m`, lines 432-436

```matlab
labels = lf_group; % 0 = LC, 1 = LF
valid_roc = ~isnan(risk_scores_all_target) & ~isnan(labels);
[roc_X, roc_Y, roc_T, roc_AUC] = perfcurve(labels(valid_roc), risk_scores_all_target(valid_roc), 1);
```

The elastic net correctly excludes competing-risk patients (lf==2) from model training (line 104-105), but the ROC analysis uses the original `lf_group` which still contains 2s. `perfcurve` with `PositiveClass=1` treats both lf==0 and lf==2 as negatives. Competing-risk patients who died before local failure could be observed are lumped with genuine local-control patients.

**Impact:** AUC is inflated because competing-risk patients may be easily separable from LF patients, creating an artificially favorable ROC curve. All derived metrics (sensitivity, specificity, Youden threshold) are affected.

**Recommendation:** Filter to `labels <= 1` before calling `perfcurve`.

---

### 11. Warning Messages Silently Swallowed Throughout Pipeline
**File:** `run_dwi_pipeline.m`, lines 270, 353, 431, 483, 515, 572, 620, 676, 716

```matlab
[warn_msg, warn_id] = lastwarn('');
if ~isempty(warn_msg) && log_fid > 0
    fprintf(log_fid, ...);
end
```

`lastwarn('')` simultaneously sets the warning to `''` and returns that empty string. So `warn_msg` is always `''`, the `if` never fires, and warnings from every pipeline step are silently discarded instead of being written to `error.log`. The correct pattern is:

```matlab
[warn_msg, warn_id] = lastwarn;   % read
lastwarn('');                       % then clear
```

**Impact:** All pipeline warnings are lost. Researchers cannot review non-fatal issues after a run completes.

---

### 12. Parameter Maps Always Show Standard DWI Type
**File:** `core/plot_parameter_maps.m`, line 32

```matlab
s = data_vectors_gtvp(j, 1, 1);  % 3rd index = 1 = Standard
```

When the pipeline runs dnCNN or IVIMnet, parameter maps still display Standard ADC data. The function is called once outside the DWI type loop in `visualize_results.m`, so it never adapts to the current pipeline type. Output appears in the dnCNN/IVIMnet subfolder but shows Standard data.

**Impact:** Misleading visualization — parameter maps labeled as dnCNN/IVIMnet actually show Standard pipeline results.

---

## Moderate Findings (Impact: May Affect Statistical Power or Interpretation)

### 13. Sub-Volume Morphological Filtering May Remove Clinically Relevant Scattered Voxels
**File:** `utils/calculate_subvolume_metrics.m`, lines 49-51

```matlab
vol_3d = imclose(imopen(vol_3d, se), se);
vol_3d = bwareaopen(vol_3d, min_cc_voxels);  % min_cc_voxels = 10
```

The morphological open-close and connected component filter removes small clusters (<10 voxels). In pancreatic tumors, treatment-resistant voxels can be scattered (heterogeneous response). Removing them biases the dose-response analysis toward contiguous resistant regions only.

**Impact:** Moderate. D95 and V50 may overestimate dose coverage to resistant tissue by excluding scattered resistant voxels.

---

### 14. f_delta Uses Absolute Change While Other Parameters Use Percent Change
**File:** `core/metrics_baseline.m`, line 373

```matlab
f_delta = (f_abs - f_abs(:,1));
```

This computes absolute change for f (perfusion fraction, range 0-1), while ADC, D, and D* use percent change. This is intentional and documented, but in the elastic net model (`metrics_stats_predictive.m` line 52), f_delta is included alongside percent-change features without explicit re-scaling. Since `lassoglm` uses `'Standardize', true`, this is handled automatically, but the `feat_names_lasso` labels it as `'f_Pct'` which is misleading:

**File:** `core/metrics_stats_predictive.m`, line 62: `'f_Pct'` should be `'f_Delta'`.

**Impact:** No analytical error (standardization handles it), but the mislabeling could lead to misinterpretation of results.

---

### 15. In-Sample KNN Imputation for Final Elastic Net Model
**File:** `core/metrics_stats_predictive.m`, line 210

```matlab
X_clean_all = knn_impute_train_test(X_impute, [], 5, id_list_impute);
```

The final model uses the entire dataset as both reference and target for KNN imputation. While the code comment explains this is intentional (LOOCV provides unbiased estimates), the final model coefficients are fitted on optimistically imputed data. If the model is deployed for new patients, the imputation would need to be re-done using only the training cohort as reference.

**Impact:** Minor for internal validation, but important for external deployment.

---

### 16. Competing Risk Patients Excluded Inconsistently Across Modules
**Files:** Multiple

- `metrics_stats_predictive.m` line 104: sets lf==2 to NaN (exclusion)
- `metrics_stats_comparisons.m` line 82: filters g <= 1 (exclusion)
- `metrics_survival.m` line 155: sets event==2 to 0 (CSH censoring)
- `metrics_baseline.m` line 162: codes as lf=2

The handling is correct for each model type (exclusion for binomial, CSH censoring for Cox), but a researcher unaware of this might not realize the effective sample sizes differ across analyses.

**Impact:** Not an error, but a documentation/interpretation concern.

---

### 17. Scan Day Assumption Creates Immortal Time Bias Risk
**File:** `core/metrics_survival.m`, lines 44-53, `utils/build_td_panel.m` line 54

The default scan days `[0, 5, 10, 15, 20, 90]` assume equally spaced 5-day intervals. In practice, scans may be irregularly spaced. The warning on line 47-52 correctly flags this, but if the researcher doesn't provide actual scan days, the counting-process intervals will be misaligned with reality.

**Impact:** Potentially introduces immortal time bias or misattributes covariates to wrong time windows.

---

## Minor Findings

### 18. CLAUDE.md Documents Wrong Default for `adc_thresh`
**File:** `CLAUDE.md` config section

CLAUDE.md shows `"adc_thresh": 0.00115` but both `config.example.json` and `parse_config.m` use `0.001`. The `0.00115` value is actually the `high_adc_thresh` default.

**Impact:** Could mislead researchers relying on CLAUDE.md for configuration guidance.

---

### 19. No Validation of `adc_thresh <= high_adc_thresh` Ordering
**File:** `utils/parse_config.m`, lines 28-33

The comments state `adc_thresh` "must be <= high_adc_thresh" but this constraint is never enforced. A misconfigured `adc_thresh > high_adc_thresh` would produce empty or overlapping sub-volume selections.

**Impact:** Silent misconfiguration could produce empty subvolumes and all-NaN dose metrics.

---

### 20. Data-Adaptive Epsilon for Percent Change Reduces Reproducibility
**File:** `core/metrics_baseline.m`, lines 352-361

The epsilon threshold for percent change denominators is `max(1e-8, 0.01 * iqr(baseline))`. Since IQR varies with cohort composition, different cohorts (or the same cohort with different outlier exclusions) will use different epsilon values, potentially filtering different patients.

**Impact:** Minor reproducibility concern across cohorts.

---

## Summary

| Severity | Count | Key Concerns |
|----------|-------|-------------|
| Critical | 12 | IPCW bias, LOOCV feature inconsistency, lambda path mismatch, decay-to-zero, zero-variance scaling, uniform decay half-life, landmark selection, DVH overwrite, dose resampling crash, ROC competing-risk inflation, lastwarn bug, parameter map DWI type |
| Moderate | 5 | Morphological filtering, f_delta mislabeling, in-sample imputation, competing risk documentation, immortal time bias |
| Minor | 3 | CLAUDE.md threshold docs, threshold ordering validation, epsilon reproducibility |

### Top Priority Fixes
1. **Fix ROC competing-risk inflation** (`metrics_stats_predictive.m:432`) — Filter `labels <= 1` before `perfcurve`
2. **Fix lastwarn bug** (`run_dwi_pipeline.m`) — Split `[w,id] = lastwarn; lastwarn('');` into two calls
3. **Guard dose resampling inputs** (`process_single_scan.m`) — Add `if b0count > 0 && ~isempty(rtdosefile)` guard; init `b0list = cell(0,1)`
4. **Fix DVH repeat indexing** (`load_dwi_data.m`) — Add `rpi` dimension to DVH metric arrays
5. **Fix parameter map DWI type** (`plot_parameter_maps.m`) — Pass dtype and index correct pipeline vectors
6. **Fix zero-variance scaling inconsistency** (`scale_td_panel.m`) — Zero out constant features in both modes
7. **Fix LOOCV feature consistency** — Apply a fixed collinearity mask across all folds
8. **Fix lambda path issue** — Recompute lambda path when feature dimensionality changes
9. **Fix decay-to-zero imputation** — Decay toward patient baseline, not zero
10. **Rename `f_Pct` to `f_Delta`** — Prevents misinterpretation of elastic net coefficients
11. **Per-parameter decay half-lives** — Allow different biological decay rates for ADC/D vs f vs D*
12. **Add threshold ordering validation** — Enforce `adc_thresh <= high_adc_thresh` in `parse_config.m`
