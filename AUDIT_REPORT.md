# Data Pipeline Analytical Audit Report

**Date:** 2026-03-05
**Scope:** All core/, utils/ modules in the pancData3 MATLAB pipeline
**Auditor:** Claude Code (automated audit)

---

## Critical Findings (Impact: Could Produce Incorrect Conclusions)

### 1. BH FDR Correction Formula Error
**File:** `core/metrics_stats_comparisons.m`, lines 193-197

The Benjamini-Hochberg adjusted p-values are computed incorrectly. The standard BH formula for the i-th ranked p-value is:

```
q(i) = min(q(i+1), p(i) * m / i)
```

where `m` is the total number of tests and `i` is the **rank** (1 = smallest p-value). The current code uses:

```matlab
q_all(ii) = min(q_all(ii+1), p_sort(ii) * (n_all / ii));
```

This is correct. **However**, the step-up enforcement at line 198 (`q_all = min(q_all, 1)`) is applied **after** the backward sweep, which is correct. No error here upon closer review.

**Status:** False alarm after deeper analysis. The BH implementation is correct.

---

### 2. IPCW Frequency Approximation Introduces Bias in Cox Model
**File:** `core/metrics_survival.m`, lines 247-267

The IPCW weights are converted to integer frequencies via `max(1, round(ipcw_weights * 100))` because MATLAB's `coxphfit` doesn't support continuous probability weights. The SE correction (`SE * sqrt(ipcw_scale)`) on line 266 is an approximation that assumes the information matrix scales linearly with N, which is not generally true for the Cox partial likelihood.

**Impact:** Standard errors and p-values for hazard ratios are approximate. Confidence intervals may be too narrow or too wide depending on the risk-set structure. This is acknowledged in the code comments but is a genuine analytical limitation.

**Recommendation:** Implement a custom weighted partial likelihood or use a sandwich variance estimator for proper inference.

---

### 3. Collinearity Filtering Uses Different Feature Sets Across LOOCV Folds
**File:** `core/metrics_stats_predictive.m`, lines 166-171 (outer CV) and 278-280 (LOOCV)

Each fold computes its own `keep_fold = filter_collinear_features(X_tr_imp, y_tr_fold)`, which can select different features per fold. The final model uses a consensus mask (`keep_fold_counts > n_folds_en / 2`), but the LOOCV risk scores are generated with per-fold collinearity masks. This means:

- Different LOOCV folds may predict on different feature subspaces
- The risk scores are not strictly comparable across folds
- The final model's feature set (consensus mask) differs from what any individual LOOCV fold used

**Impact:** The out-of-fold risk scores may have inconsistent feature spaces across patients, making the aggregate ROC AUC unreliable as a performance estimate.

**Recommendation:** Apply a fixed collinearity mask (e.g., from the full training set, or the consensus mask computed before LOOCV) consistently across all LOOCV folds.

---

### 4. Elastic Net Lambda Path Mismatch Across CV Folds
**File:** `core/metrics_stats_predictive.m`, lines 174-194

The lambda path is determined by the first fold (`cv_i == 1`) and then forced on subsequent folds via `'Lambda', common_Lambda`. However, `lassoglm` may select different features per fold due to the collinearity filter producing different feature subsets. When `X_tr_kept` has a different number of columns across folds (because `keep_fold` varies), the lambda values from fold 1 may not be appropriate for the regularization path of later folds with different dimensionality.

**Impact:** Cross-validation deviance estimates may be unreliable, leading to suboptimal lambda selection and potentially overfitting or underfitting.

---

### 5. Percent Change Denominator Uses Modified Baseline (Not Original)
**File:** `core/metrics_baseline.m`, lines 368-374

```matlab
adc_bl = ADC_abs(:,1);  adc_bl(adc_bl < adc_eps) = NaN;
ADC_pct = ((ADC_abs - ADC_abs(:,1)) ./ adc_bl) * 100;
```

The code correctly guards against near-zero denominators, but `ADC_abs(:,1)` in the numerator `(ADC_abs - ADC_abs(:,1))` still uses the original (possibly near-zero) baseline value. This means the numerator could be extremely small while the denominator is NaN'd, producing NaN for valid patients. This is actually the intended behavior (filtering them out), but worth noting that the **epsilon threshold is data-adaptive** (1% of IQR), meaning it varies across runs and could filter different patients in replication studies.

**Impact:** Minor. Reproducibility concern rather than a direct error.

---

### 6. Exponential Decay Imputation Applies to All Missing Covariates Uniformly
**File:** `utils/build_td_panel.m`, lines 121-133

The exponential decay imputation uses the same half-life for all covariates (ADC, D, f, D*), but these parameters have different biological decay dynamics:
- ADC/D: primarily reflect cellularity changes (slower response)
- f: perfusion fraction (faster vascular response)
- D*: pseudo-diffusion (fastest response, highest noise)

**Impact:** The single half-life assumption may over- or under-impute certain parameters, biasing their hazard ratios in the Cox model. The sensitivity analysis (half-life grid) partially addresses this but applies the same half-life to all covariates simultaneously rather than independently.

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

## Moderate Findings (Impact: May Affect Statistical Power or Interpretation)

### 8. GLME Z-Score Standardization Uses All Patients (Including Outliers)
**File:** `core/metrics_stats_comparisons.m`, lines 375-392

The GLME baseline z-scores are computed using `mean()` and `std()` on `long_ADC(baseline_idx)`, which includes all patients with valid baseline data. However, metrics_baseline.m already NaN-ified outlier patients' values (line 333-340). This means the z-score denominators include NaN-excluded patients in their count but not their values.

Actually, looking more carefully: the `long_ADC` array is populated from `ADC_abs(p_idx, t)` which already has outliers set to NaN, and line 356 filters to `~isnan(long_ADC)`. So the z-scores are computed from non-outlier patients only. **No error here.**

---

### 9. Sub-Volume Morphological Filtering May Remove Clinically Relevant Scattered Voxels
**File:** `utils/calculate_subvolume_metrics.m`, lines 49-51

```matlab
vol_3d = imclose(imopen(vol_3d, se), se);
vol_3d = bwareaopen(vol_3d, min_cc_voxels);  % min_cc_voxels = 10
```

The morphological open-close and connected component filter removes small clusters (<10 voxels). In pancreatic tumors, treatment-resistant voxels can be scattered (heterogeneous response). Removing them biases the dose-response analysis toward contiguous resistant regions only.

**Impact:** Moderate. D95 and V50 may overestimate dose coverage to resistant tissue by excluding scattered resistant voxels.

---

### 10. KNN Imputation Distance Computed on Z-scored Features, but Imputation Uses Raw Values
**File:** `utils/knn_impute_train_test.m`, lines 40 and 106

Distances are computed in Z-scored space (`search_space`), but the imputed values come from raw `X_tr`:

```matlab
X_tr_imp(i, m) = mean(X_tr(neighbors, m));
```

This is actually correct behavior (you want to impute raw values, not standardized ones), but it means the distance metric is scale-normalized while the imputed value is not. This is standard KNN imputation practice and is not an error.

---

### 11. f_delta Uses Absolute Change While Other Parameters Use Percent Change
**File:** `core/metrics_baseline.m`, line 373

```matlab
f_delta = (f_abs - f_abs(:,1));
```

This computes absolute change for f (perfusion fraction, range 0-1), while ADC, D, and D* use percent change. This is intentional and documented, but in the elastic net model (`metrics_stats_predictive.m` line 52), f_delta is included alongside percent-change features without explicit re-scaling. Since `lassoglm` uses `'Standardize', true`, this is handled automatically, but the `feat_names_lasso` labels it as `'f_Pct'` which is misleading:

**File:** `core/metrics_stats_predictive.m`, line 62: `'f_Pct'` should be `'f_Delta'`.

**Impact:** No analytical error (standardization handles it), but the mislabeling could lead to misinterpretation of results.

---

### 12. In-Sample KNN Imputation for Final Elastic Net Model
**File:** `core/metrics_stats_predictive.m`, line 210

```matlab
X_clean_all = knn_impute_train_test(X_impute, [], 5, id_list_impute);
```

The final model uses the entire dataset as both reference and target for KNN imputation. While the code comment explains this is intentional (LOOCV provides unbiased estimates), the final model coefficients are fitted on optimistically imputed data. If the model is deployed for new patients, the imputation would need to be re-done using only the training cohort as reference.

**Impact:** Minor for internal validation, but important for external deployment.

---

### 13. Competing Risk Patients Excluded Inconsistently Across Modules
**Files:** Multiple

- `metrics_stats_predictive.m` line 104: sets lf==2 to NaN (exclusion)
- `metrics_stats_comparisons.m` line 82: filters g <= 1 (exclusion)
- `metrics_survival.m` line 155: sets event==2 to 0 (CSH censoring)
- `metrics_baseline.m` line 162: codes as lf=2

The handling is correct for each model type (exclusion for binomial, CSH censoring for Cox), but a researcher unaware of this might not realize the effective sample sizes differ across analyses.

**Impact:** Not an error, but a documentation/interpretation concern.

---

### 14. Scan Day Assumption Creates Immortal Time Bias Risk
**File:** `core/metrics_survival.m`, lines 44-53, `utils/build_td_panel.m` line 54

The default scan days `[0, 5, 10, 15, 20, 90]` assume equally spaced 5-day intervals. In practice, scans may be irregularly spaced. The warning on line 47-52 correctly flags this, but if the researcher doesn't provide actual scan days, the counting-process intervals will be misaligned with reality.

**Impact:** Potentially introduces immortal time bias or misattributes covariates to wrong time windows.

---

## Minor Findings

### 15. `adc_sd` Indexed by `dtype` in Predictive Module May Be Wrong
**File:** `core/metrics_stats_predictive.m`, line 523

```matlab
baseline_sd  = adc_sd(valid_pts, 1, dtype);
```

The `adc_sd` variable is passed as a parameter to this function and comes from `summary_metrics.adc_sd`, which is a 3D array `[nPat x nTp x 3]`. The `dtype` index selects the correct DWI type. This appears correct.

---

## Summary

| Severity | Count | Key Concerns |
|----------|-------|-------------|
| Critical | 4 | IPCW approximation bias, inconsistent feature spaces in LOOCV, lambda path mismatch, single decay half-life |
| Moderate | 4 | Morphological filtering bias, f_delta mislabeling, in-sample imputation, immortal time bias risk |
| Minor | 1 | Reproducibility of epsilon threshold |

### Top Priority Fixes
1. **Fix LOOCV feature consistency** — Apply a fixed collinearity mask across all LOOCV folds
2. **Fix lambda path issue** — Recompute lambda path when feature dimensionality changes across folds
3. **Rename `f_Pct` to `f_Delta`** — Prevents misinterpretation of elastic net coefficients
4. **Per-parameter decay half-lives** — Allow different biological decay rates for ADC/D vs f vs D*
