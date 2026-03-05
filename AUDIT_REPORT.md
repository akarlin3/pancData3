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

## Moderate Findings (Impact: May Affect Statistical Power or Interpretation)

### 8. Sub-Volume Morphological Filtering May Remove Clinically Relevant Scattered Voxels
**File:** `utils/calculate_subvolume_metrics.m`, lines 49-51

```matlab
vol_3d = imclose(imopen(vol_3d, se), se);
vol_3d = bwareaopen(vol_3d, min_cc_voxels);  % min_cc_voxels = 10
```

The morphological open-close and connected component filter removes small clusters (<10 voxels). In pancreatic tumors, treatment-resistant voxels can be scattered (heterogeneous response). Removing them biases the dose-response analysis toward contiguous resistant regions only.

**Impact:** Moderate. D95 and V50 may overestimate dose coverage to resistant tissue by excluding scattered resistant voxels.

---

### 9. f_delta Uses Absolute Change While Other Parameters Use Percent Change
**File:** `core/metrics_baseline.m`, line 373

```matlab
f_delta = (f_abs - f_abs(:,1));
```

This computes absolute change for f (perfusion fraction, range 0-1), while ADC, D, and D* use percent change. This is intentional and documented, but in the elastic net model (`metrics_stats_predictive.m` line 52), f_delta is included alongside percent-change features without explicit re-scaling. Since `lassoglm` uses `'Standardize', true`, this is handled automatically, but the `feat_names_lasso` labels it as `'f_Pct'` which is misleading:

**File:** `core/metrics_stats_predictive.m`, line 62: `'f_Pct'` should be `'f_Delta'`.

**Impact:** No analytical error (standardization handles it), but the mislabeling could lead to misinterpretation of results.

---

### 10. In-Sample KNN Imputation for Final Elastic Net Model
**File:** `core/metrics_stats_predictive.m`, line 210

```matlab
X_clean_all = knn_impute_train_test(X_impute, [], 5, id_list_impute);
```

The final model uses the entire dataset as both reference and target for KNN imputation. While the code comment explains this is intentional (LOOCV provides unbiased estimates), the final model coefficients are fitted on optimistically imputed data. If the model is deployed for new patients, the imputation would need to be re-done using only the training cohort as reference.

**Impact:** Minor for internal validation, but important for external deployment.

---

### 11. Competing Risk Patients Excluded Inconsistently Across Modules
**Files:** Multiple

- `metrics_stats_predictive.m` line 104: sets lf==2 to NaN (exclusion)
- `metrics_stats_comparisons.m` line 82: filters g <= 1 (exclusion)
- `metrics_survival.m` line 155: sets event==2 to 0 (CSH censoring)
- `metrics_baseline.m` line 162: codes as lf=2

The handling is correct for each model type (exclusion for binomial, CSH censoring for Cox), but a researcher unaware of this might not realize the effective sample sizes differ across analyses.

**Impact:** Not an error, but a documentation/interpretation concern.

---

### 12. Scan Day Assumption Creates Immortal Time Bias Risk
**File:** `core/metrics_survival.m`, lines 44-53, `utils/build_td_panel.m` line 54

The default scan days `[0, 5, 10, 15, 20, 90]` assume equally spaced 5-day intervals. In practice, scans may be irregularly spaced. The warning on line 47-52 correctly flags this, but if the researcher doesn't provide actual scan days, the counting-process intervals will be misaligned with reality.

**Impact:** Potentially introduces immortal time bias or misattributes covariates to wrong time windows.

---

## Minor Findings

### 13. CLAUDE.md Documents Wrong Default for `adc_thresh`
**File:** `CLAUDE.md` config section

CLAUDE.md shows `"adc_thresh": 0.00115` but both `config.example.json` and `parse_config.m` use `0.001`. The `0.00115` value is actually the `high_adc_thresh` default.

**Impact:** Could mislead researchers relying on CLAUDE.md for configuration guidance.

---

### 14. No Validation of `adc_thresh <= high_adc_thresh` Ordering
**File:** `utils/parse_config.m`, lines 28-33

The comments state `adc_thresh` "must be <= high_adc_thresh" but this constraint is never enforced. A misconfigured `adc_thresh > high_adc_thresh` would produce empty or overlapping sub-volume selections.

**Impact:** Silent misconfiguration could produce empty subvolumes and all-NaN dose metrics.

---

### 15. Data-Adaptive Epsilon for Percent Change Reduces Reproducibility
**File:** `core/metrics_baseline.m`, lines 352-361

The epsilon threshold for percent change denominators is `max(1e-8, 0.01 * iqr(baseline))`. Since IQR varies with cohort composition, different cohorts (or the same cohort with different outlier exclusions) will use different epsilon values, potentially filtering different patients.

**Impact:** Minor reproducibility concern across cohorts.

---

## Summary

| Severity | Count | Key Concerns |
|----------|-------|-------------|
| Critical | 7 | IPCW bias, LOOCV feature inconsistency, lambda path mismatch, decay-to-zero, zero-variance scaling inconsistency, uniform decay half-life, landmark selection |
| Moderate | 5 | Morphological filtering, f_delta mislabeling, in-sample imputation, competing risk documentation, immortal time bias |
| Minor | 3 | CLAUDE.md threshold docs, threshold ordering validation, epsilon reproducibility |

### Top Priority Fixes
1. **Fix zero-variance scaling inconsistency** (`scale_td_panel.m`) — Zero out constant features in both modes
2. **Fix LOOCV feature consistency** — Apply a fixed collinearity mask across all LOOCV folds
3. **Fix lambda path issue** — Recompute lambda path when feature dimensionality changes across folds
4. **Fix decay-to-zero imputation** — Decay toward patient baseline, not zero
5. **Rename `f_Pct` to `f_Delta`** — Prevents misinterpretation of elastic net coefficients
6. **Per-parameter decay half-lives** — Allow different biological decay rates for ADC/D vs f vs D*
7. **Add threshold ordering validation** — Enforce `adc_thresh <= high_adc_thresh` in `parse_config.m`
