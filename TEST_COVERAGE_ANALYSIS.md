# Test Coverage Analysis — pancData3

**Date:** 2026-03-02
**Scope:** All 41 test files in `tests/`, all 17 `core/` modules, all 17 `utils/` modules

---

## Executive Summary

The test suite is strong in **data leakage prevention**, **security validation**, and **happy-path coverage**. The most impactful gaps are:

1. **Two core modules have zero dedicated tests** (`process_single_scan.m`, `fit_models.m`)
2. **One high-complexity utility has zero tests** (`apply_dir_mask_propagation.m`)
3. **Statistical correctness is never validated** — tests check that models *run*, not that outputs are *correct*
4. **Error/edge-case coverage is thin** — malformed inputs, boundary conditions, and failure modes are underrepresented
5. **`parse_config.m` lacks tests for invalid input** — malformed JSON, missing required fields

---

## 1. Modules With No Dedicated Tests

### 1a. `core/process_single_scan.m` — HIGH priority

This is the single-scan workhorse (DICOM conversion → model fitting → DIR → biomarker extraction). It is only reached indirectly through `test_dwi_pipeline.m`, which uses heavily mocked data and covers only the happy path.

**What to test:**
- DICOM conversion failure recovery (dcm2niix returns non-zero)
- Missing GTV mask graceful fallback (GTVn absent)
- DnCNN cache hit vs cache miss paths
- IVIMnet pre-computed result loading
- DIR field caching for Fx2+ scans (deformation field reuse)
- 4D volume dimension validation (mismatched b-value count)
- b-value sorting correctness
- Biomarker extraction with partially empty voxel vectors

### 1b. `core/fit_models.m` — MEDIUM priority

The dependency functions (`IVIMmodelfit`, `fit_adc_mono`) are tested individually, but the wrapper that handles 3D flattening, mask extraction, padding for even dimensions, and zero-output-to-NaN conversion is untested.

**What to test:**
- 3D-to-1D flattening and reconstruction roundtrip
- Padding to even voxel count (odd mask sum)
- Insufficient b-values above threshold (< 2) → graceful NaN return
- Zero output from IVIM fit replaced with NaN
- Non-positive signal values handled correctly in vectorized OLS
- Mask with no true voxels (empty ROI)

### 1c. `utils/apply_dir_mask_propagation.m` — HIGH priority

This is the most algorithmically complex utility (symmetric diffeomorphic Demons registration). It has zero tests despite being critical for longitudinal mask alignment.

**What to test:**
- Identity registration (identical fixed/moving → displacement ≈ 0)
- Known synthetic translation (shift by N voxels → verify displacement field)
- Mask warping correctness (Dice coefficient after roundtrip)
- Empty input handling (empty arrays → empty outputs)
- Dimension mismatch between fixed and moving images
- Robust percentile normalization (handles images with bright outlier voxels)

### 1d. `utils/init_scan_structs.m` — LOW priority

Simple struct initialization. Low complexity, but a quick test would verify field names and array dimensions.

**What to test:**
- Output has expected 27 fields
- Dimensions match requested `n_fx × n_rp`
- All fields initialize to empty `[]`

---

## 2. Thin Error-Handling and Edge-Case Coverage

### 2a. `utils/parse_config.m`

Current tests only check default-filling for missing fields. No tests for:
- **Malformed JSON** (syntax errors, truncated file)
- **Missing required fields** (`dataloc`, `dcm2nii_call`)
- **Invalid `dwi_type`** (typos, unsupported strings)
- **Empty config file** (0-byte file)
- **Non-existent file path**

### 2b. `utils/escape_shell_arg.m`

Current tests cover standard quoting but miss:
- **Non-ASCII / Unicode characters** (patient names with accents, CJK paths)
- **Empty string input**
- **Extremely long paths** (> 260 chars on Windows)
- **Strings containing newlines or null bytes**

### 2c. `utils/safe_load_mask.m`

Good coverage of unsafe type rejection, but missing:
- **Corrupted `.mat` files** (truncated, invalid header)
- **Multi-variable `.mat` files** (requested variable exists among many)
- **Very large arrays** (memory pressure scenarios)

### 2d. `core/compute_summary_metrics.m`

Tested for NaN handling and empty vectors, but missing:
- **Extreme values** (Inf in voxel data, negative ADC values)
- **Single-voxel GTV** (edge case for kurtosis/skewness — undefined for n=1)
- **Histogram bin edge cases** (all values identical → zero-width bins)
- **Octave compatibility path** (nanmean/nanstd fallbacks are untested in CI)

### 2e. `core/load_dwi_data.m`

Tests cover skip-to-reload and fallback naming, but missing:
- **Checkpoint recovery after partial completion** (some patients done, others not)
- **Patient ID normalization edge cases** (underscore/hyphen mismatches between clinical sheet and folder names)
- **Clinical data merging with categorical vs char vs cell ID types**
- **Empty clinical spreadsheet** (no matching patients)

---

## 3. Statistical Correctness Never Validated

Many tests verify that statistical functions *execute without error* but never check that the *outputs are correct*.

### 3a. `core/metrics_survival.m`

Tests confirm the Cox model runs with sufficient events and exits early with insufficient events. Never tested:
- **Known-answer Cox coefficients** (synthetic data with known HR)
- **Competing risk censoring correctness** (event code 2 treated as censored)
- **Likelihood ratio test statistic** against known chi-squared value
- **Counting process interval construction** verified against hand-computed intervals

### 3b. `core/metrics_stats_predictive.m`

Tests confirm Elastic Net produces output and leakage detection fires. Never tested:
- **Known-answer Elastic Net coefficients** (e.g., two correlated features → correct shrinkage)
- **LOOCV risk score distribution** (synthetic separable data should yield AUC ≈ 1.0)
- **ROC AUC correctness** against `perfcurve` with known inputs
- **Youden threshold calculation** verified numerically

### 3c. `core/metrics_stats_comparisons.m`

Tests exercise `perform_statistical_test` but not the full `metrics_stats_comparisons.m` orchestrator. Never tested:
- **Benjamini-Hochberg FDR correction** applied to real p-value vectors (note: `test_statistical_methods.m` tests BH inline, but not via the actual module)
- **GLME mixed-effects model** convergence and coefficient signs
- **CSV output format** (significant metrics table correctness)

### 3d. `utils/build_td_panel.m`

Tests confirm panel construction runs without error. Never tested:
- **Interval boundaries** against hand-computed start/stop times
- **Event assignment** to correct (final) interval only
- **Competing risk code propagation** through intervals
- **Exponential decay imputation** produces expected covariate values at intermediate times

---

## 4. Integration Test Gaps

### 4a. Multi-step pipeline integration

`test_dwi_pipeline.m` and `test_modularity.m` run abbreviated pipelines with mock data. Missing:
- **Load → sanity → visualize → metrics end-to-end** with a small but realistic synthetic dataset
- **Step-skipping correctness** (e.g., skip `load`, verify downstream steps use correct cached data)
- **Config mutation between runs** (`execute_all_workflows` modifies `config.json` — verify intermediate states)

### 4b. Cross-module data contract tests

No tests verify the *shape and type contracts* between modules:
- `load_dwi_data` output struct → `compute_summary_metrics` input expectations
- `metrics_baseline` output → `metrics_stats_comparisons` input expectations
- `build_td_panel` output → `coxphfit` input format

A **contract test** for each handoff point would catch silent breakage when one module changes its output format.

---

## 5. Visualization Tests Are Smoke-Only

All visualization tests (`test_visualize_smoke.m`, `test_plot_parameter_maps.m`, etc.) verify that PNG files are created but never inspect contents. This is reasonable for pixel-level validation, but missing:
- **Figure element assertions** (e.g., verify that axes labels, titles, legend entries exist)
- **Data-driven checks** (e.g., scatter plot with known data produces correct number of points)
- **Empty/degenerate input handling** (single patient, all-NaN data, zero-variance features)

---

## 6. Suggested Priority Order for New Tests

| Priority | Module / Gap | Effort | Impact |
|----------|-------------|--------|--------|
| **P0** | `process_single_scan.m` — unit tests for each subsection | High | Covers the untested workhorse; highest risk module |
| **P0** | `apply_dir_mask_propagation.m` — identity + synthetic shift tests | Medium | High-complexity algorithm with zero coverage |
| **P1** | `fit_models.m` — 3D flatten/reconstruct, edge cases | Medium | Validates the ADC/IVIM wrapper layer |
| **P1** | `parse_config.m` — malformed JSON, invalid types, missing fields | Low | Prevents silent misconfiguration |
| **P1** | `metrics_survival.m` — known-answer Cox model | Medium | Validates core statistical output |
| **P1** | `build_td_panel.m` — interval boundary correctness | Medium | Validates counting process construction |
| **P2** | `metrics_stats_predictive.m` — known-answer Elastic Net + AUC | Medium | Validates predictive modeling correctness |
| **P2** | Cross-module data contract tests | Medium | Prevents silent interface breakage |
| **P2** | `escape_shell_arg.m` — Unicode, empty, newlines | Low | Hardens security-critical utility |
| **P2** | `compute_summary_metrics.m` — extreme values, single-voxel | Low | Edge case robustness |
| **P3** | `metrics_stats_comparisons.m` — full orchestrator + FDR + GLME | Medium | Validates downstream stat pipeline |
| **P3** | `load_dwi_data.m` — checkpoint recovery, ID normalization | Medium | Robustness under real-world conditions |
| **P3** | Visualization content assertions | Low | Catches label/legend regressions |

---

## 7. Coverage Heatmap (Module × Test Dimension)

| Module | Happy Path | Error Handling | Edge Cases | Statistical Correctness | Integration |
|--------|:---:|:---:|:---:|:---:|:---:|
| load_dwi_data | Partial | Minimal | Minimal | N/A | Partial |
| process_single_scan | None | None | None | N/A | Indirect only |
| discover_patient_files | Good | Partial | Good | N/A | N/A |
| sanity_checks | Good | Good | Partial | N/A | N/A |
| compute_summary_metrics | Good | Partial | Minimal | N/A | N/A |
| visualize_results | Smoke | Minimal | Minimal | N/A | Smoke |
| fit_models | None | None | None | Indirect | None |
| convert_dicom | Good | Partial | N/A | N/A | N/A |
| metrics_baseline | Good | Minimal | Partial | N/A | N/A |
| metrics_longitudinal | Good | Partial | Good | N/A | N/A |
| metrics_dosimetry | Good | Partial | Good | N/A | N/A |
| metrics_stats_comparisons | Partial | Partial | Partial | None | None |
| metrics_stats_predictive | Good | Partial | Partial | None | None |
| metrics_survival | Good | Good | Partial | None | None |
| parse_config | Good | None | None | N/A | N/A |
| escape_shell_arg | Good | N/A | Partial | N/A | N/A |
| safe_load_mask | Good | Partial | Partial | N/A | N/A |
| build_td_panel | Good | Partial | Minimal | None | N/A |
| scale_td_panel | Good | Good | Good | N/A | N/A |
| make_grouped_folds | Good | Partial | Good | N/A | N/A |
| knn_impute_train_test | Good | Partial | Partial | None | N/A |
| filter_collinear_features | Good | Partial | Partial | None | N/A |
| apply_dir_mask_propagation | None | None | None | None | None |

**Legend:** Good = well covered | Partial = some cases | Minimal = 1-2 cases | None = untested | Smoke = runs without error only | N/A = not applicable

---

## 8. Specific Test Proposals

### Test 1: `test_process_single_scan.m`
```
- test_dicom_conversion_failure_returns_empty_result
- test_missing_gtvn_mask_returns_empty_gtvn_data
- test_dncnn_cache_hit_skips_recomputation
- test_bvalue_sorting_correctness
- test_4d_dimension_mismatch_error
- test_dir_field_caching_fx2_reuses_fx1_deformation
- test_biomarker_extraction_with_partial_voxels
```

### Test 2: `test_apply_dir_mask_propagation.m`
```
- test_identity_registration_zero_displacement
- test_known_translation_displacement_field
- test_mask_warp_dice_above_threshold
- test_empty_inputs_return_empty
- test_dimension_mismatch_error
- test_robust_normalization_ignores_outliers
```

### Test 3: `test_fit_models.m`
```
- test_3d_flatten_reconstruct_roundtrip
- test_odd_mask_sum_padding
- test_insufficient_bvalues_returns_nan
- test_zero_ivim_output_replaced_with_nan
- test_empty_mask_returns_nan
- test_nonpositive_signal_handled
```

### Test 4: `test_parse_config_errors.m`
```
- test_malformed_json_throws_error
- test_empty_file_throws_error
- test_nonexistent_file_throws_error
- test_invalid_dwi_type_uses_default
- test_missing_dataloc_still_returns_struct
```

### Test 5: `test_survival_correctness.m`
```
- test_known_answer_cox_coefficients
- test_competing_risk_censoring
- test_counting_process_intervals_match_expected
- test_likelihood_ratio_against_known_chi2
```

### Test 6: `test_build_td_panel_correctness.m`
```
- test_interval_boundaries_match_scan_schedule
- test_event_assigned_to_final_interval_only
- test_competing_risk_code_propagated
- test_all_nan_patient_skipped
- test_decay_imputation_values
```
