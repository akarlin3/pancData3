# Test Coverage Analysis — pancData3

**Date:** 2026-03-07
**Scope:** All 18 core modules, 28 utility modules, and 72 test files

---

## Executive Summary

The test suite is strong in critical areas (data leakage prevention, security, statistical correctness, integration testing) but has notable gaps in several modules that lack dedicated unit tests, and in cross-cutting concerns like error recovery, config edge cases, and the DICOM/NIfTI I/O pipeline. Below are the specific areas where coverage should be improved.

---

## 1. Modules With No Dedicated Test File

These source files have **zero** dedicated test coverage (only referenced indirectly through integration tests or `test_source_code_standards.m`):

| Module | Type | Risk Level | What's Missing |
|--------|------|------------|----------------|
| `core/metrics_stats_comparisons.m` | Core | **High** | No `test_metrics_stats_comparisons.m` exists. Only exercised indirectly via `test_dwi_pipeline.m` end-to-end test. This module performs Wilcoxon rank-sum tests and GLME mixed-effects models — both have edge cases (all-NaN groups, single-patient groups, identical values in both arms) that should be tested. |
| `core/process_single_scan.m` | Core | **High** | No dedicated tests. This is a 720-line function that orchestrates DICOM conversion, mask loading, model fitting, DIR registration, and biomarker extraction — the largest untested core module. |
| `utils/format_p_value.m` | Utils | **Medium** | Simple function but no tests for boundary values (p = 0, p = 1, p = 0.001 exactly, negative p, Inf). |
| `utils/remove_constant_columns.m` | Utils | **Medium** | No tests despite being used in survival analysis preprocessing. Edge cases: all-NaN matrix, single-column matrix, single-row matrix, mixed NaN/constant columns. |
| `utils/PipelineProgressGUI.m` | Utils | Low | GUI wrapper; hard to unit test, but the step-mapping logic and edge cases (unknown step keys, empty step list, duplicate steps) could be tested. |

---

## 2. Modules With Shallow Test Coverage

These have test files but the coverage is thin relative to the module's complexity:

### 2a. `core/plot_feature_distributions.m` — 2 test methods
- Only smoke-tests that figures are created without errors
- **Missing:** Verification of ANOVA p-value annotation text, correct group separation (LC vs LF), handling when one group has 0 patients, all-NaN parameter vectors

### 2b. `core/plot_scatter_correlations.m` — 6 test methods
- Tests basic figure generation and output files
- **Missing:** Spearman correlation correctness (verify rho/p values), empty dose data, single-patient groups, dose=0 for all patients

### 2c. `core/plot_parameter_maps.m` — 4 test methods
- Smoke tests for figure output
- **Missing:** Edge cases with empty masks, single-slice volumes, protocol deviations (missing b-values)

### 2d. `core/convert_dicom.m` — 3 test methods
- Tests basic conversion path
- **Missing:** Malformed DICOM directories, permission errors, `dcm2niix` failure modes (non-zero exit code, missing executable), directories with no `.dcm` files, mixed DICOM series in one folder

### 2e. `core/discover_patient_files.m` — 4 test methods (via `test_discover_patient_files.m`)
- Tests happy-path file discovery
- **Missing:** Broken symlinks, permission-denied directories, deeply nested patient folder structures, patients with partial data (some fractions missing)

---

## 3. Cross-Cutting Gaps

### 3a. Error Recovery and Graceful Degradation
The pipeline design uses `try/catch` with warning-and-continue for non-fatal modules. However, **no tests verify that the pipeline actually continues** when individual modules throw errors. Recommended:
- Test that `run_dwi_pipeline` completes successfully when `metrics_longitudinal` throws an error
- Test that `run_dwi_pipeline` completes when `metrics_dosimetry` throws (e.g., no dose data at all)
- Test that `metrics_survival` failure doesn't prevent output `.mat` files from being saved

### 3b. Config Validation Edge Cases
`test_parse_config.m` exists and covers basics, but is missing:
- **Invalid field values:** `dwi_type` set to an unrecognized string (e.g., `"BadType"`), negative thresholds (`adc_thresh = -1`), `core_method` set to unknown method
- **Type coercion:** Numeric fields passed as strings in JSON (e.g., `"ivim_bthr": "100"` instead of `100`)
- **Contradictory config:** `skip_to_reload = true` with no existing checkpoint files
- **Extra/unknown fields:** Config file with spurious fields (should be silently ignored, not cause errors)

### 3c. Parallel Execution Safety
- `parsave_dir_cache.m` has 6 test methods, but no tests verify behavior under **concurrent writes** (the actual `parfor` scenario)
- No tests for `parfor_progress.m` race conditions when multiple workers update progress simultaneously
- No tests for checkpoint recovery when `parfor` loop is interrupted mid-cohort

### 3d. NIfTI I/O Edge Cases
No tests exercise:
- Corrupted NIfTI headers (truncated files, wrong magic number)
- NIfTI files with unexpected dimensionality (2D instead of 3D, 5D)
- `.bval` files with non-numeric content, extra whitespace, or wrong delimiter
- Empty NIfTI volumes (all zeros)

### 3e. DIR (Deformable Image Registration) Error Paths
`test_apply_dir_mask_propagation.m` has 8 tests, but:
- No tests for `process_single_scan`'s DIR caching logic (loading cached displacement fields, cache invalidation)
- No tests for the warping path where `D_forward` exists but the mask is empty
- No tests verifying that parameter maps warped with the same displacement field remain spatially consistent

---

## 4. Statistical Method Coverage Gaps

### 4a. Survival Analysis (`metrics_survival.m`)
- 6 test methods exist but key scenarios are missing:
  - **All patients censored** (no events) — Cox model should handle gracefully
  - **Single covariate with zero variance** — should be filtered before `coxphfit`
  - **Competing risks with zero events in one cause** — cause-specific hazard model edge case
  - **Tied event times** — common in clinical data, affects concordance index

### 4b. Predictive Modeling (`metrics_stats_predictive.m`)
- 6 test methods, but missing:
  - **Feature matrix with more features than observations** (p > n) — common in DWI with many parameters
  - **Cross-validation with very small folds** (< 5 patients per fold)
  - **Collinear feature filtering interaction** with downstream model fitting
  - **Reproducibility test:** Same data should produce same AUC/C-index across runs (random seed stability)

### 4c. Longitudinal Metrics (`metrics_longitudinal.m`)
- 9 test methods, but missing:
  - **Single timepoint per patient** (no longitudinal data available) — should skip gracefully
  - **Non-monotonic timepoints** (out-of-order scan dates)
  - **Duplicate timepoints** for the same patient-fraction

---

## 5. Recommendations — Priority Order

### Priority 1: High-Risk Gaps (should add tests now)

1. **`test_metrics_stats_comparisons.m`** — New test file for rank-sum and GLME testing with edge cases (empty groups, identical distributions, single-patient groups, all-NaN features)

2. **`test_format_p_value.m`** — Small but important: verify p=0, p=1, p=0.001, p=NaN, p<0, p=Inf

3. **`test_remove_constant_columns.m`** — Verify all-NaN columns removed, single-value columns removed, mixed columns handled, empty matrix input, single-row input

4. **Config validation hardening** in `test_parse_config.m` — Add tests for invalid `dwi_type`, invalid `core_method`, negative thresholds, string-typed numeric fields

5. **Pipeline resilience tests** — Verify `run_dwi_pipeline` continues when optional modules throw errors

### Priority 2: Medium-Risk Gaps (should add before next major release)

6. **`process_single_scan.m` unit tests** — Test the helper functions (`load_nifti_mask`, `extract_biomarkers`, `compute_dncnn_fallback`) in isolation with mock data

7. **Survival analysis edge cases** — All-censored cohort, zero-variance covariates, tied event times

8. **Predictive modeling robustness** — p > n scenario, random seed reproducibility, very small folds

9. **NIfTI I/O edge cases** — Corrupted headers, dimension mismatches, malformed `.bval` files

### Priority 3: Nice-to-Have Improvements

10. **MATLAB `TestParameterization`** — Convert manual loops (e.g., iterating 11 core methods) to native parameterized tests for better test report granularity

11. **Concurrent write tests for `parsave_dir_cache.m`** — Verify file integrity under parallel `save` operations

12. **Octave CI** — Run the full test suite under GNU Octave to validate the 21 `.octave_compat/` shims

---

## 6. Current Coverage Summary

| Category | Modules | With Dedicated Tests | Coverage |
|----------|---------|---------------------|----------|
| Core modules | 18 | 15 | 83% |
| Utility modules | 28 | 23 | 82% |
| **Overall** | **46** | **38** | **83%** |

**Modules with no dedicated tests:** `metrics_stats_comparisons`, `process_single_scan`, `format_p_value`, `remove_constant_columns`, `PipelineProgressGUI`

**Test suite size:** 72 test files, ~300+ individual test methods, 27 parallel-safe test classes
