# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [1.1.0] - 2026-03-08

### Added
- **Predictive modeling** (`metrics_stats_predictive.m`): Elastic Net with nested LOOCV for pancreatic cancer treatment response prediction, integrated into pipeline orchestration
- **Centralized configuration** (`parse_config.m`): Configuration parsing with field defaults and file I/O error handling, integrated into `run_dwi_pipeline.m` and `execute_all_workflows.m`
- New tests for `process_single_scan`, `metrics_stats_comparisons`, `format_p_value`, `PipelineProgressGUI`, `plot_feature_distribution`, and `filter_collinear_features`

### Fixed
- Error message in `execute_all_workflows.m` incorrectly referenced `'gqae run'` instead of `'dnCNN run'`
- Extraneous characters in `metrics_stats_predictive.m` file header
- Missing `spectral_min_voxels` field in `config.example.json`
- Missing `patient_ids` backwards-compatible default in `parse_config.m`
- Windows `%` characters in paths not escaped for cmd.exe (`escape_shell_arg.m`)
- `compute_dice_hausdorff.m` did not validate mask dimensions before comparison
- `build_td_panel.m` did not validate feature array row counts against patient count

### Test Suite
- Upgraded `test_metrics_survival.m` assertions from `verifyTrue(true)` to actual content verification via `evalc()`
- Added n=24 LOOCV path test with signal injection to `test_metrics_stats_predictive.m`
- Added content assertions to `test_dwi_pipeline.m` end-to-end test (load `.mat` files, verify finite values)
## [1.2.0-alpha.1] - 2026-03-08

### Added
- **Python post-hoc analysis scripts** (`analysis/`): Five new scripts for automated graph extraction and cross-DWI comparison using Claude vision API
  - `batch_graph_analysis.py`: Async batch processing of pipeline graph images via Claude vision API; outputs structured CSV with axes, trends, and inflection points
  - `cross_reference_dwi.py`: Full cross-DWI comparison (Standard vs dnCNN vs IVIMnet) of trends, inflection points, and summaries
  - `cross_reference_summary.py`: Concise cross-DWI summary focusing on priority clinical graphs and trend agreement/disagreement
  - `statistical_relevance.py`: Extracts p-values and correlation coefficients; reports significant findings and cross-DWI significance
  - `statistical_by_graph_type.py`: Filters statistical findings by graph type (scatter, box, line, heatmap, bar, histogram, parameter_map)

### Changed
- Version bump to 1.2.0-alpha.1 to begin 1.2.0 development cycle
- Updated documentation (`CLAUDE.md`, `README.md`) with analysis script descriptions and requirements

## [1.1.0-rc.1] - 2026-03-08

### Added
- **Predictive modeling** (`metrics_stats_predictive.m`): Elastic Net with nested LOOCV for pancreatic cancer treatment response prediction, integrated into pipeline orchestration
- **Centralized configuration** (`parse_config.m`): Configuration parsing with field defaults and file I/O error handling, integrated into `run_dwi_pipeline.m` and `execute_all_workflows.m`
- New tests for `process_single_scan`, `metrics_stats_comparisons`, `format_p_value`, `PipelineProgressGUI`, `plot_feature_distribution`, and `filter_collinear_features`

### Fixed
- Error message in `execute_all_workflows.m` incorrectly referenced `'gqae run'` instead of `'dnCNN run'`
- Extraneous characters in `metrics_stats_predictive.m` file header
- Missing `spectral_min_voxels` field in `config.example.json`

## [1.1.0-beta.1] - 2026-03-07

### Added
- **Tumor core method comparison** (`compare_core_methods.m`): Pairwise comparison of all 11 tumor core delineation methods with Dice coefficient and Hausdorff distance metrics, integrated as an optional `compare_cores` pipeline step
- **Patient data check** (`patient_data_check.m`): Pre-pipeline data integrity scanner that validates directory structure, DICOM availability, GTV masks, and RT dose folders
- **Pipeline progress GUI** (`PipelineProgressGUI.m`): Pipeline-aware progress bar wrapper that maps step keys to display names
- **Professional progress bar** (`ProgressGUI.m`): Custom-figure progress bar for MATLAB pipelines and test suite runs
- **Longitudinal summary metrics** (`compute_summary_metrics.m`): Voxel-to-summary-metric aggregation for tumor and dosimetry analysis
- **Single scan processing** (`process_single_scan.m`): Per-scan DICOM conversion, model fitting, and caching module
- **Cross-DWI subvolume comparison** (`plot_cross_dwi_subvolume_comparison.m`): Visualization for cross-DWI-type ADC subvolume comparison
- **Dice/Hausdorff utility** (`compute_dice_hausdorff.m`): Dice coefficient and Hausdorff distance computation between 3D binary masks
- **JSON field utility** (`json_set_field.m`): Targeted regex replacement of field values in raw JSON strings
- New config fields: `core_method`, `core_percentile`, `core_n_clusters`, `fdm_parameter`, `fdm_thresh`
- Documentation maintenance guidelines in `CLAUDE.md`

### Changed
- Professional progress bar GUI replaces basic waitbar for test suite and pipeline steps

## [1.0.0-beta.2] - 2026-03-06

### Fixed
- **Baseline ref_col selection**: `ref_col` used wrong column instead of earliest fraction as baseline (`metrics_baseline.m`)
- **Iteration limit warning**: `lassoglm`/`fitglm` calls hit iteration cap with `MaxIter` 1e5; increased to 1e7 (`metrics_stats_predictive.m`)
- **Config restore after execute_all_workflows**: `config.json` was not restored after multi-type pipeline execution (`execute_all_workflows.m`)
- **Config .bak overwrite**: `.bak` file restore silently overwrote user config edits between pipeline runs (`execute_all_workflows.m`)
- **Config formatting loss**: `config.json` formatting was not preserved during pipeline execution
- **DnCNN mat2gray normalization**: `mat2gray` normalization was missing for DnCNN-denoised DWI loading (`load_dwi_data.m`)
- **Warning log pollution**: `KFoldMissingGrp` and `defaultScanDays` warnings leaked to pipeline error log

### Changed
- Removed `.bak` file mechanism from `execute_all_workflows` in favor of in-memory config preservation

### Test Suite
- Fix `testDefaultOutputFolder` deleting `execute_all_workflows` output folder during concurrent runs
- Add test for empty string field preservation (`cause_of_death_column`)

## [1.0.0-beta.1] - 2026-03-06

### Added
- `clear_cache` config option to delete all cached pipeline files before execution
- `CauseOfDeath` column support for competing risk classification in survival analysis
- Detailed analytical comments across all 39 non-dependency MATLAB files explaining medical physics and statistical rationale
- Protocol deviation check in `plot_parameter_maps.m` that skips patients with non-standard b-values
- Config backwards-compatibility rule in `CLAUDE.md` covering both field addition and removal

### Fixed

#### Analytical / Logical Errors
- **Immortal time bias**: Scan days were derived from folder names instead of DICOM headers, inflating early-timepoint survival estimates (`metrics_survival.m`)
- **Anti-conservative LRT p-value**: IPCW frequency-weight scaling inflated both the log-likelihood and degrees of freedom, producing artificially small p-values in survival analysis (`metrics_survival.m`)
- **IPCW censoring model**: Non-terminal intervals were misclassified as censoring events, biasing inverse-probability weights (`metrics_survival.m`)
- **Competing risk event counting**: Events were miscounted in survival analysis, and dosimetry used incorrect indexing (`metrics_survival.m`, `metrics_dosimetry.m`)
- **NaN-unsafe column removal**: `remove_constant_columns` in survival dropped columns containing NaN instead of ignoring them, and imputation used decay-to-zero instead of proper fill (`metrics_survival.m`)
- **Repeatability bias**: Failed IVIM fits (f=0, D*=0/NaN) were included in summary metrics instead of being filtered to NaN, biasing voxel-level statistics (`compute_summary_metrics.m`)
- **Simpson's paradox risk**: KS test was applied with k≤1 preventing meaningful distribution comparison; added guard requiring k>1 (`compute_summary_metrics.m`)
- **ROC inflation**: Audit identified inflated ROC curves from data leakage in predictive modeling (`metrics_stats_predictive.m`)
- **Baseline scaling after landmark subsetting**: Hardcoded `t_start==0` failed after landmark subsetting removed time-zero rows; replaced with `min(t_start)` (`scale_td_panel.m`)
- **scale_td_panel zeroing features**: Features with zero standard deviation were zeroed out instead of being clamped to σ=1 (`scale_td_panel.m`)
- **Silent patient dropout**: Missing quote stripping in clinical spreadsheet matching caused patients with quoted IDs to silently drop from analysis (`metrics_baseline.m`)
- **wCV dimension mismatch crash**: Compound masking caused reshape errors during IVIMnet within-coefficient-of-variation computation (`metrics_baseline.m`)
- **Reshape crash in wCV for IVIMnet**: Used self-referential reshape dimensions instead of deriving from potentially stale arrays (`metrics_baseline.m`)
- **Hardcoded DWI type index**: ADC SD heterogeneity check used a hardcoded pipeline index instead of the current DWI type (`metrics_stats_predictive.m`)
- **NaN risk scores from failed LOOCV**: Risk score array initialized as logical (cannot hold NaN); changed to double (`metrics_stats_predictive.m`)
- **Missing dosimetry univariate analysis**: Univariate dose-response analysis was skipped silently (`metrics_dosimetry.m`)
- **GTVn DVH spatial mismatch**: Dose-volume histogram computed on misaligned GTVn masks (`metrics_dosimetry.m`)
- **LF rate denominator error**: Local failure rate denominator was incorrect (`metrics_dosimetry.m`)
- **DVH running-mean accumulation error**: DVH values accumulated incorrectly for patients with 3+ repeat scans (`metrics_dosimetry.m`)
- **Dosimetry 3D mask mismatch**: DIR-warped fraction masks had spatial dimension mismatches (`metrics_dosimetry.m`)
- **Division-by-zero in metrics**: Unguarded division produced Inf/NaN in edge cases (`metrics_stats_comparisons.m`)
- **Silent GLME failure**: Generalized linear mixed-effects model failures were silently swallowed (`metrics_stats_comparisons.m`)
- **Swapped ADC threshold defaults**: `adc_thresh` and `high_adc_thresh` config defaults were inverted (`parse_config.m`)
- **Missing IVIMnet repeat ADC**: IVIMnet repeat ADC vectors were not loaded (`load_dwi_data.m`)
- **Fraction folder substring matching**: `'Fx1'` incorrectly matched `'Fx10'` due to `strfind` substring semantics; switched to exact matching (`discover_patient_files.m`, `load_dwi_data.m`)
- **Stale checkpoint dimensions**: Cached summary metrics with mismatched cohort dimensions were silently loaded instead of recomputed (`compute_summary_metrics.m`)
- **IVIM squeeze vs reshape**: `squeeze()` on IVIM output could collapse wrong dimensions; replaced with explicit `reshape()` (`fit_models.m`)
- **Scatter plot column indexing**: Baseline feature scatter plots indexed wrong columns (`plot_scatter_correlations.m`)
- **Silent dtype clamping and trend line extrapolation**: Plot data was silently clamped to display range, and trend lines extrapolated beyond data bounds (`visualize_results.m`)
- **Colorbar range error**: Parameter map colorbar limits were set incorrectly (`plot_parameter_maps.m`)
- **Competing risk fold counting**: Competing risk patients (y=2) were miscounted as events during fold stratification (`make_grouped_folds.m`)
- **Collinear feature scan-day leakage**: Scan day columns leaked into feature matrix (`filter_collinear_features.m`)
- **wCV denominator floor guard**: Applied floor guard before division to prevent transient Inf values (`metrics_baseline.m`)
- **Hardcoded constants in load path**: Fixed hardcoded b-value and indexing constants (`load_dwi_data.m`)
- **CauseOfDeath case-insensitive mismatch**: Column lookup failed silently when spreadsheet header casing differed across MATLAB versions (`metrics_survival.m`)
- **NaN-to-logical conversion**: Risk score stratification crashed when assigning NaN to a logical array (`metrics_stats_predictive.m`)
- **Missing summary_metrics in sanity-skip fallback**: Pipeline crashed when both load and sanity steps were skipped (`run_dwi_pipeline.m`)
- **Diary restart after load step**: Orchestrator diary was not restarted after load module overrode it (`run_dwi_pipeline.m`)
- **skip_to_reload fallback for typed DWI runs**: First-run of dnCNN/IVIMnet failed because skip_to_reload was true but no prior .mat files existed (`load_dwi_data.m`)

#### Octave Compatibility
- `hist()` bin-count vs bin-edge semantics difference producing wrong numerical results (`plot_feature_distribution.m`)
- SEM computation inconsistency across MATLAB/Octave (`metrics_baseline.m`)

#### Test Suite
- Boxplot crash when only one outcome group is present (`plot_feature_distribution.m`)
- Swapped `adc_thresh`/`high_adc_thresh` in 7 test files
- Test cleanup deleting orchestrator output folder during pipeline runs (`run_all_tests.m`, test files)
- `testVis_DeviationExcludesPatient` false regex match on comment text
- `testCauseOfDeathNoWarningWhenPresent` failures from `readtable` dropping empty columns and `verifyWarningFree` conflicts
- `test_modularity` RepoRoot path computation bug
- 11 additional test fixes across pipeline modules for consistent error IDs and threshold values

#### UI
- Waitbar text centering and layout in test suite GUI
- Test progress bar window sizing and axes centering

### Changed
- Version references updated to 1.0.0 pre-release scheme
- Parallel pool now shuts down at end of standalone test suite runs

## [1.0.0-alpha.1] - 2026-03-03

### Added
- Time-dependent Cox PH model for survival analysis with IPCW weighting
- DICOM-to-NIfTI conversion module (`convert_dicom.m`) using `dcm2niix`
- Per-scan processing module (`process_single_scan.m`) for modular scan handling
- Scan data structure initialization utility (`init_scan_structs.m`)
- Octave compatibility guards across core, utils, and tests
- Octave-compatible test runner (`run_all_tests_octave.m`) and `matlab.unittest` shim layer
- Comprehensive unit and integration tests for DWI data loading, model fitting, mask propagation, and pipeline metrics
- Tests for safe mask loading, statistical methods, and grouped folds
- Core scripts for DWI/IVIM metrics baseline, statistical comparisons, and predictive analysis
- `CONTRIBUTING.md` with contribution guidelines
- `CITATION.cff` for machine-readable academic citation
- `SECURITY.md` with vulnerability reporting policy
- `CHANGELOG.md` for tracking version history
- `.editorconfig` for consistent cross-editor formatting
- `.gitattributes` for consistent line endings and diff behavior
- GitHub Actions CI workflow for automated testing
- GitHub issue and pull request templates
- Improved `README.md` with badges, table of contents, and structured documentation

### Fixed
- Octave compatibility issues across core, utils, and test modules

## [0.0.0] - 2026-01-01

### Added
- Master orchestrator pipeline (`run_dwi_pipeline.m`) with modular step execution
- Multi-type sequential runner (`execute_all_workflows.m`) for Standard, dnCNN, and IVIMnet
- IVIM and ADC model fitting (`fit_models.m`)
- Deep learning denoising integration (DnCNN, IVIMnet)
- Sanity checks with convergence validation and spatial alignment
- Visualization suite: parameter maps, feature distributions, scatter correlations
- Metrics pipeline: baseline, longitudinal, dosimetry, statistical comparisons, predictive modeling
- Survival analysis with competing risks (Cause-Specific Hazards) and IPCW weighting
- Comprehensive test suite with code coverage
- Octave compatibility shims
- Checkpointing system for large cohort recovery
- Security utilities: `safe_load_mask`, `escape_shell_arg`
- Data leakage prevention: patient-stratified folds, temporal KNN bounds, DL provenance tracking
- Configuration-driven pipeline via `config.json`
