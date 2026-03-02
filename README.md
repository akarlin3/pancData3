# pancData

Analysis pipeline for pancreatic DWI (Diffusion-Weighted Imaging) data, including IVIM model fitting, ADC computation, spatial Deep Learning denoising, and correlation with RT dose maps.

## Agent Collaboration & Safety Rules
This repository is maintained with the assistance of two AI agents:
- **Antigravity (Local)**: Handles core physics modeling, MRI calibration logic, and specialized scripts. 
- **Jules (Background)**: Delegated for background tasks including unit testing, documentation, and code styling.
> [!CAUTION]
> **Safety Rule**: Never send patient data or sensitive CSVs to the Jules cloud. Only logic and code structures are permitted. Do not modify files in the `dependencies/` folder.

## Workflows
- **`/run_data`**: A structured workflow that executes the DWI pipeline sequentially for all DWI types (`Standard`, `dnCNN`, and `IVIMnet`). It iteratively modifies `config.json` and evaluates the master orchestrator, finishing by running the full test suite.

## Running the Pipeline Manually

To run the full pipeline sequentially for all DWI types (`Standard`, `dnCNN`, and `IVIMnet`), open MATLAB, navigate to the repository root, and execute the wrapper script in the Command Window:

```matlab
execute_all_workflows
```

To run a single specific phase or test a configuration, you can call the main orchestrator directly:

```matlab
% 1. Add required paths
addpath('core', 'utils', 'dependencies');

% 2. Run the pipeline (e.g., only the 'load' step)
run_dwi_pipeline('config.json', {'load'});
```

## Repository Structure & Scripts

### Root Directory
- **`run_dwi_pipeline.m`**: The main orchestrator of the pipeline. It handles environment setup, loads configurations, and calls modular processing steps. Supports sequential processing of multiple DWI types based on `config.json` configuration.
- **`execute_all_workflows.m`**: Wrapper that sets up a parallel pool (max 2 workers) and runs all three DWI types (`Standard`, `dnCNN`, `IVIMnet`) sequentially, modifying `config.json` between runs and skipping reload after the first.
- **`run_all_tests.m`**: The master test runner for the MATLAB DWI pipeline. Discovers and executes all tests in the `tests/` directory and outputs a coverage report.
- **`test_octave.m`**: Compatibility test script for GNU Octave.
- **`config.json`** & **`config.example.json`**: Configuration files used to customize parameters (data paths, `dwi_type` selection, temporal groupings) without modifying source code.

### `core/`
Contains the primary logic blocks for the DWI analysis (17 files):
- **`load_dwi_data.m`**: Loads DICOM DWI data, applies Deep Learning denoising, fits models, extracts dose metrics, and builds baseline clinical features. Includes a robust **Checkpointing** mechanism and `parfor` loop to safely process large cohorts iteratively and recover from interruptions.
- **`process_single_scan.m`**: Processes a single fraction-by-repeat scan for a patient, handling DICOM conversion, mask saving, dose resampling, volume loading, model fitting, DIR registration, DL loading, and biomarker extraction.
- **`sanity_checks.m`**: Validates data by checking convergence (identifying `NaN`, `Inf`, or negative outputs), summarizing missingness, outlining outliers, and confirming spatial alignment between Dose and DWI maps.
- **`compute_summary_metrics.m`**: Aggregates voxel-level data into patient-level summary metrics.
- **`metrics_baseline.m`**: Computes baseline measures, outlier cleaning, and percent delta calculations.
- **`metrics_longitudinal.m`**: Longitudinal change analysis across treatment timepoints.
- **`metrics_dosimetry.m`**: Dose-related metrics computation (D95, V50) within diffusion-defined subvolumes.
- **`metrics_stats_comparisons.m`**: Statistical group comparisons between clinical outcome groups.
- **`metrics_stats_predictive.m`**: Predictive modeling, feature selection, and cross-validated model evaluation.
- **`metrics_survival.m`**: Survival analysis, competing risks modeling (Cause-Specific Hazards), and IPCW weighting.
- **`discover_patient_files.m`**: File system navigation for discovering and organizing patient cohort data.
- **`convert_dicom.m`**: DICOM-to-NIfTI conversion via `dcm2niix`.
- **`fit_models.m`**: ADC mono-exponential and IVIM model fitting with 3D-to-1D flattening and padding.
- **`visualize_results.m`**: Generates visual outputs, including parameter maps overlaid on anatomical MRI, feature distributions by outcome, and longitudinal trajectories.
- **`plot_parameter_maps.m`**: Visualization helpers for parameter map overlays on anatomical images.
- **`plot_feature_distributions.m`**: Feature histogram and boxplot rendering split by clinical outcome.
- **`plot_scatter_correlations.m`**: Correlation scatter plots between dose and diffusion metrics.

### `utils/`
Helper scripts for specific data processing and modeling tasks (17 files):
- **`parse_config.m`**: Parses `config.json` into MATLAB structs.
- **`build_td_panel.m`**: Structures patients' longitudinal data into time-dependent panels.
- **`scale_td_panel.m`**: Handles timepoint-specific feature scaling across varying fractions.
- **`make_grouped_folds.m`**: Generates patient-stratified folds for grouped cross-validation, preventing intra-patient data leakage.
- **`knn_impute_train_test.m`**: A specialized KNN imputer enforcing strict bounds against temporal leakage.
- **`filter_collinear_features.m`**: Prunes highly collinear features before model fitting.
- **`apply_dir_mask_propagation.m`**: Handles spatial deformations and rigid alignment.
- **`safe_load_mask.m`**: Securely loads a mask variable from a `.mat` file by inspecting variable type before loading, rejecting unsafe or non-numeric classes to prevent arbitrary code execution.
- **`calculate_subvolume_metrics.m`**: Computes dose coverage metrics (D95, V50) for a diffusion-defined sub-volume within the GTV, including morphological cleanup via opening/closing.
- **`load_dl_provenance.m`**: Loads the deep learning training provenance manifest (dnCNN and IVIMnet training IDs) to guard against data leakage between model training and patient analysis sets.
- **`perform_statistical_test.m`**: Encapsulates statistical hypothesis testing (Wilcoxon rank-sum) with NaN-safe group extraction and safe handling of insufficient data.
- **`parsave_dir_cache.m`**: Wrapper enabling `save` calls inside `parfor` loops for caching deformable image registration (DIR) results.
- **`escape_shell_arg.m`**: Cross-platform shell argument escaping (Windows and Unix) for safe file path use in `system()` calls.
- **`discover_gtv_file.m`**: Low-level helper that locates a GTV mask file in a folder using flexible naming patterns and a repeat scan index, with a single-mask fallback.
- **`find_gtv_files.m`**: Locates primary (GTVp) and nodal (GTVn) GTV mask files, handling multiple naming conventions for patients with complex tumor anatomy.
- **`init_scan_structs.m`**: Creates pre-populated empty GTVp and GTVn struct arrays with a consistent field ordering to prevent struct concatenation errors in `parfor` loops.
- **`plot_feature_distribution.m`**: Renders histogram or boxplot visualizations of a feature split by clinical outcome (Local Control vs. Local Failure), with one-way ANOVA p-value annotation for boxplots.

### `tests/`
Contains an extensive suite of MATLAB tests (45 files) simulating scenarios and asserting correctness of statistical models and data transformations:
- **`test_dwi_pipeline.m`**: Comprehensive tests covering imputation bounds, leakage, IPCW weighting, competing risks analysis, and feature space matching.
- **`test_statistical_methods.m`**: Pure algorithmic tests for statistical and numerical routines including BH/Holm corrections, correlation filtering, LOOCV, FDR, imputation, and time-dependent panel construction.
- **`test_source_code_standards.m`**: Static code analysis tests verifying source-code patterns, naming conventions, and architectural invariants without executing pipeline logic.
- **`test_corr_filter.m`**: Validates correlation pruning logic.
- **`test_grouped_folds.m`**: Asserts no data leakage occurs across generated folds.
- **`test_knn_temporal_leakage.m`**: Confirms KNN imputation enforces infinite distances for rows originating from the same patient.
- **`test_mask_loading.m`**: Checks target-constrained normalization logic and boundary padding.
- **`test_safe_load_mask.m`**: Verifies that `safe_load_mask` correctly rejects unsafe variable classes and handles missing files or variables.
- **`test_escape_shell_arg.m`**: Tests shell argument escaping across Windows and Unix styles, including edge cases with embedded quotes and trailing backslashes.
- **`test_perform_statistical_test.m`**: Validates rank-sum test behavior for various group configurations and degenerate inputs.
- **`test_scale_td_panel.m`**: Asserts correct timepoint-specific scaling without cross-timepoint leakage.
- **`test_build_td_panel.m`**: Tests longitudinal panel construction for correctness of patient ordering and data alignment.
- **`test_parse_config.m`**: Validates config parsing for required fields and type correctness.
- **`test_fit_adc_mono.m`**: Unit tests for mono-exponential ADC model fitting.
- **`test_fit_models.m`**: Tests the `fit_models` wrapper: 3D-to-1D flattening roundtrips, padding, graceful NaN returns for insufficient b-values, and empty mask handling.
- **`test_IVIMmodelfit.m`**: Tests IVIM model fitting routines.
- **`test_dvh.m`**: Tests dose-volume histogram computation functions.
- **`test_load_dwi_data.m`**: Integration tests for the data loading and checkpointing logic.
- **`test_metrics_baseline.m`**: Tests baseline statistical metric computations.
- **`test_metrics_longitudinal.m`**: Smoke tests for longitudinal analysis including figure creation, single-patient handling, and all-NaN data resilience.
- **`test_metrics_dosimetry.m`**: Tests dosimetry metric computation including output dimensions, NaN passthrough, dtype dispatch, and absent GTV mask handling.
- **`test_metrics_stats_predictive.m`**: Tests predictive modeling including no-op paths, all-NaN features, correct output types, and DL provenance leakage detection.
- **`test_metrics_survival.m`**: Tests survival analysis including insufficient events early return, competing-risk events, and valid_pts logical mask subsetting.
- **`test_sanity_checks.m`**: Validates convergence and spatial alignment checks in `sanity_checks.m`.
- **`test_compute_summary_metrics.m`**: Asserts correctness of summary metric aggregation across patients.
- **`test_visualize_smoke.m`**: Smoke tests for `visualize_results`, verifying that output figures (`Parameter_Maps`, `Feature_Histograms`, `Feature_BoxPlots`, `Dose_vs_Diffusion`) are produced under valid data, missing b-value, and protocol deviation conditions.
- **`test_visualize_refactor.m`**: Regression tests ensuring visualization refactoring preserves prior output behavior.
- **`test_plot_parameter_maps.m`**: Smoke tests for `plot_parameter_maps` covering no-valid-patient paths, protocol deviations, and valid data producing PNGs.
- **`test_plot_feature_distributions.m`**: Tests feature distribution plot generation with temporary output directories.
- **`test_plot_scatter_correlations.m`**: Smoke tests for `plot_scatter_correlations` covering valid data, insufficient points, dtype index dispatch, and all-NaN dose vectors.
- **`test_modularity.m`**: Checks that pipeline modules can be called independently without side effects.
- **`test_apply_dir_mask_propagation.m`**: Tests deformable image registration: identity registration, mask warping, empty input handling, dimension mismatch rejection, and outlier normalization.
- **`test_calculate_subvolume_metrics.m`**: Tests D95 and V50 dose-metric computation for diffusion-defined resistant sub-volumes.
- **`test_convert_dicom.m`**: Tests DICOM-to-NIfTI conversion using mock `dcm2niix` scripts on both platforms.
- **`test_discover_gtv_file.m`**: Tests GTV mask file discovery logic with temporary directories and dummy files.
- **`test_discover_patient_files.m`**: Tests patient file system navigation using mock data directories and mock `dicominfo`.
- **`test_find_gtv_files.m`**: Tests GTVp/GTVn mask file location for single-GTV and dual-GTV patients.
- **`test_init_scan_structs.m`**: Tests struct array creation with expected field names, dimensions, and empty initialization.
- **`test_load_dl_provenance.m`**: Tests loading of deep learning provenance manifests from `.mat` files.
- **`test_parsave_dir_cache.m`**: Tests parallel-safe save wrapper for file creation, variable persistence, and overwriting.
- **`test_fix_verify.m`**: Reproduces and verifies a path construction bug fix related to GTV mask file loading.
- **`test_perf_knn.m`**: Performance benchmark for `knn_impute_train_test` on a 500-row, 50-feature matrix.
- **`test_normalization_logic.m`**, **`test_landmark_cindex.m`**, **`test_landmark_cindex_mock.m`**: Specialized unit tests for normalization, concordance index computation, and mock-based landmark validation.

#### `tests/benchmarks/`
Performance benchmarks comparing algorithm implementations.
- **`benchmark_filter_collinear.m`**: Benchmarks collinearity filtering speed.
- **`benchmark_metrics_opt.m`**: Benchmarks metric computation optimizations.
- **`benchmark_make_grouped_folds.m`**: Benchmarks patient-stratified fold generation.
- **`benchmark_scale_td.m`**: Benchmarks temporal panel scaling (original vs. optimized).
- **`test_opt.m`**: Compares patient-level outcome aggregation strategies (loop vs. `accumarray`).
- **`test_opt2.m`**: Compares fold-assignment optimization approaches.
- **`test_perf.m`**: Benchmarks GTV file search pattern matching.

#### `tests/diagnostics/`
Lightweight diagnostic tests for spot-checking during development.
- **`run_manual_test.m`**: Manually exercises `compute_summary_metrics` with inline mock data.
- **`run_test_local.m`**: Quick local harness that calls `test_compute_summary_metrics` directly.
- **`test_accumarray.m`**: Spot-check of MATLAB's `accumarray` for patient-level aggregation.
- **`test_make_grouped_folds.m`**: Minimal sanity check for `make_grouped_folds`.
- **`test_vis_script.m`**: Smoke runner for `visualize_results` without the full test framework.

### Dependencies
- MATLAB R2021a or later (with Image Processing Toolbox and Statistics and Machine Learning Toolbox).
- MRIcroGL’s `dcm2niix` executable or equivalent added to the system path for DICOM-to-NIFTI conversion.
- Script implementations for Intravoxel Incoherent Motion (IVIM) fitting (`IVIMmodelfit.m`, `IVIM_seg.m`, `IVIM_bayes.m`) and mono-exponential ADC (`fit_adc_mono.m`).
- **`apply_dncnn_symmetric.m`**: DnCNN deep learning denoising applied to MRI images within a GTV bounding box.
- **`clean_dir_command.m`**: Wraps MATLAB’s `dir()` to filter out `.`, `..`, and non-folder entries.
- **`halfSampleMode.m`**: Computes the half-sample mode (robust location estimator) for each row of a matrix.
- **`im2Y.m`**: Transforms a 3D/4D functional image array into a 2D data matrix, optionally masked to a region of interest.
- Functions for DVH processing (`dvh.m`, `sample_rtdose_on_image.m`).
