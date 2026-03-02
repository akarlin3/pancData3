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
Alternatively, you can use the `execute_pipeline.m` script which sets up the parallel pool and environment for a single run.

## Repository Structure & Scripts

### Root Directory
- **`run_dwi_pipeline.m`**: The main orchestrator of the pipeline. It handles environment setup, loads configurations, and calls modular processing steps. Supports sequential processing of multiple DWI types based on `config.json` configuration.
- **`run_all_tests.m`**: The master test runner for the MATLAB DWI pipeline. Discovers and executes all tests in the `tests/` directory and outputs a coverage report.
- **`config.json`** & **`config.example.json`**: Configuration files used to customize parameters (data paths, `dwi_type` selection, temporal groupings) without modifying source code.

### `core/`
Contains the primary logic blocks for the DWI analysis:
- **`load_dwi_data.m`**: Loads DICOM DWI data, applies Deep Learning denoising, fits models, extracts dose metrics, and builds baseline clinical features. Includes a robust **Checkpointing** mechanism and `parfor` loop to safely process large cohorts iteratively and recover from interruptions.
- **`sanity_checks.m`**: Validates data by checking convergence (identifying `NaN`, `Inf`, or negative outputs), summarizing missingness, outlining outliers, and confirming spatial alignment between Dose and DWI maps.
- **`metrics.m`**: The main statistical engine. It evaluates repeatability, builds longitudinal data, handles temporal grouping, applies stabilized IPCW techniques for survival analysis, performs feature selection, and evaluates Competing Risks models (Cause-Specific Hazards) or Logistic Regression.
- **`visualize_results.m`**: Generates visual outputs, including parameter maps overlaid on anatomical MRI, feature distributions by outcome, and longitudinal trajectories.

### `utils/`
Helper scripts for specific data processing and modeling tasks:
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
- **`plot_feature_distribution.m`**: Renders histogram or boxplot visualizations of a feature split by clinical outcome (Local Control vs. Local Failure), with one-way ANOVA p-value annotation for boxplots.

### `tests/`
Contains an extensive suite of MATLAB tests simulating scenarios and asserting correctness of statistical models and data transformations:
- **`test_dwi_pipeline.m`**: Comprehensive tests covering imputation bounds, leakage, IPCW weighting, competing risks analysis, and feature space matching.
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
- **`test_IVIMmodelfit.m`**: Tests IVIM model fitting routines.
- **`test_dvh.m`**: Tests dose-volume histogram computation functions.
- **`test_load_dwi_data.m`**: Integration tests for the data loading and checkpointing logic.
- **`test_metrics_baseline.m`**: Tests baseline statistical metric computations in `metrics.m`.
- **`test_sanity_checks.m`**: Validates convergence and spatial alignment checks in `sanity_checks.m`.
- **`test_compute_summary_metrics.m`**: Asserts correctness of summary metric aggregation across patients.
- **`test_visualize_smoke.m`**: Smoke tests for `visualize_results`, verifying that output figures (`Parameter_Maps`, `Feature_Histograms`, `Feature_BoxPlots`, `Dose_vs_Diffusion`) are produced under valid data, missing b-value, and protocol deviation conditions.
- **`test_visualize_refactor.m`**: Regression tests ensuring visualization refactoring preserves prior output behavior.
- **`test_modularity.m`**: Checks that pipeline modules can be called independently without side effects.
- **`test_normalization_logic.m`**, **`test_landmark_cindex.m`**, **`test_landmark_cindex_mock.m`**: Specialized unit tests for normalization, concordance index computation, and mock-based landmark validation.
- **`benchmark_filter_collinear.m`**, **`benchmark_metrics_opt.m`**: Performance benchmarks for collinearity filtering and metrics optimization.

### Dependencies
- MATLAB R2021a or later (with Image Processing Toolbox and Statistics and Machine Learning Toolbox).
- MRIcroGL’s `dcm2niix` executable or equivalent added to the system path for DICOM-to-NIFTI conversion.
- Script implementations for Intravoxel Incoherent Motion (IVIM) fitting (`IVIMmodelfit.m`, `IVIM_seg.m`, `IVIM_bayes.m`) and mono-exponential ADC (`fit_adc_mono.m`).
- Deep learning architectures appropriately configured for inference (e.g., `dncnn_model.mat`).
- Functions for DVH processing (`dvh.m`, `sample_rtdose_on_image.m`).
