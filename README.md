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

### `tests/`
Contains an extensive suite of MATLAB tests simulating scenarios and asserting correctness of statistical models and data transformations:
- **`test_dwi_pipeline.m`**: Comprehensive tests covering imputation bounds, leakage, IPCW weighting, competing risks analysis, and feature space matching.
- **`test_corr_filter.m`**: Validates correlation pruning logic.
- **`test_grouped_folds.m`**: Asserts no data leakage occurs across generated folds.
- **`test_knn_temporal_leakage.m`**: Confirms KNN imputation enforces infinite distances for rows originating from the same patient.
- **`test_mask_loading.m`**: Checks target-constrained normalization logic and boundary padding.
- **`test_normalization_logic.m`**, **`test_landmark_cindex.m`**, etc.: Specialized unit tests for various mathematical operations in the codebase.

### Dependencies
- MATLAB R2021a or later (with Image Processing Toolbox and Statistics and Machine Learning Toolbox).
- MRIcroGL’s `dcm2niix` executable or equivalent added to the system path for DICOM-to-NIFTI conversion.
- Script implementations for Intravoxel Incoherent Motion (IVIM) fitting (`IVIMmodelfit.m`, `IVIM_seg.m`, `IVIM_bayes.m`) and mono-exponential ADC (`fit_adc_mono.m`).
- Deep learning architectures appropriately configured for inference (e.g., `dncnn_model.mat`).
- Functions for DVH processing (`dvh.m`, `sample_rtdose_on_image.m`).
