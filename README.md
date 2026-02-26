# pancData

Analysis pipeline for pancreatic DWI (Diffusion-Weighted Imaging) data, including IVIM model fitting, ADC computation, spatial Deep Learning denoising, and correlation with RT dose maps.

## Repository Structure & Scripts

### Root Directory
- **`run_dwi_pipeline.m`**: The main orchestrator of the pipeline. It handles environment setup (adding paths, checking toolboxes), loads configurations, and subsequently calls modular processing steps.
- **`run_all_tests.m`**: The master test runner for the MATLAB DWI pipeline. It automatically discovers and executes all class-based and function-based tests located in the `tests/` directory, outputting a formatted coverage report and asserting success for CI/CD environments.
- **`config.json`** & **`config.example.json`**: Configuration files used by the pipeline to customize parameters such as data paths, feature selection criteria, and temporal groupings without modifying source code.

### `core/`
Contains the primary logic blocks for the DWI analysis:
- **`load_dwi_data.m`**: Loads DICOM DWI data, applies Deep Learning denoising, fits ADC and IVIM models, extracts dose metrics, and builds the baseline clinical features.
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
- **`TestDwiPipeline.m`**: Comprehensive tests covering imputation bounds, leakage, IPCW weighting, competing risks analysis, and feature space matching.
- **`test_corr_filter.m`**: Validates correlation pruning logic.
- **`test_grouped_folds.m`**: Asserts no data leakage occurs across generated folds.
- **`test_knn_temporal_leakage.m`**: Confirms KNN imputation enforces infinite distances for rows originating from the same patient.
- **`test_mask_loading.m`**: Checks target-constrained normalization logic and boundary padding.
- **`test_normalization_logic.m`**, **`test_landmark_cindex.m`**, etc.: Specialized unit tests for various mathematical operations in the codebase.

### `dependencies/`
External resources and third-party logic used directly by the pipeline:
- NIfTI reading/writing scripts (`load_untouch_nii`, `load_nii_ext`, etc.).
- Custom deep learning routines like `apply_dncnn_symmetric.m` for spatial denoising.
- Script implementations for Intravoxel Incoherent Motion (IVIM) fitting (`IVIMmodelfit.m`, `IVIM_seg.m`, `IVIM_bayes.m`) and mono-exponential ADC (`fit_adc_mono.m`).
- Functions for DVH processing (`dvh.m`, `sample_rtdose_on_image.m`).
