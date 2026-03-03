# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [1.1.0] - 2026-03-03

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

## [1.0.0] - 2026-01-01

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
