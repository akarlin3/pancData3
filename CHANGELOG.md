# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- `CONTRIBUTING.md` with contribution guidelines
- `CITATION.cff` for machine-readable academic citation
- `SECURITY.md` with vulnerability reporting policy
- `CHANGELOG.md` for tracking version history
- `.editorconfig` for consistent cross-editor formatting
- `.gitattributes` for consistent line endings and diff behavior
- GitHub Actions CI workflow for automated testing
- GitHub issue and pull request templates
- Improved `README.md` with badges, table of contents, and structured documentation

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
- Comprehensive test suite (47 test files) with code coverage
- Octave compatibility shims
- Checkpointing system for large cohort recovery
- Security utilities: `safe_load_mask`, `escape_shell_arg`
- Data leakage prevention: patient-stratified folds, temporal KNN bounds, DL provenance tracking
- Configuration-driven pipeline via `config.json`
