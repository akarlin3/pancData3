# CLAUDE.md — AI Assistant Guide for pancData3

This file provides essential context for AI assistants (Claude, Antigravity, etc.) working in this repository.

For detailed module tables, test file lists, utility descriptions, and analysis script references, see [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md).

---

## Project Overview

**pancData3** is a MATLAB-based analysis pipeline for pancreatic DWI (Diffusion-Weighted Imaging) research at Memorial Sloan Kettering Cancer Center. It processes MRI data to:

- Fit IVIM (Intravoxel Incoherent Motion) and ADC (Apparent Diffusion Coefficient) models
- Apply Deep Learning denoising (DnCNN, IVIMnet)
- Correlate findings with RT (Radiotherapy) dose maps
- Perform survival analysis, competing risks modeling, and treatment response prediction

**Language:** MATLAB (R2021a+)
**License:** AGPL-3.0 (Copyright 2026 Avery Karlin)
**Domain:** Medical Physics / Oncology Research
**Platforms:** Windows 10/11, macOS 13+, Linux (Ubuntu 22.04+) — CI-tested on all three

---

## Multi-Agent Collaboration Model

This repository uses a three-agent architecture:

| Agent | Role | Scope |
|---|---|---|
| **Claude Code** (interactive) | Feature implementation, pipeline enhancements, debugging, code review | Runs locally with full repository access |
| **Antigravity** (local) | Core physics modeling, MRI calibration, specialized scripts | Runs locally with access to patient data |

### Critical Safety Rules

> **NEVER send patient data, sensitive CSVs, or PHI to any cloud agent or external service.**
> Only code logic and structure may be shared with cloud-based AI.

- Do **not** modify files in the `pipeline/dependencies/` folder under any circumstances.

---

## Repository Structure

```
pancData3/
├── config.json                          # Active configuration (not committed)
├── config.example.json                  # Configuration template (committed)
├── pipeline/                            # MATLAB pipeline
│   ├── run_dwi_pipeline.m              # Master orchestrator — main entry point
│   ├── execute_all_workflows.m         # Runs all 3 DWI types sequentially
│   ├── patient_data_check.m            # Pre-pipeline data integrity scanner
│   ├── core/                           # Primary pipeline modules (18 files)
│   ├── utils/                          # Helper utilities (61 files)
│   ├── .octave_compat/                 # Octave compatibility shims (21 files)
│   ├── tests/                          # Full test suite (100 test files)
│   │   ├── run_all_tests.m             # MATLAB unittest test runner
│   │   ├── benchmarks/                 # Performance benchmarks (7 files)
│   │   └── diagnostics/                # Diagnostic spot-check scripts (5 files)
│   └── dependencies/                   # Third-party scripts — DO NOT MODIFY
├── analysis/                            # Python post-hoc analysis scripts
│   ├── run_analysis.py                 # Orchestrator entry point
│   ├── shared.py                       # Shared utilities and config loading
│   ├── analysis_config.json            # Centralised configuration
│   ├── parsers/                        # Log/CSV/MAT/vision parsing (4 files)
│   ├── cross_reference/                # Cross-DWI comparison scripts (5 files)
│   ├── report/                         # HTML+PDF report generation
│   │   ├── generate_report.py          # Report orchestrator
│   │   ├── report_formatters.py        # Formatting utilities
│   │   ├── report_constants.py         # CSS, JS, references, templates
│   │   └── sections/                   # Section builders (17 files)
│   └── tests/                          # Python test suite — 30 test files, 1186 tests (pytest)
├── .agents/
│   ├── rules/physics_rules.md          # Agent safety and delegation rules
│   └── workflows/run_data.md           # Structured /run_data workflow definition
├── README.md                            # Human-facing documentation
├── CLAUDE.md                            # This file (essentials for every session)
└── CLAUDE_REFERENCE.md                  # Detailed module/test/utility reference tables
```

---

## Configuration

Copy `config.example.json` to `config.json` and fill in site-specific paths before running the pipeline.

Key fields:

```json
{
  "dataloc": "/path/to/pancreas_dwi/",
  "dcm2nii_call": "/path/to/dcm2niix",
  "clinical_data_sheet": "MASTER_pancreas_DWIanalysis.xlsx",
  "clinical_sheet_name": "Clin List_MR",
  "skip_to_reload": false,
  "skip_tests": false,
  "ivim_bthr": 100,
  "adc_thresh": 0.001,
  "high_adc_thresh": 0.00115,
  "d_thresh": 0.001,
  "f_thresh": 0.1,
  "dstar_thresh": 0.01,
  "min_vox_hist": 100,
  "adc_max": 0.003,
  "use_checkpoints": true,
  "patient_ids": [],
  "dwi_type": "Standard",
  "td_scan_days": [],
  "cause_of_death_column": "",
  "clear_cache": true,
  "core_method": "adc_threshold",
  "core_percentile": 25,
  "core_n_clusters": 2,
  "fdm_parameter": "adc",
  "fdm_thresh": 0.0004,
  "spectral_min_voxels": 20,
  "run_compare_cores": false,
  "run_all_core_methods": false,
  "store_core_masks": false,
  "use_firth_refit": true,
  "compute_fine_gray": true,
  "exclude_motion_volumes": false,
  "use_texture_features": false,
  "use_gpu": false,
  "gpu_device": 1
}
```

`dwi_type` must be one of: `"Standard"`, `"dnCNN"`, or `"IVIMnet"`.

`core_method` selects the tumor core delineation algorithm. Options: `"adc_threshold"` (default), `"d_threshold"`, `"df_intersection"`, `"otsu"`, `"gmm"`, `"kmeans"`, `"region_growing"`, `"active_contours"`, `"percentile"`, `"spectral"`, `"fdm"`. The `percentile`, `spectral`, and `fdm` methods use a unified core mask across all parameters (replacing individual D/f/D* thresholds).

### Config Backwards Compatibility (mandatory)

When adding a new config field, you **must** add a corresponding default in `pipeline/utils/parse_config.m` using the existing `isfield` + fallback pattern so that config files without the new field continue to work unchanged. Do **not** modify `config.json` to add a new field without this default in place.

When removing a config field, you **must** ensure that all code referencing the field is updated so that existing config files still containing the removed field do not cause errors (e.g., the field is simply ignored). Do **not** remove a field from `config.example.json` or `pipeline/utils/parse_config.m` without verifying that every consumer of that field has been cleaned up.

If a change (addition or removal) truly cannot be made backwards-compatible, you **must** ask the user for explicit permission before proceeding.

---

## Running the Pipeline

### Full multi-type run (all 3 DWI types)

```matlab
cd pipeline
execute_all_workflows
```

This creates a timestamped output folder (`saved_files_YYYYMMDD_HHMMSS/`), cleans up stale parallel jobs, sets up a parallel pool (max 2 workers), runs Standard → dnCNN → IVIMnet sequentially, modifying `config.json` between runs and skipping reload phases after the first. All console output and logs are saved to the timestamped folder.

### Single step / targeted run

```matlab
addpath('pipeline/core', 'pipeline/utils', 'pipeline/dependencies');

% Run only the 'load' step
run_dwi_pipeline('config.json', {'load'});

% Run load + sanity checks
run_dwi_pipeline('config.json', {'load', 'sanity'});

% Run with a pre-existing output folder (used by execute_all_workflows)
run_dwi_pipeline('config.json', {'load', 'visualize'}, 'path/to/saved_files_folder');
```

The optional 3rd argument specifies a parent output folder. If omitted, `run_dwi_pipeline` creates its own timestamped folder. When called from `execute_all_workflows`, the same timestamped folder is reused across all three DWI type runs.

### Available pipeline steps (in order)

1. `load` — DICOM conversion, model fitting, checkpointing
2. `sanity` — Convergence validation, missingness, spatial alignment
3. `visualize` — Parameter maps, distributions, trajectories
4. `metrics` — Summary metrics, baseline, longitudinal, dosimetry, stats, survival
5. `compare_cores` — Pairwise comparison of all 11 core methods (Dice, Hausdorff, volumes). Not in default steps; invoke explicitly: `run_dwi_pipeline('config.json', {'compare_cores'})`

### Agent workflow (`/run_data`)

The structured `/run_data` workflow in `.agents/workflows/run_data.md` does:
1. Verifies `config.json` exists
2. Runs `pipeline/tests/test_dwi_pipeline.m` for stability check
3. Executes `matlab -batch "cd pipeline; execute_all_workflows"` as a background command
4. Waits for completion and verifies output MAT files and figures

---

## Running Tests

```matlab
run('pipeline/tests/run_all_tests.m')
```

- Uses MATLAB's built-in `unittest` framework
- Auto-discovers all test files in `pipeline/tests/` (including subfolders)
- Generates code coverage report for `pipeline/core/` and `pipeline/utils/`
- Compatible with MATLAB R2017b+
- Returns non-zero exit code on failure (CI-safe)

### Test categories

| Location | Description | Run by CI? |
|---|---|---|
| `pipeline/tests/` | Full functional and integration tests | Yes |
| `pipeline/tests/benchmarks/` | Performance benchmarks | Yes |
| `pipeline/tests/diagnostics/` | Diagnostic spot-check tests | Yes |

### Key test files

| File | What it covers |
|---|---|
| `test_dwi_pipeline.m` | Integration: imputation, IPCW, competing risks, leakage |
| `test_load_dwi_data.m` | Data loading and checkpointing |
| `test_sanity_checks.m` | Convergence and spatial alignment |
| `test_metrics_baseline.m` | Baseline metric computation |
| `test_grouped_folds.m` | Cross-validation fold generation (no leakage) |
| `test_knn_temporal_leakage.m` | KNN imputation temporal safety |
| `test_safe_load_mask.m` | Secure mask file loading |
| `test_escape_shell_arg.m` | Cross-platform shell argument escaping |
| `test_visualize_smoke.m` | Visualization output smoke tests |
| `test_fit_models.m` | IVIM + ADC model fitting (dimensions, padding, b-value validation, known value recovery) |
| `test_core_methods.m` | Tumor core extraction method validation (all 11 methods) |
| `test_new_core_methods.m` | Detailed tests for percentile, spectral, and fDM core methods |
| `test_compute_dice_hausdorff.m` | Dice/Hausdorff distance computation tests |
| `test_json_set_field.m` | JSON field replacement utility tests |
| `test_cross_dwi_subvolume.m` | Cross-DWI-type subvolume comparison tests |
| `test_landmark_cindex.m` | Landmark concordance index tests |
| `test_source_code_standards.m` | Source code standards enforcement |
| `test_modularity.m` | Module independence and interface tests |
| `test_statistical_methods.m` | Statistical methods validation |
| `test_compare_core_methods.m` | Core method pairwise comparison validation |
| `test_multi_core_methods.m` | Multi-method core integration and backward compatibility |
| `test_process_single_scan.m` | Per-scan pipeline processing (init, NaN defaults, struct layout) |
| `test_metrics_stats_comparisons.m` | Wilcoxon rank-sum, BH FDR correction, GLME mixed-effects |
| `test_format_p_value.m` | P-value formatting for biomedical reporting conventions |
| `test_pipeline_progress_gui.m` | Pipeline progress bar wrapper (step mapping, lifecycle) |
| `test_plot_feature_distribution.m` | Feature distribution visualization (histogram/boxplot modes) |
| `test_filter_collinear_features.m` | Collinearity pruning, AUC tie-breaking, time-stratification |
| `test_metrics_dosimetry.m` | Dosimetry metric computation (D95, V50, DVH) |
| `test_metrics_longitudinal.m` | Longitudinal change analysis |
| `test_metrics_survival.m` | Survival analysis, competing risks, IPCW |
| `test_metrics_stats_predictive.m` | Predictive modeling, elastic net, LOOCV |
| `test_compute_summary_metrics.m` | Voxel-to-summary-metric aggregation |
| `test_build_td_panel.m` | Time-dependent panel construction |
| `test_scale_td_panel.m` | Timepoint-specific feature scaling |
| `test_parse_config.m` | Config loading, defaults, backwards compatibility |
| `test_discover_patient_files.m` | Patient file system navigation |
| `test_apply_dir_mask_propagation.m` | Deformable image registration and mask alignment |
| `test_calculate_subvolume_metrics.m` | Dose coverage within diffusion-defined subvolumes |
| `test_init_scan_structs.m` | Scan data structure initialization |
| `test_compute_scan_days_from_dates.m` | DICOM-derived scan day computation |
| `test_perform_statistical_test.m` | Wilcoxon rank-sum with NaN-safe extraction |
| `test_remove_constant_columns.m` | Zero-variance and all-NaN column removal |
| `test_execute_pipeline_step.m` | Non-fatal step executor: success, error, warning, diary |
| `test_load_data_from_disk.m` | DWI vector loading with legacy fallback |
| `test_compute_adc_metrics.m` | ADC metric computation (volume, sub-volume, histogram, KS) |
| `test_compute_ivim_metrics.m` | IVIM metric computation (failed-fit filtering, unified methods) |
| `test_benjamini_hochberg_fdr.m` | BH FDR correction: step-up procedure, q-value capping, monotonicity, order preservation |
| `test_compute_ipcw_weights.m` | IPCW weights: censoring model, mean normalization, competing event exclusion, floor truncation |
| `test_assemble_predictive_features.m` | Feature matrix assembly: 22-column layout, post-treatment dose exclusion, NaN column removal |
| `test_compute_multi_core_metrics.m` | Multi-method core metrics: all 11 methods, unified mask sharing, fDM volume fractions |
| `test_compute_spatial_repeatability.m` | Spatial repeatability: Dice/Hausdorff across Fx1 repeats, 12-output validation |
| `test_gpu_available.m` | GPU detection utility: availability check, graceful fallback, invalid device handling |
| `test_clear_pipeline_cache.m` | Cache clearing: deletion, protection, sentinel, once-per-session guard |
| `test_setup_output_folders.m` | Output folder creation: explicit reuse, timestamped auto-creation, sentinel |
| `test_load_baseline_from_disk.m` | Baseline loading: field access, missing file error |
| `test_resolve_scan_days.m` | Scan day resolution: DICOM preferred, config fallback, empty fallback |
| `test_schoenfeld_residuals.m` | Schoenfeld residuals: PH holds, time-varying effect, few events, diagnostic figure |
| `test_fine_gray.m` | Fine-Gray model: competing events, no competing events, all censored, CIF plot |
| `test_calibration_metrics.m` | Model calibration: perfect/miscalibrated, single class, Brier decomposition |
| `test_bootstrap_ci.m` | Bootstrap CI: normal mean, ordering, degenerate input, NaN resilience, median |
| `test_compute_texture_features.m` | Texture features: checkerboard, uniform, field count, 3D, empty mask |
| `test_compute_registration_quality.m` | Registration quality: identity transform, known shift, mutual information |
| `test_detect_motion_artifacts.m` | Motion artifacts: clean DWI, signal dropout, single slice, output structure |

---

## Core Modules (`pipeline/core/`)

| File | Purpose |
|---|---|
| `load_dwi_data.m` | DICOM loading, DL denoising, model fitting, checkpointing via `parfor` |
| `sanity_checks.m` | NaN/Inf/negative detection, outlier summary, spatial alignment |
| `visualize_results.m` | Parameter map overlays, feature distributions, longitudinal trajectories |
| `compute_summary_metrics.m` | Voxel-to-summary-metric aggregation |
| `metrics_baseline.m` | Baseline measures, outlier cleaning, percent delta |
| `metrics_longitudinal.m` | Longitudinal change analysis |
| `metrics_dosimetry.m` | Dose-related metrics (D95, V50) |
| `metrics_stats_comparisons.m` | Statistical group comparisons |
| `metrics_stats_predictive.m` | Predictive modeling, feature selection |
| `metrics_survival.m` | Survival analysis, competing risks (Cause-Specific Hazards) |
| `discover_patient_files.m` | File system navigation for patient cohort |
| `convert_dicom.m` | DICOM-to-NIFTI via `dcm2niix` |
| `fit_models.m` | ADC mono-exponential and IVIM model fitting |
| `plot_parameter_maps.m` | Visualization helpers for parameter overlays |
| `plot_feature_distributions.m` | Feature histogram/boxplot rendering |
| `plot_scatter_correlations.m` | Correlation scatter plots |
| `process_single_scan.m` | Per-scan DICOM conversion, model fitting, and caching |
| `compare_core_methods.m` | Pairwise comparison of all 11 tumor core methods (Dice, Hausdorff, volume) |

---

## Utility Modules (`pipeline/utils/`)

| File | Purpose |
|---|---|
| `parse_config.m` | Loads and validates `config.json` |
| `build_td_panel.m` | Constructs longitudinal time-dependent data panels |
| `scale_td_panel.m` | Timepoint-specific feature scaling (prevents cross-timepoint leakage) |
| `make_grouped_folds.m` | Patient-stratified cross-validation folds |
| `knn_impute_train_test.m` | KNN imputation with strict temporal leakage bounds |
| `filter_collinear_features.m` | Prunes collinear features before model fitting |
| `apply_dir_mask_propagation.m` | Deformable image registration and rigid alignment |
| `safe_load_mask.m` | Securely loads `.mat` mask files (rejects unsafe variable classes) |
| `calculate_subvolume_metrics.m` | Dose coverage metrics within diffusion-defined GTV subvolume |
| `load_dl_provenance.m` | Loads DL training provenance to guard against data leakage |
| `nanmean_safe.m` | Octave-compatible NaN-ignoring mean |
| `nanstd_safe.m` | Octave-compatible NaN-ignoring standard deviation |
| `perform_statistical_test.m` | Wilcoxon rank-sum testing with NaN-safe group extraction |
| `parsave_dir_cache.m` | Parallel-safe `save` wrapper for `parfor` caching |
| `escape_shell_arg.m` | Cross-platform shell argument escaping (Windows and Unix) |
| `discover_gtv_file.m` | Locates GTV mask file with flexible naming patterns |
| `find_gtv_files.m` | Locates GTVp and GTVn masks for complex tumor anatomy |
| `plot_feature_distribution.m` | Histogram/boxplot with ANOVA p-value annotation |
| `init_scan_structs.m` | Initializes scan data structures for pipeline processing |
| `compute_scan_days_from_dates.m` | Derives scan days from DICOM acquisition dates |
| `detect_baseline_outliers.m` | Outcome-blinded IQR outlier detection for baseline metrics |
| `format_p_value.m` | Formats p-values for display with appropriate precision |
| `remove_constant_columns.m` | Removes constant/NaN-only columns from feature matrices |
| `parfor_progress.m` | Parallel loop progress reporting |
| `text_progress_bar.m` | Text-based progress bar display |
| `extract_tumor_core.m` | Configurable tumor core delineation (11 methods) |
| `PipelineProgressGUI.m` | Pipeline-aware progress bar wrapper (maps step keys to display names) |
| `ProgressGUI.m` | Professional custom-figure progress bar for MATLAB pipelines |
| `compute_dice_hausdorff.m` | Dice coefficient and Hausdorff distance between 3D binary masks |
| `compute_histogram_laplace.m` | Laplace-smoothed histogram probability distribution |
| `json_set_field.m` | Targeted regex replacement of a field value in raw JSON strings |
| `plot_cross_dwi_subvolume_comparison.m` | Cross-DWI-type ADC subvolume comparison visualization |
| `compute_adc_metrics.m` | ADC summary metrics for a single patient/timepoint/DWI-type (extracted from compute_summary_metrics) |
| `compute_ivim_metrics.m` | IVIM (D/f/D*) summary metrics for a single patient/timepoint/DWI-type (extracted from compute_summary_metrics) |
| `compute_kurt_skew.m` | Kurtosis/skewness computation with minimum sample guard |
| `compute_spatial_repeatability.m` | Dice and Hausdorff spatial repeatability between Fx1 repeat sub-volumes |
| `compute_multi_core_metrics.m` | Multi-method (11 core methods) sub-volume metrics per patient/timepoint |
| `compute_percent_deltas.m` | Treatment-induced percent/absolute changes from baseline |
| `assemble_predictive_features.m` | Builds 22-column feature matrix for elastic net (extracted from metrics_stats_predictive) |
| `run_elastic_net_cv.m` | 5-fold elastic net CV + final model fitting (extracted from metrics_stats_predictive) |
| `run_loocv_risk_scores.m` | Nested LOOCV for unbiased out-of-fold risk scores (extracted from metrics_stats_predictive) |
| `plot_predictive_diagnostics.m` | ROC curve, sanity check panels, and 2D scatter plots (extracted from metrics_stats_predictive) |
| `execute_pipeline_step.m` | Generic non-fatal pipeline step executor with try-catch, diary, GUI, warning logging (extracted from run_dwi_pipeline) |
| `initialize_pipeline.m` | Pipeline initialization: path setup, pre-flight tests, toolbox license checks (extracted from run_dwi_pipeline) |
| `load_data_from_disk.m` | Load DWI vectors and summary metrics from disk with legacy fallback (extracted from run_dwi_pipeline) |
| `normalize_patient_ids.m` | Octave-compatible patient ID normalization for spreadsheet/folder matching |
| `select_dwi_vectors.m` | Extract ADC/D/f/D* voxel vectors by DWI processing type (Standard/dnCNN/IVIMnet) |
| `write_sentinel_file.m` | Write pipeline step completion sentinel files |
| `benjamini_hochberg_fdr.m` | Benjamini-Hochberg FDR correction for multiple hypothesis testing |
| `compute_ipcw_weights.m` | Inverse probability of censoring weights for Cox PH survival models |
| `gpu_available.m` | GPU availability detection with graceful CPU fallback |
| `clear_pipeline_cache.m` | Remove pipeline-generated .mat cache files (once-per-session guard, protected files, sentinel checks) |
| `setup_output_folders.m` | Create or reuse the master pipeline output folder (timestamped auto-creation with sentinel) |
| `load_baseline_from_disk.m` | Load persisted metrics_baseline outputs from .mat file |
| `resolve_scan_days.m` | Three-level scan day resolution for survival analysis (DICOM dates -> config -> defaults) |
| `bootstrap_ci.m` | BCa bootstrap confidence intervals for arbitrary scalar metric functions |
| `compute_schoenfeld_residuals.m` | Scaled Schoenfeld residuals and PH assumption testing via Spearman correlation |
| `compute_calibration_metrics.m` | Calibration assessment: Brier score, Hosmer-Lemeshow test, calibration slope/intercept |
| `compute_texture_features.m` | First-order and GLCM texture feature extraction from parameter maps |
| `compute_registration_quality.m` | Registration quality metrics: Jacobian determinant, NCC, mutual information |
| `detect_motion_artifacts.m` | DWI volume quality assessment: CV, NMI, signal dropout detection |

### Octave Compatibility (`pipeline/.octave_compat/`)

Contains 21 shim files for GNU Octave compatibility, including:

- `@table/` class implementation (`table.m`, `subsasgn.m`, `subsref.m`, `display.m`)
- `+matlab/+unittest/` namespace shims (`TestSuite.m`, `TestCase.m`, `TestRunner.m`)
- `+matlab/+unittest/+fixtures/` shim (`PathFixture.m`)
- `+matlab/+unittest/+plugins/` shim (`CodeCoveragePlugin.m`)
- Standard function replacements: `cvpartition.m`, `nanmean.m`, `nanstd.m`, `categorical.m`, `niftiread.m`, `niftiwrite.m`, `niftiinfo.m`, `fitglme.m`, `contains.m`, `sgtitle.m`, `yline.m`, `spectralcluster.m`

### Analysis Scripts (`analysis/`)

Python scripts for post-hoc analysis of pipeline outputs, organized into subpackages. The suite includes vision-based graph analysis (via Google Gemini API), direct log/CSV parsing, cross-DWI comparison, and automated HTML/PDF report generation.

**Requirements:** Python 3.12+, `google-genai`, `pydantic`, `tqdm`, `weasyprint` (install via `pip install -r analysis/requirements.txt`). Vision analysis requires `GEMINI_API_KEY` environment variable; PDF generation requires `weasyprint`; all other scripts work without these optional dependencies. All scripts display `tqdm` progress bars during processing.

**Configuration:** All analysis scripts share a centralised config loaded by `shared.load_analysis_config()`. Defaults are built into `shared.py`; overrides come from `analysis/analysis_config.json` (committed) and optionally from the MATLAB `config.json` (for `dwi_type`). The `run_analysis.py` orchestrator also accepts `--gemini-model`, `--concurrency`, `--config`, `--skip-checks`, and `--interactive` CLI flags. By default, the orchestrator verifies that all `requirements.txt` packages are installed and runs the full pytest suite before starting the analysis pipeline; `--skip-checks` bypasses these pre-flight checks.

| File | Purpose |
|---|---|
| `run_analysis.py` | Orchestrator: runs the full analysis workflow with `--folder`, `--skip-vision`, `--report-only`, `--no-pdf`, `--html`, `--skip-checks`, `--interactive` flags; verifies requirements and runs tests before starting |
| `analysis_config.json` | Centralised configuration: vision model, concurrency, statistical thresholds, priority graphs |
| `shared.py` | Shared utilities: folder discovery, DWI type parsing, p-value/correlation regex extraction, config loading |
| `parsers/batch_graph_analysis.py` | Async batch processing of all graph images via Google Gemini vision API; outputs structured CSV with axes, trends, inflection points, statistical tests, outliers, reference lines, clinical relevance, and metadata |
| `parsers/parse_log_metrics.py` | Direct parsing of MATLAB log files: Wilcoxon p-values, AUC, hazard ratios, GLME interaction terms, sanity check convergence/alignment |
| `parsers/parse_csv_results.py` | Direct parsing of pipeline CSV exports (Significant_LF_Metrics.csv, FDR_Sig_Global.csv) with cross-DWI comparison |
| `parsers/parse_mat_metrics.py` | Parses MATLAB `.mat` output files (core comparison, dosimetry, summary metrics) into JSON for downstream analysis |
| `report/generate_report.py` | HTML+PDF report orchestrator: data loading, section assembly, CLI entry point for `analysis_report.html` and `analysis_report.pdf` |
| `report/report_formatters.py` | Formatting utilities for the HTML report (escaping, badges, nav bar, stat cards, forest plot cells, effect size helpers, table/figure numbering, figure captions, citation system, manuscript sentence helpers) |
| `report/report_constants.py` | Large constants extracted from report_formatters (CSS stylesheet, JavaScript, publication references with BibTeX, HTML template) |
| `report/generate_interactive_report.py` | Interactive HTML report with client-side filtering, Chart.js visualisations, patient drill-down, sortable tables, and DWI/core-method comparison |
| `report/interactive_constants.py` | CSS and JavaScript constants for the interactive report (sidebar, tabs, chart rendering, filter logic) |
| `report/sections/` | Section builder package for the HTML report, split into 16 submodules: `metadata.py`, `main_results.py`, `statistical_reporting.py`, `manuscript.py`, `enrollment.py`, `supplemental.py`, `gallery.py`, `graph_overview.py`, `cross_dwi.py`, `correlations.py`, `effect_sizes.py`, `model_diagnostics.py`, `power_analysis.py`, `discussion.py`, `publication.py`, `_helpers.py` (shared utility functions). Legacy shims (`analysis_sections.py`, `statistics.py`, `data_sections.py`) re-export for backward compatibility. |
| `cross_reference/cross_reference_dwi.py` | Full cross-DWI comparison (Standard vs dnCNN vs IVIMnet) of trends, inflection points, and summaries |
| `cross_reference/cross_reference_summary.py` | Concise cross-DWI summary focusing on priority clinical graphs and trend agreement/disagreement |
| `cross_reference/statistical_relevance.py` | Extracts p-values and correlation coefficients; reports significant findings, notable correlations, and cross-DWI significance |
| `cross_reference/statistical_by_graph_type.py` | Filters statistical findings by graph type (scatter, box, line, heatmap, bar, histogram, parameter_map) |
| `cross_reference/cross_dwi_agreement.py` | Bland-Altman, Lin's CCC, and ICC agreement analysis between DWI types |
| `report/sections/forest_plot.py` | Forest plot section builder: HR extraction, matplotlib forest plot, report integration |

**Python Test Suite (pytest):** 30 test files with 1186 tests in `analysis/tests/`. Run with `cd analysis/tests && python -m pytest -v`.

| File | What it covers |
|---|---|
| `conftest.py` | Shared fixtures: synthetic saved_files directories, graph CSVs, log files, pipeline CSV exports |
| `test_shared.py` | DWI type parsing, p-value/correlation extraction, CSV loading, folder resolution |
| `test_parse_log_metrics.py` | GLME, ROC/AUC, survival, baseline, sanity check regex parsing; integration with log files |
| `test_parse_csv_results.py` | CSV reading, cross-DWI significance consistency analysis |
| `test_batch_graph_analysis.py` | Image collection, base64 encoding, MIME types, Pydantic schemas (Axis, Trend, InflectionPoint, StatisticalTest, Outlier, ReferenceLine, GraphAnalysis), CSV flattening |
| `test_generate_report_helpers.py` | Formatting helpers: series normalization, significance tags, section headers, effect sizes, copy buttons, figure captions |
| `test_generate_report_integration.py` | Full HTML report generation, data quality, Cox PH direction, correlations, new sections integration |
| `test_generate_report_manuscript.py` | Manuscript findings, reporting checklist, table/figure index, BibTeX export, results draft, journal guide |
| `test_generate_report_sections.py` | Section builders: forest plots, data completeness, feature overlap, power analysis, patient flow, sensitivity, appendix, figure gallery |
| `test_interactive_report.py` | Interactive report: HTML escaping, DWI badges, significance classes, trend tags, patient extraction, core method extraction, section builders (overview, patient explorer, visualisations, significance, graph explorer, core comparison, dosimetry), Chart.js integration, JSON data blob, sidebar filters, sortable tables |
| `test_treatment_plan.py` | Suggested treatment plan: core recommendations, survival/predictive integration, timing guidance, backward compatibility |
| `test_script_outputs.py` | stdout-based tests for cross_reference, statistical, and run_analysis scripts |
| `test_report_formatters.py` | HTML escaping, significance markers, DWI badges, trend tags, effect sizes, consensus, figure captions, nav sections |
| `test_parse_mat_metrics.py` | MAT file parsing, dosimetry/core/longitudinal extraction, scipy graceful degradation |
| `test_analysis_config.py` | Config loading, deep merge, layered overrides, MATLAB config integration, caching |
| `test_report_sections_helpers.py` | Report helper functions: JSON loading, series normalization, cohort size, AUC finding, trend agreement |
| `test_report_sections_metadata.py` | Report metadata sections: cover page, part breaks, TOC, publication header, data availability |
| `test_report_sections_main_results.py` | Main results sections: executive summary, hypothesis, statistical significance, treatment response |
| `test_report_sections_data_sections.py` | Data sections: cohort_overview (cohort, patient flow, completeness), supplemental (MAT data), gallery (appendix, figure gallery) |
| `test_report_sections_analysis.py` | Analysis sections: graph_overview (overview, issues, stats by type), cross_dwi (comparison, feature overlap), correlations |
| `test_report_sections_statistics.py` | Statistics sections: effect_sizes (HR effects, multiple comparisons), model_diagnostics (diagnostics, sensitivity), power_analysis |
| `test_report_sections_discussion.py` | Discussion sections: methods, limitations, conclusions, reporting checklist, journal guide |
| `test_cross_reference_dwi.py` | Cross-DWI comparison output: graph matching, trend display, stat test formatting, truncation, JSON robustness |
| `test_cross_reference_summary.py` | Cross-DWI summary: trend agreement/disagreement, priority ordering, parameter maps, inflection points |
| `test_statistical_relevance.py` | Statistical findings: p-value extraction, Bonferroni correction, significance markers, correlations, cross-DWI comparison |
| `test_statistical_by_graph_type.py` | Per-graph-type analysis: grouping, trend directions, top-5 non-sig, density/comparison aggregation, summary table |
| `test_xref_unit.py` | Cross-reference correctness: safe_text, p-value/correlation edge cases, trend agreement logic, significance markers, Bonferroni, direction classification, priority ordering |
| `test_integration.py` | End-to-end analysis pipeline integration: runs run_analysis.py on synthetic data, verifies HTML output sections and tables |
| `test_api_connection.py` | Gemini API connection smoke test (skipped without API key) |
| `test_cross_dwi_agreement.py` | Bland-Altman, Lin's CCC, ICC agreement analysis tests |
| `test_forest_plot.py` | HR data extraction and forest plot generation tests |
For the full list of 100 MATLAB test files and 30 Python test files with descriptions, see [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md#key-matlab-test-files).

---

## Code Conventions

### Architecture

- **Orchestrator pattern**: `pipeline/run_dwi_pipeline.m` sequences all steps; modules are independently callable.
- Explicit parameter passing between modules — no shared global workspace state.
- Outputs written to a timestamped folder; path passed as argument, not hardcoded.

### Parallelization

- `parfor` loops process patient cohorts in parallel.
- `parsave_dir_cache.m` enables safe `save` inside `parfor`.
- Checkpointing allows recovery from interruptions mid-cohort.
- Parallel pool capped at 2 workers (`pipeline/execute_all_workflows.m`).

### Cross-Platform Compatibility

- `escape_shell_arg.m` auto-detects `ispc()` for Windows (double-quote) vs Unix (single-quote) shell escaping.
- All file paths use `fullfile()`, `filesep`, and `pathsep` — never hardcoded separators.
- Python analysis scripts use `pathlib.Path` throughout and reconfigure `stdout` to UTF-8 on Windows for emoji support.
- CI runs the full MATLAB and Python test suites on Linux, macOS, and Windows.

### Security

- `safe_load_mask.m` inspects variable type before loading; rejects non-numeric classes to prevent arbitrary code execution from `.mat` files.
- `escape_shell_arg.m` must be used for all paths passed to `system()`.
- Never pass unsanitized user strings to shell commands.

### File Deletion Safety

The pipeline must **never delete a file or directory it did not create**. All deletion sites use provenance verification:

- **Sentinel files**: Pipeline-created directories (`saved_files_*`, `processed_patients/`) contain a `.pipeline_created` sentinel file. Before `rmdir`, code must verify this sentinel exists.
- **Cache clearing** (`pipeline/run_dwi_pipeline.m`): Only deletes `.mat` files matching pipeline-generated patterns (`dwi_vectors_*.mat`, `summary_metrics_*.mat`, `adc_vectors.mat`). Manually curated files are protected via the `protected_files` list.
- **Lock files**: Only deleted when orphaned (stale from a crashed worker) or after successful checkpoint completion.
- **Diary/log files**: Only deleted immediately before being recreated by the same module.
- **Test cleanup**: Must only remove artifacts the test itself created — use pre/post directory snapshots or sentinel checks.
- When writing new code that deletes files or directories, always add a provenance check: verify a `.pipeline_created` sentinel or confirm the file was created in the same function scope.

### Data Leakage Prevention

- `make_grouped_folds.m` — patient-stratified folds prevent intra-patient leakage across CV splits.
- `knn_impute_train_test.m` — enforces infinite distances for rows from the same patient.
- `scale_td_panel.m` — scaling is strictly timepoint-specific.
- `load_dl_provenance.m` — guards against overlap between DL training and analysis patient sets.

### Logging Style

Progress output uses emoji prefixes for clarity:
- `🚀` — pipeline start
- `⚙️` — processing step
- `✅` — success
- `❌` — fatal error / pipeline halt
- `⚠️` — non-fatal warning (module failed but pipeline continues)
- `📁` — file output
- `💡` — informational note

### Diary / Console Logging Architecture

All console output is captured to log files via MATLAB's `diary` command. MATLAB only supports one active diary at a time, so the architecture uses a restart pattern:

1. **`execute_all_workflows.m`** opens a master diary (`execute_all_workflows.log`) in the timestamped output folder.
2. **`run_dwi_pipeline.m`** opens a per-DWI-type diary (`pipeline_log_{type}.txt`) in the type subfolder.
3. Each **core module** (sanity_checks, visualize_results, metrics_*) opens its own diary, overriding the orchestrator's.
4. After each module returns, the orchestrator **restarts** its diary to resume capturing.

For the full output folder structure, see [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md#diary--console-logging--output-folder-structure).

**Important for tests:** Any test that exercises a core module must call `diary off;` in its `TestMethodTeardown` before calling `rmdir` on temp directories, because the module's diary file will still be open (locked on Windows).

### Error Handling

- Validate MATLAB toolbox licenses at startup (Statistics + Image Processing).
- Halt with a clear message on any unrecoverable error.
- Non-fatal module failures (compare_core_methods, metrics_longitudinal, dosimetry, stats, survival) log a `⚠️` warning and continue the pipeline.
- Graceful halting is preferred over silent continuation with bad data.

---

## Dependencies

### Required MATLAB Toolboxes

- Statistics and Machine Learning Toolbox
- Image Processing Toolbox

### External Tools

- `dcm2niix` (MRIcroGL) — DICOM-to-NIFTI conversion; must be on system path or specified in `config.json` as `dcm2nii_call`.

### `pipeline/dependencies/` folder (DO NOT MODIFY)

Contains third-party scripts. Treat as read-only. For the full file listing, see [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md#pipelinedependencies-contents-do-not-modify).

---

## Module Reference

For detailed tables of all core modules (18 files), utility modules (48 files), Octave compatibility shims (21 files), analysis scripts, and Python test files, see [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md).

---

## Documentation Maintenance (mandatory)

After **every feature implementation** (adding a new file, adding a config field, changing module signatures, adding tests, etc.), you **must** update all affected markdown documentation files before considering the task complete. This is not optional — stale documentation is a recurring problem.

### Files to check and update

| File | What to update |
|---|---|
| `CLAUDE.md` | File counts in Repository Structure, config example block |
| `CLAUDE_REFERENCE.md` | Module tables (Core/Utils/Analysis), key test files, Octave compat listing |
| `README.md` | File counts (test badge, Repository Structure tree, utils/tests/analysis counts), config field table if a user-facing field was added |
| `MEMORY.md` (auto-memory) | File signatures, new patterns, any architectural decisions made during the feature |

### Checklist (run mentally after each feature)

1. **New `.m` file in `pipeline/core/` or `pipeline/utils/`?** → Add to the corresponding CLAUDE_REFERENCE.md module table with a one-line purpose. Update the file count in Repository Structure (CLAUDE.md and README.md).
2. **New test file?** → Update test file count in Repository Structure (CLAUDE.md and README.md) and the README test badge number. If it covers a notable area, add to CLAUDE_REFERENCE.md Key test files table.
3. **New config field?** → Add to the CLAUDE.md config JSON block. If user-facing, add to README.md config field table. Ensure `pipeline/utils/parse_config.m` default exists (per Config Backwards Compatibility rules).
4. **New `pipeline/.octave_compat/` shim?** → Add to the CLAUDE_REFERENCE.md Octave Compatibility listing and update the file count.
5. **New top-level `.m` file in `pipeline/`?** → Add to CLAUDE.md Repository Structure tree.
6. **Changed module signature?** → Update MEMORY.md File Signatures section.
7. **New Python script in `analysis/`?** → Add to CLAUDE_REFERENCE.md Analysis Scripts table. Update the file count in Repository Structure (CLAUDE.md and README.md).

---

## Git Workflow

- Development happens on feature branches; never commit directly to `main`.
- Branch naming follows: `claude/<description>-<session-id>` for Claude branches.
- `.gitignore` excludes: MATLAB autosave files (`*.asv`, `*.m~`, `*.mex*`), imaging data (`*.nii`, `*.dcm`, `*.h5`), CSVs, and clinical spreadsheets (`*.xlsx`) to prevent accidental PHI commits.
- After pushing a branch, Claude should create a pull request targeting `main` using `gh pr create` **only if `gh` is available** (check with `which gh` first). Include a summary of changes and a test plan in the PR body. If `gh` is not installed, skip PR creation silently (do not notify the user).
- CI (`.github/workflows/ci.yml`) runs MATLAB tests and Python analysis tests on all three platforms (Linux, macOS, Windows) plus an Octave smoke test on Linux.

---

## What AI Assistants Should and Should Not Do

### Do
- Read and understand existing code before suggesting changes.
- Use `safe_load_mask` and `escape_shell_arg` when handling file I/O.
- Follow the orchestrator pattern — keep pipeline steps modular and independently callable.
- Preserve checkpointing logic in `load_dwi_data.m`; it is critical for large cohort recovery.
- Consult `config.example.json` before adding new configuration fields.
- **Update documentation after every feature implementation** — see [Documentation Maintenance](#documentation-maintenance-mandatory) above.

### Do Not
- Modify anything in `pipeline/dependencies/`.
- Send patient data, clinical CSVs, or any PHI to cloud services.
- Delete files or directories the pipeline did not create — always verify a `.pipeline_created` sentinel or known pipeline-generated naming pattern before deletion.
- Introduce global variables or persistent workspace state between pipeline steps.
- Bypass the temporal leakage safeguards in imputation or cross-validation.
- Use unsanitized strings in `system()` calls.
- Hard-code file paths — all paths must flow through `config.json`.
- Run the pipeline (`execute_all_workflows`, `run_dwi_pipeline`) without explicit researcher approval. Pipeline runs are initiated by the researcher, not by AI assistants. Running the test suite (`pipeline/tests/run_all_tests.m`) for verification is encouraged.
- Add, remove, or rename fields in `config.json` / `config.example.json` without ensuring **backwards compatibility**: new fields must have a default in `pipeline/utils/parse_config.m` (via `isfield` + fallback), and removed fields must be cleaned up from all consumers so that configs still containing them do not cause errors. If a change truly cannot be made backwards-compatible, you **must** ask the user for explicit permission before proceeding.
