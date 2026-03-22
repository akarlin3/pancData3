# CLAUDE.md — AI Assistant Guide for pancData3

This file provides essential context for AI assistants (Claude, Antigravity, etc.) working in this repository.

For detailed module tables, test file lists, utility descriptions, and analysis script references, see [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md).

For running the pipeline, running tests, git workflow, documentation maintenance rules, and the improvement loop, see [CLAUDE_WORKFLOWS.md](CLAUDE_WORKFLOWS.md).

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

This repository uses a multi-agent architecture:

| Agent | Role | Scope |
|---|---|---|
| **Claude Code** (interactive) | Feature implementation, pipeline enhancements, debugging, code review | Runs locally with full repository access |
| **Antigravity** (local) | Core physics modeling, MRI calibration, specialized scripts | Runs locally with access to patient data |
| **Improvement Loop Agents** | Autonomous audit → implement → review → merge (RAG-enhanced) | API-driven, no patient data access |

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
├── analysis_config.json                 # Active analysis configuration (not committed)
├── analysis_config.example.json         # Analysis configuration template (committed)
├── .matlab_version                      # Expected MCR version for Docker builds
├── Dockerfile                           # Multi-stage Docker build (MATLAB Runtime + Python)
├── docker-compose.yml                   # Pipeline + analysis Docker services
├── .dockerignore                        # Docker build exclusions
├── docker/                              # Docker support files
│   └── entrypoint.sh                   #   Container entrypoint (mode dispatch + validation)
├── docs/                                # Additional documentation
│   └── DOCKER.md                       #   Docker usage guide
├── pipeline/                            # MATLAB pipeline
│   ├── run_dwi_pipeline.m              # Master orchestrator — main entry point
│   ├── execute_all_workflows.m         # Runs all 3 DWI types sequentially
│   ├── patient_data_check.m            # Pre-pipeline data integrity scanner
│   ├── core/                           # Primary pipeline modules (18 files)
│   ├── utils/                          # Helper utilities (72 files)
│   ├── .octave_compat/                 # Octave compatibility shims (21 files)
│   ├── tests/                          # Full test suite (122 test files)
│   │   ├── run_all_tests.m             # MATLAB unittest test runner
│   │   ├── benchmarks/                 # Performance benchmarks (7 files)
│   │   └── diagnostics/                # Diagnostic spot-check scripts (5 files)
│   └── dependencies/                   # Third-party scripts — DO NOT MODIFY
├── analysis/                            # Python post-hoc analysis scripts
│   ├── run_analysis.py                 # Orchestrator entry point
│   ├── shared.py                       # Shared utilities and config loading
│   ├── parsers/                        # Log/CSV/MAT/vision parsing (4 files)
│   ├── cross_reference/                # Cross-DWI comparison scripts (5 files)
│   ├── report/                         # HTML+PDF report generation
│   │   ├── generate_report.py          # Report orchestrator
│   │   ├── report_formatters.py        # Formatting utilities
│   │   ├── report_constants.py         # CSS, JS, references, templates
│   │   └── sections/                   # Section builders (37 files)
│   └── tests/                          # Python test suite — 45 test files, 1795 tests (pytest)
├── improvement_loop/                    # Automated audit/fix loop
│   ├── orchestrator_v2.py              #   Pipeline orchestrator (audit → implement → review → merge)
│   ├── orchestrator_v1.py              #   Legacy single-pass orchestrator
│   ├── evaluator.py                    #   Finding schema + audit scoring + exit logic
│   ├── loop_tracker.py                 #   Iteration logging + context generation
│   ├── loop_config.py                  #   Centralised config (LoopConfig dataclass)
│   ├── git_utils.py                    #   Subprocess-based git operations
│   ├── agents/                         #   Agent modules
│   │   ├── auditor.py                 #     Code audit agent (RAG-enhanced)
│   │   ├── implementer.py            #     Fix implementation agent (RAG-enhanced)
│   │   ├── reviewer.py               #     Code review quality gate (RAG-enhanced)
│   │   └── _api.py                   #     Shared API retry helper
│   └── rag/                           #   Retrieval-Augmented Generation
│       ├── chunker.py                 #     Semantic code chunker (MATLAB/Python/MD/JSON)
│       ├── indexer.py                 #     ChromaDB vector index (build/update/query)
│       └── retriever.py               #     Query interface + agent context builders
├── improvement_loop_config.example.json # Improvement loop config template (committed)
├── .agents/
│   ├── rules/physics_rules.md          # Agent safety and delegation rules
│   └── workflows/run_data.md           # Structured /run_data workflow definition
├── README.md                            # Human-facing documentation
├── CLAUDE.md                            # This file (project info and reference)
├── CLAUDE_WORKFLOWS.md                  # Workflow and process instructions
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
  "texture_quantization_method": "fixed_bin_number",
  "use_gpu": false,
  "gpu_device": 1,
  "run_imputation_sensitivity": false,
  "fit_time_varying_cox": true,
  "export_validation_model": false,
  "external_validation_data": "",
  "auxiliary_biomarker_csv": "",
  "use_auxiliary_biomarkers": false,
  "run_trajectory_plots": true
}
```

`dwi_type` must be one of: `"Standard"`, `"dnCNN"`, or `"IVIMnet"`.

`core_method` selects the tumor core delineation algorithm. Options: `"adc_threshold"` (default), `"d_threshold"`, `"df_intersection"`, `"otsu"`, `"gmm"`, `"kmeans"`, `"region_growing"`, `"active_contours"`, `"percentile"`, `"spectral"`, `"fdm"`. The `percentile`, `spectral`, and `fdm` methods use a unified core mask across all parameters (replacing individual D/f/D* thresholds).

### Config Backwards Compatibility (mandatory)

When adding a new config field, you **must** add a corresponding default in `pipeline/utils/parse_config.m` using the existing `isfield` + fallback pattern so that config files without the new field continue to work unchanged. Do **not** modify `config.json` to add a new field without this default in place.

When removing a config field, you **must** ensure that all code referencing the field is updated so that existing config files still containing the removed field do not cause errors (e.g., the field is simply ignored). Do **not** remove a field from `config.example.json` or `pipeline/utils/parse_config.m` without verifying that every consumer of that field has been cleaned up.

If a change (addition or removal) truly cannot be made backwards-compatible, you **must** ask the user for explicit permission before proceeding.

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
| `compute_texture_features.m` | First-order, GLCM, GLRLM, shape, and uniformity texture feature extraction from parameter maps (24 features) |
| `compute_registration_quality.m` | Registration quality metrics: Jacobian determinant, NCC, mutual information |
| `detect_motion_artifacts.m` | DWI volume quality assessment: CV, NMI, signal dropout detection |
| `imputation_sensitivity.m` | Imputation sensitivity analysis (KNN vs LOCF vs Mean vs Linear Interp) |
| `fit_time_varying_cox.m` | Stratified and extended Cox models for PH violations |
| `decision_curve_analysis.m` | Decision curve analysis for clinical utility |
| `compute_nri.m` | Net reclassification improvement (NRI, cNRI, IDI) |
| `prepare_external_validation.m` | Export trained model for external validation |
| `apply_external_validation.m` | Apply saved model to external dataset |
| `load_auxiliary_biomarkers.m` | Load non-DWI biomarker data from CSV |
| `suppress_core_warnings.m` | Suppress expected warnings from `extract_tumor_core` during batch operations |
| `dispatch_load_and_sanity.m` | Extracted dispatch logic for load and sanity check pipeline steps |
| `dispatch_pipeline_steps.m` | Extracted dispatch logic for metrics, visualization, and comparison pipeline steps |
| `prepare_pipeline_session.m` | Pipeline session initialization with try-catch error handling |

### Octave Compatibility (`pipeline/.octave_compat/`)

Contains 21 shim files for GNU Octave compatibility, including:

- `@table/` class implementation (`table.m`, `subsasgn.m`, `subsref.m`, `display.m`)
- `+matlab/+unittest/` namespace shims (`TestSuite.m`, `TestCase.m`, `TestRunner.m`)
- `+matlab/+unittest/+fixtures/` shim (`PathFixture.m`)
- `+matlab/+unittest/+plugins/` shim (`CodeCoveragePlugin.m`)
- Standard function replacements: `cvpartition.m`, `nanmean.m`, `nanstd.m`, `categorical.m`, `niftiread.m`, `niftiwrite.m`, `niftiinfo.m`, `fitglme.m`, `contains.m`, `sgtitle.m`, `yline.m`, `spectralcluster.m`

### Analysis Scripts (`analysis/`)

Python scripts for post-hoc analysis of pipeline outputs, organized into subpackages. The suite includes vision-based graph analysis (via Google Gemini and/or Anthropic Claude APIs), direct log/CSV parsing, cross-DWI comparison, and automated HTML/PDF report generation.

**Requirements:** Python 3.12+, `anthropic`, `google-genai`, `pydantic`, `tqdm`, `weasyprint` (install via `pip install -r analysis/requirements.txt`). Vision analysis requires `GEMINI_API_KEY` and/or `ANTHROPIC_API_KEY` environment variables depending on the selected provider (`--provider gemini|claude|both`); PDF generation requires `weasyprint`; all other scripts work without these optional dependencies. All scripts display `tqdm` progress bars during processing.

**Configuration:** All analysis scripts share a centralised config loaded by `shared.load_analysis_config()`. Defaults are built into `shared.py`; overrides come from `analysis_config.json` at the repository root (copy from `analysis_config.example.json`) and optionally from the MATLAB `config.json` (for `dwi_type`). The `run_analysis.py` orchestrator also accepts `--provider`, `--gemini-model`, `--claude-model`, `--concurrency`, `--config`, `--skip-checks`, and `--interactive` CLI flags. By default, the orchestrator verifies that all `requirements.txt` packages are installed and runs the full pytest suite before starting the analysis pipeline; `--skip-checks` bypasses these pre-flight checks.

| File | Purpose |
|---|---|
| `run_analysis.py` | Orchestrator: runs the full analysis workflow with `--folder`, `--skip-vision`, `--report-only`, `--no-pdf`, `--html`, `--skip-checks`, `--interactive`, `--provider` flags; verifies requirements and runs tests before starting |
| `shared.py` | Shared utilities: folder discovery, DWI type parsing, p-value/correlation regex extraction, config loading |
| `parsers/batch_graph_analysis.py` | Async batch processing of all graph images via Google Gemini and/or Anthropic Claude vision APIs; supports `--provider gemini\|claude\|both` for dual-provider comparison; outputs structured CSV with axes, trends, inflection points, statistical tests, outliers, reference lines, clinical relevance, and metadata |
| `parsers/parse_log_metrics.py` | Direct parsing of MATLAB log files: Wilcoxon p-values, AUC, hazard ratios, GLME interaction terms, sanity check convergence/alignment |
| `parsers/parse_csv_results.py` | Direct parsing of pipeline CSV exports (Significant_LF_Metrics.csv, FDR_Sig_Global.csv) with cross-DWI comparison |
| `parsers/parse_mat_metrics.py` | Parses MATLAB `.mat` output files (core comparison, dosimetry, summary metrics) into JSON for downstream analysis |
| `report/generate_report.py` | HTML+PDF report orchestrator: data loading, section assembly, CLI entry point for `analysis_report.html` and `analysis_report.pdf` |
| `report/report_formatters.py` | Formatting utilities for the HTML report (escaping, badges, nav bar, stat cards, forest plot cells, effect size helpers, table/figure numbering, figure captions, citation system, manuscript sentence helpers) |
| `report/report_constants.py` | Large constants extracted from report_formatters (CSS stylesheet, JavaScript, publication references with BibTeX, HTML template) |
| `report/generate_interactive_report.py` | Interactive HTML report with client-side filtering, Chart.js visualisations, patient drill-down, sortable tables, and DWI/core-method comparison |
| `report/interactive_constants.py` | CSS and JavaScript constants for the interactive report (sidebar, tabs, chart rendering, filter logic) |
| `report/sections/` | Section builder package for the HTML report (36 submodules). Core modules: `_helpers.py`, `metadata.py`, `enrollment.py`, `publication.py`, `discussion.py`. Main results: `main_results.py`, `main_results_summary.py`, `main_results_hypothesis.py`, `main_results_trends.py`. Manuscript: `manuscript.py`, `manuscript_findings.py`, `manuscript_performance.py`, `manuscript_results.py`. Data: `data_overview.py`, `data_quality.py`, `data_supplemental.py`, `supplemental.py`, `gallery.py`. Analysis: `graph_overview.py`, `analysis_graphs.py`, `cross_dwi.py`, `analysis_cross_dwi.py`, `analysis_features.py`, `correlations.py`. Statistics: `statistical_reporting.py`, `effect_sizes.py`, `statistics_effects.py`, `statistics_diagnostics.py`, `statistics_robustness.py`, `model_diagnostics.py`, `model_robustness.py`, `power_analysis.py`, `forest_plot.py`. Legacy shims (`analysis_sections.py`, `statistics.py`, `data_sections.py`) re-export for backward compatibility. |
| `cross_reference/cross_reference_dwi.py` | Full cross-DWI comparison (Standard vs dnCNN vs IVIMnet) of trends, inflection points, and summaries |
| `cross_reference/cross_reference_summary.py` | Concise cross-DWI summary focusing on priority clinical graphs and trend agreement/disagreement |
| `cross_reference/statistical_relevance.py` | Extracts p-values and correlation coefficients; reports significant findings, notable correlations, and cross-DWI significance |
| `cross_reference/statistical_by_graph_type.py` | Filters statistical findings by graph type (scatter, box, line, heatmap, bar, histogram, parameter_map) |
| `cross_reference/cross_dwi_agreement.py` | Bland-Altman, Lin's CCC, and ICC agreement analysis between DWI types |
| `report/sections/forest_plot.py` | Forest plot section builder: HR extraction, matplotlib forest plot, report integration |

**Python Test Suite (pytest):** 45 test files with 1795 tests in `analysis/tests/`. Run with `cd analysis/tests && python -m pytest -v`.

| File | What it covers |
|---|---|
| `conftest.py` | Shared fixtures: synthetic saved_files directories, graph CSVs, log files, pipeline CSV exports |
| `test_shared.py` | DWI type parsing, p-value/correlation extraction, CSV loading, folder resolution |
| `test_parse_log_metrics.py` | GLME, ROC/AUC, survival, baseline, sanity check regex parsing; integration with log files |
| `test_parse_csv_results.py` | CSV reading, cross-DWI significance consistency analysis |
| `test_batch_graph_analysis.py` | Image collection, base64 encoding, MIME types, Pydantic schemas (Axis, Trend, InflectionPoint, StatisticalTest, Outlier, ReferenceLine, GraphAnalysis), CSV flattening, Claude rate-limit detection, provider selection, dual-provider comparison, shared response parsing |
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
| `test_parse_imputation_and_tv_cox.py` | Imputation sensitivity AUC parsing and time-varying Cox HR extraction tests |
| `test_report_sections_robustness.py` | Model robustness report section: imputation comparison table, time-varying Cox summary |
| `test_evaluator_finding.py` | Improvement loop evaluator: Finding Pydantic model, audit scoring, exit condition logic, diminishing returns detection, loop config loading |
| `test_git_utils.py` | Improvement loop git utilities: branch operations, test runners, commit helpers |
| `test_loop_tracker.py` | Improvement loop tracker: iteration logging, context generation, score drift detection |
| `test_orchestrator.py` | Improvement loop orchestrator: audit/fix/evaluate cycle, self-healing protocol, JSON escape sanitization |
| `test_new_report_sections.py` | New report section builders: data overview, data quality, manuscript sub-sections, analysis features, statistics sub-sections |
| `test_auditor_agent.py` | Improvement loop auditor agent: RAG-enhanced code audit, source file collection, finding parsing |
| `test_implementer_agent.py` | Improvement loop implementer agent: branch creation, code fix generation, syntax checking, dry-run mode |
| `test_reviewer_agent.py` | Improvement loop reviewer agent: diff generation, quality gate verdicts, risk flag enforcement |
| `test_orchestrator_v2.py` | Improvement loop v2 orchestrator: four-agent pipeline (audit → implement → review → merge), RAG integration |
| `test_chunker.py` | RAG semantic code chunker: MATLAB/Python/Markdown/JSON splitting, oversized chunk handling, metadata extraction |
| `test_retriever.py` | RAG retriever: ChromaDB query interface, agent-specific context builders, deduplication, relevance filtering |
| `test_rag_integration.py` | RAG integration: end-to-end chunk → index → retrieve pipeline, incremental updates, history indexing |
| `test_indexer.py` | RAG indexer: ChromaDB build/update, incremental indexing, duplicate chunk ID handling, improvement history |
For the full list of 105 MATLAB test files and 45 Python test files with descriptions, see [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md#key-matlab-test-files).

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

For detailed tables of all core modules (18 files), utility modules (72 files), Octave compatibility shims (21 files), analysis scripts (37 report section files), and Python test files (45 files), see [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md).
