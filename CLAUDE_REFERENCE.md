# CLAUDE_REFERENCE.md — Detailed Module & Test Reference

This file contains detailed reference tables for modules, utilities, tests, and analysis scripts. It is split from `CLAUDE.md` to reduce context window usage — read this file on demand when you need to look up a specific module, test, or utility.

For project overview, safety rules, configuration, conventions, and workflow instructions, see [CLAUDE.md](CLAUDE.md).

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
| `analyze_core_method_outcomes.m` | Per-method univariable Cox PH + KM analysis of dose coverage vs local control |
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
| `perform_statistical_test.m` | Wilcoxon rank-sum testing with NaN-safe group extraction |
| `parsave_dir_cache.m` | Parallel-safe `save` wrapper for `parfor` caching |
| `escape_shell_arg.m` | Cross-platform shell argument escaping (Windows and Unix) |
| `discover_gtv_file.m` | Locates GTV mask file with flexible naming patterns |
| `find_gtv_files.m` | Locates GTVp and GTVn masks for complex tumor anatomy |
| `plot_feature_distribution.m` | Histogram/boxplot with ANOVA p-value annotation |
| `init_scan_structs.m` | Initializes scan data structures for pipeline processing |
| `compute_scan_days_from_dates.m` | Derives scan days from DICOM acquisition dates |
| `format_p_value.m` | Formats p-values for display with appropriate precision |
| `remove_constant_columns.m` | Removes constant/NaN-only columns from feature matrices |
| `parfor_progress.m` | Parallel loop progress reporting |
| `text_progress_bar.m` | Text-based progress bar display |
| `extract_tumor_core.m` | Configurable tumor core delineation (11 methods) |
| `PipelineProgressGUI.m` | Pipeline-aware progress bar wrapper (maps step keys to display names) |
| `ProgressGUI.m` | Professional custom-figure progress bar for MATLAB pipelines |
| `compute_core_failure_rates.m` | Aggregate failure rate breakdown (fallback/empty/insufficient/all-NaN) for all 11 core methods x 3 pipelines |
| `filter_core_methods.m` | Prune core methods by failure rate threshold, manual exclusion, and minimum voxel count |
| `compute_cross_pipeline_dice.m` | Cross-pipeline (Standard/DnCNN/IVIMNet) Dice coefficients for all 11 core methods at Fx1 |
| `compute_dice_hausdorff.m` | Dice coefficient and Hausdorff distance between 3D binary masks |
| `json_set_field.m` | Targeted regex replacement of a field value in raw JSON strings |
| `plot_cross_dwi_subvolume_comparison.m` | Cross-DWI-type ADC subvolume comparison visualization |
| `compute_adc_metrics.m` | ADC summary metrics for a single patient/timepoint/DWI-type (extracted from compute_summary_metrics) |
| `compute_ivim_metrics.m` | IVIM (D/f/D*) summary metrics for a single patient/timepoint/DWI-type (extracted from compute_summary_metrics) |
| `compute_spatial_repeatability.m` | Dice and Hausdorff spatial repeatability between Fx1 repeat sub-volumes |
| `compute_multi_core_metrics.m` | Multi-method (11 core methods) sub-volume metrics per patient/timepoint |
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
| `compute_histogram_laplace.m` | Laplace-smoothed histogram probability distribution |
| `compute_kurt_skew.m` | Kurtosis/skewness computation with minimum sample guard |
| `detect_baseline_outliers.m` | Outcome-blinded IQR outlier detection for baseline metrics |
| `nanmean_safe.m` | Octave-compatible NaN-ignoring mean |
| `nanstd_safe.m` | Octave-compatible NaN-ignoring standard deviation |
| `compute_percent_deltas.m` | Treatment-induced percent/absolute changes from baseline |
| `dispatch_load_and_sanity.m` | Extracted dispatch logic for load and sanity check pipeline steps |
| `dispatch_pipeline_steps.m` | Extracted dispatch logic for metrics, visualization, and comparison pipeline steps |
| `prepare_pipeline_session.m` | Pipeline session initialization with try-catch error handling |
| `get_system_memory.m` | Cross-platform physical memory query (total and available GB); returns `[NaN, NaN]` on unsupported platforms |

---

## Octave Compatibility (`pipeline/.octave_compat/`)

Contains 24 shim files for GNU Octave compatibility, including:

- `@table/` class implementation (`table.m`, `subsasgn.m`, `subsref.m`, `display.m`)
- `+matlab/+unittest/` namespace shims (`TestSuite.m`, `TestCase.m`, `TestRunner.m`)
- `+matlab/+unittest/+fixtures/` shim (`PathFixture.m`)
- `+matlab/+unittest/+plugins/` shim (`CodeCoveragePlugin.m`)
- Standard function replacements: `cvpartition.m`, `nanmean.m`, `nanmedian.m`, `nanstd.m`, `categorical.m`, `iscategorical.m`, `niftiread.m`, `niftiwrite.m`, `niftiinfo.m`, `fitglme.m`, `contains.m`, `sgtitle.m`, `yline.m`, `spectralcluster.m`, `string.m`

---

## Key MATLAB Test Files

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
| `test_compute_core_failure_rates.m` | Core failure rate aggregation tests (struct fields, dimensions, rates range, NaN patient, adc_threshold no fallback) |
| `test_filter_core_methods.m` | Core method pruning tests (no pruning, failure rate, manual exclusion, adc_threshold safety, min voxels, combined) |
| `test_extract_tumor_core_fit_info.m` | extract_tumor_core fit_info output tests (fields, backward compat, success, fallback, all-NaN, empty mask, insufficient voxels) |
| `test_compute_cross_pipeline_dice.m` | Cross-pipeline Dice coefficient tests (struct fields, dimensions, range, identical pipelines, edge cases) |
| `test_compute_dice_hausdorff.m` | Dice/Hausdorff distance computation tests |
| `test_json_set_field.m` | JSON field replacement utility tests |
| `test_cross_dwi_subvolume.m` | Cross-DWI-type subvolume comparison tests |
| `test_landmark_cindex.m` | Landmark concordance index tests |
| `test_source_code_standards.m` | Source code standards enforcement |
| `test_modularity.m` | Module independence and interface tests |
| `test_statistical_methods.m` | Statistical methods validation |
| `test_analyze_core_method_outcomes.m` | Core method outcome analysis tests (struct fields, HR positivity, KM fields, ranking, significance, edge cases) |
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
| `test_imputation_sensitivity.m` | Imputation sensitivity: KNN vs LOCF vs Mean vs Linear Interp comparison |
| `test_time_varying_cox.m` | Time-varying Cox: stratified models, extended Cox, PH violation handling |
| `test_decision_curve_analysis.m` | Decision curve analysis: net benefit, treat-all/none comparison, clinical utility |
| `test_compute_nri.m` | NRI computation: reclassification tables, continuous NRI, IDI |
| `test_prepare_external_validation.m` | External validation: model export, external dataset application, portability |
| `test_load_auxiliary_biomarkers.m` | Auxiliary biomarkers: CSV loading, missing file handling, column validation |
| `test_IVIMmodelfit.m` | IVIM model fitting dependency validation |
| `test_convert_dicom.m` | DICOM-to-NIfTI conversion via dcm2niix |
| `test_corr_filter.m` | Correlation-based feature filtering |
| `test_data_integrity_check.m` | Pre-pipeline data integrity validation |
| `test_discover_gtv_file.m` | GTV mask file discovery with flexible naming |
| `test_dispatch_pipeline_steps.m` | Pipeline step dispatch logic |
| `test_dvh.m` | Dose-volume histogram computation |
| `test_escape_shell_arg.m` | Cross-platform shell argument escaping |
| `test_execute_all_workflows.m` | Multi-DWI-type sequential workflow |
| `test_external_validation_wiring.m` | External validation model export/import wiring |
| `test_find_gtv_files.m` | GTVp/GTVn mask file discovery |
| `test_fit_adc_mono.m` | Mono-exponential ADC fitting |
| `test_fix_verify.m` | Post-fix regression verification |
| `diagnostics/diag_landmark_cindex_mock.m` | Landmark concordance index with mocked data |
| `test_load_dl_provenance.m` | DL training provenance loading and leakage detection |
| `test_mask_loading.m` | Secure mask loading from .mat files |
| `test_normalize_patient_ids.m` | Patient ID normalization for spreadsheet/folder matching |
| `test_normalization_logic.m` | Feature normalization logic validation |
| `test_octave.m` | GNU Octave compatibility smoke tests |
| `test_octave_shims.m` | Octave compatibility shim function validation |
| `test_parfor_progress.m` | Parallel loop progress reporting |
| `test_parsave_dir_cache.m` | Parallel-safe save wrapper for parfor |
| `test_perf_knn.m` | KNN imputation performance benchmarks |
| `test_plot_feature_distributions.m` | Feature distribution visualization (multi-plot) |
| `test_plot_parameter_maps.m` | Parameter map overlay visualization |
| `test_plot_predictive_diagnostics.m` | ROC curve and predictive diagnostic panels |
| `test_plot_scatter_correlations.m` | Correlation scatter plot rendering |
| `test_prepare_pipeline_session.m` | Pipeline session initialization |
| `test_progress_gui.m` | Custom figure progress bar lifecycle |
| `test_run_elastic_net_cv.m` | Elastic net cross-validation fitting |
| `test_run_loocv_risk_scores.m` | Nested LOOCV risk score computation |
| `test_select_dwi_vectors.m` | DWI vector extraction by processing type |
| `test_text_progress_bar.m` | Text-based progress bar display |
| `test_trajectory_visualizations.m` | Waterfall, swimmer, and spider plot rendering |
| `test_visualize_refactor.m` | Visualization module refactoring validation |
| `test_write_sentinel_file.m` | Pipeline sentinel file writing |

---

## Analysis Scripts (`analysis/`)

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
| `report/sections/` | Section builder package for the HTML report (36 submodules). Core: `_helpers.py`, `metadata.py`, `enrollment.py`, `publication.py`, `discussion.py`. Main results: `main_results.py`, `main_results_summary.py`, `main_results_hypothesis.py`, `main_results_trends.py`. Manuscript: `manuscript.py`, `manuscript_findings.py`, `manuscript_performance.py`, `manuscript_results.py`. Data: `data_overview.py`, `data_quality.py`, `data_supplemental.py`, `supplemental.py`, `gallery.py`. Analysis: `graph_overview.py`, `analysis_graphs.py`, `cross_dwi.py`, `analysis_cross_dwi.py`, `analysis_features.py`, `correlations.py`. Statistics: `statistical_reporting.py`, `effect_sizes.py`, `statistics_effects.py`, `statistics_diagnostics.py`, `statistics_robustness.py`, `model_diagnostics.py`, `model_robustness.py`, `power_analysis.py`, `forest_plot.py`. Legacy shims: `analysis_sections.py`, `statistics.py`, `data_sections.py`. |
| `cross_reference/cross_reference_dwi.py` | Full cross-DWI comparison (Standard vs dnCNN vs IVIMnet) of trends, inflection points, and summaries |
| `cross_reference/cross_reference_summary.py` | Concise cross-DWI summary focusing on priority clinical graphs and trend agreement/disagreement |
| `cross_reference/statistical_relevance.py` | Extracts p-values and correlation coefficients; reports significant findings, notable correlations, and cross-DWI significance |
| `cross_reference/statistical_by_graph_type.py` | Filters statistical findings by graph type (scatter, box, line, heatmap, bar, histogram, parameter_map) |
| `cross_reference/cross_dwi_agreement.py` | Bland-Altman, Lin's CCC, and ICC agreement analysis between DWI types |
| `report/sections/cross_pipeline_dice.py` | Cross-pipeline Dice section builder: per-method Dice table across Standard/DnCNN/IVIMNet pipelines |
| `report/sections/failure_rates.py` | Core method failure rate section builder: color-coded table of fallback/empty/insufficient/NaN rates |
| `report/sections/pruning_results.py` | Core method pruning results section builder: pruned/retained method tables with reasons |
| `report/sections/core_method_outcomes.py` | Core method outcome analysis section builder: Cox PH ranking table with HR, CI, p-values |
| `report/sections/forest_plot.py` | Forest plot section builder: HR extraction, matplotlib forest plot, report integration |
| `report/sections/model_robustness.py` | Model robustness section: imputation sensitivity AUC comparison, time-varying Cox HR summary |

---

## Python Test Suite (`analysis/tests/`)

42 test files. Run with `cd analysis/tests && python -m pytest -v`. (Improvement loop tests are in the [code-improvement-loop](https://github.com/akarlin3/improvementLoop) package.)

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
| `test_analysis_config.py` | Config loading, deep merge, layered overrides, MATLAB config integration, caching, Claude/provider config |
| `test_report_sections_helpers.py` | Report helper functions: JSON loading, series normalization, cohort size, AUC finding, trend agreement |
| `test_report_sections_metadata.py` | Report metadata sections: cover page, part breaks, TOC, publication header, data availability |
| `test_report_sections_main_results.py` | Main results sections: executive summary, hypothesis, statistical significance, treatment response |
| `test_report_sections_data_sections.py` | Data sections: cohort overview, patient flow, data completeness, MAT data, appendix, figure gallery |
| `test_report_sections_analysis.py` | Analysis sections: graph overview/issues, stats by type, cross-DWI comparison, correlations, feature overlap |
| `test_report_sections_statistics.py` | Statistics sections: effect sizes, multiple comparisons, model diagnostics, sensitivity, power analysis |
| `test_report_sections_discussion.py` | Discussion sections: methods, limitations, conclusions, reporting checklist, journal guide |
| `test_cross_reference_dwi.py` | Cross-DWI comparison output: graph matching, trend display, stat test formatting, truncation, JSON robustness |
| `test_cross_reference_summary.py` | Cross-DWI summary: trend agreement/disagreement, priority ordering, parameter maps, inflection points |
| `test_statistical_relevance.py` | Statistical findings: p-value extraction, Bonferroni correction, significance markers, correlations, cross-DWI comparison |
| `test_statistical_by_graph_type.py` | Per-graph-type analysis: grouping, trend directions, top-5 non-sig, density/comparison aggregation, summary table |
| `test_xref_unit.py` | Cross-reference correctness: safe_text, p-value/correlation edge cases, trend agreement logic, significance markers, Bonferroni, direction classification, priority ordering |
| `test_integration.py` | End-to-end analysis pipeline integration: runs run_analysis.py on synthetic data, verifies HTML output sections and tables |
| `test_api_connection.py` | Gemini and Claude API connection smoke tests (skipped without API keys) |
| `test_cross_dwi_agreement.py` | Bland-Altman, Lin's CCC, ICC agreement analysis tests |
| `test_cross_pipeline_dice_section.py` | Cross-pipeline Dice report section: synthetic data, empty/missing data graceful handling, table structure |
| `test_failure_rates_section.py` | Core method failure rates report section: synthetic data, color coding, sorting, empty data handling |
| `test_pruning_section.py` | Core method pruning results report section: pruned/retained methods, reasons, edge cases |
| `test_core_method_outcomes_section.py` | Core method outcomes report section: significant/non-significant results, ranking table, HR display |
| `test_forest_plot.py` | HR data extraction and forest plot generation tests |
| `test_parse_imputation_and_tv_cox.py` | Imputation sensitivity AUC parsing and time-varying Cox HR extraction tests |
| `test_report_sections_robustness.py` | Model robustness report section: imputation comparison table, time-varying Cox summary |
| `test_new_report_sections.py` | New report section builders: data overview, data quality, manuscript sub-sections, analysis features, statistics sub-sections |
| `test_implementer_agent.py` | `ImplementResult` dataclass, `implement()` dry-run / file-not-found / branch-exists paths, `_generate_fix()` API call contract |
| `test_orchestrator_v2.py` | `FindingState`/`IterationState` dataclasses, `run_loop()` dry-run and exit-condition behavior, rejected-finding non-merge guarantee, `_print_agent_summary()` output |
| `test_reviewer_agent.py` | `ReviewVerdict` dataclass, `_generate_diff()`, `_parse_review_verdict()` (valid/invalid JSON, fenced, preamble), critical-flag override, parse-failure fallback |

---

## `pipeline/dependencies/` Contents (DO NOT MODIFY)

Contains third-party scripts. Treat as read-only:

| File | Purpose |
|---|---|
| `IVIMmodelfit.m`, `IVIM_seg.m`, `IVIM_bayes.m` | IVIM model fitting implementations |
| `fit_adc_mono.m` | Mono-exponential ADC fitting |
| `apply_dncnn_symmetric.m` | DnCNN deep learning denoising |
| `dvh.m`, `sample_rtdose_on_image.m` | Dose-volume histogram processing |
| `clean_dir_command.m` | Directory cleaning helper |
| `halfSampleMode.m` | Half-sample mode statistical estimator |
| `im2Y.m` | Image-to-luminance conversion |

See `pipeline/dependencies/README_DEPENDENCIES.md` for licenses and attribution.

---

## Diary / Console Logging — Output Folder Structure

```
saved_files_YYYYMMDD_HHMMSS/
├── execute_all_workflows.log        # Top-level workflow log
├── test_suite_output.log            # Full test suite output
├── preflight_tests_output.log       # Pre-flight test output (run_dwi_pipeline)
├── error.log                        # Error/warning log
├── Standard/                        # DWI type subfolder
│   ├── pipeline_log_Standard.txt    # Orchestrator log
│   ├── sanity_checks_output.txt
│   ├── visualize_results_output.txt
│   ├── metrics_baseline_output_Standard.txt
│   ├── compare_core_methods_output_Standard.txt
│   ├── compare_core_results_Standard.mat
│   ├── metrics_longitudinal_output_Standard.txt
│   ├── metrics_dosimetry_output.txt
│   ├── metrics_stats_comparisons_output_Standard.txt
│   ├── metrics_stats_predictive_output_Standard.txt
│   └── metrics_survival_output_Standard.txt
├── dnCNN/                           # Same structure as Standard/
└── IVIMnet/                         # Same structure as Standard/
```
