# CLAUDE_WORKFLOWS.md — Workflow & Process Instructions

This file contains all actionable workflow instructions, rules, and processes for AI assistants working in this repository.

For project info, module tables, config reference, and code conventions, see [CLAUDE.md](CLAUDE.md).
For detailed module/test/utility reference tables, see [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md).

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
| `test_imputation_sensitivity.m` | Imputation sensitivity: KNN vs LOCF vs Mean vs Linear Interp comparison |
| `test_time_varying_cox.m` | Time-varying Cox: stratified models, extended Cox, PH violation handling |
| `test_decision_curve_analysis.m` | Decision curve analysis: net benefit, treat-all/none comparison, clinical utility |
| `test_compute_nri.m` | NRI computation: reclassification tables, continuous NRI, IDI |
| `test_compute_median_followup.m` | Median follow-up: all-censored, mixed cohort, competing risk, valid_pts subsetting, reverse KM |
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

---

# Pipeline Improvement Loop

> The improvement loop is now an **external package**: [`code-improvement-loop`](https://github.com/akarlin3/improvementLoop).
> The package name (`improvement_loop`) is unchanged — it is now installed from an external source via `pip`.

## Installation

```bash
pip install -r analysis/requirements.txt
```

This installs `code-improvement-loop` from GitHub along with all other dependencies.

## Configuration

Copy `project_config.example.yaml` to `project_config.yaml` at the repo root and edit as needed:

```bash
cp project_config.example.yaml project_config.yaml
```

The config file defines:
- Project metadata and git default branch
- Source directories and test commands
- Key files for audit context
- System prompts for audit, review, and judge agents
- Safety critical flags (`LEAKAGE_RISK`, `PHI_RISK`)
- RAG collection name

Runtime tuning (API models, token limits, exit strategy thresholds, RAG settings) is still configured via `improvement_loop_config.json` at the repo root.

## Running

```bash
# Full run
python -m improvement_loop.orchestrator_v2 [--max-iterations N] [--dry-run] [--single-iteration]

# Legacy single-pass fallback
python -m improvement_loop.orchestrator_v1 [--max-iterations N] [--dry-run] [--single-iteration]
```

**v2 pipeline flow:** Each iteration runs four agent phases in sequence:
1. **Audit** — `auditor.audit()` scans the codebase and returns `Finding` objects
2. **Implement** — `implementer.implement()` creates a branch, generates a fix, commits
3. **Review** — `reviewer.review()` diffs original vs new content, returns approve/request_changes/reject
4. **Test & Merge** — runs test suite, merges approved branches, post-merge sanity check

**Requirement:** `ANTHROPIC_API_KEY` must be set.

## RAG Index Management

The RAG index is automatically built/updated on the first `orchestrator_v2` run when `rag_enabled` is `True` (the default). Manual management:

```bash
# View index statistics
python -m improvement_loop.rag.indexer --stats

# Force full rebuild (drops and recreates)
python -m improvement_loop.rag.indexer --force-rebuild --stats
```

RAG can be disabled by setting `"rag_enabled": false` in `improvement_loop_config.json`.

## Exit Condition

The loop exits when `should_continue_loop()` returns False, which requires ALL of:
1. Every finding in the current audit has importance STRICTLY LESS THAN 2 out of 10
2. The evaluator's coverage score is at least 6/10 (audit was thorough enough to trust)
3. The evaluator raised no flags (other than EVALUATION_FAILED)

Additionally, you may ONLY stop when:
4. All tests pass on the dev branch
5. You have explicitly listed every finding from the current audit with its score
6. You have written a one-sentence justification for why each finding is below 2/10

If ANY finding scores 2/10 or higher, you MUST continue the loop.

## Completion
After the loop exits, print the full history:
```python
from improvement_loop.loop_tracker import print_full_summary
print_full_summary()
```

Report:
- Total iterations run
- All improvements made, with their scores
- Evaluator score trend across iterations
- Final state assessment
