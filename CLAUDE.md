# CLAUDE.md — AI Assistant Guide for pancData3

This file provides essential context for AI assistants (Claude, Jules, Antigravity, etc.) working in this repository.

---

## Project Overview

**pancData3** is a MATLAB-based analysis pipeline for pancreatic DWI (Diffusion-Weighted Imaging) research at Memorial Sloan Kettering Cancer Center. It processes MRI data to:

- Fit IVIM (Intravoxel Incoherent Motion) and ADC (Apparent Diffusion Coefficient) models
- Apply Deep Learning denoising (DnCNN, IVIMnet)
- Correlate findings with RT (Radiotherapy) dose maps
- Perform survival analysis, competing risks modeling, and treatment response prediction

**Language:** MATLAB (R2021a+)
**License:** MIT (Copyright 2026 Avery Karlin)
**Domain:** Medical Physics / Oncology Research

---

## Multi-Agent Collaboration Model

This repository uses a three-agent architecture:

| Agent | Role | Scope |
|---|---|---|
| **Claude Code** (interactive) | Feature implementation, pipeline enhancements, debugging, code review | Runs locally with full repository access |
| **Antigravity** (local) | Core physics modeling, MRI calibration, specialized scripts | Runs locally with access to patient data |
| **Jules** (cloud/background) | Unit testing, documentation, code styling | Cloud — never receives patient data |

### Critical Safety Rules

> **NEVER send patient data, sensitive CSVs, or PHI to any cloud agent or external service.**
> Only code logic and structure may be shared with Jules or any cloud-based AI.

- Do **not** modify files in the `dependencies/` folder under any circumstances.
- When Jules submits a PR, notify the researcher for joint review before merging.
- All cloud agent work (Jules) is limited to: unit tests, documentation, and PEP 8/style fixes.

---

## Repository Structure

```
pancData3/
├── run_dwi_pipeline.m          # Master orchestrator — main entry point
├── execute_all_workflows.m     # Runs all 3 DWI types sequentially
├── config.json                 # Active configuration (not committed)
├── config.example.json         # Configuration template (committed)
├── core/                       # Primary pipeline modules (17 files)
├── utils/                      # Helper utilities (17 files + octave_compat/)
│   └── octave_compat/          # Octave compatibility shims (20 files)
├── tests/                      # Full test suite (63 test files)
│   ├── run_all_tests.m         # MATLAB unittest test runner
│   ├── benchmarks/             # Performance benchmarks (7 files)
│   └── diagnostics/            # Diagnostic spot-check scripts (5 files)
├── dependencies/               # Third-party scripts — DO NOT MODIFY
├── .agents/
│   ├── rules/physics_rules.md  # Agent safety and delegation rules
│   └── workflows/run_data.md   # Structured /run_data workflow definition
├── README.md                   # Human-facing documentation
└── CLAUDE.md                   # This file
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
  "adc_thresh": 0.00115,
  "dwi_type": "Standard"
}
```

`dwi_type` must be one of: `"Standard"`, `"dnCNN"`, or `"IVIMnet"`.

---

## Running the Pipeline

### Full multi-type run (all 3 DWI types)

```matlab
execute_all_workflows
```

This creates a timestamped output folder (`saved_files_YYYYMMDD_HHMMSS/`), cleans up stale parallel jobs, sets up a parallel pool (max 2 workers), runs Standard → dnCNN → IVIMnet sequentially, modifying `config.json` between runs and skipping reload phases after the first. All console output and logs are saved to the timestamped folder.

### Single step / targeted run

```matlab
addpath('core', 'utils', 'dependencies');

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

### Agent workflow (`/run_data`)

The structured `/run_data` workflow in `.agents/workflows/run_data.md` does:
1. Verifies `config.json` exists
2. Runs `tests/test_dwi_pipeline.m` for stability check
3. Executes `matlab -batch "execute_all_workflows"` as a background command
4. Waits for completion and verifies output MAT files and figures

---

## Running Tests

```matlab
run('tests/run_all_tests.m')
```

- Uses MATLAB's built-in `unittest` framework
- Auto-discovers all test files in `tests/` (including subfolders)
- Generates code coverage report for `core/` and `utils/`
- Compatible with MATLAB R2017b+
- Returns non-zero exit code on failure (CI-safe)

### Test categories

| Location | Description | Run by CI? |
|---|---|---|
| `tests/` | Full functional and integration tests | Yes |
| `tests/benchmarks/` | Performance benchmarks | Yes |
| `tests/diagnostics/` | Diagnostic spot-check tests | Yes |

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

---

## Core Modules (`core/`)

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

---

## Utility Modules (`utils/`)

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

### Octave Compatibility (`utils/octave_compat/`)

Contains 20 shim files for GNU Octave compatibility, including:

- `@table/` class implementation (`table.m`, `subsasgn.m`, `subsref.m`, `display.m`)
- `+matlab/+unittest/` namespace shims (`TestSuite.m`, `TestCase.m`, `TestRunner.m`)
- Standard function replacements: `cvpartition.m`, `nanmean.m`, `nanstd.m`, `categorical.m`, `niftiread.m`, `niftiwrite.m`, `niftiinfo.m`, `fitglme.m`, `contains.m`, `sgtitle.m`, `yline.m`

---

## Code Conventions

### Architecture

- **Orchestrator pattern**: `run_dwi_pipeline.m` sequences all steps; modules are independently callable.
- Explicit parameter passing between modules — no shared global workspace state.
- Outputs written to a timestamped folder; path passed as argument, not hardcoded.

### Parallelization

- `parfor` loops process patient cohorts in parallel.
- `parsave_dir_cache.m` enables safe `save` inside `parfor`.
- Checkpointing allows recovery from interruptions mid-cohort.
- Parallel pool capped at 2 workers (`execute_all_workflows.m`).

### Security

- `safe_load_mask.m` inspects variable type before loading; rejects non-numeric classes to prevent arbitrary code execution from `.mat` files.
- `escape_shell_arg.m` must be used for all paths passed to `system()`.
- Never pass unsanitized user strings to shell commands.

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

**Output folder structure:**

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
│   ├── metrics_longitudinal_output_Standard.txt
│   ├── metrics_dosimetry_output.txt
│   ├── metrics_stats_comparisons_output_Standard.txt
│   ├── metrics_stats_predictive_output_Standard.txt
│   └── metrics_survival_output_Standard.txt
├── dnCNN/                           # Same structure as Standard/
└── IVIMnet/                         # Same structure as Standard/
```

**Important for tests:** Any test that exercises a core module must call `diary off;` in its `TestMethodTeardown` before calling `rmdir` on temp directories, because the module's diary file will still be open (locked on Windows).

### Error Handling

- Validate MATLAB toolbox licenses at startup (Statistics + Image Processing).
- Halt with a clear message on any unrecoverable error.
- Non-fatal module failures (metrics_longitudinal, dosimetry, stats, survival) log a `⚠️` warning and continue the pipeline.
- Graceful halting is preferred over silent continuation with bad data.

---

## Dependencies

### Required MATLAB Toolboxes

- Statistics and Machine Learning Toolbox
- Image Processing Toolbox

### External Tools

- `dcm2niix` (MRIcroGL) — DICOM-to-NIFTI conversion; must be on system path or specified in `config.json` as `dcm2nii_call`.

### `dependencies/` folder (DO NOT MODIFY)

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

See `dependencies/README_DEPENDENCIES.md` for licenses and attribution.

---

## Git Workflow

- Development happens on feature branches; never commit directly to `main`.
- Branch naming follows: `claude/<description>-<session-id>` for Claude branches, `jules-<description>` for Jules branches.
- PRs from Jules must be reviewed before merging.
- `.gitignore` excludes: MATLAB autosave files (`*.asv`, `*.m~`, `*.mex*`), imaging data (`*.nii`, `*.dcm`, `*.h5`), CSVs, and clinical spreadsheets (`*.xlsx`) to prevent accidental PHI commits.
- After pushing a branch, Claude should create a pull request targeting `main` using `gh pr create` **only if `gh` is available** (check with `which gh` first). Include a summary of changes and a test plan in the PR body. If `gh` is not installed, skip PR creation silently (do not notify the user).

---

## What AI Assistants Should and Should Not Do

### Do
- Read and understand existing code before suggesting changes.
- Use `safe_load_mask` and `escape_shell_arg` when handling file I/O.
- Follow the orchestrator pattern — keep pipeline steps modular and independently callable.
- Preserve checkpointing logic in `load_dwi_data.m`; it is critical for large cohort recovery.
- Consult `config.example.json` before adding new configuration fields.

### Do Not
- Modify anything in `dependencies/`.
- Send patient data, clinical CSVs, or any PHI to cloud services.
- Introduce global variables or persistent workspace state between pipeline steps.
- Bypass the temporal leakage safeguards in imputation or cross-validation.
- Use unsanitized strings in `system()` calls.
- Hard-code file paths — all paths must flow through `config.json`.
- Run the test suite (`run_all_tests.m`) or the pipeline (`execute_all_workflows`, `run_dwi_pipeline`) as a verification step. Tests and pipeline runs are initiated by the researcher, not by AI assistants.
