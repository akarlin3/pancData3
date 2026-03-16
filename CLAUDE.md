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
│   ├── utils/                          # Helper utilities (48 files)
│   ├── .octave_compat/                 # Octave compatibility shims (21 files)
│   ├── tests/                          # Full test suite (92 test files)
│   │   ├── run_all_tests.m             # MATLAB unittest test runner
│   │   ├── benchmarks/                 # Performance benchmarks (7 files)
│   │   └── diagnostics/                # Diagnostic spot-check scripts (5 files)
│   └── dependencies/                   # Third-party scripts — DO NOT MODIFY
├── analysis/                            # Python post-hoc analysis scripts
│   ├── run_analysis.py                 # Orchestrator entry point
│   ├── shared.py                       # Shared utilities and config loading
│   ├── analysis_config.json            # Centralised configuration
│   ├── parsers/                        # Log/CSV/MAT/vision parsing (4 files)
│   ├── cross_reference/                # Cross-DWI comparison scripts (4 files)
│   ├── report/                         # HTML+PDF report generation
│   │   ├── generate_report.py          # Report orchestrator
│   │   ├── report_formatters.py        # Formatting utilities
│   │   ├── report_constants.py         # CSS, JS, references, templates
│   │   └── sections/                   # Section builders (8 files)
│   └── tests/                          # Python test suite — 23 test files, 720 tests (pytest)
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

For the full list of 92 MATLAB test files and 23 Python test files with descriptions, see [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md#key-matlab-test-files).

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
