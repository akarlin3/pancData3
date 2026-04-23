# CLAUDE.md — AI Assistant Guide for pancData3

Essential context for AI assistants working in this repo. Detailed tables (modules, utilities, tests, analysis scripts, dependencies, diary folder structure) live in [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md). Pipeline runs, tests, git workflow, and AveryLoop in [CLAUDE_WORKFLOWS.md](CLAUDE_WORKFLOWS.md).

---

## Project Overview

**pancData3** is a MATLAB pipeline for pancreatic DWI (Diffusion-Weighted Imaging) research at Memorial Sloan Kettering. It fits IVIM/ADC models, applies DL denoising (DnCNN, IVIMnet), correlates with RT dose maps, and runs survival / competing-risks / treatment-response analysis.

- **Language:** MATLAB (R2021a+), with Python 3.12+ for post-hoc analysis
- **Version:** 2.3.2 · **License:** AGPL-3.0 (Copyright 2026 Avery Karlin)
- **Platforms:** Windows 10/11, macOS 13+, Linux (Ubuntu 22.04+) — CI-tested on all three

---

## Multi-Agent Collaboration

| Agent | Role | Scope |
|---|---|---|
| **Claude Code** (interactive) | Feature implementation, debugging, review | Local, full repo access |
| **Antigravity** (local) | Core physics, MRI calibration, specialized scripts | Local, with patient data |
| **AveryLoop** (external) | Autonomous audit → implement → review → merge (RAG) | External [averyloop](https://github.com/akarlin3/averyLoop), configured via `project_config.yaml` |

### Critical Safety Rules

> **NEVER send patient data, sensitive CSVs, or PHI to any cloud agent or external service.** Only code logic and structure may leave the machine.

- Do **not** modify files in `pipeline/dependencies/` (third-party, read-only).

---

## Repository Structure

```
pancData3/
├── config.json / config.example.json            # Pipeline config (active / template)
├── analysis_config.json / .example.json         # Analysis config (active / template)
├── Dockerfile · docker-compose.yml · docker/    # Docker build (see docs/DOCKER.md)
├── pipeline/                                    # MATLAB pipeline
│   ├── run_dwi_pipeline.m                       #   Master orchestrator (entry point)
│   ├── execute_all_workflows.m                  #   Runs all 3 DWI types sequentially
│   ├── patient_data_check.m                     #   Pre-pipeline data integrity scanner
│   ├── generate_patient_exclusion_report.m      #   Patient exclusion report
│   ├── core/                                    #   Primary modules (19 files)
│   ├── utils/                                   #   Helper utilities (85 files)
│   ├── .octave_compat/                          #   Octave shims (24 files)
│   ├── tests/                                   #   Test suite (136 files)
│   └── dependencies/                            #   Third-party — DO NOT MODIFY
├── analysis/                                    # Python post-hoc analysis
│   ├── run_analysis.py · shared.py              #   Orchestrator + shared utilities
│   ├── parsers/ · cross_reference/              #   Log/CSV/MAT/vision parsers; cross-DWI
│   ├── report/sections/                         #   HTML+PDF report builders (51 files)
│   └── tests/                                   #   Python pytest suite (52 files)
├── project_config.yaml / .example.yaml          # AveryLoop project config
├── averyloop_config.example.json                # AveryLoop runtime config template
├── .agents/rules · .agents/workflows            # Agent safety rules & /run_data workflow
└── CLAUDE.md · CLAUDE_WORKFLOWS.md · CLAUDE_REFERENCE.md
```

See [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md) for per-file tables.

---

## Configuration

Copy `config.example.json` to `config.json` and fill in site-specific paths. The example file is the source of truth for field names and defaults — consult it rather than duplicating a list here.

Key enumerated fields:
- `dwi_type`: `"Standard"` | `"dnCNN"` | `"IVIMnet"`
- `core_method`: `"adc_threshold"` (default) | `"d_threshold"` | `"df_intersection"` | `"otsu"` | `"gmm"` | `"kmeans"` | `"region_growing"` | `"active_contours"` | `"percentile"` | `"spectral"` | `"fdm"`. The `percentile`, `spectral`, and `fdm` methods use a unified core mask across all parameters (replacing individual D/f/D* thresholds).

### Config Backwards Compatibility (mandatory)

- **Adding a field:** add a corresponding default in `pipeline/utils/parse_config.m` using the existing `isfield` + fallback pattern. Config files without the new field must continue to work.
- **Removing a field:** update every consumer so existing configs that still contain the field are silently ignored. Do not remove from `config.example.json` or `parse_config.m` until consumers are clean.
- If a change truly cannot be backwards-compatible, **ask the user first**.

---

## Code Conventions

### Architecture
- **Orchestrator pattern:** `pipeline/run_dwi_pipeline.m` sequences steps; modules are independently callable.
- Explicit parameter passing between modules — no shared global workspace state.
- Outputs written to a timestamped folder; the path is passed as an argument, never hardcoded.

### Parallelization
- `parfor` processes patient cohorts in parallel; `parsave_dir_cache.m` enables safe `save` inside `parfor`.
- Checkpointing allows recovery from mid-cohort interruptions.
- Parallel pool capped at 2 workers in `pipeline/execute_all_workflows.m`.

### Cross-Platform Compatibility
- `escape_shell_arg.m` auto-detects `ispc()` for Windows (double-quote) vs Unix (single-quote) shell escaping.
- All file paths use `fullfile()`, `filesep`, `pathsep` — never hardcoded separators.
- Python scripts use `pathlib.Path` and reconfigure stdout to UTF-8 on Windows.
- CI runs the full MATLAB and Python test suites on Linux, macOS, and Windows.

### Security
- `safe_load_mask.m` inspects variable types and rejects non-numeric classes to prevent arbitrary code execution from `.mat` files.
- Use `escape_shell_arg.m` for every path passed to `system()`. Never pass unsanitized user strings to shell commands.

### File Deletion Safety
The pipeline must **never delete a file or directory it did not create**. All deletion sites use provenance verification:

- **`dwi_vectors*.mat` is the user's nameset — off-limits.** The pipeline must never read, write, or delete any file matching `dwi_vectors*.mat` (including variants like `dwi_vectors.mat`, `dwi_vectors_Standard.mat`, `dwi_vectors_ea.mat`, and date-stamped backups). That nameset belongs to the user's curated files. The pipeline's own upstream voxel-extraction checkpoint lives under `pipeline_voxels*.mat` instead. `clear_pipeline_cache.m` has a defensive `startsWith('dwi_vectors')` guard; new code must not add any I/O path targeting these files.
- **Sentinel files:** Pipeline-created directories (`saved_files_*`, `processed_patients/`) contain a `.pipeline_created` sentinel. Verify it exists before `rmdir`. Legacy `processed_patients/` directories that predate the sentinel convention are healed at the start of each `load_dwi_data` run via `pipeline/utils/backfill_checkpoint_sentinel.m`, which writes the sentinel iff the directory contains at least one `patient_NNN_*.mat` file and no foreign entries (a content-based provenance check that refuses to claim unrelated user folders).
- **Cache clearing** (`pipeline/utils/clear_pipeline_cache.m`): deletes only pipeline-owned caches (`pipeline_voxels_*.mat`, `summary_metrics_*.mat`, `adc_vectors.mat`). `dwi_vectors*.mat` is excluded by design.
- **Lock files:** deleted only when orphaned (stale crashed worker) or after successful checkpoint completion.
- **Diary/log files:** deleted only immediately before recreation by the same module.
- **Test cleanup:** only remove artifacts the test itself created — use pre/post directory snapshots or sentinel checks.
- New code that deletes must add a provenance check: sentinel file, or confirmation the file was created in the same function scope.

### Data Leakage Prevention
- `make_grouped_folds.m` — patient-stratified folds prevent intra-patient CV leakage.
- `knn_impute_train_test.m` — infinite distances for rows from the same patient.
- `scale_td_panel.m` — scaling is strictly timepoint-specific.
- `load_dl_provenance.m` — guards against overlap between DL training and analysis patient sets.

### Logging & Diaries
Emoji prefixes: 🚀 start · ⚙️ step · ✅ success · ❌ fatal · ⚠️ non-fatal · 📁 file · 💡 note.

MATLAB allows only one active diary at a time, so the architecture uses a restart pattern:
1. `execute_all_workflows.m` opens a master diary (`execute_all_workflows.log`).
2. `run_dwi_pipeline.m` opens a per-DWI-type diary (`pipeline_log_{type}.txt`).
3. Each core module (sanity_checks, visualize_results, metrics_*) opens its own diary, overriding the orchestrator.
4. After each module returns, the orchestrator **restarts** its diary.

**Important for tests:** any test that exercises a core module must call `diary off;` in `TestMethodTeardown` before `rmdir` on temp dirs — on Windows the diary file stays locked. Full folder layout: [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md#diary--console-logging--output-folder-structure).

### Error Handling
- Validate MATLAB toolbox licenses at startup (Statistics + Image Processing).
- Halt with a clear message on unrecoverable errors; graceful halting over silent bad data.
- Non-fatal module failures (`compare_core_methods`, `metrics_longitudinal`, dosimetry, stats, survival) log a ⚠️ and continue.

---

## Dependencies

- **MATLAB toolboxes:** Statistics and Machine Learning · Image Processing
- **External tools:** `dcm2niix` (MRIcroGL) — on PATH or set `dcm2nii_call` in `config.json`
- **Python (analysis):** 3.12+; `pip install -r analysis/requirements.txt` (`anthropic`, `google-genai`, `pydantic`, `tqdm`, `weasyprint`). Vision analysis needs `GEMINI_API_KEY` and/or `ANTHROPIC_API_KEY`; PDF needs `weasyprint`.
- `pipeline/dependencies/` — third-party, read-only. File listing in [CLAUDE_REFERENCE.md](CLAUDE_REFERENCE.md#pipelinedependencies-contents-do-not-modify).
