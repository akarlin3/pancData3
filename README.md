# pancData3

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-blue?logo=mathworks)](https://www.mathworks.com/products/matlab.html)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-2.1.0--beta.1-blue)](#citation)
[![Tests](https://img.shields.io/badge/tests-120%20MATLAB%20%2B%2036%20Python%20files-brightgreen)](#running-tests)

**A MATLAB-based analysis pipeline for pancreatic DWI (Diffusion-Weighted Imaging) research.**

Developed at [Memorial Sloan Kettering Cancer Center](https://www.mskcc.org/), this pipeline processes MRI data to fit IVIM and ADC diffusion models, apply deep learning denoising, correlate findings with radiotherapy dose maps, and perform survival analysis for treatment response prediction.

**Current version:** 2.1.0-beta.1 — see [CHANGELOG.md](CHANGELOG.md) for details.

---

## Table of Contents

- [Supported Platforms](#supported-platforms)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
  - [Docker](#docker)
- [Configuration](#configuration)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Running Tests](#running-tests)
- [Post-Hoc Analysis Scripts](#post-hoc-analysis-scripts)
- [Improvement Loop](#improvement-loop)
- [Repository Structure](#repository-structure)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

---

## Features

- **IVIM & ADC Model Fitting** -- Intravoxel Incoherent Motion (segmented, Bayesian) and mono-exponential ADC fitting
- **Deep Learning Denoising** -- DnCNN and IVIMnet integration for spatial denoising of DWI data
- **Radiotherapy Dose Correlation** -- Dose-volume histogram analysis with diffusion-defined GTV subvolumes
- **Survival Analysis** -- Competing risks modeling (Cause-Specific Hazards), IPCW weighting, and Kaplan-Meier estimation
- **Predictive Modeling** -- Cross-validated feature selection with strict data leakage prevention
- **Robust Checkpointing** -- Parallel-safe `parfor` processing with recovery from interruptions
- **Cache Management** -- `clear_cache` option to delete all cached pipeline files before re-execution

---

## Supported Platforms

| Platform | MATLAB Pipeline | Python Analysis | CI Tested |
|---|---|---|---|
| **Linux** (Ubuntu 22.04+) | Yes | Yes | Yes |
| **macOS** (13 Ventura+) | Yes | Yes | Yes |
| **Windows** (10/11) | Yes | Yes | Yes |

The codebase uses platform-aware code paths throughout:

- **Shell escaping** -- `escape_shell_arg.m` detects `ispc()` and applies Windows (double-quote) or Unix (single-quote) escaping for all `system()` calls
- **Path handling** -- All paths use `fullfile()`, `filesep`, and `pathsep` (MATLAB) or `pathlib.Path` (Python) instead of hardcoded separators
- **Console encoding** -- Python analysis scripts reconfigure `stdout` to UTF-8 on Windows to handle emoji log prefixes
- **Mock test scripts** -- Tests generate `.bat` files on Windows and shell scripts on Unix

---

## Requirements

### MATLAB Toolboxes

| Toolbox | Required |
|---|---|
| MATLAB R2021a+ | Yes |
| Statistics and Machine Learning Toolbox | Yes |
| Image Processing Toolbox | Yes |

### External Tools

| Tool | Purpose | Installation |
|---|---|---|
| [dcm2niix](https://github.com/rordenlab/dcm2niix) (MRIcroGL) | DICOM-to-NIfTI conversion | [Download](https://github.com/rordenlab/dcm2niix/releases) for Windows, macOS, or Linux. Must be on your system PATH or specified in `config.json` as `dcm2nii_call`. |

---

## Installation

```bash
git clone https://github.com/akarlin3/pancData3.git
cd pancData3
```

Then in MATLAB:

```matlab
addpath('pipeline/core', 'pipeline/utils', 'pipeline/dependencies');
```

### Docker

A Docker image provides fully reproducible execution without installing MATLAB toolboxes, `dcm2niix`, or Python dependencies on the host. For the full reference (MCR versioning, troubleshooting, data safety), see [`docs/DOCKER.md`](docs/DOCKER.md).

#### Prerequisites

| Platform | Requirement |
|---|---|
| **Windows** | [Docker Desktop](https://www.docker.com/products/docker-desktop/) with WSL 2 backend enabled |
| **macOS** | [Docker Desktop](https://www.docker.com/products/docker-desktop/) |
| **Linux** | Docker Engine 20.10+ and Docker Compose v2.0+ |

Verify your installation:

```bash
docker --version          # Docker version 20.10+
docker compose version    # Docker Compose version v2.0+
```

#### Building the image

```bash
docker build -t pancdata3:latest .
```

The first build takes approximately 10--15 minutes (compiles `dcm2niix` from source and installs the MATLAB Runtime). Subsequent builds use Docker layer caching and complete in seconds unless `pipeline/` or `analysis/` code has changed.

#### Config setup

Copy the example config and set **container-internal** paths (not host paths):

```bash
cp config.example.json config.json
```

```json
{
  "dataloc": "/opt/pancData3/data/",
  "dcm2nii_call": "dcm2niix",
  "skip_to_reload": false,
  "dwi_type": "Standard"
}
```

> **Note:** Inside the container, patient data is mounted at `/opt/pancData3/data/` and `dcm2niix` is on the system PATH, so `dcm2nii_call` should be `"dcm2niix"` (not a host path).

#### Environment setup (`.env` file)

Docker Compose automatically reads a `.env` file in the project root. Copy the template and fill in your values:

```bash
cp .env.example .env
```

Then edit `.env` with your host paths and (optionally) API keys:

```dotenv
# Required
DATA_DIR=/path/to/patient_dwi_data

# Optional (defaults shown)
OUTPUT_DIR=./output
CONFIG_FILE=./config.json

# Optional: vision-based graph analysis API keys
GEMINI_API_KEY=your-gemini-key
ANTHROPIC_API_KEY=your-anthropic-key
```

> **Security:** The `.env` file is gitignored and should never be committed, as it may contain API keys. Only `.env.example` (which contains no real values) is tracked.

Alternatively, you can export the variables in your shell instead of using a `.env` file.

#### Running with Docker Compose

Use Compose targets to run the pipeline, analysis, or both:

```bash
# Pipeline + analysis (default)
docker compose up

# Pipeline only
docker compose up pipeline

# Analysis only (requires existing pipeline output in OUTPUT_DIR)
docker compose up analysis
```

#### Running with `docker run`

```bash
# Full pipeline (all DWI types)
docker run --rm \
  -v /path/to/patient_data:/opt/pancData3/data:ro \
  -v /path/to/output:/opt/pancData3/output \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  pancdata3:latest pipeline

# Analysis only
docker run --rm \
  -v /path/to/output:/opt/pancData3/output \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  pancdata3:latest analysis

# Both pipeline and analysis
docker run --rm \
  -v /path/to/patient_data:/opt/pancData3/data:ro \
  -v /path/to/output:/opt/pancData3/output \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  pancdata3:latest all
```

#### Dry run

Validate your setup (config, volume mounts, MATLAB Runtime, Python) without executing the pipeline:

```bash
docker run --rm \
  -v /path/to/patient_data:/opt/pancData3/data:ro \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  pancdata3:latest --dry-run pipeline
```

#### Optional: Vision analysis API keys

To enable vision-based graph analysis, pass API key(s) for your chosen provider:

```bash
# Gemini (default provider)
docker run --rm \
  -e GEMINI_API_KEY=your_key_here \
  -v /path/to/output:/opt/pancData3/output \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  pancdata3:latest analysis --skip-checks

# Claude
docker run --rm \
  -e ANTHROPIC_API_KEY=your_key_here \
  -v /path/to/output:/opt/pancData3/output \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  pancdata3:latest analysis --skip-checks --provider claude

# Both (comparison mode)
docker run --rm \
  -e GEMINI_API_KEY=your_gemini_key \
  -e ANTHROPIC_API_KEY=your_anthropic_key \
  -v /path/to/output:/opt/pancData3/output \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  pancdata3:latest analysis --skip-checks --provider both
```

#### Output

Results are written to the `/opt/pancData3/output` mount point inside the container, which maps to your host `OUTPUT_DIR` (or the `-v .../output` path). The pipeline creates a timestamped `saved_files_YYYYMMDD_HHMMSS/` folder containing all logs, figures, and MAT files.

#### Rebuilding after code changes

After modifying pipeline or analysis code, re-run the build command. Docker layer caching ensures only changed layers are rebuilt:

```bash
docker build -t pancdata3:latest .
```

---

## Configuration

Copy the example configuration files and update paths for your environment:

```bash
cp config.example.json config.json
cp improvement_loop_config.example.json improvement_loop_config.json
```

Edit `config.json` with your local paths:

```json
{
  "dataloc": "/path/to/pancreas_dwi/",
  "dcm2nii_call": "/path/to/dcm2niix",
  "clinical_data_sheet": "MASTER_pancreas_DWIanalysis.xlsx",
  "clinical_sheet_name": "Clin List_MR",
  "skip_to_reload": false,
  "dwi_type": "Standard"
}
```

| Field | Description |
|---|---|
| `dataloc` | Root directory containing patient DWI data |
| `dcm2nii_call` | Path to `dcm2niix` executable |
| `dwi_type` | One of: `"Standard"`, `"dnCNN"`, `"IVIMnet"` |
| `skip_to_reload` | Skip data loading and use cached results |
| `use_checkpoints` | Enable/disable checkpoint recovery |
| `clear_cache` | Delete all cached pipeline files before execution |
| `cause_of_death_column` | Column name for competing risk classification (default: `"CauseOfDeath"`) |
| `core_method` | Method for defining the tumor core. Options: `"adc_threshold"`, `"d_threshold"`, `"df_intersection"`, `"otsu"`, `"gmm"`, `"kmeans"`, `"region_growing"`, `"active_contours"`, `"percentile"`, `"spectral"`, `"fdm"` (default: `"adc_threshold"`) |
| `core_percentile` | Percentile cutoff for the `"percentile"` core method (default: `25`) |
| `core_n_clusters` | Number of tissue clusters for the `"spectral"` core method (default: `2`) |
| `fdm_parameter` | Diffusion parameter for the `"fdm"` core method: `"adc"` or `"d"` (default: `"adc"`) |
| `fdm_thresh` | Fallback fDM significance threshold in mm²/s when no repeat-scan data is available (default: `0.0004`) |
| `run_compare_cores` | Auto-include `compare_cores` step in the default pipeline (default: `false`) |
| `run_all_core_methods` | Compute sub-volume metrics for all 11 core methods per patient/timepoint (default: `false`) |
| `spectral_min_voxels` | Minimum valid voxels required for the `"spectral"` core method (default: `20`) |
| `store_core_masks` | Store per-method 1D core masks for reuse by `compare_core_methods` (default: `false`) |
| `use_firth_refit` | Refit predictive models with Firth penalized logistic regression after elastic net feature selection to handle perfect separation (default: `true`) |
| `compute_fine_gray` | Compute Fine-Gray subdistribution hazard model alongside cause-specific hazard (default: `true`) |
| `exclude_motion_volumes` | Exclude motion-corrupted DWI volumes before model fitting (default: `false`) |
| `use_texture_features` | Compute GLCM and first-order texture features from parameter maps (default: `false`) |
| `texture_quantization_method` | IBSI quantization for texture features: `"fixed_bin_number"` (default) or `"fixed_bin_width"` |
| `use_gpu` | Offload ADC WLS fitting and DnCNN inference to a CUDA GPU via `gpuArray`. Requires Parallel Computing Toolbox and a CUDA-capable GPU. Falls back to CPU when unavailable (default: `false`) |
| `gpu_device` | 1-based index of the CUDA GPU device to use when `use_gpu` is true (default: `1`) |
| `run_imputation_sensitivity` | Compare KNN imputation against LOCF, Mean, and Linear Interpolation alternatives (default: `false`) |
| `fit_time_varying_cox` | Fit stratified and extended Cox models when proportional hazards assumption is violated (default: `true`) |
| `export_validation_model` | Export trained predictive model for external validation (default: `false`) |
| `external_validation_data` | Path to external validation data folder (default: `""`) |
| `auxiliary_biomarker_csv` | Path to CSV file containing non-DWI auxiliary biomarkers (default: `""`) |
| `use_auxiliary_biomarkers` | Include auxiliary biomarkers in the predictive modeling feature matrix (default: `false`) |

See [`config.example.json`](config.example.json) for all available fields and threshold parameters.

### Improvement Loop Configuration

The automated improvement loop has its own config file. Copy the example and customize:

```bash
cp improvement_loop_config.example.json improvement_loop_config.json
```

Key fields include `exit_strategy` (`"classic"`, `"diminishing_returns"`, or `"both"`), diminishing returns thresholds (`dr_window`, `dr_max_merge_rate`, `dr_max_avg_importance`, `dr_min_file_repeats`, `dr_max_audit_score`), API settings (`anthropic_api_key`, `audit_model`, `fix_model`, `judge_model`), and orchestrator knobs (`max_api_retries`, `retry_base_delay`, `max_file_chars`). All fields have sensible defaults — the file is optional. See [`improvement_loop_config.example.json`](improvement_loop_config.example.json) for the full template.

---

## Usage

### Full Pipeline (all DWI types)

Runs Standard, dnCNN, and IVIMnet sequentially with a parallel pool (max 2 workers). Creates a timestamped output folder (`saved_files_YYYYMMDD_HHMMSS/`) containing all logs, figures, and results:

```matlab
run('pipeline/execute_all_workflows.m')
```

### Single Pipeline Run

```matlab
addpath('pipeline/core', 'pipeline/utils', 'pipeline/dependencies');

% Run all steps
run_dwi_pipeline('config.json');

% Run specific steps only
run_dwi_pipeline('config.json', {'load', 'sanity'});

% Run with a pre-existing output folder
run_dwi_pipeline('config.json', {'load', 'visualize'}, 'path/to/saved_files_folder');
```

The optional 3rd argument specifies a parent output folder. If omitted, a new timestamped folder is created automatically.

### Output Structure

All pipeline output is organized into a single timestamped folder:

```
saved_files_YYYYMMDD_HHMMSS/
├── .pipeline_created                # Provenance sentinel (safe to delete)
├── execute_all_workflows.log        # Master workflow log
├── test_suite_output.log            # Test suite output
├── Standard/                        # Per-DWI-type results
│   ├── pipeline_log_Standard.txt    # Pipeline console log
│   ├── *.txt                        # Per-module logs
│   └── *.fig / *.png                # Figures
├── dnCNN/
└── IVIMnet/
```

> **File deletion safety:** The pipeline never deletes files or directories it did not create. Pipeline-created directories contain a `.pipeline_created` sentinel file; all cleanup code verifies this sentinel before removing a directory. Manually curated files (e.g., `dwi_vectors_ea.mat`) are explicitly protected from cache clearing.

### Pre-Pipeline Data Check

Before running the pipeline, you can verify that your patient data directory is properly structured:

```matlab
addpath('pipeline/core', 'pipeline/utils', 'pipeline/dependencies');
report = enrollment_check('config.json');
```

This scans the file system and reports:

- Missing or inaccessible data directories
- `dcm2niix` availability
- Clinical spreadsheet accessibility
- Per-patient: missing fraction folders, DWI DICOM files, GTV masks, and RT dose folders

Issues are classified as errors (must fix before running), warnings (may affect results), or informational. The returned `report` struct contains a summary table for programmatic use.

---

## Pipeline Steps

The pipeline executes the following steps in order:

| Step | Module | Description |
|---|---|---|
| 1. **Load** | `load_dwi_data` | DICOM conversion, model fitting, DL denoising, checkpointing |
| 2. **Sanity** | `sanity_checks` | Convergence validation, missingness, spatial alignment |
| 3. **Visualize** | `visualize_results` | Parameter maps, distributions, longitudinal trajectories |
| 4. **Metrics** | `compute_summary_metrics` | Summary metric aggregation, then sub-steps below |
| 4a. Baseline | `metrics_baseline` | Summary metrics, outlier cleaning, percent delta |
| 4b. Longitudinal | `metrics_longitudinal` | Longitudinal change analysis |
| 4c. Dosimetry | `metrics_dosimetry` | Dose metrics (D95, V50) within diffusion-defined subvolumes |
| 4d. Comparisons | `metrics_stats_comparisons` | Statistical group comparisons |
| 4e. Predictive | `metrics_stats_predictive` | Feature selection, cross-validated predictive modeling |
| 4f. Survival | `metrics_survival` | Competing risks, Cause-Specific Hazards analysis |
| 5. **Compare Cores** | `compare_core_methods` | Pairwise comparison of all 11 tumor core methods (optional, invoke explicitly) |

### Comparing Core Delineation Methods

The `compare_cores` step runs all 11 tumor core delineation methods on every patient and timepoint, then computes pairwise spatial agreement metrics. This is not part of the default pipeline — invoke it explicitly:

```matlab
addpath('pipeline/core', 'pipeline/utils', 'pipeline/dependencies');

% Run after data has been loaded (requires dwi_vectors and summary_metrics on disk)
run_dwi_pipeline('config.json', {'compare_cores'});
```

**Output files** (saved to the DWI-type subfolder):

| File | Description |
|---|---|
| `compare_core_results_{type}.mat` | Full comparison struct: pairwise Dice and HD95 matrices, per-method volume fractions, fallback flags |
| `core_method_dice_heatmap_{type}.png` | 11x11 heatmap of mean pairwise Dice coefficients |
| `core_method_hd95_heatmap_{type}.png` | 11x11 heatmap of mean pairwise 95th-percentile Hausdorff distance (only when 3D GTV masks are available) |
| `core_method_volume_comparison_{type}.png` | Bar chart of mean core volume (% of GTV) per method with error bars |
| `core_method_fallbacks_{type}.png` | Bar chart showing how often each method fell back to `adc_threshold` (e.g., fDM at baseline, spatial methods without 3D masks) |
| `compare_core_methods_output_{type}.txt` | Console log |

The MAT file contains a `compare_results` struct with fields: `method_names`, `mean_dice_matrix`, `mean_hd95_matrix`, `volume_fractions`, `fallback_flags`, `all_dice` (per-patient Dice matrices), and `all_hd95` (per-patient HD95 matrices).

---

## Running Tests

```matlab
run('pipeline/tests/run_all_tests.m')
```

The test suite includes 120 test files covering:

- **Integration tests** -- End-to-end pipeline validation
- **Unit tests** -- Individual module correctness
- **Leakage tests** -- Data leakage prevention in CV and imputation
- **Security tests** -- Safe file loading, shell argument escaping, file deletion provenance
- **Smoke tests** -- Visualization output verification
- **Benchmarks** -- Performance comparisons for optimized algorithms
- **Static analysis** -- Source code standards and naming conventions

Tests generate a code coverage report for `pipeline/core/` and `pipeline/utils/`.

---

## Post-Hoc Analysis Scripts

The `analysis/` folder contains a comprehensive Python analysis suite for automated extraction, cross-referencing, and reporting of pipeline results.

### Requirements

- Python 3.12+
- `pip install -r analysis/requirements.txt`
- Vision API key(s) — see [API Key Setup](#api-key-setup) below

#### API Key Setup

Vision-based graph analysis requires API key(s) depending on your chosen provider:

| Provider | Environment Variable | Get a Key |
|---|---|---|
| Gemini (default) | `GEMINI_API_KEY` | [Google AI Studio](https://aistudio.google.com/) |
| Claude | `ANTHROPIC_API_KEY` | [Anthropic Console](https://console.anthropic.com/) |
| Both (comparison) | Both of the above | Both links above |

Set the key(s) as environment variables before running:

```bash
# Option 1: Export in your shell
export GEMINI_API_KEY="your-gemini-key-here"
export ANTHROPIC_API_KEY="your-claude-key-here"

# Option 2: Inline with the command
GEMINI_API_KEY="..." python analysis/run_analysis.py

# Option 3: Interactive prompt — if a required key is missing,
# run_analysis.py will prompt you to paste it at startup
python analysis/run_analysis.py --provider claude
# → ANTHROPIC_API_KEY is not set. Paste your key: ▌
```

Provider selection: `--provider gemini` (default), `--provider claude`, or `--provider both` (runs both APIs and writes a comparison CSV noting per-image differences).

You can also set defaults in `analysis/analysis_config.json` or via environment variables:

```bash
export PANCDATA3_VISION_PROVIDER="both"
export PANCDATA3_GEMINI_MODEL="gemini-2.5-flash"
export PANCDATA3_CLAUDE_MODEL="claude-opus-4-6"
```

#### WeasyPrint System Dependencies (PDF report generation)

WeasyPrint requires **Cairo**, **Pango**, and **GDK-PixBuf** installed at the system level. These are *not* installed by `pip` and must be set up separately. If you only need HTML reports (use `--no-pdf`), you can skip this step.

**Ubuntu / Debian:**

```bash
sudo apt-get install -y libcairo2-dev libpango1.0-dev libgdk-pixbuf-2.0-dev libffi-dev
```

**macOS (Homebrew):**

```bash
brew install cairo pango gdk-pixbuf libffi
```

**Windows:**

Install the GTK3 runtime, which bundles Cairo, Pango, and GDK-PixBuf:

1. Download the latest GTK3 installer from <https://github.com/nickvdp/gtk3-windows-installer/releases> or install via MSYS2: `pacman -S mingw-w64-x86_64-gtk3`
2. Add the GTK3 `bin` directory to your system `PATH`
3. Restart your terminal before running `pip install weasyprint`

See the [WeasyPrint installation docs](https://doc.courtbouillon.org/weasyprint/stable/first_steps.html) for the latest platform-specific instructions.

### Usage

```bash
# Full analysis workflow (auto-detects latest output folder)
python analysis/run_analysis.py

# Skip vision API calls (use existing CSV or only direct parsing)
python analysis/run_analysis.py --skip-vision

# Only generate the HTML report from existing data
python analysis/run_analysis.py --report-only

# Also keep the HTML report on disk (default: PDF only)
python analysis/run_analysis.py --html

# Specify a particular output folder
python analysis/run_analysis.py --folder saved_files_20260308_010713

# Skip pre-flight checks (requirements verification and test suite)
python analysis/run_analysis.py --skip-checks

# Generate interactive report with client-side filtering and Chart.js charts
python analysis/run_analysis.py --interactive

# Use Claude API instead of Gemini for vision analysis
python analysis/run_analysis.py --provider claude

# Run both Gemini and Claude, compare results (writes graph_analysis_comparison.csv)
python analysis/run_analysis.py --provider both

# Individual scripts can still be run standalone
python analysis/parsers/batch_graph_analysis.py
python analysis/parsers/parse_log_metrics.py [saved_files_path]
python analysis/parsers/statistical_relevance.py [saved_files_path]
```

| Script | Description |
|---|---|
| `run_analysis.py` | Orchestrator: runs the full analysis workflow with CLI flags |
| `shared.py` | Shared utilities: folder discovery, DWI type parsing, regex extractors |
| **`parsers/`** | **Data extraction subpackage** |
| `parsers/batch_graph_analysis.py` | Sends all pipeline graph images to Google Gemini vision API; extracts axes, trends, inflection points into a structured CSV |
| `parsers/parse_log_metrics.py` | Direct parsing of MATLAB log files for Wilcoxon p-values, AUC, hazard ratios, GLME results |
| `parsers/parse_csv_results.py` | Direct parsing of pipeline CSV exports with cross-DWI significance comparison |
| `parsers/parse_mat_metrics.py` | Parses MATLAB `.mat` output files (core comparison, dosimetry, summary metrics) into JSON |
| `parsers/statistical_relevance.py` | Extracts p-values and correlation coefficients; reports significant findings |
| `parsers/statistical_by_graph_type.py` | Filters statistical findings by graph type (scatter, box, line, etc.) |
| **`cross_reference/`** | **Cross-DWI comparison subpackage** |
| `cross_reference/cross_reference_dwi.py` | Full side-by-side comparison of Standard vs dnCNN vs IVIMnet results |
| `cross_reference/cross_reference_summary.py` | Concise summary of trend agreement/disagreement across DWI types |
| **`report/`** | **Report generation subpackage** |
| `report/generate_report.py` | HTML+PDF report generator combining all data sources into `analysis_report.html` and `analysis_report.pdf` |
| `report/report_formatters.py` | Formatting utilities for the HTML report (escaping, badges, nav bar, stat cards, etc.) |
| `report/report_constants.py` | Large constants (CSS stylesheet, JavaScript, publication references, HTML template) |
| `report/sections/` | Section builder modules for the HTML report (17 submodules: metadata, main_results, statistical_reporting, manuscript, enrollment, supplemental, gallery, graph_overview, cross_dwi, correlations, effect_sizes, model_diagnostics, model_robustness, power_analysis, discussion, publication, _helpers) |

### Report Features

The generated HTML/PDF report includes:

- **Cover page** — print-only title page with run timestamp, DWI types, and graph count
- **Table of Contents** — grouped two-column TOC visible on screen; occupies its own page in the PDF
- **Part breaks** — 6 logical divisions (Overview / Data / Statistics / Outcomes / Discussion / Appendix) force page breaks in the PDF
- **PDF page numbers** — WeasyPrint-native footer: `Page N of M`, suppressed on the cover page
- **Clinical context** — RTOG D95/V50 benchmarks, evidence hierarchy (High/Moderate/Exploratory), Research Use Only disclaimer
- **Improved parsers** — NaN/Inf → JSON `null`, timepoint normalization, parse-failure warnings, IPCW range validation

### Analysis Test Suite

The analysis scripts have a comprehensive Python test suite (1482 tests across 32 files) using pytest:

```bash
cd analysis/tests && python -m pytest -v
```

---

## Improvement Loop

The `improvement_loop/` package provides an automated audit-fix-evaluate cycle that uses the Claude API to iteratively improve the codebase. It audits source files, parses structured findings, applies fixes on isolated git branches, runs tests, and logs each iteration to a persistent JSON log.

### Requirements

- `ANTHROPIC_API_KEY` environment variable set
- Python packages: `anthropic`, `pydantic` (already in `analysis/requirements.txt`)

### Usage

All commands must be run from the repository root (`pancData3/`), since the package uses relative imports:

```bash
# Live run — audits codebase, applies fixes on branches, runs tests
python -m improvement_loop.orchestrator_v1

# Dry run — no API calls, no code changes, validates plumbing
python -m improvement_loop.orchestrator_v1 --dry-run

# Limit to 3 audit/fix cycles
python -m improvement_loop.orchestrator_v1 --max-iterations 3

# View iteration history
python -m improvement_loop.loop_tracker summary

# Print context for the next iteration (what's already been done)
python -m improvement_loop.loop_tracker context
```

### How It Works

Each iteration of the loop:

1. **Audit** — Collects key pipeline and analysis source files, sends them to Claude Sonnet with context from prior iterations, and receives a JSON array of structured findings (dimension, file, description, fix, importance, branch name).
2. **Parse** — Validates each finding against a Pydantic schema (`Finding` model in `evaluator.py`). Malformed findings are skipped with a warning.
3. **Fix** — For each finding, creates a git branch (`improvement/<slug>`), calls Claude to generate the updated file, commits, and runs the Python test suite. Findings are tagged as `implemented` (tests pass) or `pending` (tests fail).
4. **Evaluate** — Sends the raw audit output to a judge model (`evaluator.score_audit`) that scores it on 6 dimensions (specificity, accuracy, coverage, prioritization, domain appropriateness, overall). Flags like `LEAKAGE_RISK` or `PHI_RISK` are surfaced.
5. **Exit check** — The loop stops when no findings have importance >= 2, audit coverage is adequate (>= 6/10), and no critical flags are raised.

All iterations are logged to `improvement_loop_log.json` (gitignored). The log tracks audit scores, findings with unique IDs, branch status, and exit conditions.

### Package Structure

| File | Purpose |
|---|---|
| `orchestrator_v1.py` | Main loop: audit → parse → fix → evaluate → exit check |
| `evaluator.py` | `Finding` schema, Claude-based audit scoring, exit condition logic |
| `loop_tracker.py` | Persistent JSON logging, iteration context generation, CLI interface |
| `git_utils.py` | Branch management, test runners, commit helpers |

---

## Repository Structure

```
pancData3/
├── config.example.json             # Configuration template
├── Dockerfile                      # Multi-stage Docker build
├── docker-compose.yml              # Pipeline + analysis services
├── .dockerignore                   # Docker build exclusions
├── docker/                         # Docker support files
│   └── entrypoint.sh              #   Container entrypoint script
├── docs/                           # Additional documentation
│   └── DOCKER.md                  #   Docker usage guide
├── pipeline/                       # MATLAB pipeline
│   ├── run_dwi_pipeline.m          #   Main orchestrator entry point
│   ├── execute_all_workflows.m     #   Sequential multi-type runner
│   ├── enrollment_check.m        #   Pre-pipeline data validation
│   ├── core/                       #   Pipeline modules (18 files)
│   │   ├── load_dwi_data.m         #     Data loading & model fitting
│   │   ├── sanity_checks.m         #     Data validation
│   │   ├── visualize_results.m     #     Visualization generation
│   │   ├── process_single_scan.m   #     Per-scan DICOM/model processing
│   │   ├── metrics_baseline.m      #     Baseline metric computation
│   │   ├── metrics_survival.m      #     Survival analysis
│   │   └── ...
│   ├── utils/                      #   Helper utilities (71 files)
│   │   ├── parse_config.m          #     Configuration parser
│   │   ├── safe_load_mask.m        #     Secure .mat loading
│   │   ├── escape_shell_arg.m      #     Shell argument escaping
│   │   ├── init_scan_structs.m     #     Scan data structure initialization
│   │   └── ...
│   ├── tests/                      #   Test suite (120 test files)
│   │   ├── run_all_tests.m         #     Master test runner
│   │   ├── benchmarks/             #     Performance benchmarks (7 files)
│   │   └── diagnostics/            #     Diagnostic spot-checks (5 files)
│   ├── dependencies/               #   Third-party scripts (read-only)
│   └── .octave_compat/             #   GNU Octave compatibility shims (21 files)
├── analysis/                       # Python post-hoc analysis suite
│   ├── run_analysis.py             #   Orchestrator (full workflow runner)
│   ├── shared.py                   #   Shared utilities
│   ├── parsers/                    #   Data extraction subpackage
│   │   ├── batch_graph_analysis.py #     Vision API batch graph extraction
│   │   ├── parse_log_metrics.py    #     Direct MATLAB log parsing
│   │   ├── parse_csv_results.py    #     Direct CSV export parsing
│   │   ├── parse_mat_metrics.py    #     MATLAB .mat file parser → JSON
│   │   ├── statistical_relevance.py #    Statistical significance extraction
│   │   └── statistical_by_graph_type.py # Stats filtered by graph type
│   ├── cross_reference/            #   Cross-DWI comparison subpackage
│   │   ├── cross_reference_dwi.py  #     Full cross-DWI type comparison
│   │   └── cross_reference_summary.py #  Concise cross-DWI summary
│   ├── report/                     #   Report generation subpackage
│   │   ├── generate_report.py      #     HTML+PDF report generator
│   │   ├── report_formatters.py    #     Formatting utilities
│   │   ├── report_constants.py     #     CSS, JS, references, HTML template
│   │   ├── generate_interactive_report.py  # Interactive HTML report with filtering
│   │   ├── interactive_constants.py #     CSS/JS for interactive report
│   │   └── sections/              #     Section builder modules
│   └── tests/                      #   Python test suite (37 test files, 1576 tests)
├── improvement_loop/               # Automated audit-fix-evaluate loop
│   ├── orchestrator_v1.py          #   Main loop orchestrator
│   ├── evaluator.py                #   Finding schema & audit scoring
│   ├── loop_tracker.py             #   Persistent JSON logging & CLI
│   ├── loop_config.py              #   Centralised config (LoopConfig dataclass)
│   └── git_utils.py                #   Branch management & test runners
├── improvement_loop_config.example.json  # Loop config template
└── .agents/                        # AI agent configuration
    ├── rules/                      #   Agent safety rules
    └── workflows/                  #   Structured workflows
```

---

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

---

## Citation

If you use this software in your research, please cite it:

```bibtex
@software{karlin2026pancdata3,
  author    = {Karlin, Avery},
  title     = {pancData3: Pancreatic DWI Analysis Pipeline},
  year      = {2026},
  version   = {2.1.0-alpha.1},
  url       = {https://github.com/akarlin3/pancData3},
  license   = {AGPL-3.0}
}
```

See [CITATION.cff](CITATION.cff) for a machine-readable citation file.

---

## License

This project is licensed under the GNU Affero General Public License v3.0. See [LICENSE](LICENSE) for details.

Copyright (c) 2026 Avery Karlin
