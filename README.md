# pancData3

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-blue?logo=mathworks)](https://www.mathworks.com/products/matlab.html)
[![Python](https://img.shields.io/badge/Python-3.12%2B-blue?logo=python)](https://www.python.org/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](https://github.com/akarlin3/pancData3/blob/main/LICENSE)
[![Version](https://img.shields.io/badge/version-2.3.2-blue)](#citation)
[![Tests](https://img.shields.io/badge/tests-133%20MATLAB%20%2B%2047%20Python%20files-brightgreen)](#running-tests)

**An open-source MATLAB and Python pipeline for MRI-guided pancreatic cancer radiotherapy research.**

Developed at [Memorial Sloan Kettering Cancer Center](https://www.mskcc.org/), this pipeline processes diffusion-weighted MRI data to fit biophysical models, apply deep learning denoising, quantify tumor–dose relationships, and predict treatment outcomes using survival analysis — all from a single configuration file.

**Current version:** 2.3.2 — see [CHANGELOG.md](CHANGELOG.md) for details.

---

## Why This Exists

Pancreatic cancer has some of the lowest survival rates of any solid tumor. MRI-guided adaptive radiotherapy is emerging as a promising treatment approach, but analyzing diffusion-weighted imaging (DWI) data to predict which patients will respond is a complex, multi-step computational problem. This pipeline automates the entire workflow: from raw DICOM images through diffusion model fitting, deep learning denoising, dose–response quantification, and survival modeling — enabling researchers to go from scan to clinical prediction in a single reproducible run.

---

## Key Technical Highlights

| Area | What It Does |
|------|-------------|
| **Deep Learning** | DnCNN and IVIMnet models for MRI denoising and diffusion parameter estimation |
| **Biophysical Modeling** | IVIM (segmented, Bayesian) and mono-exponential ADC fitting on DWI data |
| **Survival Analysis** | Competing risks (Cause-Specific Hazards), IPCW weighting, Kaplan-Meier estimation |
| **Predictive Modeling** | Cross-validated elastic net with strict train/test data leakage prevention and optional Firth penalized refit |
| **Tumor Segmentation** | 11 core delineation methods with pairwise spatial agreement benchmarking (Dice, HD95) |
| **Dose–Response** | Dose-volume histogram analysis correlating diffusion-defined GTV subvolumes with radiotherapy dose maps |
| **Reproducibility** | Docker support, parallel-safe checkpointing, config-driven execution, 180-file test suite |
| **Cross-Platform** | Runs on Linux, macOS, and Windows with platform-aware shell escaping and path handling |

---

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                         config.json                                 │
│            (data paths, DWI type, model parameters)                 │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────────────┐
│  1. LOAD                                                             │
│  DICOM → NIfTI (dcm2niix) → IVIM/ADC model fitting                  │
│  Optional: DnCNN or IVIMnet deep learning denoising                  │
│  Parallel-safe checkpointing with parfor                             │
└──────────────────────────┬───────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────────────┐
│  2. SANITY CHECKS                                                    │
│  Convergence validation · Missingness audit · Spatial alignment      │
└──────────────────────────┬───────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────────────┐
│  3. VISUALIZE                                                        │
│  Parameter maps · Distributions · Longitudinal trajectories          │
└──────────────────────────┬───────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────────────┐
│  4. METRICS                                                          │
│  ┌────────────┐ ┌──────────────┐ ┌────────────┐ ┌───────────────┐   │
│  │ Baseline   │ │ Longitudinal │ │ Dosimetry  │ │ Comparisons   │   │
│  │ metrics,   │ │ change       │ │ D95, V50   │ │ group stats   │   │
│  │ outliers   │ │ analysis     │ │ in GTV     │ │               │   │
│  └────────────┘ └──────────────┘ └────────────┘ └───────────────┘   │
│  ┌─────────────────────┐ ┌──────────────────────────────────────┐   │
│  │ Predictive modeling │ │ Survival analysis                    │   │
│  │ elastic net + CV    │ │ competing risks, IPCW, Kaplan-Meier │   │
│  └─────────────────────┘ └──────────────────────────────────────┘   │
└──────────────────────────┬───────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────────────┐
│  5. COMPARE CORES (optional)                                         │
│  11 tumor core delineation methods · Dice/HD95 agreement matrices    │
└──────────────────────────┬───────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────────────┐
│  POST-HOC ANALYSIS (Python)                                          │
│  Log parsing · Graph extraction (Gemini vision) · Statistical        │
│  relevance scoring · Cross-DWI-type comparison · HTML/PDF report     │
└──────────────────────────────────────────────────────────────────────┘
```

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
- [AveryLoop](#averyloop)
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
|----------|:-:|:-:|:-:|
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
|---------|:-:|
| MATLAB R2021a+ | Yes |
| Statistics and Machine Learning Toolbox | Yes |
| Image Processing Toolbox | Yes |

### External Tools

| Tool | Purpose | Installation |
|------|---------|-------------|
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
|----------|------------|
| **Windows** | [Docker Desktop](https://www.docker.com/products/docker-desktop/) with WSL 2 backend enabled |
| **macOS** | [Docker Desktop](https://www.docker.com/products/docker-desktop/) |
| **Linux** | Docker Engine 20.10+ and Docker Compose v2.0+ |

Verify your installation:

```bash
docker --version
docker compose version
```

#### Quick Start

```bash
# Build the image
docker compose build

# Run the full pipeline (mount your data directory)
docker compose run --rm pipeline

# Run tests
docker compose run --rm test
```

See [`docs/DOCKER.md`](docs/DOCKER.md) for data mounting, GPU passthrough, and troubleshooting.

---

## Configuration

Copy the example configuration and update paths for your environment:

```bash
cp config.example.json config.json
cp averyloop_config.example.json averyloop_config.json
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
|-------|------------|
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

See [`config.example.json`](config.example.json) for all available fields and threshold parameters.

### AveryLoop Configuration

AveryLoop has its own config file. Copy the example and customize:

```bash
cp averyloop_config.example.json averyloop_config.json
```

Key fields include `exit_strategy` (`"classic"`, `"diminishing_returns"`, or `"both"`), diminishing returns thresholds (`dr_window`, `dr_max_merge_rate`, `dr_max_avg_importance`, `dr_min_file_repeats`, `dr_max_audit_score`), API settings (`anthropic_api_key`, `audit_model`, `fix_model`, `judge_model`, `review_model`), orchestrator knobs (`max_api_retries`, `retry_base_delay`, `max_file_chars`), and RAG settings (`rag_enabled`, `rag_db_path`, `rag_top_k`, `rag_min_relevance`). All fields have sensible defaults — the file is optional. See [`averyloop_config.example.json`](averyloop_config.example.json) for the full template.

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

> **File deletion safety:** The pipeline never deletes files or directories it did not create. Pipeline-created directories contain a `.pipeline_created` sentinel file; all cleanup code verifies this sentinel before removing a directory. The `dwi_vectors*.mat` nameset is reserved for the user's curated files — the pipeline never reads, writes, or deletes files matching that pattern. The pipeline's own voxel cache lives under `pipeline_voxels*.mat`.

### Pre-Pipeline Data Check

Before running the pipeline, you can verify that your patient data directory is properly structured:

```matlab
addpath('pipeline/core', 'pipeline/utils', 'pipeline/dependencies');
report = patient_data_check('config.json');
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
|------|--------|------------|
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

% Run after data has been loaded (requires pipeline_voxels and summary_metrics on disk)
run_dwi_pipeline('config.json', {'compare_cores'});
```

**Output files** (saved to the DWI-type subfolder):

| File | Description |
|------|------------|
| `compare_core_results_{type}.mat` | Full comparison struct: pairwise Dice and HD95 matrices, per-method volume fractions, fallback flags |
| `core_method_dice_heatmap_{type}.png` | 11×11 heatmap of mean pairwise Dice coefficients |
| `core_method_hd95_heatmap_{type}.png` | 11×11 heatmap of mean pairwise 95th-percentile Hausdorff distance (only when 3D GTV masks are available) |
| `core_method_volume_comparison_{type}.png` | Bar chart of mean core volume (% of GTV) per method with error bars |
| `core_method_fallbacks_{type}.png` | Bar chart showing how often each method fell back to `adc_threshold` (e.g., fDM at baseline, spatial methods without 3D masks) |
| `compare_core_methods_output_{type}.txt` | Console log |

The MAT file contains a `compare_results` struct with fields: `method_names`, `mean_dice_matrix`, `mean_hd95_matrix`, `volume_fractions`, `fallback_flags`, `all_dice` (per-patient Dice matrices), and `all_hd95` (per-patient HD95 matrices).

---

## Running Tests

```matlab
run('pipeline/tests/run_all_tests.m')
```

The test suite includes 133 MATLAB test files and 47 Python test files covering:

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
- `GEMINI_API_KEY` environment variable (only needed for vision-based graph analysis)

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

# Individual scripts can still be run standalone
python analysis/parsers/batch_graph_analysis.py
python analysis/parsers/parse_log_metrics.py [saved_files_path]
python analysis/parsers/statistical_relevance.py [saved_files_path]
```

| Script | Description |
|--------|------------|
| `run_analysis.py` | Orchestrator: runs the full analysis workflow with CLI flags |
| `shared.py` | Shared utilities: folder discovery, DWI type parsing, regex extractors |
| **`parsers/`** | **Data extraction subpackage** |
| `parsers/batch_graph_analysis.py` | Sends all pipeline graph images to Google Gemini vision API; extracts axes, trends, inflection points into a structured CSV |
| `parsers/parse_log_metrics.py` | Direct parsing of MATLAB log files for Wilcoxon p-values, AUC, hazard ratios, GLME results |
| `parsers/parse_csv_results.py` | Direct parsing of pipeline CSV exports with cross-DWI significance comparison |
| `parsers/parse_mat_metrics.py` | Parses MATLAB `.mat` output files (core comparison, dosimetry, summary metrics) into JSON |
| `parsers/statistical_relevance.py` | Scores and ranks findings by clinical significance, effect size, and novelty |
| **`reports/`** | **Report generation subpackage** |
| `reports/generate_report.py` | Builds final HTML/PDF report with templated sections, embedded figures, and ranked findings |

---

## AveryLoop

AveryLoop is an **external package**: [`averyloop`](https://github.com/akarlin3/averyLoop). It provides an automated audit-fix-evaluate cycle that uses the Claude API to iteratively improve the codebase via a four-agent pipeline (audit → implement → review → merge) with RAG-enhanced context retrieval.

### Installation

```bash
pip install -r analysis/requirements.txt
```

### Configuration

Copy `project_config.example.yaml` to `project_config.yaml` and edit as needed. Runtime tuning (API models, token limits, RAG settings) uses `averyloop_config.json`.

### Usage

```bash
# v2 pipeline (recommended) — four-agent RAG-enhanced loop
python -m averyloop.orchestrator_v2 [--max-iterations N] [--dry-run] [--single-iteration]

# v1 fallback — single-pass orchestrator
python -m averyloop.orchestrator_v1 [--max-iterations N] [--dry-run] [--single-iteration]

# RAG index management
python -m averyloop.rag.indexer --stats
python -m averyloop.rag.indexer --force-rebuild --stats
```

**Requirement:** `ANTHROPIC_API_KEY` environment variable must be set.

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
│   ├── utils/                      #   Helper utilities (76 files)
│   │   ├── parse_config.m          #     Configuration parser
│   │   ├── safe_load_mask.m        #     Secure .mat loading
│   │   ├── escape_shell_arg.m      #     Shell argument escaping
│   │   ├── init_scan_structs.m     #     Scan data structure initialization
│   │   └── ...
│   ├── tests/                      #   Test suite (127 test files)
│   │   ├── run_all_tests.m         #     Master test runner
│   │   ├── benchmarks/             #     Performance benchmarks (7 files)
│   │   └── diagnostics/            #     Diagnostic spot-checks (6 files)
│   ├── dependencies/               #   Third-party scripts (read-only)
│   └── .octave_compat/             #   GNU Octave compatibility shims (24 files)
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
│   └── tests/                      #   Python test suite (43 test files)
├── project_config.example.yaml     # AveryLoop project config template
├── averyloop_config.example.json   # AveryLoop runtime config template
└── .agents/                        # AI agent configuration
    ├── rules/                      #   Agent safety rules
    └── workflows/                  #   Structured workflows
```

---

## Contributing

Contributions are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

---

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{pancData3,
  author = {Karlin, Avery and Aliotta, Eric},
  title = {pancData3: A Pipeline for Pancreatic DWI Analysis},
  url = {https://github.com/akarlin3/pancData3},
  version = {2.3.2},
  year = {2026}
}
```

See also [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

---

## License

This project is licensed under the [GNU Affero General Public License v3.0](LICENSE).
