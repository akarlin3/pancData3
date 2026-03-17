# pancData3

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-blue?logo=mathworks)](https://www.mathworks.com/products/matlab.html)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-2.0.1-blue)](#citation)
[![Tests](https://img.shields.io/badge/tests-106%20MATLAB%20%2B%2030%20Python%20files-brightgreen)](#running-tests)

**A MATLAB-based analysis pipeline for pancreatic DWI (Diffusion-Weighted Imaging) research.**

Developed at [Memorial Sloan Kettering Cancer Center](https://www.mskcc.org/), this pipeline processes MRI data to fit IVIM and ADC diffusion models, apply deep learning denoising, correlate findings with radiotherapy dose maps, and perform survival analysis for treatment response prediction.

**Current version:** 2.0.1 — see [CHANGELOG.md](CHANGELOG.md) for details.

---

## Table of Contents

- [Supported Platforms](#supported-platforms)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Running Tests](#running-tests)
- [Post-Hoc Analysis Scripts](#post-hoc-analysis-scripts)
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

For reproducible containerized execution, see the [Docker usage guide](docs/DOCKER.md).

```bash
# Build the image
docker build -t pancdata3:latest .

# Run the full pipeline
docker run --rm \
  -v /path/to/patient_data:/opt/pancData3/data:ro \
  -v /path/to/output:/opt/pancData3/output \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  pancdata3:latest pipeline
```

---

## Configuration

Copy the example configuration and update paths for your environment:

```bash
cp config.example.json config.json
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
| `texture_3d` | When true, compute GLRLM texture features in 13 3D directions (axis-aligned, face-diagonal, body-diagonal) instead of 4 in-plane 2D directions. Falls back to 2D for single-slice inputs (default: `true`) |
| `use_gpu` | Offload ADC WLS fitting and DnCNN inference to a CUDA GPU via `gpuArray`. Requires Parallel Computing Toolbox and a CUDA-capable GPU. Falls back to CPU when unavailable (default: `false`) |
| `gpu_device` | 1-based index of the CUDA GPU device to use when `use_gpu` is true (default: `1`) |
| `run_imputation_sensitivity` | Compare KNN imputation against LOCF, Mean, and Linear Interpolation alternatives (default: `false`) |
| `fit_time_varying_cox` | Fit stratified and extended Cox models when proportional hazards assumption is violated (default: `true`) |
| `export_validation_model` | Export trained predictive model for external validation (default: `false`) |
| `external_validation_data` | Path to external validation data folder (default: `""`) |
| `auxiliary_biomarker_csv` | Path to CSV file containing non-DWI auxiliary biomarkers (default: `""`) |
| `use_auxiliary_biomarkers` | Include auxiliary biomarkers in the predictive modeling feature matrix (default: `false`) |

See [`config.example.json`](config.example.json) for all available fields and threshold parameters.

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

The test suite includes 106 test files covering:

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
| `report/sections/` | Section builder modules for the HTML report (16 submodules: metadata, main_results, statistical_reporting, manuscript, enrollment, supplemental, gallery, graph_overview, cross_dwi, correlations, effect_sizes, model_diagnostics, power_analysis, discussion, publication, _helpers) |

### Report Features (v2.0.0)

The generated HTML/PDF report includes:

- **Cover page** — print-only title page with run timestamp, DWI types, and graph count
- **Table of Contents** — grouped two-column TOC visible on screen; occupies its own page in the PDF
- **Part breaks** — 6 logical divisions (Overview / Data / Statistics / Outcomes / Discussion / Appendix) force page breaks in the PDF
- **PDF page numbers** — WeasyPrint-native footer: `Page N of M`, suppressed on the cover page
- **Clinical context** — RTOG D95/V50 benchmarks, evidence hierarchy (High/Moderate/Exploratory), Research Use Only disclaimer
- **Improved parsers** — NaN/Inf → JSON `null`, timepoint normalization, parse-failure warnings, IPCW range validation

### Analysis Test Suite

The analysis scripts have a comprehensive Python test suite (1164 tests across 28 files) using pytest:

```bash
cd analysis/tests && python -m pytest -v
```

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
│   ├── utils/                      #   Helper utilities (68 files)
│   │   ├── parse_config.m          #     Configuration parser
│   │   ├── safe_load_mask.m        #     Secure .mat loading
│   │   ├── escape_shell_arg.m      #     Shell argument escaping
│   │   ├── init_scan_structs.m     #     Scan data structure initialization
│   │   └── ...
│   ├── tests/                      #   Test suite (106 test files)
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
│   └── tests/                      #   Python test suite (30 test files, 1186 tests)
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
  version   = {2.0.1},
  url       = {https://github.com/akarlin3/pancData3},
  license   = {AGPL-3.0}
}
```

See [CITATION.cff](CITATION.cff) for a machine-readable citation file.

---

## License

This project is licensed under the GNU Affero General Public License v3.0. See [LICENSE](LICENSE) for details.

Copyright (c) 2026 Avery Karlin
