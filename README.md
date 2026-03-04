# pancData3

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-blue?logo=mathworks)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-63%20files-brightgreen)](#running-tests)
[![Octave Compatible](https://img.shields.io/badge/Octave-compatible-orange?logo=gnu)](https://www.gnu.org/software/octave/)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)

**A MATLAB-based analysis pipeline for pancreatic DWI (Diffusion-Weighted Imaging) research.**

Developed at [Memorial Sloan Kettering Cancer Center](https://www.mskcc.org/), this pipeline processes MRI data to fit IVIM and ADC diffusion models, apply deep learning denoising, correlate findings with radiotherapy dose maps, and perform survival analysis for treatment response prediction.

---

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Running Tests](#running-tests)
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
- **Octave Compatibility** -- Shim layer for running on GNU Octave environments

---

## Requirements

### MATLAB Toolboxes

| Toolbox | Required |
|---|---|
| MATLAB R2021a+ | Yes |
| Statistics and Machine Learning Toolbox | Yes |
| Image Processing Toolbox | Yes |

### External Tools

| Tool | Purpose |
|---|---|
| [dcm2niix](https://github.com/rordenlab/dcm2niix) (MRIcroGL) | DICOM-to-NIfTI conversion |

---

## Installation

```bash
git clone https://github.com/akarlin3/pancData3.git
cd pancData3
```

Then in MATLAB:

```matlab
addpath('core', 'utils', 'dependencies');
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

See [`config.example.json`](config.example.json) for all available fields.

---

## Usage

### Full Pipeline (all DWI types)

Runs Standard, dnCNN, and IVIMnet sequentially with a parallel pool (max 2 workers). Creates a timestamped output folder (`saved_files_YYYYMMDD_HHMMSS/`) containing all logs, figures, and results:

```matlab
execute_all_workflows
```

### Single Pipeline Run

```matlab
addpath('core', 'utils', 'dependencies');

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
├── execute_all_workflows.log        # Master workflow log
├── test_suite_output.log            # Test suite output
├── Standard/                        # Per-DWI-type results
│   ├── pipeline_log_Standard.txt    # Pipeline console log
│   ├── *.txt                        # Per-module logs
│   └── *.fig / *.png                # Figures
├── dnCNN/
└── IVIMnet/
```

### Pre-Pipeline Data Check

Before running the pipeline, you can verify that your patient data directory is properly structured:

```matlab
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

---

## Running Tests

```matlab
run('tests/run_all_tests.m')
```

The test suite includes 63 test files covering:

- **Integration tests** -- End-to-end pipeline validation
- **Unit tests** -- Individual module correctness
- **Leakage tests** -- Data leakage prevention in CV and imputation
- **Security tests** -- Safe file loading, shell argument escaping
- **Smoke tests** -- Visualization output verification
- **Benchmarks** -- Performance comparisons for optimized algorithms
- **Static analysis** -- Source code standards and naming conventions

Tests generate a code coverage report for `core/` and `utils/`.

---

## Repository Structure

```
pancData3/
├── run_dwi_pipeline.m          # Main orchestrator entry point
├── execute_all_workflows.m     # Sequential multi-type runner
├── patient_data_check.m       # Pre-pipeline data validation
├── config.example.json         # Configuration template
├── core/                       # Pipeline modules (17 files)
│   ├── load_dwi_data.m         #   Data loading & model fitting
│   ├── sanity_checks.m         #   Data validation
│   ├── visualize_results.m     #   Visualization generation
│   ├── process_single_scan.m   #   Per-scan DICOM/model processing
│   ├── metrics_baseline.m      #   Baseline metric computation
│   ├── metrics_survival.m      #   Survival analysis
│   └── ...
├── utils/                      # Helper utilities (19 files)
│   ├── parse_config.m          #   Configuration parser
│   ├── safe_load_mask.m        #   Secure .mat loading
│   ├── escape_shell_arg.m      #   Shell argument escaping
│   ├── init_scan_structs.m     #   Scan data structure initialization
│   ├── parfor_progress.m       #   Parallel loop progress reporting
│   ├── text_progress_bar.m     #   Text-based progress bar display
│   ├── octave_compat/          #   Octave compatibility shims (20 files)
│   └── ...
├── tests/                      # Test suite (63 test files)
│   ├── run_all_tests.m         #   Master test runner
│   ├── benchmarks/             #   Performance benchmarks (7 files)
│   └── diagnostics/            #   Diagnostic spot-checks (5 files)
├── dependencies/               # Third-party scripts (read-only)
└── .agents/                    # AI agent configuration
    ├── rules/                  #   Agent safety rules
    └── workflows/              #   Structured workflows
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
  url       = {https://github.com/akarlin3/pancData3},
  license   = {MIT}
}
```

See [CITATION.cff](CITATION.cff) for a machine-readable citation file.

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

Copyright (c) 2026 Avery Karlin
