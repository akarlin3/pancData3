# pancData3 — Onboarding

A practical guide to getting the pipeline running and understanding what it produces. For the project's decision rationale, open issues, and successor pickup checklist, see the separate `HANDOFF_NOTES.md` (kept outside the public repo).

**Repository:** https://github.com/akarlin3/pancData3 (AGPL-3.0)
**License:** AGPL-3.0
**Working branch:** `v2.1-dev`

> **Version note:** README and CHANGELOG advertise `v2.1.0`. The actual implemented state is `v2.4` (with v2.3 cross-pipeline / pruning / dose-coverage workstreams + v2.4 audit-driven improvements applied). Treat the working tree as authoritative until the version string is bumped.

---

## 1. Running the pipeline cold

### 1.1 Required environment

| Layer | Requirement | Notes |
|---|---|---|
| MATLAB | R2021a+ | Statistics and Machine Learning Toolbox + Image Processing Toolbox required. Pre-flight check is in `pipeline/utils/initialize_pipeline.m`. |
| Python | 3.12 | Use the `py -3.12` launcher on Windows to avoid hitting an older interpreter. |
| dcm2niix | latest | Path goes in `config.json` under `dcm2nii_call`. |
| MATLAB Runtime | r2024a (default) | Only needed for Docker-based runs. |
| OS | Windows / macOS / Linux | All three are CI-tested. |

MATLAB itself is accessed via remote desktop into MSK servers — you do not need a local MATLAB install for pipeline execution, only for development.

### 1.2 First-run checklist

```bash
# 1. Clone
git clone https://github.com/akarlin3/pancData3.git
cd pancData3
git checkout v2.1-dev

# 2. Configure
cp config.example.json config.json
cp analysis_config.example.json analysis_config.json
cp improvement_loop_config.example.json improvement_loop_config.json
# Then edit config.json: set dataloc, dcm2nii_call, clinical_data_sheet

# 3. Python deps (analysis layer)
py -3.12 -m venv .venv
.\.venv\Scripts\activate     # PowerShell
pip install -r analysis/requirements-lock.txt

# 4. API key for the improvement loop (PowerShell)
$env:ANTHROPIC_API_KEY = "sk-..."           # session only
# Persistent:
[System.Environment]::SetEnvironmentVariable("ANTHROPIC_API_KEY", "sk-...", "User")
```

In MATLAB (from `pipeline/`):

```matlab
% Pre-flight tests
run('pipeline/tests/run_all_tests.m')

% Single DWI type
run_dwi_pipeline('config.json')

% All three DWI types in sequence
execute_all_workflows
```

### 1.3 Pipeline anatomy at a glance

- **Input:** DICOM DWI series + clinical spreadsheet (`MASTER_pancreas_DWIanalysis.xlsx` by default).
- **Processing strategies:** `Standard`, `dnCNN` (denoised), `IVIMnet` (NN-fitted) — selected via `dwi_type` in `config.json`.
- **Steps (in order):** `load → sanity → metrics_baseline → metrics_longitudinal → metrics_dosimetry → metrics_stats_comparisons → metrics_stats_predictive → metrics_survival → visualize`. The `compare_cores` step is opt-in.
- **Output:** Timestamped `saved_files_YYYYMMDD_HHMMSS/` folder with logs, figures, and `.mat` results. Routed to Google Drive via `rclone`.
- **Analysis layer:** `analysis/run_analysis.py` consumes those outputs and produces an HTML+PDF report.

### 1.4 Windows / PowerShell quirks

These trip up developers coming from a Unix background:

- Replace `&&` with separate commands, or use `;` if you do not care about exit-code chaining.
- Replace `tail -n 50 file` with `Get-Content file | Select-Object -Last 50`.
- `Select-String` requires named parameters, not positional (e.g., `Select-String -Path file.txt -Pattern "AUC"`).
- ChromaDB (used by the improvement loop's RAG layer) **must not** live on a OneDrive-synced path — OneDrive corrupts HNSW binary files. Keep the persistence directory on a non-synced local disk.
- The improvement loop log determines its own location at import time via `os.getcwd()`. Always launch it from the pancData3 repo root, not from a subdirectory.

### 1.5 Key files to know your way around

| Purpose | File |
|---|---|
| Master orchestrator | `pipeline/run_dwi_pipeline.m` |
| All-DWI sweep | `pipeline/execute_all_workflows.m` |
| Config parsing + defaults | `pipeline/utils/parse_config.m` |
| Tumor core extraction (all 11 methods) | `pipeline/utils/extract_tumor_core.m` |
| Cross-method comparison | `pipeline/core/compare_core_methods.m` |
| Summary metric aggregation | `pipeline/core/compute_summary_metrics.m` |
| Dosimetry | `pipeline/core/metrics_dosimetry.m` |
| Survival modeling | `pipeline/core/metrics_survival.m` |
| Predictive modeling | `pipeline/core/metrics_stats_predictive.m` |
| Python report orchestrator | `analysis/report/generate_report.py` |
| Improvement loop driver | `improvement_loop/orchestrator_v1.py` |

---

## 2. Scientific state

### 2.1 The cohort

42 pancreatic cancer patients treated on the Elekta Unity 1.5T MR-Linac with a 33/50 Gy SIB schedule in 5 fractions. Diffusion-weighted imaging is acquired at each fraction. Outcomes are local control vs. local failure, with competing-risk handling for non-cancer death.

### 2.2 The three DWI processing strategies

Every downstream analysis runs three times — once per strategy. This is the central methodological design of the pipeline:

- **Standard** — segmented IVIM + mono-exponential ADC fitting on raw DWI.
- **dnCNN** — DnCNN-denoised DWI, then standard fitting.
- **IVIMnet** — neural-network IVIM fitting.

The point of running all three is to test whether the biomarker findings are robust to the choice of denoising / fitting strategy. The cross-DWI agreement metric in the report tells you how often the three strategies tell the same story.

### 2.3 The 11 tumor core methods

Implemented in `extract_tumor_core.m`. Each method tries to identify the treatment-resistant sub-volume of the GTV. Categories:

- **Threshold-based (3):** `adc_threshold`, `d_threshold`, `df_intersection`
- **Statistical clustering (3):** `otsu`, `gmm`, `kmeans`
- **Spatial (2):** `region_growing`, `active_contours`
- **Distributional (1):** `percentile`
- **Spectral (1):** `spectral`
- **Temporal/longitudinal (1):** `fdm` (functional diffusion map)

`adc_threshold` is the default and the universal fallback when other methods fail.

### 2.4 Headline results (v2.4)

**Cross-pipeline Dice at Fx1 (WS1)**
Mask agreement across the 11 core methods at fraction 1, computed for each of the three DWI strategies. Used as input to the pruning decision below.

**Failure rates and pruning (WS2)**
- Failure rate is computed using `min_core_voxels=10` (not the older `min_vox_hist=100`, which inflated failure rates).
- Pruning cutoff: `max_core_failure_rate = 0.25`. **6 of 11 methods survive.**
- `adc_threshold` is retained regardless of cutoff (safety guard hardcoded in pipeline).

**Dose coverage vs. local control (WS3, the headline finding)**
- D95 (ADC-defined resistant sub-volume) ≈ **36.9 Gy** — well below the 45 Gy adequacy threshold.
- V50 (ADC sub-volume) ≈ **69%** — below the 90% radical-intent target.
- Spearman rs ≈ **−0.37** (D95 vs. local-failure proxy).
- GTV-volume confounding: **<5% HR change** after adjustment for any surviving method, i.e., the signal is not driven by tumor size.

**Clinical interpretation:** ADC-defined resistant regions are receiving roughly **26% less dose than the radical-intent threshold**. This motivates daily adaptive replanning targeted at the diffusion-defined sub-volume.

### 2.5 What the figures show

The PDF report (`analysis/report/generate_report.py`) is the canonical figure source. The most important figures, in order of weight:

1. Per-DWI D95 ADC pass/fail against the 45 Gy threshold (data_supplemental section).
2. Pairwise Dice heatmap across 11 core methods (cross_dwi section).
3. Dose-vs-diffusion scatter, LC vs LF colored, with per-group trend lines (`plot_scatter_correlations.m`).
4. Feature distribution histograms + boxplots, LC vs LF, with Wilcoxon p-values (`plot_feature_distribution.m`).
5. Survival curves (Cox PH, cause-specific hazards, IPCW-weighted).

---

## 3. Where to look next

- `CLAUDE.md` — architectural conventions and safety rules
- `CLAUDE_REFERENCE.md` — detailed module and test reference
- `CLAUDE_WORKFLOWS.md` — improvement-loop workflow
- `CONTRIBUTING.md` — review criteria, the "do not modify" list, dependency lock-file regeneration
- `CHANGELOG.md` / `CHANGELOG_ARCHIVE.md` — version history
- `docs/DOCKER.md` — containerized execution
