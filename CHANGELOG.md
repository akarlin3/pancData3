# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [2.0.1] - 2026-03-16

### Added
- **`clear_pipeline_cache.m`** (`pipeline/utils/`): Extracted from `run_dwi_pipeline.m` — removes pipeline-generated `.mat` cache files with once-per-session guard, protected files list, and sentinel checks
- **`setup_output_folders.m`** (`pipeline/utils/`): Extracted from `run_dwi_pipeline.m` — creates or reuses the master pipeline output folder with timestamped auto-creation and sentinel file
- **`load_baseline_from_disk.m`** (`pipeline/utils/`): Extracted from `run_dwi_pipeline.m` — loads persisted `metrics_baseline` outputs from `.mat` file
- **`resolve_scan_days.m`** (`pipeline/utils/`): Extracted from `run_dwi_pipeline.m` — three-level scan day resolution (DICOM dates → config fallback → defaults)
- **Report section submodules**: Split oversized report section modules into focused submodules for maintainability:
  - `gallery.py` — figure gallery and appendix sections
  - `manuscript.py` — manuscript-ready findings, predictive performance, results draft
  - `publication.py` — reporting checklist and journal guide
  - `statistical_reporting.py` — statistical significance and broad statistical overview
- **Pinned Python dependency lock file** (`analysis/requirements-lock.txt`)
- `pytest` and `scipy` added to `analysis/requirements.txt`
- 4 new Python test files: `test_cross_reference_dwi.py`, `test_cross_reference_summary.py`, `test_statistical_by_graph_type.py`, `test_statistical_relevance.py` (Python test suite now 22 files)
- New MATLAB unit tests for core pipeline modules enhancing test coverage

### Changed
- **Orchestrator refactoring**: `run_dwi_pipeline.m` reduced in complexity by extracting 4 helper functions to `pipeline/utils/`
- **Report sections refactoring**: `data_sections.py`, `main_results.py`, and `discussion.py` significantly reduced by moving code to new focused submodules; backward compatibility maintained via `__init__.py` re-exports
- CI workflow removed (continuous integration disabled)

### Removed
- WIN-53O1VVN1FV6 duplicate directories (Windows artifact cleanup) — removed duplicate copies of `pipeline/core/`, `pipeline/utils/`, `pipeline/.octave_compat/`, and `analysis/report/sections/` that were inadvertently created on a Windows machine

---

## [2.0.0] - 2026-03-14

### Changed
- **License**: Changed from MIT to GNU Affero General Public License v3.0 (AGPL-3.0)

### Fixed
- **HD95 heatmap readability**: Replaced `hot` colormap with `parula` to fix yellow-on-yellow and dark-red-on-dark-red text contrast issues; added axis labels and colorbar label
- **Dice heatmap readability**: Added missing axis labels (`Core Method`), colorbar label (`Dice Coefficient`), and increased font sizes
- **Parameter maps resolution**: Increased figure size, font sizes, and switched to 150 DPI output (`print -dpng -r150`) to fix small/blurry labels
- **Feature histogram empty-group annotation**: Made LF n=0 warning more prominent with bold text, background box, and border
- **Feature boxplot single-group annotation**: Added visible warning when only one outcome group is present (e.g., no LF events)
- **Dose vs Diffusion scatter small-group warning**: Added on-plot annotation when LF group has n<3, warning that inference is unreliable
- **Metric_Set figure improvements**: Added per-group sample sizes (LC/LF counts) to subplot titles and boxplot labels; improved "Insufficient Data" display with visible in-axes message

---

For older versions (v2.0.0-rc.1 and earlier), see [CHANGELOG_ARCHIVE.md](CHANGELOG_ARCHIVE.md).
