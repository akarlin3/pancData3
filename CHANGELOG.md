# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

#### Sub-Volume Stability Analysis
- **`pipeline/utils/compute_subvolume_stability.m`**: New pipeline step — computes Dice between each fraction's core mask and Fx1 baseline for all active core methods. Generates line plot of mean Dice vs fractions. Gated on `run_subvolume_stability` config flag.
- **`analysis/report/sections/subvolume_stability.py`**: New report section builder — method × timepoint table of mean Dice vs Fx1 with color coding (green ≥0.7, yellow 0.5–0.7, red <0.5).
- **`pipeline/tests/test_compute_subvolume_stability.m`**: 4 MATLAB unit tests (struct fields, Fx1 Dice=1.0, Dice range, dimensions).
- **`analysis/tests/test_subvolume_stability_section.py`**: 6 Python tests for the report section builder.

#### Dose-Response ROC Analysis
- **`pipeline/utils/compute_dose_response_roc.m`**: New pipeline step — ROC analysis on D95/V50 of sub-volumes to find optimal dose cutoff separating LC from LF. Bootstrap AUC CI, Youden index threshold. Gated on `run_dose_response_roc` config flag.
- **`analysis/report/sections/dose_response_roc.py`**: New report section builder — AUC, CI, optimal threshold, sensitivity/specificity table with clinical guidance text.
- **`pipeline/tests/test_compute_dose_response_roc.m`**: 4 MATLAB unit tests (struct fields, AUC above chance, threshold in range, AUC range).
- **`analysis/tests/test_dose_response_roc_section.py`**: 7 Python tests for the report section builder.

#### Risk-Dose Concordance
- **`pipeline/utils/compute_risk_dose_concordance.m`**: New pipeline step — compares elastic net risk model classifications against D95-based stratification using Cohen's kappa, concordance %, and combined AUC. Gated on `run_risk_dose_concordance` config flag.
- **`analysis/report/sections/risk_dose_concordance.py`**: New report section builder — kappa with interpretation, concordance %, complementary patients, combined AUC.
- **`pipeline/tests/test_compute_risk_dose_concordance.m`**: 4 MATLAB unit tests (struct fields, concordant/complementary scenarios, kappa range, missing data).
- **`analysis/tests/test_risk_dose_concordance_section.py`**: 7 Python tests for the report section builder.

#### Per-Method CoR from Fx1 Repeats
- **`pipeline/utils/compute_per_method_cor.m`**: New pipeline step — extends wCV/CoR repeatability analysis to all active core methods using Fx1 repeat scans. Generates bar chart with reproducibility thresholds. Gated on `run_per_method_cor` config flag.
- **`analysis/report/sections/per_method_cor.py`**: New report section builder — method × wCV/CoR table with color coding (green <15%, yellow 15–30%, red >30%).
- **`pipeline/tests/test_compute_per_method_cor.m`**: 5 MATLAB unit tests (struct fields, CoR positivity, identical repeats near-zero, dimensions, repeat count).
- **`analysis/tests/test_per_method_cor_section.py`**: 6 Python tests for the report section builder.

#### GTV Volume Confounding Check
- **`pipeline/utils/compute_gtv_confounding.m`**: New pipeline step — tests whether D95-outcome association is confounded by GTV shrinkage via Spearman correlation, Wilcoxon rank-sum, and GTV-adjusted Cox PH. Gated on `run_gtv_confounding` config flag.
- **`analysis/report/sections/gtv_confounding.py`**: New report section builder — D95-GTV correlation, unadjusted/adjusted HR, % change, confounding flag with warning box.
- **`pipeline/tests/test_compute_gtv_confounding.m`**: 4 MATLAB unit tests (struct fields, confounded/independent scenarios, missing GTV graceful handling).
- **`analysis/tests/test_gtv_confounding_section.py`**: 6 Python tests for the report section builder.

#### Pipeline Integration
- **`pipeline/utils/parse_config.m`**: 5 new config flags: `run_subvolume_stability`, `run_dose_response_roc`, `run_risk_dose_concordance`, `run_per_method_cor`, `run_gtv_confounding` (all default `false`). Added to `logical_fields` for type validation.
- **`pipeline/utils/dispatch_pipeline_steps.m`**: 5 new step dispatchers with wrapper functions, gated on config flags and data availability.
- **`pipeline/utils/prepare_pipeline_session.m`**: 5 new step injection entries with dependency-aware ordering.
- **`analysis/parsers/parse_mat_metrics.py`**: Extended to parse 5 new MAT file types (per_method_cor, subvolume_stability, dose_response_roc, gtv_confounding, risk_dose_concordance).
- **`analysis/report/generate_report.py`**: 5 new section builders wired into the statistics section pipeline.
- **`analysis/report/sections/__init__.py`**: 5 new section builder exports.

#### Documentation
- Updated CLAUDE.md: 5 new utility modules (76→81), 5 new report sections (40→45), 5 new MATLAB tests (127→132), 5 new Python tests (42→47), 5 new config fields.
- Updated CLAUDE_REFERENCE.md: matching test and module table entries.
- Updated README.md: Python test count 42→47.

## [2.3.1] - 2026-03-29

### Added

#### Core Method Failure Rates
- **`pipeline/utils/extract_tumor_core.m`**: Added optional second output `fit_info` struct with failure mode details (method, success, fallback, empty_mask, insufficient_voxels, all_nan_input, n_core_voxels, error_msg). Non-breaking change — all existing single-output calls continue to work.
- **`pipeline/utils/compute_core_failure_rates.m`**: New pipeline step — aggregates failure rates for all 11 core methods x 3 DWI pipelines x all patients x all timepoints. Tracks four failure categories: fallback, empty mask, insufficient voxels, and all-NaN input. Generates stacked bar chart and console summary table.
- **`analysis/report/sections/failure_rates.py`**: New report section builder — color-coded HTML table of failure rates sorted by total failure descending (green <10%, yellow 10-25%, red >25%).
- **`analysis/parsers/parse_mat_metrics.py`**: Extended to parse `core_failure_rates_{type}.mat` files.
- **`pipeline/utils/parse_config.m`**: New config field `run_core_failure_rates` (default `false`).
- **`pipeline/tests/test_extract_tumor_core_fit_info.m`**: 9 MATLAB unit tests for the fit_info output.
- **`pipeline/tests/test_compute_core_failure_rates.m`**: 5 MATLAB unit tests for the failure rate aggregation.
- **`analysis/tests/test_failure_rates_section.py`**: 8 Python tests for the report section builder.

#### Core Method Outcome Analysis
- **`pipeline/core/analyze_core_method_outcomes.m`**: New pipeline step — for each surviving core method, runs univariable Cox PH on 4 dose metrics (D95/V50 of ADC and D sub-volumes) against local control, performs KM + log-rank for the best metric, ranks methods by p-value, and generates KM plots and summary bar charts.
- **`analysis/report/sections/core_method_outcomes.py`**: New report section builder — ranking table with HR (95% CI), p-values, and interpretation paragraph identifying the best method.
- **`analysis/parsers/parse_mat_metrics.py`**: Extended to parse `core_method_outcomes_{type}.mat` files with full method_results/univariable/km structure.
- **`pipeline/utils/parse_config.m`**: New config field `run_core_method_outcomes` (default `false`).
- **`pipeline/tests/test_analyze_core_method_outcomes.m`**: 9 MATLAB unit tests covering struct fields, HR positivity, KM fields, ranking, significance, insufficient events, empty dosimetry, and missing methods.
- **`analysis/tests/test_core_method_outcomes_section.py`**: 7 Python tests for the report section builder.

#### Core Method Pruning
- **`pipeline/utils/filter_core_methods.m`**: New function — prunes core methods exceeding failure rate threshold (`max_core_failure_rate`), manually excluded methods (`excluded_core_methods`), or methods with insufficient median voxel counts (`min_core_voxels`). Safety guard ensures `adc_threshold` is never pruned.
- **`pipeline/utils/dispatch_pipeline_steps.m`**: Added core method pruning block that runs after failure rates and stores `active_core_methods` in `config_struct` for downstream steps.
- **`pipeline/core/compare_core_methods.m`**, **`pipeline/core/compute_summary_metrics.m`**, **`pipeline/core/metrics_dosimetry.m`**, **`pipeline/utils/compute_multi_core_metrics.m`**: Modified to check `config_struct.active_core_methods` and use it instead of hardcoded 11-method list, with backward-compatible fallback.
- **`pipeline/utils/parse_config.m`**: Three new config fields: `max_core_failure_rate` (default 1.0, disabled), `excluded_core_methods` (default empty), `min_core_voxels` (default 0, disabled).
- **`analysis/report/sections/pruning_results.py`**: New report section builder — shows pruned methods with reasons and retained methods list.
- **`analysis/parsers/parse_mat_metrics.py`**: Extended to parse `core_pruning_{type}.mat` files.
- **`pipeline/tests/test_filter_core_methods.m`**: 8 MATLAB unit tests covering all pruning modes and safety guards.
- **`analysis/tests/test_pruning_section.py`**: 7 Python tests for the report section builder.

#### Cross-Pipeline Dice Analysis
- **`pipeline/utils/compute_cross_pipeline_dice.m`**: New pipeline step — computes pairwise Dice coefficients between Standard, DnCNN, and IVIMNet pipelines for each of the 11 tumor core delineation methods at fraction 1. Measures whether DWI processing pipeline choice affects the spatial location of the treatment-resistant sub-volume.
- **`analysis/report/sections/cross_pipeline_dice.py`**: New report section builder — HTML table of mean Dice ± SD per core method and pipeline pair, with red highlighting for low agreement (Dice < 0.5).
- **`analysis/parsers/parse_mat_metrics.py`**: Extended to parse `cross_pipeline_dice_{type}.mat` files into the `mat_data` dict for report generation.
- **`pipeline/utils/parse_config.m`**: New config field `run_cross_pipeline_dice` (default `false`) with backwards-compatible default.
- **`pipeline/utils/dispatch_pipeline_steps.m`**: New `cross_pipeline_dice` step block with heatmap figure generation.
- **`pipeline/tests/test_compute_cross_pipeline_dice.m`**: 6 MATLAB unit tests covering struct fields, dimensions, range, identical pipelines, single patient, and empty pipeline edge cases.
- **`analysis/tests/test_cross_pipeline_dice_section.py`**: 8 Python tests for the report section builder.

### Documentation
- Updated version to 2.3.1 across CLAUDE.md, README.md, CHANGELOG.md, CITATION.cff, SECURITY.md
- Test file counts: MATLAB 124 → 127, Python 41 → 42 (43 including conftest.py)
- README utility count: 72 → 76, Octave shims: 21 → 24

---

## [2.2.0] - 2026-03-26

### Added

#### Octave Compatibility
- **`pipeline/.octave_compat/iscategorical.m`**: New shim — `iscategorical()` always returns `false` in Octave (no native categorical type)
- **`pipeline/.octave_compat/nanmedian.m`**: New shim — NaN-ignoring median for Octave
- **`pipeline/.octave_compat/string.m`**: New shim — `string()` type conversion for Octave
- **`TestCase.m`**: Added `verifyNotEqual`, `verifySize`, `assumeTrue`, and `assumeFail` assert methods to the Octave unittest shim
- **`TestRunner.m`**: Full `TestMethodSetup`/`TestMethodTeardown` discovery (parsed from source); per-test verification-failure accumulation and reset; robust teardown on failure
- **`TestSuite.m`**: Recursive `test_*.m` discovery via `**` glob; benchmark file inclusion; de-duplication by full path; `parseTestMethods` integration

#### Pipeline Utilities
- **`pipeline/utils/get_system_memory.m`**: New utility — cross-platform physical memory query (total and available GB); used by `load_dwi_data.m` to log memory headroom before parallel patient loading; returns `[NaN, NaN]` on unsupported platforms

#### Python Tests (improvement_loop v2 API)
- **`analysis/tests/test_implementer_agent.py`**: 10 tests for `ImplementResult` dataclass, `implement()` dry-run, file-not-found and branch-exists failure paths, and `_generate_fix()` API call contract
- **`analysis/tests/test_orchestrator_v2.py`**: 13 tests for `FindingState`/`IterationState` dataclasses, `run_loop()` dry-run and exit-condition behavior, rejected-finding non-merge guarantee, and `_print_agent_summary()` output
- **`analysis/tests/test_reviewer_agent.py`**: 23 tests for `ReviewVerdict` dataclass, `_generate_diff()`, `_parse_review_verdict()` (valid JSON, fenced, preamble, all invalid cases), critical-flag override (`LEAKAGE_RISK`/`PHI_RISK` force reject), and parse-failure fallback

### Fixed

#### Security / Error Handling
- **`pipeline/utils/safe_load_mask.m`**: Duplicate-variable security violations (unsafe class, oversized entry) now raise `error` instead of `warning` + `return`; eliminates the risk of silently continuing after a tampered `.mat` file is detected
- **`pipeline/utils/write_sentinel_file.m`**: Guard against `fopen` throwing on Windows when the parent directory doesn't exist; added `isfolder` pre-check and `try/catch` around `fopen`; also added confirmation log line on successful sentinel write

#### MATLAB Test Suite (R2025a compatibility — 10 test failures fixed)
- **Multiple test files**: Corrected `coxphfit` `Stratification` parameter usage, updated sentinel save message assertions, fixed struct field mismatches in `fit_time_varying_cox` interaction_models output (`penalized`, `stable_period_used`, `left_truncated` fields now included)
- **`pipeline/tests/ProgressBarPlugin.m`**: Added `DiaryFile` property and diary-restart logic; core modules that hijack `diary` during tests no longer silently drop the test log

#### Pipeline Logic
- **`pipeline/utils/fit_time_varying_cox.m`**: Sparse-period detection now skips the stable-period restriction when the identified stable region is < 20 % of the total time range (too short to be useful); reports an informational note instead of a warning
- **`pipeline/utils/escape_shell_arg.m`**: Extended control-character filtering (comprehensive Unicode sanitization)
- **`pipeline/utils/filter_collinear_features.m`**: Strict duplicate-variable validation
- **`pipeline/utils/run_loocv_risk_scores.m`**: Robust fallback-cleanup handling; additional edge case guards

### Changed

#### improvement_loop Package — v2 Structured API
- **`improvement_loop.agents.implementer`**: Added `ImplementResult` dataclass (`success`, `original_content`, `new_content`, `error`); `_generate_fix(finding, original_content)` API wrapper; `implement(finding, base_branch, dry_run)` function that replaces direct `apply_fix` calls in the orchestrator
- **`improvement_loop.agents.reviewer`**: Added `ReviewVerdict` dataclass (`verdict`, `reasoning`, `risk_flags`); `_generate_diff(original, new, filename)` unified-diff helper; `_parse_review_verdict(raw)` JSON parser (handles markdown fences, preamble, all validation failures); new `review(finding, old_content, new_content, dry_run)` replaces diff-based signature; `LEAKAGE_RISK`/`PHI_RISK` flags force rejection regardless of LLM verdict; parse failure returns `request_changes` + `EVALUATION_FAILED` flag
- **`improvement_loop.orchestrator_v2`**: `FindingState` updated — added `implement_result: Optional[ImplementResult]`, `review_verdict: Optional[ReviewVerdict]`, `error` defaults to `""`; `IterationState` updated — added `context: str`, `findings: List[FindingState]`, `dry_run` now optional; module-level `_audit(iteration, context, dry_run)`, `_implement(finding, base_branch, dry_run)`, `_review(finding, old_content, new_content, dry_run)` wrappers (monkeypatchable for tests); `_print_agent_summary(state)` per-iteration console summary; `run_loop` rewritten to use new agent wrappers

### Documentation
- Updated version to 2.2.0 across CLAUDE.md, README.md, CHANGELOG.md
- Octave shim count: 21 → 24; utility count: 72 → 73; Python test file count: 34 → 37
- Added new utility, shim, and test entries to CLAUDE.md, CLAUDE_REFERENCE.md

---

## [2.2.0-rc.2] - 2026-03-22

### Changed

#### Improvement Loop Extracted to External Package
- **`improvement_loop/`** directory removed — now installed as [`code-improvement-loop`](https://github.com/akarlin3/improvementLoop) via `pip install -r analysis/requirements.txt`
- **`project_config.yaml`** (gitignored) + **`project_config.example.yaml`** (committed) replace inline agent prompts and key file lists; contains full audit/review/judge system prompts, safety flags, RAG collection config
- **`improvement_loop_config.json`** retained for runtime tuning (API models, token limits, exit strategy)
- 12 improvement loop test files moved to the external package; pancData3 Python test suite: 46 → 34 files

#### Bug Fixes & Stability
- **`knn_impute_train_test.m`**: Restored after truncation; NaN handling in BH FDR correction
- **`load_dwi_data.m`**: Restored after truncation
- **`test_statistical_methods.m`**: Restored after truncation
- **`test_octave_shims.m`**: Moved to serial execution to avoid parallel conflicts
- **`test_source_code_standards.m`**: Updated to match actual implementations
- **`run_all_tests.m`**: Persistent variable scope fix for parallel capability caching
- Three improvement-loop regressions repaired (16 tests fixed)

### Documentation
- Updated CLAUDE.md, CLAUDE_WORKFLOWS.md, CLAUDE_REFERENCE.md, README.md for external improvement loop
- Python test suite count: 46 → 34 files (loop tests moved to external package)
- Repository structure updated to reflect `project_config.yaml` / `project_config.example.yaml`

---

## [2.2.0-rc.1] - 2026-03-22

### Changed

#### Test Quality & Cross-Platform Hardening (13 improvement loop iterations)
- **`test_benjamini_hochberg_fdr.m`**: Cross-platform RNG determinism — replaced `rng(123)` with platform-independent test data; added edge case tests; tests now call the production `benjamini_hochberg_fdr` function instead of reimplementing inline; assertions upgraded to `verifyEqual` with `AbsTol`
- **`test_parsave_dir_cache.m`**: Improved `parfor` test assertions
- **`test_text_progress_bar.m`**: Progress bar test assertions tightened
- **`test_statistical_methods.m`**: Added `k > available neighbors` edge case for KNN imputation; unsorted b-value test data fix
- **Octave compatibility**: Robust `strjoin` fallback, non-fatal `verify*` assertions with scalar expansion
- **`run_all_tests.m`**: Persistent variable scope fix for parallel capability caching
- **`initialize_pipeline.m`**: Preflight mode now logs active status
- **Vision skip assertions**: Tightened `test_script_outputs.py` skip condition

### Documentation
- Updated all documentation for v2.2.0 release
- Corrected MATLAB test suite count: 119 files (106 root + 7 benchmarks + 6 diagnostics); Python test suite: 46 files (1795 tests)

---

## [2.2.0-beta.1] - 2026-03-22

### Added

#### Configurable IVIM Optimizer Parameters
- **`fit_models.m`**: IVIM `lsqnonlin` optimizer settings (tolerances, iteration limits) are now configurable via `opts` struct fields (`optim_tol`, `func_tol`, `step_tol`, `max_iterations`, `max_func_evals`) with backward-compatible defaults; extracted `getfield_default()` local helper for clean fallback logic

#### Test Coverage
- **`test_implementer_agent.py`**: 5 new tests for `implement()` dry-run mode — happy path with real file, nonexistent file, no filesystem side effects, no git calls, and no API calls
- **`test_statistical_methods.m`**: Added explicit b-value count mismatch assertion in extra-value validation test

### Fixed
- **`implementer.py`**: Fixed `NameError` — `cfg` was referenced in `implement()` at line 110 but only defined inside `_generate_fix()`; added `cfg = _get_loop_config()` to `implement()` scope
- **`fit_models.m`**: Empty mask guard — when GTV mask contains no valid voxels, returns NaN-initialized parameter maps with a warning instead of crashing on empty array indexing

### Changed
- **`run_all_tests.m`**: Parallel capability check cached in a persistent variable to avoid repeated metaclass introspection; serial fallback now runs the complete test suite (parallel + serial partitions) instead of silently dropping the parallel-safe partition
- **`test_text_progress_bar.m`**: Octave skip converted from bare `return` to `assumeTrue()` for proper test framework skip reporting

### Documentation
- Python test suite: 1790 → 1795 tests

---

## [2.2.0-alpha.1] - 2026-03-21

### Added

#### Multi-Agent Improvement Loop (v2)
- **`improvement_loop/orchestrator_v2.py`**: Four-agent pipeline orchestrator (audit → implement → review → test & merge) replacing the single-pass v1 architecture
- **`improvement_loop/agents/` subpackage**: Extracted agent modules from the orchestrator:
  - **`auditor.py`** — RAG-enhanced code audit agent returning structured `Finding` objects
  - **`implementer.py`** — Fix implementation agent with branch creation, code generation, and syntax checking
  - **`reviewer.py`** — Code review quality gate with approve/request_changes/reject verdicts and critical risk flag enforcement (`LEAKAGE_RISK`, `PHI_RISK`)
  - **`_api.py`** — Shared API retry helper used by all agents

#### RAG (Retrieval-Augmented Generation)
- **`improvement_loop/rag/` subpackage**: Semantic code search via ChromaDB vector store:
  - **`chunker.py`** — Splits MATLAB, Python, Markdown, and JSON files into semantically meaningful chunks (functions, classes, sections) via regex-based parsing
  - **`indexer.py`** — ChromaDB persistent index with full/incremental build, improvement history indexing, and `--force-rebuild` / `--stats` CLI interface
  - **`retriever.py`** — Query interface with type/language/file filtering, relevance scoring, and agent-specific context builders (`get_context_for_audit`, `get_context_for_fix`, `get_context_for_review`)
- RAG-enhanced agents receive semantically relevant code context instead of hardcoded file lists; can be disabled via `rag_enabled: false` in config

#### New Report Section Builders (16 files)
- `data_overview.py`, `data_quality.py`, `data_supplemental.py` — Data inspection sections
- `main_results_summary.py`, `main_results_hypothesis.py`, `main_results_trends.py` — Granular main results
- `manuscript_findings.py`, `manuscript_performance.py`, `manuscript_results.py` — Manuscript-ready sub-sections
- `analysis_cross_dwi.py`, `analysis_graphs.py`, `analysis_features.py` — Analysis sub-sections
- `statistics_diagnostics.py`, `statistics_effects.py`, `statistics_robustness.py` — Statistics sub-sections
- `model_diagnostics.py` — Model diagnostic section builder

#### New Test Files (9 Python)
- `test_auditor_agent.py`, `test_implementer_agent.py`, `test_reviewer_agent.py` — Agent module tests
- `test_orchestrator_v2.py` — v2 orchestrator pipeline tests
- `test_chunker.py`, `test_indexer.py`, `test_retriever.py`, `test_rag_integration.py` — RAG subsystem tests
- `test_new_report_sections.py` — New report section builder tests

#### Configuration
- **`improvement_loop_config.example.json`**: Added `review_model`, `review_max_tokens`, `rag_enabled`, `rag_db_path`, `rag_top_k`, `rag_min_relevance` fields
- **`loop_config.py`**: `LoopConfig` dataclass extended with review agent and RAG settings (all with sensible defaults)

### Fixed
- **`improvement_loop/rag/indexer.py`**: Fixed duplicate chunk ID crash when a file has multiple functions with the same name (e.g., `imputation_sensitivity.m` has four local functions named `impute`). Chunk IDs now include start line: `file_path::name::L{start_line}`

### Changed
- `orchestrator_v2.py` is now the canonical loop driver; `orchestrator_v1.py` retained as legacy fallback
- Session cleanup and encapsulation improvements from improvement loop iterations
- Safe load variable name checking hardened
- JSON comment stripping edge cases fixed
- Windows shell escape ampersand handling fixed
- Survival result variable initialization guards added

### Documentation
- Updated all documentation files (CLAUDE.md, CLAUDE_REFERENCE.md, CLAUDE_WORKFLOWS.md, README.md) with accurate file counts
- MATLAB test suite: 121 → 122 files; Python test suite: 37 → 45 files (1576 → 1790 tests)
- Report sections: 18 → 37 files (16 new sub-section builders)
- Added 37 previously undocumented MATLAB test files to reference tables
- Expanded improvement loop documentation with v2 pipeline, RAG subsystem, and agent architecture

---

## [2.1.0] - 2026-03-20

### Added
- CI/CD pipeline with GitHub Actions — Python tests on Ubuntu, Windows, macOS; MATLAB file integrity check
- Auxiliary biomarker wiring — non-DWI biomarker CSV data now flows through elastic net and Cox survival models when `use_auxiliary_biomarkers` is enabled
- External validation round-trip — `apply_external_validation.m` now called from orchestrator when `external_validation_data` config field is set
- Report sections for decision curve analysis, NRI/IDI, texture features, and registration quality
- Longitudinal trajectory visualizations — waterfall, swimmer, and spider plots
- Programmatic improvement loop orchestrator (`improvement_loop/orchestrator_v1.py`)
- Structured finding schema with Pydantic validation
- Diminishing returns detection in loop exit condition
- Parallel implementation phase in orchestrator (configurable `--max-workers`)

### Fixed
- IPCW-weighted Cox standard error correction now uses sandwich variance ratio instead of naive `sqrt(mean(weights))` scaling (`metrics_survival.m`)
- Command injection vulnerability in `escape_shell_arg.m` Windows Unicode short-path lookup — unescaped argument was passed to `system()` during 8.3 name resolution
- Immortal time bias in competing risk censoring — competing event patients now use event time instead of total follow-up time (`metrics_survival.m`)
- Landmark day selection in time-dependent Cox models — largest gap heuristic replaced with clinically meaningful landmark (`metrics_survival.m`)
- Global `warning('off', 'all')` before `coxphfit` replaced with scoped warning suppression via `onCleanup` (`metrics_survival.m`)
- IVIM fitting with insufficient b-values (< 4) now skips fitting entirely instead of producing unreliable parameters (`fit_models.m`)
- Non-monotonic b-value validation relaxed to allow duplicate values, which are valid for averaging (`fit_models.m`)
- Memory optimization in DWI volume reshaping — reduced peak memory via chunked processing (`fit_models.m`)
- Diary file leak in `metrics_stats_predictive.m` — added `onCleanup` to ensure diary is turned off on error
- API key masking regex in `run_analysis.py` expanded to cover additional secret patterns

### Changed
- `improvement_loop_log.json` moved to `.gitignore` as a runtime artifact
- CLAUDE_WORKFLOWS.md updated — `orchestrator_v1.py` is now canonical loop driver

---

## [2.1.0-beta.2] - 2026-03-20

### Added

#### Improvement Loop Configuration System
- **`improvement_loop/loop_config.py`**: Centralised configuration for the improvement loop via `improvement_loop_config.json` — all tuneable knobs in a single `LoopConfig` dataclass with built-in defaults (same pattern as `parse_config.m` for the MATLAB pipeline)
- **`improvement_loop_config.example.json`**: Committed template with all configurable fields
- **Configurable exit strategy**: `exit_strategy` field selects `"classic"` (threshold-only), `"diminishing_returns"` (staleness detector), or `"both"` (default) for loop termination
- **Configurable API settings**: `anthropic_api_key`, `audit_model`, `fix_model`, `judge_model` fields — API key can live in the config file instead of requiring an environment variable
- **Configurable thresholds**: All classic exit thresholds (`importance_threshold`, `min_coverage_score`) and diminishing returns parameters (`dr_window`, `dr_max_merge_rate`, `dr_max_avg_importance`, `dr_min_file_repeats`, `dr_max_audit_score`) are tuneable via config

#### Diminishing Returns Detection
- **`evaluator.py` — `check_diminishing_returns()`**: Detects stale loops by checking four simultaneous conditions over the last N iterations (configurable via `dr_window`): low merge rate, low average importance, repeated file targeting, and no high audit scores
- **`should_continue_loop()`** updated to call the diminishing returns detector after classic checks, controlled by `exit_strategy` config

#### JSON Parsing Robustness
- **`orchestrator_v1.py` — `_sanitize_json_escapes()`**: Pre-sanitises raw API responses by escaping invalid backslash sequences (e.g. Windows paths like `C:\Users`) before JSON parsing
- **`_parse_findings()`**: Three-layer fallback chain — sanitised `json.loads()`, then `unicode_escape` decoding, then regex extraction of the JSON array

### Changed

#### Model Upgrade to Claude Opus 4.6
- **Improvement loop**: Default `audit_model`, `fix_model`, and `judge_model` switched from `claude-sonnet-4-20250514` to `claude-opus-4-6`
- **Vision analysis**: Default `claude_model` switched from `claude-sonnet-4-6` to `claude-opus-4-6` in `analysis_config.json`, `shared.py`, `batch_graph_analysis.py`, and `run_analysis.py`
- All hardcoded model strings replaced with config-driven values

#### Evaluator & Orchestrator Refactoring
- **`evaluator.py`**: `score_audit()` and `should_continue_loop()` now read model names and thresholds from `LoopConfig` instead of hardcoded constants; Anthropic client created lazily via `_get_client()` with config API key support
- **`orchestrator_v1.py`**: `_api_call_with_retry()`, `_run_audit()`, `_apply_single_fix()`, and `_collect_source_files()` all use `LoopConfig` for retries, delays, models, and file size limits; client created via `_get_client()`

#### Documentation
- **`README.md`**: Configuration section updated to mention both `config.json` and `improvement_loop_config.json`; new "Improvement Loop Configuration" subsection; Claude model example updated to Opus 4.6
- **`CLAUDE.md`**, **`CLAUDE_REFERENCE.md`**: Updated file counts, module tables, and test suite numbers

### Fixed
- `improvement_loop_config.json` added to `.gitignore` (may contain API keys)

### Test Suite
- Python test suite: 36 → 37 files (1559 → 1576 tests)
- New tests in `test_evaluator_finding.py`: `TestCheckDiminishingReturns` (5 unit tests), `TestShouldContinueLoopDiminishingReturns` (5 integration tests), `TestLoopConfig` (6 config tests)
- New test in `test_orchestrator.py`: `test_parse_findings_handles_windows_paths_in_json`

---

## [2.1.0-beta.1] - 2026-03-19

### Added

#### Automated Improvement Loop
- **`improvement_loop/` package**: Programmatic audit → fix → test → commit cycle with Claude API integration
  - **`orchestrator_v1.py`**: Main loop driver with configurable max iterations and dry-run mode
  - **`evaluator.py`**: Finding schema (Pydantic) and audit quality scoring with exit condition logic
  - **`loop_tracker.py`**: Iteration logging, context generation for subsequent iterations, score drift detection
  - **`git_utils.py`**: Subprocess-based git operations (branch management, test runners, commit helpers)
- **Self-healing protocol**: When a fix causes test regressions, the orchestrator automatically attempts to fix the regression up to 2 times before reverting (`MAX_SELF_HEAL_ATTEMPTS`)
- **Run summary**: End-of-loop summary with per-iteration stats, dimension breakdown, and convergence status

#### New Report Sections
- **`analysis/report/sections/analysis_cross_dwi.py`**: Cross-DWI comparison section builder
- **`analysis/report/sections/analysis_features.py`**: Feature overlap analysis section builder

#### New Test Files (4 Python)
- `test_evaluator_finding.py`, `test_git_utils.py`, `test_loop_tracker.py`, `test_orchestrator.py` for the improvement loop package
- Python test suite: 32 → 36 files (~1559 tests)

### Changed

#### Correctness Fixes (Improvement Loop — 5 iterations)
- **`metrics_survival.m`**: Replaced invalid MATLAB ternary syntax (`? :`) with proper if/else blocks
- **`compute_calibration_metrics.m`**: Replaced 3 invalid MATLAB ternary expressions with if/else blocks
- **`fit_time_varying_cox.m`**: Replaced 2 invalid MATLAB ternaries; guard empty finite values before median; fixed entry weight denominator; changed log(0) minimum from 0.5 to 1 day
- **`extract_tumor_core.m`**: NaN/zero `max_adc` guard in active contours; small-sample fallback in region growing
- **`detect_motion_artifacts.m`**: Zero noise floor fallback to 1% of masked signal median; NaN NMI flagging
- **`bootstrap_ci.m`**: Consolidated BCa fallback; clamped quantile indices; removed redundant degenerate check
- **`apply_external_validation.m`**: Center zero-variance features; warn on padded missing features
- **`build_td_panel.m`**: Warning when scan_days exceeds available timepoints
- **`metrics_stats_comparisons.m`**: Guard before ranksum when fewer than 2 groups
- **`metrics_longitudinal.m`**: Fixed Octave SEM computation for N=1 case
- **`compute_texture_features.m`**: Constant-image guard for 3D GLRLM; existence check before `compute_glrlm_3d`
- **`make_grouped_folds.m`**: Warning on k≤1 degenerate fallback
- **`imputation_sensitivity.m`**: Case-insensitive patient ID matching (`strcmpi`)
- **`compute_nri.m`**: Minimum 3 events/non-events threshold; `isfinite` guards on variance and z-score
- **`fit_models.m`**: `rcond` check to skip warm start on singular matrices
- **`process_single_scan.m`**: Check `mkdir` return code with error on failure
- **`decision_curve_analysis.m`**: Floating-point safety for threshold comparison
- **`compute_ipcw_weights.m`**: Guard against division by zero risk sum

#### Analysis Pipeline
- **`batch_graph_analysis.py`**: Moved `queue.task_done()` into `finally:` block to prevent asyncio deadlock on interruption
- **`cross_dwi_agreement.py`**: Fisher z-transform requires n>3; replaced manual exp formula with `math.tanh()` for overflow safety
- **`run_analysis.py`**: `TeeWriter.fileno()` wrapped in try/except for captured streams (pytest compatibility)

#### Documentation
- **`CLAUDE_WORKFLOWS.md`**: Added self-healing protocol documentation to the Validation Protocol section
- **`CLAUDE.md`**, **`CLAUDE_REFERENCE.md`**: Updated file counts and module tables

### Fixed
- 5 instances of invalid MATLAB ternary syntax (`? :`) that would crash at runtime
- Asyncio queue deadlock when `batch_graph_analysis.py` is interrupted via KeyboardInterrupt
- CCC confidence interval overflow when sample size ≤ 3
- IPCW weight computation NaN propagation from exp() underflow

---

## [2.1.0-alpha.2] - 2026-03-17

### Added

#### Pipeline Infrastructure
- **`suppress_core_warnings.m`** (`pipeline/utils/`): Extracted duplicate warning suppression into a reusable utility for cleaner core module code
- **`CLAUDE_WORKFLOWS.md`**: Extracted workflow and process documentation from `CLAUDE.md` into a dedicated file for better separation of concerns

#### Input Validation & Error Handling
- **`extract_tumor_core.m`**: Input validation for mask dimensions and method parameters
- **`filter_collinear_features.m`**: Input validation for feature matrix and threshold arguments
- **`load_auxiliary_biomarkers.m`**: Column dimension validation and data loss warnings
- **`fit_models.m`**: Early b-value validation moved to function entry for fail-fast behavior
- **`load_dwi_data.m`**: Guard against empty DICOM arrays
- **`compute_texture_features.m`**: GLCM error logging with message IDs
- **`discover_patient_files.m`**: Error instead of warning when clinical spreadsheet has no patient column

### Changed
- **`compute_summary_metrics.m`**: Cache all GTV masks using `containers.Map` for improved performance on large cohorts
- **`scale_td_panel.m`**: Added `omitnan` flag to baseline scaling path for consistency with other scaling operations
- **`build_td_panel.m`**: Clarified NaN propagation behavior in decay imputation documentation
- **`compute_schoenfeld_residuals.m`**: Improved LOWESS smoother numerical stability
- **`assemble_predictive_features.m`**, **`compute_percent_deltas.m`**: Pre-allocated arrays to avoid repeated memory reallocation in loops
- **`dispatch_pipeline_steps.m`**: Simplified by extracting warning suppression to `suppress_core_warnings.m`
- **`batch_graph_analysis.py`**: Re-raise `KeyboardInterrupt` and `CancelledError` in batch worker instead of silently swallowing; added error logging to silent exceptions
- **CLAUDE.md**: Split workflow documentation into `CLAUDE_WORKFLOWS.md`

### Fixed
- **`dispatch_pipeline_steps.m`**: Fixed critical bugs — undefined variable reference and empty struct initialization
- **`compare_core_methods.m`**: Fixed incorrect variable name `nTp` in bootstrap CI loop
- **`detect_motion_artifacts.m`**: Flag motion artifacts when NMI is `NaN` due to insufficient voxels
- **`process_single_scan.m`**, **`metrics_dosimetry.m`**, **`metrics_longitudinal.m`**: Fixed `fopen` return value checks and division-by-zero in standard error computation
- **`batch_graph_analysis.py`**: Fixed regex deduplication overlap in error handling
- **`test_convert_dicom.m`**: Escaped shell arguments in `chmod` calls for path safety

---

## [2.1.0-alpha.1] - 2026-03-17

### Added

#### Advanced Modeling Framework
- **Imputation sensitivity analysis** (`imputation_sensitivity.m`): Compare KNN imputation against LOCF, Mean, and Linear Interpolation alternatives; controlled via `run_imputation_sensitivity` config field
- **Time-varying Cox models** (`fit_time_varying_cox.m`): Stratified and extended Cox models for when proportional hazards assumption is violated; controlled via `fit_time_varying_cox` config field
- **Decision curve analysis** (`decision_curve_analysis.m`): Net benefit calculation and treat-all/none comparison for clinical utility assessment
- **Net reclassification improvement** (`compute_nri.m`): NRI, continuous NRI, and IDI for comparing predictive models
- **External validation** (`prepare_external_validation.m`, `apply_external_validation.m`): Export trained models and apply to external datasets; controlled via `export_validation_model` and `external_validation_data` config fields
- **Auxiliary biomarker integration** (`load_auxiliary_biomarkers.m`): Load non-DWI biomarker data from CSV for multi-modal predictive modeling; controlled via `auxiliary_biomarker_csv` and `use_auxiliary_biomarkers` config fields

#### Survival & Predictive Enhancements
- **Schoenfeld residuals** (`compute_schoenfeld_residuals.m`): Scaled Schoenfeld residuals and PH assumption testing via Spearman correlation with diagnostic figures
- **Fine-Gray subdistribution hazard model**: Proper subdistribution weights replacing ad-hoc 0.5 multiplier
- **Calibration metrics** (`compute_calibration_metrics.m`): Brier score, Hosmer-Lemeshow test, calibration slope/intercept
- **Bootstrap confidence intervals** (`bootstrap_ci.m`): BCa bootstrap CIs for arbitrary scalar metric functions; vectorized resampling with parfor support
- **Forest plot section** (`analysis/report/sections/forest_plot.py`): Hazard ratio extraction and matplotlib forest plot generation for the HTML report

#### Radiomics & Image Quality
- **Texture features** (`compute_texture_features.m`): First-order, GLCM, GLRLM (3D with 13 directions), shape, and uniformity features (24 total); IBSI-compliant quantization; controlled via `use_texture_features` and `texture_quantization_method` config fields
- **Registration quality metrics** (`compute_registration_quality.m`): Jacobian determinant, NCC, mutual information; configurable `voxel_spacing` parameter
- **Motion artifact detection** (`detect_motion_artifacts.m`): DWI volume quality assessment — CV, NMI, signal dropout detection; controlled via `exclude_motion_volumes` config field

#### GPU Acceleration
- **GPU-accelerated fitting** (`gpu_available.m`): Offload ADC WLS fitting and DnCNN inference to CUDA GPUs via `gpuArray`; graceful CPU fallback; controlled via `use_gpu` and `gpu_device` config fields
- **GPU memory safety**: Automatic memory checks and fallback in `fit_models.m` when GPU memory is insufficient

#### Docker Support
- **Dockerfile**: Multi-stage Docker build with MATLAB Runtime and Python environment
- **docker-compose.yml**: Pipeline and analysis service definitions
- **`docker/entrypoint.sh`**: Container entrypoint with mode dispatch, validation, and `--dry-run` flag
- **`docs/DOCKER.md`**: Comprehensive Docker usage guide
- Configurable MCR version build arg (`.matlab_version`) with runtime version check
- MATLAB/Python pre-flight checks in Docker entrypoint
- Relaxed apt-get version pins to prevent stale build failures

#### Vision Analysis
- **Claude API provider** (`analysis/parsers/batch_graph_analysis.py`): Anthropic Claude vision analysis alongside Gemini; `--provider gemini|claude|both` for single or dual-provider comparison with per-image difference CSV
- **Local Gemini fallback**: Automatic retry with local Gemini when API calls fail
- **Cross-DWI agreement analysis** (`analysis/cross_reference/cross_dwi_agreement.py`): Bland-Altman, Lin's CCC, and ICC agreement metrics between DWI types

#### Interactive Reporting
- **Interactive HTML report** (`analysis/report/generate_interactive_report.py`): Client-side filtering, Chart.js visualizations, patient drill-down, sortable tables, and DWI/core-method comparison
- **Interactive report constants** (`analysis/report/interactive_constants.py`): CSS and JavaScript for sidebar, tabs, chart rendering, and filter logic

#### Pipeline Infrastructure
- **`prepare_pipeline_session.m`**: Pipeline session initialization with try-catch error handling
- **`dispatch_load_and_sanity.m`**, **`dispatch_pipeline_steps.m`**: Extracted pipeline step dispatch logic
- **`compute_percent_deltas.m`**: Treatment-induced percent/absolute changes from baseline (extracted from metrics_baseline)
- **`compute_histogram_laplace.m`**: Laplace-smoothed histogram probability distribution
- **`compute_kurt_skew.m`**: Kurtosis/skewness computation with minimum sample guard
- **`detect_baseline_outliers.m`**: Outcome-blinded IQR outlier detection for baseline metrics
- **`nanmean_safe.m`**, **`nanstd_safe.m`**: Octave-compatible NaN-ignoring mean/std

#### New Test Files (14 MATLAB, 10 Python)
- MATLAB: `test_gpu_available.m`, `test_imputation_sensitivity.m`, `test_initialize_pipeline.m`, `test_load_auxiliary_biomarkers.m`, `test_normalize_patient_ids.m`, `test_octave_shims.m`, `test_parfor_progress.m`, `test_plot_predictive_diagnostics.m`, `test_prepare_external_validation.m`, `test_prepare_pipeline_session.m`, `test_progress_gui.m`, `test_run_elastic_net_cv.m`, `test_run_loocv_risk_scores.m`, `test_schoenfeld_residuals.m`, `test_select_dwi_vectors.m`, `test_time_varying_cox.m`, `test_write_sentinel_file.m`
- Python: `test_cross_dwi_agreement.py`, `test_forest_plot.py`, `test_generate_report_helpers.py`, `test_generate_report_integration.py`, `test_generate_report_manuscript.py`, `test_generate_report_sections.py`, `test_integration.py`, `test_parse_imputation_and_tv_cox.py`, `test_report_sections_robustness.py`, `test_interactive_report.py`
- MATLAB test suite: 106 → 120 files; Python test suite: 22 → 32 files (~1482 tests)

### Changed
- **`metrics_baseline.m`**: Refactored to return a single struct instead of 29 positional outputs
- **Pipeline step functions**: Refactored to use `baseline_results` and `session` structs for cleaner parameter passing
- **`run_dwi_pipeline.m`**: Further reduced from 60KB to 7KB by extracting orchestrator logic into `dispatch_load_and_sanity.m`, `dispatch_pipeline_steps.m`, and `prepare_pipeline_session.m`
- **`bootstrap_ci.m`**: Optimized with vectorized resampling, pre-generated indices, and parfor support
- **NMI computation**: Replaced hand-rolled NMI with `histcounts2` and added constant-signal guard
- **`imputation_sensitivity.m`**: Refactored `evaluate_imputed` to delegate to `run_elastic_net_cv` and `run_loocv_risk_scores`
- **Report section submodules**: Further split into finer-grained submodules for analysis, statistics, data, and discussion sections
- **`analysis/__init__.py`**: Updated to import directly from sub-modules instead of wrappers
- **Per-method dosimetry**: Extracted from disk loading when loading dosimetry results
- **UTF-8 handling**: Replaced raw UTF-8 byte escape sequences with native Unicode emoji characters
- **Python dependencies**: Added `numpy`, `matplotlib`; regenerated `requirements-lock.txt`; added upper version bounds
- **CLAUDE.md**: Split into two files (`CLAUDE.md` + `CLAUDE_REFERENCE.md`) to reduce context window usage
- **Changelog**: Archived pre-v2.0.0 entries to `CHANGELOG_ARCHIVE.md`

### Fixed
- **Fine-Gray weights**: Replaced ad-hoc 0.5 multiplier with proper subdistribution weights for competing risk modeling
- **Dockerfile**: Relaxed apt-get version pins that caused stale build failures

---

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
