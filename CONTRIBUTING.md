# Contributing to pancData3

Thank you for your interest in contributing to pancData3! This document provides guidelines and information for contributors.

## Table of Contents

- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Making Changes](#making-changes)
- [Testing](#testing)
- [Code Style](#code-style)
- [Pull Request Process](#pull-request-process)
- [Important Constraints](#important-constraints)

---

## Getting Started

1. Fork the repository
2. Clone your fork locally
3. Create a feature branch from `main`
4. Make your changes
5. Run the test suite
6. Submit a pull request

## Development Setup

### Prerequisites

- MATLAB R2021a or later
- Statistics and Machine Learning Toolbox
- Image Processing Toolbox
- [dcm2niix](https://github.com/rordenlab/dcm2niix) for DICOM conversion

### Setup

```bash
git clone https://github.com/<your-username>/pancData3.git
cd pancData3
cp config.example.json config.json
# Edit config.json with your local paths
```

In MATLAB:

```matlab
addpath('core', 'utils', 'dependencies');
```

## Making Changes

### Branch Naming

- Feature branches: `feature/<description>`
- Bug fixes: `fix/<description>`
- Documentation: `docs/<description>`

### Architecture Guidelines

- **Orchestrator pattern**: All pipeline steps are modular and independently callable via `run_dwi_pipeline.m`.
- **Explicit parameter passing**: No global variables or shared workspace state between pipeline steps.
- **Checkpoint-safe**: Preserve checkpointing logic in `load_dwi_data.m`; it is critical for large cohort recovery.
- **Security-first**: Use `safe_load_mask` for `.mat` file loading and `escape_shell_arg` for all `system()` calls.

## Testing

Run the full test suite before submitting any changes:

```matlab
run('tests/run_all_tests.m')
```

### Test Requirements

- All existing tests must pass
- New functionality should include corresponding tests
- Tests go in the `tests/` directory
- Use MATLAB's `unittest` framework (class-based tests preferred)
- Test file names must start with `test_`

### Test Categories

| Location | Purpose |
|---|---|
| `tests/` | Functional and integration tests |
| `tests/benchmarks/` | Performance benchmarks |
| `tests/diagnostics/` | Diagnostic spot-check scripts |

### Octave Compatibility

The pipeline supports GNU Octave via shims in `.octave_compat/`. When adding new functionality:

- Avoid MATLAB-only syntax not covered by the existing shim layer
- Run `tests/run_all_tests_octave.m` to verify Octave compatibility if available
- Compatibility shims are in `.octave_compat/` -- add new shims there if needed

## Code Style

### General

- Use descriptive variable and function names
- Functions should have a help block (H1 line + description)
- Progress output uses emoji prefixes: `🚀` start, `⚙️` processing, `✅` success, `❌` failure, `📁` file output, `💡` info

### Data Leakage Prevention

This is a research pipeline with strict leakage prevention requirements. When modifying cross-validation, imputation, or scaling code:

- Patient-stratified folds: no intra-patient leakage across CV splits
- KNN imputation: infinite distances for same-patient rows
- Scaling: strictly timepoint-specific
- DL provenance: guard against training/analysis overlap

## Pull Request Process

1. Ensure all tests pass (`run('tests/run_all_tests.m')`)
2. Update documentation if adding new features or configuration fields
3. Keep PRs focused -- one feature or fix per PR
4. Write a clear PR description explaining the "why" behind changes
5. Request review from a maintainer

### Review Criteria

- Code correctness and test coverage
- No data leakage introduced
- No modification of files in `dependencies/`
- No hardcoded file paths (uses `config.json`)
- No unsanitized strings passed to `system()`
- Octave compatibility preserved (no MATLAB-only syntax without shims in `.octave_compat/`)

## Python Dependencies

When modifying `analysis/requirements.txt`, you **must** regenerate the lock file to keep pinned versions in sync:

```bash
python3.12 -m venv /tmp/venv
/tmp/venv/bin/pip install -r analysis/requirements.txt
/tmp/venv/bin/pip freeze > analysis/requirements-lock.txt
rm -rf /tmp/venv
```

Commit both `requirements.txt` and `requirements-lock.txt` together. The lock file ensures reproducible installs across environments.

---

## Important Constraints

> **Do not modify files in the `dependencies/` folder.** These are third-party scripts maintained under their own licenses.

> **Never commit patient data, clinical CSVs, or PHI.** The `.gitignore` is configured to prevent this, but always verify before committing.

> **Do not introduce global variables** or persistent workspace state between pipeline steps.

---

## Questions?

Open an issue on the [GitHub Issues](https://github.com/akarlin3/pancData3/issues) page for questions or discussion.
