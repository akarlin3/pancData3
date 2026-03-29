# Security Policy

## Reporting a Vulnerability

If you discover a security vulnerability in this project, please report it responsibly.

**Do not open a public issue for security vulnerabilities.**

Instead, please email the maintainer directly:

- **Contact**: Open a private security advisory via [GitHub Security Advisories](https://github.com/akarlin3/pancData3/security/advisories/new)

### What to Include

- A description of the vulnerability
- Steps to reproduce the issue
- Potential impact assessment
- Suggested fix (if any)

### Response Timeline

- **Acknowledgment**: Within 48 hours
- **Assessment**: Within 1 week
- **Fix**: As soon as reasonably possible, depending on severity

## Security Considerations

This project processes medical imaging data. The following security measures are built into the codebase:

### Protected Health Information (PHI)

- Patient data, clinical CSVs, and PHI must **never** be committed to the repository
- `.gitignore` is configured to exclude common PHI-containing file types (`.xlsx`, `.csv`, `.dcm`, `.nii`)
- Cloud-based AI agents are restricted from accessing patient data

### Safe File Loading

- `safe_load_mask.m` validates variable types before loading `.mat` files, preventing arbitrary code execution from maliciously crafted files

### Shell Injection Prevention

- `escape_shell_arg.m` must be used for all paths passed to `system()` calls
- Unsanitized user strings must never be passed to shell commands

### Data Integrity

- Cross-validation uses patient-stratified folds to prevent data leakage
- KNN imputation enforces temporal leakage bounds
- DL provenance tracking prevents training/analysis overlap

## Supported Versions

| Version | Supported |
|---|---|
| 2.3.1 (latest) | Yes |
| 2.2.0 | Yes |
| 2.1.x and earlier | No |
