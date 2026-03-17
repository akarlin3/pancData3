# Docker Usage Guide

This guide covers running the pancData3 pipeline and analysis suite in Docker containers for reproducible execution.

---

## Prerequisites

- [Docker](https://docs.docker.com/get-docker/) 20.10+
- [Docker Compose](https://docs.docker.com/compose/install/) v2.0+ (included with Docker Desktop)
- A valid `config.json` (copy from `config.example.json` and configure)
- Patient DWI data directory

---

## Building the Image

From the repository root:

```bash
docker build -t pancdata3:latest .
```

The multi-stage build:
1. **Stage 1** — Compiles `dcm2niix` from source in a slim Debian image
2. **Stage 2** — Copies the compiled binary into the MATLAB Runtime image, installs Python 3.12 and all analysis dependencies

Build time is approximately 10–15 minutes on the first run (cached thereafter).

---

## Preparing Your Configuration

Copy the example config and edit it for Docker paths:

```bash
cp config.example.json config.json
```

Inside the container, patient data is mounted at `/opt/pancData3/data` and output goes to `/opt/pancData3/output`. Set your `config.json` accordingly:

```json
{
  "dataloc": "/opt/pancData3/data/",
  "dcm2nii_call": "dcm2niix",
  "skip_to_reload": false,
  "dwi_type": "Standard"
}
```

> **Important:** `dcm2nii_call` should be set to `"dcm2niix"` (it is on the container PATH). `dataloc` must point to `/opt/pancData3/data/` where your host data directory is mounted.

---

## Running with Docker Compose

### Full Pipeline + Analysis

```bash
# Set required environment variables
export DATA_DIR=/path/to/patient_dwi_data
export OUTPUT_DIR=/path/to/output
export CONFIG_FILE=/path/to/config.json

# Run pipeline followed by analysis
docker compose up
```

### Pipeline Only

```bash
docker compose up pipeline
```

### Analysis Only

If you already have pipeline output and only want to run the Python analysis:

```bash
docker compose up analysis
```

> **Note:** By default, the `analysis` service depends on `pipeline` completing successfully. To run analysis independently, use `docker run` directly (see below).

---

## Running with Docker Run

### Full Pipeline (All DWI Types)

```bash
docker run --rm \
  -v /path/to/patient_data:/opt/pancData3/data:ro \
  -v /path/to/output:/opt/pancData3/output \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  pancdata3:latest pipeline
```

### Analysis Only

```bash
docker run --rm \
  -v /path/to/output:/opt/pancData3/output \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  pancdata3:latest analysis
```

### Full Pipeline + Analysis

```bash
docker run --rm \
  -v /path/to/patient_data:/opt/pancData3/data:ro \
  -v /path/to/output:/opt/pancData3/output \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  pancdata3:latest all
```

### Analysis with Vision (Gemini API)

To enable vision-based graph analysis, pass the `GEMINI_API_KEY` environment variable and omit the `--skip-vision` flag:

```bash
docker run --rm \
  -e GEMINI_API_KEY=your_api_key_here \
  -v /path/to/output:/opt/pancData3/output \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  pancdata3:latest analysis --skip-checks
```

---

## Volume Mounts

| Host Path | Container Path | Mode | Purpose |
|---|---|---|---|
| Patient DWI data directory | `/opt/pancData3/data` | `ro` (read-only) | Input DICOM/NIfTI patient files |
| Output directory | `/opt/pancData3/output` | `rw` (read-write) | Pipeline results, figures, logs |
| `config.json` | `/opt/pancData3/config.json` | `ro` (read-only) | Pipeline configuration |

---

## Environment Variables

| Variable | Required | Default | Description |
|---|---|---|---|
| `DATA_DIR` | Yes (compose) | — | Host path to patient DWI data directory |
| `OUTPUT_DIR` | No (compose) | `./output` | Host path for pipeline output |
| `CONFIG_FILE` | No (compose) | `./config.json` | Host path to `config.json` |
| `GEMINI_API_KEY` | No | — | Google Gemini API key for vision-based graph analysis |
| `MCR_ROOT` | No | `/opt/matlabruntime/${MCR_VERSION}` | MATLAB Runtime installation path (set automatically) |
| `MCR_VERSION` | No | `r2024a` | MATLAB Runtime version (set from build arg; see [MATLAB Runtime Version](#matlab-runtime-version)) |

---

## Entrypoint Modes

The container entrypoint (`docker/entrypoint.sh`) accepts these modes:

| Mode | Command | Description |
|---|---|---|
| `pipeline` | `docker run ... pancdata3 pipeline` | Run the MATLAB DWI pipeline only |
| `analysis` | `docker run ... pancdata3 analysis` | Run the Python analysis suite only |
| `all` | `docker run ... pancdata3 all` | Run pipeline, then analysis (default) |

### Flags

| Flag | Command | Description |
|---|---|---|
| `--dry-run` | `docker run ... pancdata3 --dry-run pipeline` | Run only pre-flight checks without executing the pipeline or analysis |

The `--dry-run` flag must appear **before** the mode argument. It validates the environment (config, data volumes, MATLAB Runtime, Python) and exits without running the actual pipeline or analysis. This is useful for verifying that a Docker setup is correct before committing to a full run.

### Pre-flight Checks

The entrypoint validates that:
1. `config.json` exists at `/opt/pancData3/config.json`
2. `/opt/pancData3/data` is mounted and non-empty (pipeline mode only)
3. `/opt/pancData3/output` directory exists (creates it if missing)
4. **MATLAB Runtime** is functional — runs `matlab -batch "disp('OK')"` (or verifies `MCR_ROOT` and the compiled binary exist)
5. **Python 3.12** is available and executable

If a pre-flight check fails, a clear error message is printed with troubleshooting guidance.

---

## MATLAB Runtime Version

The Docker image includes a specific version of the MATLAB Compiler Runtime (MCR). **The MCR version must exactly match the MATLAB version used to compile the pipeline binary.** A version mismatch will cause runtime errors when executing the compiled pipeline.

### Default version

The default MCR version is `r2024a`, defined by the `MCR_VERSION` build argument in the `Dockerfile`. The expected version is also recorded in the `.matlab_version` file at the repository root.

### Changing the MCR version

If the pipeline binary was compiled with a different MATLAB version, override the build argument:

```bash
docker build --build-arg MCR_VERSION=r2023b -t pancdata3:latest .
```

Then update `.matlab_version` in the repository root to match:

```bash
echo "r2023b" > .matlab_version
```

Available MCR versions correspond to tags on the `mathworks/matlab-runtime` Docker Hub image (e.g., `r2023a`, `r2023b`, `r2024a`).

### Runtime version check

On startup, the entrypoint script prints the installed MCR version and compares it against the `.matlab_version` file baked into the image. If there is a mismatch, a warning is displayed with the correct `docker build` command to fix it. This catches cases where the image was built with one MCR version but the repository expects another (e.g., after upgrading the pipeline compiler).

---

## Data Safety

The Docker image **never** contains:
- Patient data, clinical spreadsheets, or any PHI
- `config.json` (must be mounted at runtime)
- `.csv`, `.xlsx`, `.nii`, `.dcm`, or `.h5` files

Patient data volumes should always be mounted read-only (`:ro`) to prevent accidental modification.

---

## Troubleshooting

### "config.json not found"

Mount your config file:
```bash
-v /path/to/config.json:/opt/pancData3/config.json:ro
```

### "/data is empty"

Ensure your data directory is correctly mounted and contains patient folders:
```bash
-v /path/to/patient_data:/opt/pancData3/data:ro
```

### MATLAB Runtime pre-flight fails

Run with `--dry-run` to isolate the issue:
```bash
docker run --rm \
  -v /path/to/config.json:/opt/pancData3/config.json:ro \
  -v /path/to/patient_data:/opt/pancData3/data:ro \
  pancdata3:latest --dry-run pipeline
```

Common causes:
- **License missing/expired**: Ensure valid MATLAB license is configured
- **MCR_ROOT not set**: The `MCR_ROOT` env var must point to the MATLAB Runtime installation
- **Corrupted installation**: Rebuild the Docker image from scratch

### MATLAB Runtime errors

The container uses MATLAB Runtime at the version specified by the `MCR_VERSION` build argument (default: `r2024a`). Ensure your compiled pipeline is compatible with this version — see [MATLAB Runtime Version](#matlab-runtime-version). If running uncompiled `.m` files, a full MATLAB installation is required (not included in the Runtime image).

### WeasyPrint PDF generation fails

WeasyPrint system dependencies (Cairo, Pango, GDK-PixBuf) are pre-installed in the image. If PDF generation still fails, try running with `--no-pdf` to generate HTML reports only.
