#!/usr/bin/env bash
# =============================================================================
# pancData3 Docker entrypoint
# Usage: entrypoint.sh [--dry-run] [pipeline|analysis|all]
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Parse flags
# ---------------------------------------------------------------------------
DRY_RUN=0
while [[ "${1:-}" == --* ]]; do
    case "$1" in
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        *)
            echo "❌ Unknown flag: $1"
            exit 1
            ;;
    esac
done

MODE="${1:-all}"

# ---------------------------------------------------------------------------
# MCR version check
# ---------------------------------------------------------------------------
check_mcr_version() {
    local installed_version="${MCR_VERSION:-unknown}"
    echo "💡 Installed MATLAB Runtime version: ${installed_version}"

    local version_file="/opt/pancData3/.matlab_version"
    if [ -f "${version_file}" ]; then
        local expected_version
        expected_version="$(tr -d '[:space:]' < "${version_file}")"
        if [ "${installed_version}" != "${expected_version}" ]; then
            echo "⚠️  WARNING: MCR version mismatch!"
            echo "   Installed: ${installed_version}"
            echo "   Expected (from .matlab_version): ${expected_version}"
            echo "   The MCR version must match the MATLAB version used to compile the pipeline binary."
            echo "   Rebuild with: docker build --build-arg MCR_VERSION=${expected_version} -t pancdata3:latest ."
        else
            echo "✅ MCR version matches .matlab_version"
        fi
    else
        echo "💡 No .matlab_version file found — skipping version verification"
    fi
}

check_mcr_version

# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------
validate_config() {
    if [ ! -f /opt/pancData3/config.json ]; then
        echo "❌ ERROR: config.json not found at /opt/pancData3/config.json"
        echo "   Mount your config file: -v /path/to/config.json:/opt/pancData3/config.json:ro"
        exit 1
    fi
    echo "✅ config.json found"
}

validate_data() {
    if [ ! -d /opt/pancData3/data ]; then
        echo "❌ ERROR: /opt/pancData3/data directory not found"
        echo "   Mount your data directory: -v /path/to/patient_data:/opt/pancData3/data:ro"
        exit 1
    fi
    if [ -z "$(ls -A /opt/pancData3/data 2>/dev/null)" ]; then
        echo "❌ ERROR: /opt/pancData3/data is empty"
        echo "   Mount a directory containing patient DWI files"
        exit 1
    fi
    echo "✅ /data volume mounted and non-empty"
}

validate_output() {
    if [ ! -d /opt/pancData3/output ]; then
        mkdir -p /opt/pancData3/output
    fi
    echo "✅ /output volume ready"
}

preflight_matlab() {
    echo "⚙️  Running MATLAB Runtime pre-flight check..."
    if command -v matlab &>/dev/null; then
        if matlab -batch "disp('OK')" &>/dev/null; then
            echo "✅ MATLAB pre-flight passed"
        else
            echo "❌ ERROR: MATLAB pre-flight failed (matlab -batch returned non-zero)"
            echo "   Possible causes:"
            echo "   - MATLAB Runtime license is missing or expired"
            echo "   - MATLAB installation is corrupted"
            echo "   - Required toolboxes are not installed"
            echo "   Check your MATLAB Runtime installation and license configuration."
            exit 1
        fi
    elif [ -f /opt/pancData3/pipeline/run_execute_all_workflows.sh ]; then
        if [ -z "${MCR_ROOT:-}" ]; then
            echo "❌ ERROR: MCR_ROOT environment variable is not set"
            echo "   Set MCR_ROOT to the MATLAB Runtime installation path"
            echo "   (e.g., /opt/matlabruntime/R2024a)"
            exit 1
        fi
        if [ ! -d "${MCR_ROOT}" ]; then
            echo "❌ ERROR: MCR_ROOT directory does not exist: ${MCR_ROOT}"
            echo "   Verify the MATLAB Runtime is installed at the specified path"
            exit 1
        fi
        # Test the compiled binary with a quick invocation
        if /opt/pancData3/pipeline/run_execute_all_workflows.sh "${MCR_ROOT}" --help &>/dev/null; then
            echo "✅ Compiled MATLAB binary pre-flight passed"
        else
            echo "⚠️  Compiled binary --help check returned non-zero (may be expected)"
            echo "   MCR_ROOT=${MCR_ROOT} exists and compiled binary is present"
            echo "✅ MATLAB Runtime pre-flight passed (basic checks)"
        fi
    else
        echo "❌ ERROR: No MATLAB or compiled pipeline binary found"
        echo "   Ensure the Docker image includes MATLAB Runtime or a compiled binary"
        echo "   Expected one of:"
        echo "   - 'matlab' command on PATH"
        echo "   - /opt/pancData3/pipeline/run_execute_all_workflows.sh + MCR_ROOT"
        exit 1
    fi
}

preflight_python() {
    echo "⚙️  Running Python pre-flight check..."
    if ! command -v python3.12 &>/dev/null; then
        echo "❌ ERROR: python3.12 not found on PATH"
        exit 1
    fi
    if ! python3.12 -c "import sys; print(f'Python {sys.version}')" &>/dev/null; then
        echo "❌ ERROR: python3.12 failed to execute"
        exit 1
    fi
    echo "✅ Python pre-flight passed"
}

# ---------------------------------------------------------------------------
# Run modes
# ---------------------------------------------------------------------------
run_pipeline() {
    echo "🚀 Starting MATLAB pipeline (execute_all_workflows)..."
    validate_config
    validate_data
    validate_output
    preflight_matlab

    if [ "${DRY_RUN}" -eq 1 ]; then
        echo "💡 --dry-run: pre-flight checks passed, skipping pipeline execution"
        return 0
    fi

    cd /opt/pancData3
    # Run via MATLAB Runtime compiled application
    # The MATLAB Runtime expects the compiled binary; fall back to matlab -batch
    if command -v matlab &>/dev/null; then
        matlab -batch "addpath('pipeline/core','pipeline/utils','pipeline/dependencies'); cd pipeline; execute_all_workflows"
    elif [ -f /opt/pancData3/pipeline/run_execute_all_workflows.sh ]; then
        /opt/pancData3/pipeline/run_execute_all_workflows.sh "${MCR_ROOT}"
    else
        echo "❌ ERROR: No MATLAB or compiled pipeline binary found"
        echo "   Ensure the MATLAB Runtime image is correctly configured"
        exit 1
    fi
    echo "✅ Pipeline completed"
}

run_analysis() {
    echo "🚀 Starting Python analysis suite..."
    validate_config
    validate_output
    preflight_python

    if [ "${DRY_RUN}" -eq 1 ]; then
        echo "💡 --dry-run: pre-flight checks passed, skipping analysis execution"
        return 0
    fi

    cd /opt/pancData3
    python3.12 analysis/run_analysis.py \
        --folder /opt/pancData3/output \
        --skip-vision \
        --skip-checks \
        "$@"
    echo "✅ Analysis completed"
}

# ---------------------------------------------------------------------------
# Main dispatch
# ---------------------------------------------------------------------------
case "${MODE}" in
    pipeline)
        run_pipeline
        ;;
    analysis)
        shift  # Remove "analysis" from args, pass rest to run_analysis
        run_analysis "$@"
        ;;
    all)
        run_pipeline
        shift 2>/dev/null || true
        run_analysis "$@"
        ;;
    *)
        echo "Usage: entrypoint.sh [--dry-run] [pipeline|analysis|all]"
        echo ""
        echo "Modes:"
        echo "  pipeline  — Run the MATLAB DWI analysis pipeline"
        echo "  analysis  — Run the Python post-hoc analysis suite"
        echo "  all       — Run pipeline followed by analysis (default)"
        echo ""
        echo "Flags:"
        echo "  --dry-run — Run only pre-flight checks (config, data, MATLAB/Python)"
        echo "              without executing the pipeline or analysis"
        exit 1
        ;;
esac
