#!/usr/bin/env bash
# =============================================================================
# pancData3 Docker entrypoint
# Usage: entrypoint.sh [pipeline|analysis|all]
# =============================================================================
set -euo pipefail

MODE="${1:-all}"

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

# ---------------------------------------------------------------------------
# Run modes
# ---------------------------------------------------------------------------
run_pipeline() {
    echo "🚀 Starting MATLAB pipeline (execute_all_workflows)..."
    validate_config
    validate_data
    validate_output

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
        echo "Usage: entrypoint.sh [pipeline|analysis|all]"
        echo ""
        echo "Modes:"
        echo "  pipeline  — Run the MATLAB DWI analysis pipeline"
        echo "  analysis  — Run the Python post-hoc analysis suite"
        echo "  all       — Run pipeline followed by analysis (default)"
        exit 1
        ;;
esac
