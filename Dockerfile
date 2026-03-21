# =============================================================================
# pancData3 — Multi-stage Docker build
# Stage 1: Build dcm2niix from source
# Stage 2: Runtime with MATLAB Runtime + Python 3.12
# =============================================================================

# ---------------------------------------------------------------------------
# Stage 1 — Build dcm2niix from latest GitHub release
# ---------------------------------------------------------------------------
FROM debian:bookworm-20240211-slim AS builder

# NOTE: No exact version pins on system packages — exact pins cause apt-get
# failures when the base image updates its package index (the pinned version
# is removed and only the newer version is available).
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        ca-certificates \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone and build dcm2niix (pinned to stable release tag for reproducibility)
RUN git clone --branch v1.0.20250506 --depth 1 \
        https://github.com/rordenlab/dcm2niix.git /tmp/dcm2niix \
    && mkdir /tmp/dcm2niix/build \
    && cd /tmp/dcm2niix/build \
    && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/dcm2niix .. \
    && make -j"$(nproc)" \
    && make install

# ---------------------------------------------------------------------------
# Stage 2 — Runtime image with MATLAB Runtime and Python 3.12
# ---------------------------------------------------------------------------
ARG MCR_VERSION=r2024a
FROM mathworks/matlab-runtime:${MCR_VERSION}

LABEL maintainer="Avery Karlin <akarlin3>" \
      description="pancData3: Pancreatic DWI Analysis Pipeline" \
      license="AGPL-3.0"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install Python 3.12, WeasyPrint system dependencies, and runtime libraries.
# NOTE: No exact version pins — exact pins cause apt-get failures when the
# base image updates its package index (the pinned version is removed and
# only the newer version is available).
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3.12 \
        python3.12-venv \
        python3-pip \
        libcairo2 \
        libpango-1.0-0 \
        libpangocairo-1.0-0 \
        libgdk-pixbuf-2.0-0 \
        libffi-dev \
        shared-mime-info \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

# Copy dcm2niix from builder stage
COPY --from=builder /opt/dcm2niix/bin/dcm2niix /usr/local/bin/dcm2niix

# Set up application directory
WORKDIR /opt/pancData3

# Install Python dependencies first (layer caching)
COPY analysis/requirements.txt /opt/pancData3/analysis/requirements.txt
RUN python3.12 -m pip install --no-cache-dir --break-system-packages \
        -r /opt/pancData3/analysis/requirements.txt

# Copy pipeline and analysis code (no patient data or config)
COPY pipeline/ /opt/pancData3/pipeline/
COPY analysis/ /opt/pancData3/analysis/

# Copy MCR version file for runtime verification
COPY .matlab_version /opt/pancData3/.matlab_version

# Copy entrypoint script
COPY docker/entrypoint.sh /opt/pancData3/entrypoint.sh
RUN chmod +x /opt/pancData3/entrypoint.sh

# Create mount points
RUN mkdir -p /opt/pancData3/data /opt/pancData3/output

# Re-declare ARG after FROM so it's available in this stage
ARG MCR_VERSION=r2024a

# Set PATH to include MATLAB Runtime, dcm2niix, and Python
ENV PATH="/usr/local/bin:/opt/pancData3:${PATH}"
# MATLAB Runtime paths (location in the mathworks image, version-dependent)
ENV LD_LIBRARY_PATH="/opt/matlabruntime/${MCR_VERSION}/runtime/glnxa64:/opt/matlabruntime/${MCR_VERSION}/bin/glnxa64:/opt/matlabruntime/${MCR_VERSION}/sys/os/glnxa64:${LD_LIBRARY_PATH}"
ENV MCR_ROOT="/opt/matlabruntime/${MCR_VERSION}"
ENV MCR_VERSION="${MCR_VERSION}"

# Volumes for patient data and results
VOLUME ["/opt/pancData3/data", "/opt/pancData3/output"]

ENTRYPOINT ["/opt/pancData3/entrypoint.sh"]
CMD ["all"]
