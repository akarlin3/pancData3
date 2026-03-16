# =============================================================================
# pancData3 — Multi-stage Docker build
# Stage 1: Build dcm2niix from source
# Stage 2: Runtime with MATLAB Runtime + Python 3.12
# =============================================================================

# ---------------------------------------------------------------------------
# Stage 1 — Build dcm2niix from latest GitHub release
# ---------------------------------------------------------------------------
FROM debian:bookworm-20240211-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential=12.9 \
        cmake=3.25.1-1 \
        git=1:2.39.5-0+deb12u1 \
        ca-certificates=20230311 \
        zlib1g-dev=1:1.2.13.dfsg-1 \
    && rm -rf /var/lib/apt/lists/*

# Clone and build dcm2niix (pinned to latest stable tag)
RUN git clone --branch v1.0.20240202 --depth 1 \
        https://github.com/rordenlab/dcm2niix.git /tmp/dcm2niix \
    && mkdir /tmp/dcm2niix/build \
    && cd /tmp/dcm2niix/build \
    && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/dcm2niix .. \
    && make -j"$(nproc)" \
    && make install

# ---------------------------------------------------------------------------
# Stage 2 — Runtime image with MATLAB Runtime and Python 3.12
# ---------------------------------------------------------------------------
FROM mathworks/matlab-runtime:r2024a

LABEL maintainer="Avery Karlin <akarlin3>" \
      description="pancData3: Pancreatic DWI Analysis Pipeline" \
      license="AGPL-3.0"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install Python 3.12, WeasyPrint system dependencies, and runtime libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3.12=3.12.1-2+b1 \
        python3.12-venv=3.12.1-2+b1 \
        python3-pip=23.0.1+dfsg-1 \
        libcairo2=1.16.0-7 \
        libpango-1.0-0=1.50.12+ds-1 \
        libpangocairo-1.0-0=1.50.12+ds-1 \
        libgdk-pixbuf-2.0-0=2.42.10+dfsg-1+b1 \
        libffi-dev=3.4.4-1 \
        shared-mime-info=2.2-1 \
        zlib1g=1:1.2.13.dfsg-1 \
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

# Copy entrypoint script
COPY docker/entrypoint.sh /opt/pancData3/entrypoint.sh
RUN chmod +x /opt/pancData3/entrypoint.sh

# Create mount points
RUN mkdir -p /opt/pancData3/data /opt/pancData3/output

# Set PATH to include MATLAB Runtime, dcm2niix, and Python
ENV PATH="/usr/local/bin:/opt/pancData3:${PATH}"
# MATLAB Runtime paths (r2024a default location in the mathworks image)
ENV LD_LIBRARY_PATH="/opt/matlabruntime/R2024a/runtime/glnxa64:/opt/matlabruntime/R2024a/bin/glnxa64:/opt/matlabruntime/R2024a/sys/os/glnxa64:${LD_LIBRARY_PATH}"
ENV MCR_ROOT="/opt/matlabruntime/R2024a"

# Volumes for patient data and results
VOLUME ["/opt/pancData3/data", "/opt/pancData3/output"]

ENTRYPOINT ["/opt/pancData3/entrypoint.sh"]
CMD ["all"]
