#!/bin/bash

MATLAB_CACHE_DIR="/workspace/.cache/matlab_R2025a"
if [ -f "$MATLAB_CACHE_DIR/bin/matlab" ]; then
    echo "MATLAB is already installed. Skipping installation."
else
    echo "Installing MATLAB..."
    mkdir -p $MATLAB_CACHE_DIR
    wget -q https://www.mathworks.com/mpm/get/linux -O mpm
    chmod +x mpm
    sudo ./mpm install --release=R2025a --destination=$MATLAB_CACHE_DIR --products MATLAB
    rm mpm
    echo "Installation complete."
fi

export PATH="$MATLAB_CACHE_DIR/bin:$PATH"