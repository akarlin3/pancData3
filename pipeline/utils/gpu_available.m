function [available, gpu_dev] = gpu_available(device_index)
% GPU_AVAILABLE Check whether a usable CUDA GPU is available in MATLAB.
%
%   [available, gpu_dev] = gpu_available()
%   [available, gpu_dev] = gpu_available(device_index)
%
%   Returns available=true and the gpuDevice object when:
%     1. The Parallel Computing Toolbox is installed
%     2. canUseGPU() returns true
%     3. A GPU with sufficient compute capability (>= 3.5) exists
%
%   device_index (optional): 1-based GPU index to select (default 1).
%
%   When no GPU is available, returns available=false and gpu_dev=[].
%   Never errors — designed for graceful fallback to CPU.

    if nargin < 1 || isempty(device_index)
        device_index = 1;
    end

    available = false;
    gpu_dev = [];

    % Guard 1: Parallel Computing Toolbox must be installed
    if ~exist('gpuDevice', 'file')
        return;
    end

    % Guard 2: Octave does not support gpuDevice
    if exist('OCTAVE_VERSION', 'builtin')
        return;
    end

    % Guard 3: Quick MATLAB built-in check (R2022a+)
    if exist('canUseGPU', 'file') && ~canUseGPU()
        return;
    end

    try
        gpu_dev = gpuDevice(device_index);

        % Guard 4: Minimum compute capability for modern CUDA kernels
        if str2double(gpu_dev.ComputeCapability) < 3.5
            fprintf('  [GPU] Device "%s" compute capability %s < 3.5; falling back to CPU.\n', ...
                gpu_dev.Name, gpu_dev.ComputeCapability);
            gpu_dev = [];
            return;
        end

        available = true;
        fprintf('  [GPU] Using device %d: %s (%.1f GB, compute %s)\n', ...
            device_index, gpu_dev.Name, gpu_dev.AvailableMemory / 1e9, ...
            gpu_dev.ComputeCapability);
    catch
        % gpuDevice() can throw if CUDA driver is missing/incompatible
        gpu_dev = [];
    end
end
