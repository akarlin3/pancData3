function [data_gtvp, data_gtvn, summary_metrics_out] = load_data_from_disk( ...
    dwi_vectors_file, fallback_dwi_vectors_file, summary_metrics_file, current_dtype, current_name)
% LOAD_DATA_FROM_DISK  Load DWI vectors and summary metrics from disk with legacy fallback.
%
%   [data_gtvp, data_gtvn, summary_metrics_out] = load_data_from_disk( ...
%       dwi_vectors_file, fallback_dwi_vectors_file, summary_metrics_file, ...
%       current_dtype, current_name)
%
%   Loads previously computed voxel-wise parameter vectors and summary
%   metrics from disk.  Implements the legacy fallback logic: if the
%   type-specific dwi_vectors file is not found, falls back to the
%   fallback_dwi_vectors_file (legacy dwi_vectors.mat without type suffix)
%   ONLY for Standard (current_dtype == 1) to prevent cross-contamination.
%
%   Inputs:
%     dwi_vectors_file          - Full path to type-specific dwi_vectors .mat file
%     fallback_dwi_vectors_file - Full path to legacy dwi_vectors.mat (no type suffix)
%     summary_metrics_file      - Full path to summary_metrics .mat file
%     current_dtype             - Numeric DWI type index (1=Standard, 2=dnCNN, 3=IVIMnet)
%     current_name              - Human-readable DWI type name (e.g., 'Standard')
%
%   Outputs:
%     data_gtvp             - Loaded data_vectors_gtvp struct
%     data_gtvn             - Loaded data_vectors_gtvn struct
%     summary_metrics_out   - Loaded summary_metrics struct (empty [] if
%                             summary_metrics_file not found)
%
%   Errors:
%     Throws an error if neither the type-specific nor the legacy fallback
%     DWI vectors file can be found.

    % Determine which file to load, with legacy fallback
    if exist(dwi_vectors_file, 'file')
        target_dwi_file = dwi_vectors_file;
    elseif exist(fallback_dwi_vectors_file, 'file')
        % Legacy file without DWI-type suffix exists.  Only allow
        % fallback for Standard (type 1) to prevent silently loading
        % data from a different DWI processing method.
        if current_dtype == 1
            fprintf('  💡 Using legacy dwi_vectors.mat (no type suffix) for Standard.\n');
            target_dwi_file = fallback_dwi_vectors_file;
        else
            fprintf('  ❌ Type-specific file %s not found and legacy dwi_vectors.mat cannot be used for %s (risk of cross-contamination).\n', ...
                dwi_vectors_file, current_name);
            target_dwi_file = '';
        end
    else
        target_dwi_file = '';
    end

    if isempty(target_dwi_file)
        error('LoadDataFromDisk:NotFound', 'Data vectors not found. Please run "load" step first.');
    end

    % Load the DWI vectors.
    % data_vectors_gtvp: struct array (patients x timepoints x repeats) with
    %   per-voxel parameter vectors (ADC, D, f, D*) for the primary GTV (GTVp).
    % data_vectors_gtvn: same structure for nodal GTV (GTVn), if available.
    % Using explicit variable names in load() prevents loading the entire
    % .mat file into memory (which may contain large intermediate variables).
    tmp_vectors = load(target_dwi_file, 'data_vectors_gtvp', 'data_vectors_gtvn');
    data_gtvp = tmp_vectors.data_vectors_gtvp;
    data_gtvn = tmp_vectors.data_vectors_gtvn;
    fprintf('      💾 Loaded data from disk (%s).\n', target_dwi_file);

    % Load summary metrics (patient-level aggregated statistics) if available.
    % Summary metrics are computed by compute_summary_metrics.m and saved
    % separately. They may not exist if only the 'load' step was run
    % previously without the 'metrics' step.
    summary_metrics_out = [];
    if ~isempty(summary_metrics_file) && exist(summary_metrics_file, 'file')
        tmp_metrics = load(summary_metrics_file, 'summary_metrics');
        summary_metrics_out = tmp_metrics.summary_metrics;
        fprintf('      💾 Loaded summary_metrics from disk (%s).\n', summary_metrics_file);
    end
end
