function [data_gtvp, data_gtvn, summary_metrics_out] = load_data_from_disk( ...
    voxel_cache_file, voxel_cache_fallback_file, summary_metrics_file, current_dtype, current_name)
% LOAD_DATA_FROM_DISK  Load voxel cache and summary metrics from disk with fallback.
%
%   [data_gtvp, data_gtvn, summary_metrics_out] = load_data_from_disk( ...
%       voxel_cache_file, voxel_cache_fallback_file, summary_metrics_file, ...
%       current_dtype, current_name)
%
%   Loads previously computed voxel-wise parameter vectors and summary
%   metrics from disk.  Implements the fallback logic: if the
%   type-specific pipeline_voxels file is not found, falls back to the
%   un-typed pipeline_voxels.mat ONLY for Standard (current_dtype == 1)
%   to prevent cross-contamination between DWI processing variants.
%
%   Inputs:
%     voxel_cache_file          - Full path to type-specific pipeline_voxels_<type>.mat
%     voxel_cache_fallback_file - Full path to un-typed pipeline_voxels.mat fallback
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
%     Throws an error if neither the type-specific nor the fallback
%     voxel-cache file can be found.

    % Determine which file to load, with fallback
    if exist(voxel_cache_file, 'file')
        target_cache_file = voxel_cache_file;
    elseif exist(voxel_cache_fallback_file, 'file')
        % Un-typed file exists.  Only allow fallback for Standard
        % (type 1) to prevent silently loading data from a different
        % DWI processing method.
        if current_dtype == 1
            fprintf('  💡 Using un-typed pipeline_voxels.mat for Standard.\n');
            target_cache_file = voxel_cache_fallback_file;
        else
            fprintf('  ❌ Type-specific file %s not found and un-typed pipeline_voxels.mat cannot be used for %s (risk of cross-contamination).\n', ...
                voxel_cache_file, current_name);
            target_cache_file = '';
        end
    else
        target_cache_file = '';
    end

    if isempty(target_cache_file)
        error('LoadDataFromDisk:NotFound', 'Voxel cache not found. Please run "load" step first.');
    end

    % Load the voxel-level vectors.
    % data_vectors_gtvp: struct array (patients x timepoints x repeats) with
    %   per-voxel parameter vectors (ADC, D, f, D*) for the primary GTV (GTVp).
    % data_vectors_gtvn: same structure for nodal GTV (GTVn), if available.
    % Using explicit variable names in load() prevents loading the entire
    % .mat file into memory (which may contain large intermediate variables).
    tmp_vectors = load(target_cache_file, 'data_vectors_gtvp', 'data_vectors_gtvn');
    data_gtvp = tmp_vectors.data_vectors_gtvp;
    data_gtvn = tmp_vectors.data_vectors_gtvn;
    fprintf('      💾 Loaded data from disk (%s).\n', target_cache_file);

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
