function [dataloc, data_vectors_gtvn, data_vectors_gtvp, lf, immuno, ...
    mrn_list, id_list, fx_dates, dwi_locations, rtdose_locations, ...
    gtv_locations, gtvn_locations, dmean_gtvp, dmean_gtvn, d95_gtvp, ...
    d95_gtvn, v50gy_gtvp, v50gy_gtvn, bad_dwi_locations, bad_dwi_count] = ...
    load_dwi_data_reload(config_struct)
    % LOAD_DWI_DATA_RELOAD Section 4 reload path of load_dwi_data.m.
    %   Extracted verbatim from the skip_to_reload (else) branch. Reloads the
    %   pre-processed pipeline_voxels[.type].mat checkpoint so downstream
    %   analysis (Section 5) can skip the expensive Section 1-3 processing.

    %% ========================================================================
    fprintf('\n--- SECTION 4: Reload Saved Data ---\n');
    %  SECTION 4 — RELOAD SAVED DATA
    %  [ENTRY POINT]: If skip_to_reload=true, execution begins here.
    %  Loads the pre-processed 'pipeline_voxels.mat' containing voxel-level data.

    % Set data path from configuration
    dataloc = config_struct.dataloc;

if isfield(config_struct, 'dwi_type_name')
    file_prefix = ['_' config_struct.dwi_type_name];
else
    file_prefix = '';
end
datasave = fullfile(dataloc, ['pipeline_voxels' file_prefix '.mat']);
if ~exist(datasave, 'file') && ~isempty(file_prefix)
    % Fallback to the default (un-typed) file when the variant-specific
    % file has not been created yet (e.g. first run before per-type saves).
    % NOTE: run_dwi_pipeline.m (reload branch, ~line 288) restricts this
    % fallback to Standard (dtype==1) only, preventing cross-type
    % contamination.  This path is only reached during the initial 'load'
    % step or direct calls outside the orchestrator.
    datasave_fallback = fullfile(dataloc, 'pipeline_voxels.mat');
    if exist(datasave_fallback, 'file')
        fprintf('  %s not found — falling back to %s\n', ...
            ['pipeline_voxels' file_prefix '.mat'], 'pipeline_voxels.mat');
        datasave = datasave_fallback;
    end
end
if ~exist(datasave, 'file')
    type_label = '';
    if isfield(config_struct, 'dwi_type_name')
        type_label = config_struct.dwi_type_name;
    end
    error('load_dwi_data:fileNotFound', ...
        'Required data file ''%s'' not found. Run the load step for DWI type ''%s'' before reloading.', ...
        datasave, type_label);
end
tmp_data = load(datasave);
data_vectors_gtvn = tmp_data.data_vectors_gtvn; data_vectors_gtvp = tmp_data.data_vectors_gtvp; lf = tmp_data.lf;
immuno = tmp_data.immuno; mrn_list = tmp_data.mrn_list; id_list = tmp_data.id_list; fx_dates = tmp_data.fx_dates;
dwi_locations = tmp_data.dwi_locations; rtdose_locations = tmp_data.rtdose_locations; gtv_locations = tmp_data.gtv_locations;
gtvn_locations = tmp_data.gtvn_locations; dmean_gtvp = tmp_data.dmean_gtvp; dmean_gtvn = tmp_data.dmean_gtvn;
d95_gtvp = tmp_data.d95_gtvp; d95_gtvn = tmp_data.d95_gtvn; v50gy_gtvp = tmp_data.v50gy_gtvp; v50gy_gtvn = tmp_data.v50gy_gtvn;
bad_dwi_locations = tmp_data.bad_dwi_locations; bad_dwi_count = tmp_data.bad_dwi_count;

if exist('OCTAVE_VERSION', 'builtin') && ~exist('id_list', 'var')
    warning('id_list not loaded from save file. This may occur during mock tests. Proceeding with dummy data.');
    id_list = {}; mrn_list = {}; lf = []; immuno = {}; gtv_locations = []; dwi_locations = []; dmean_gtvp = []; d95_gtvp = []; v50gy_gtvp = []; data_vectors_gtvp = []; data_vectors_gtvn = [];
end
end
