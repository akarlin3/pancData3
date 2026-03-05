function [d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub] = metrics_dosimetry(m_id_list, id_list, nTp, config_struct, m_data_vectors_gtvp, gtv_locations)
% METRICS_DOSIMETRY — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 3/5 of the metrics step. Computes dose metrics (D95, V50) for resistant sub-volumes.
%
% ANALYTICAL OVERVIEW:
%   This module answers a key question in adaptive radiotherapy: "Is the
%   treatment-resistant tumour sub-region receiving adequate radiation dose?"
%
%   The concept of "resistant sub-volumes" is central to this analysis:
%   - Low ADC/D regions within the GTV may indicate areas of high cellularity
%     (densely packed tumour cells) that are more resistant to radiation.
%   - Low f regions may indicate hypoxic areas with poor blood supply, which
%     are known to be more radioresistant (oxygen fixation hypothesis).
%   - If these resistant sub-volumes receive inadequate dose coverage (low
%     D95 or V50), the tumour is more likely to recur locally.
%
%   Sub-volume definition: Voxels within the GTV whose diffusion parameter
%   falls below a configured threshold are classified as "resistant."  The
%   thresholds (adc_thresh, d_thresh, f_thresh, dstar_thresh) are set in
%   config.json based on prior literature or institutional data.
%
%   Dose metrics computed for each sub-volume:
%     D95 = minimum dose (Gy) delivered to 95% of the sub-volume.
%           A high D95 means nearly all resistant tissue received adequate dose.
%     V50 = fraction of the sub-volume receiving >= 50 Gy.
%           50 Gy is a typical prescription dose for conventional fractionation;
%           for SBRT, this threshold may need adjustment.
%
%   These metrics are computed at each timepoint because the resistant
%   sub-volume may shrink during treatment (as radiation kills sensitive
%   cells, the remaining resistant population may become more concentrated).
%
% Inputs:
%   m_id_list         - Cell array of valid patient identifiers
%   id_list           - Original cell array of all patient identifiers
%   nTp               - Number of timepoints
%   config_struct     - Configuration settings (thresholds)
%   m_data_vectors_gtvp- Struct array containing parametric vectors and dose vectors
%   gtv_locations     - Cell array showing paths to 3D mask arrays
%
% Outputs:
%   d95_*_sub         - D95 dose (Gy) delivered to the resistant sub-volume defined by *
%   v50_*_sub         - V50Gy dose coverage (%) delivered to the resistant sub-volume defined by *
%

fprintf('  --- SECTION 6: Target Coverage (Sub-Volume Dose Metrics) ---\n');

% Diary: capture console output to output_folder
if isfield(config_struct, 'output_folder')
    diary_file = fullfile(config_struct.output_folder, 'metrics_dosimetry_output.txt');
    if exist(diary_file, 'file'), delete(diary_file); end
    diary(diary_file);
end

% Thresholds define the boundary between "resistant" and "non-resistant"
% voxels within the GTV for each diffusion parameter.  Voxels BELOW
% the threshold for ADC, D, and D* (or ABOVE for f in some protocols)
% are classified as resistant.  These are typically derived from
% literature values or ROC analysis on a training cohort.
adc_thresh   = config_struct.adc_thresh;    % mm^2/s — low ADC = high cellularity
d_thresh     = config_struct.d_thresh;      % mm^2/s — low D = restricted diffusion
f_thresh     = config_struct.f_thresh;      % fraction — low f = poor perfusion (hypoxic)
dstar_thresh = config_struct.dstar_thresh;  % mm^2/s — low D* = slow microcirculation  

if iscell(nTp)
    nTp = nTp{1};
end

% Initialize output arrays with NaN (patients x timepoints).
% NaN indicates missing data (no dose map available, no GTV segmentation,
% or sub-volume too small to compute meaningful dose statistics).
% Each diffusion parameter defines its own sub-volume, so dose coverage
% metrics are computed independently for each parameter's resistant region.
d95_adc_sub = nan(length(m_id_list), nTp);
v50_adc_sub = nan(length(m_id_list), nTp);

d95_d_sub = nan(length(m_id_list), nTp);
v50_d_sub = nan(length(m_id_list), nTp);

d95_f_sub = nan(length(m_id_list), nTp);
v50_f_sub = nan(length(m_id_list), nTp);

d95_dstar_sub = nan(length(m_id_list), nTp);
v50_dstar_sub = nan(length(m_id_list), nTp);

% Cache the last loaded GTV mask to avoid redundant disk I/O.
% Multiple timepoints for the same patient often share the same GTV mask
% (when the contour was not re-drawn at each fraction), so caching
% provides significant speedup for large cohorts.
last_gtv_mat = '';
last_gtv_mask_3d = [];

n_pat_dosimetry = length(m_id_list);
for j = 1:n_pat_dosimetry
    text_progress_bar(j, n_pat_dosimetry, 'Computing dosimetry metrics');
    j_orig = find(strcmp(id_list, m_id_list{j}), 1, 'first');
    for k = 1:nTp
        if isfield(config_struct, 'dwi_types_to_run') && isscalar(config_struct.dwi_types_to_run)
            dtype_idx = config_struct.dwi_types_to_run;
        else
            dtype_idx = 1;
        end
        
        % Select the voxel-level parameter vectors for the current DWI
        % processing pipeline.  Each vector contains one value per GTV voxel,
        % enabling voxel-wise thresholding to define resistant sub-volumes.
        % The voxel ordering matches the GTV mask and dose vector, so
        % element-wise comparisons correctly identify spatial sub-regions.
        switch dtype_idx
            case 1
                % Standard pipeline: conventional IVIM fitting without denoising
                adc_vec = m_data_vectors_gtvp(j,k,1).adc_vector;
                d_vec   = m_data_vectors_gtvp(j,k,1).d_vector;
                f_vec   = m_data_vectors_gtvp(j,k,1).f_vector;
                dstar_vec = m_data_vectors_gtvp(j,k,1).dstar_vector;
            case 2
                % dnCNN pipeline: deep learning denoised signal before IVIM fit
                adc_vec = m_data_vectors_gtvp(j,k,1).adc_vector_dncnn;
                d_vec   = m_data_vectors_gtvp(j,k,1).d_vector_dncnn;
                f_vec   = m_data_vectors_gtvp(j,k,1).f_vector_dncnn;
                dstar_vec = m_data_vectors_gtvp(j,k,1).dstar_vector_dncnn;
            case 3
                % IVIMnet pipeline: neural network directly estimates IVIM parameters
                % ADC is always from the standard mono-exponential fit (IVIMnet
                % does not produce its own ADC estimate because ADC is model-free)
                adc_vec = m_data_vectors_gtvp(j,k,1).adc_vector;
                d_vec   = m_data_vectors_gtvp(j,k,1).d_vector_ivimnet;
                f_vec   = m_data_vectors_gtvp(j,k,1).f_vector_ivimnet;
                dstar_vec = m_data_vectors_gtvp(j,k,1).dstar_vector_ivimnet;
        end
        % The dose vector contains the radiation dose (Gy) at each GTV voxel,
        % sampled from the RT dose grid onto the DWI image space.  This
        % spatial co-registration enables voxel-by-voxel dose-response analysis.
        dose_vec  = m_data_vectors_gtvp(j,k,1).dose_vector;

        % Both dose and diffusion data must be available for sub-volume analysis
        if ~isempty(dose_vec) && ~isempty(adc_vec)
            gtv_mat = gtv_locations{j_orig, k, 1};
            has_3d = false;
            gtv_mask_3d = [];
            if ~isempty(gtv_mat)
                path_parts = strsplit(gtv_mat, {'/', '\'});
                gtv_mat = fullfile(path_parts{:});
                if isunix && ~startsWith(gtv_mat, filesep) && isempty(path_parts{1})
                    gtv_mat = [filesep gtv_mat];
                end

                if exist(gtv_mat, 'file')
                    if strcmp(gtv_mat, last_gtv_mat)
                        gtv_mask_3d = last_gtv_mask_3d;
                    else
                        % SECURITY: Use safe_load_mask to prevent unsafe deserialization of untrusted .mat files
                        gtv_mask_3d = safe_load_mask(gtv_mat, 'Stvol3d');
                        last_gtv_mat = gtv_mat;
                        last_gtv_mask_3d = gtv_mask_3d;
                    end

                    if ~isempty(gtv_mask_3d) && sum(gtv_mask_3d(:) == 1) == length(adc_vec)
                        has_3d = true;
                    elseif ~isempty(gtv_mask_3d) && k > 1
                        % Native fraction mask doesn't match vector length.
                        % When DIR warping is active, vectors are in Fx1 space
                        % but gtv_locations{j_orig, k, 1} points to the native
                        % fraction mask. Fall back to the Fx1 reference mask.
                        fx1_mat = gtv_locations{j_orig, 1, 1};
                        if ~isempty(fx1_mat)
                            fx1_path_parts = strsplit(fx1_mat, {'/', '\'});
                            fx1_mat = fullfile(fx1_path_parts{:});
                            if isunix && ~startsWith(fx1_mat, filesep) && isempty(fx1_path_parts{1})
                                fx1_mat = [filesep fx1_mat];
                            end
                            if exist(fx1_mat, 'file')
                                fx1_mask_3d = safe_load_mask(fx1_mat, 'Stvol3d');
                                if ~isempty(fx1_mask_3d) && sum(fx1_mask_3d(:) == 1) == length(adc_vec)
                                    gtv_mask_3d = fx1_mask_3d;
                                    has_3d = true;
                                    warning('metrics_dosimetry:fx1MaskFallback', ...
                                        'Patient %s Fx%d: native mask size mismatch (%d vs %d voxels). Using Fx1 reference mask for 3D sub-volume cleanup.', ...
                                        m_id_list{j}, k, sum(last_gtv_mask_3d(:) == 1), length(adc_vec));
                                end
                            end
                        end
                    end
                end
            end
            
            % Compute D95 and V50 for each diffusion-defined resistant sub-volume.
            % calculate_subvolume_metrics identifies voxels below threshold,
            % optionally applies 3D morphological cleanup (erosion to remove
            % single-voxel noise), and computes dose-volume histogram statistics
            % on the resulting sub-region.  If the sub-volume is empty (no
            % voxels below threshold), NaN is returned — indicating the entire
            % GTV appears non-resistant by that metric.
            [d95_adc_sub(j,k), v50_adc_sub(j,k)] = calculate_subvolume_metrics(adc_vec, adc_thresh, dose_vec, has_3d, gtv_mask_3d);
            [d95_d_sub(j,k), v50_d_sub(j,k)]     = calculate_subvolume_metrics(d_vec, d_thresh, dose_vec, has_3d, gtv_mask_3d);
            [d95_f_sub(j,k), v50_f_sub(j,k)]     = calculate_subvolume_metrics(f_vec, f_thresh, dose_vec, has_3d, gtv_mask_3d);
            [d95_dstar_sub(j,k), v50_dstar_sub(j,k)] = calculate_subvolume_metrics(dstar_vec, dstar_thresh, dose_vec, has_3d, gtv_mask_3d);
            
        end
    end
end

% Close diary for this module
if isfield(config_struct, 'output_folder')
    diary off;
end
end
