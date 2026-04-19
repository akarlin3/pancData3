function [dice_adc, hd_max_adc, hd95_adc, dice_d, hd_max_d, hd95_d, ...
    dice_f, hd_max_f, hd95_f, dice_dstar, hd_max_dstar, hd95_dstar] = ...
    compute_spatial_repeatability(data_vectors_gtvp, j, dwi_type, ...
    gtv_locations, adc_thresh, d_thresh, f_thresh, dstar_thresh, ...
    morph_se, morph_min_cc, last_rpt_gtv_mat, last_rpt_gtv_mask)
% COMPUTE_SPATIAL_REPEATABILITY — Computes Dice and Hausdorff between Fx1 repeat sub-volumes
%
% Loads 3D GTV masks for each valid Fx1 repeat scan, constructs
% threshold-defined binary sub-volumes, applies morphological cleanup,
% and computes pairwise Dice/Hausdorff metrics averaged across all pairs.
%
% Inputs:
%   data_vectors_gtvp - Full data vectors struct array
%   j                 - Patient index
%   dwi_type          - DWI pipeline type index (1=Standard, 2=DnCNN, 3=IVIMnet)
%   gtv_locations     - Cell array of GTV path locations (patients x timepoints x repeats)
%   adc_thresh        - ADC threshold for restricted diffusion sub-volume
%   d_thresh          - D threshold for restricted diffusion sub-volume
%   f_thresh          - f threshold for low-perfusion sub-volume
%   dstar_thresh      - D* threshold for sub-volume
%   morph_se          - Morphological structuring element for cleanup
%   morph_min_cc      - Minimum connected component size (voxels)
%   last_rpt_gtv_mat  - Path of last cached GTV mask file (for I/O reduction)
%   last_rpt_gtv_mask - Cached 3D GTV mask corresponding to last_rpt_gtv_mat
%
% Outputs:
%   dice_adc, hd_max_adc, hd95_adc       - ADC sub-volume spatial repeatability
%   dice_d, hd_max_d, hd95_d             - D sub-volume spatial repeatability
%   dice_f, hd_max_f, hd95_f             - f sub-volume spatial repeatability
%   dice_dstar, hd_max_dstar, hd95_dstar - D* sub-volume spatial repeatability
%   (all NaN if insufficient repeats or missing 3D masks)

% Initialize all outputs to NaN
dice_adc = NaN; hd_max_adc = NaN; hd95_adc = NaN;
dice_d = NaN; hd_max_d = NaN; hd95_d = NaN;
dice_f = NaN; hd_max_f = NaN; hd95_f = NaN;
dice_dstar = NaN; hd_max_dstar = NaN; hd95_dstar = NaN;

% Determine voxel dimensions for physical Hausdorff distances (mm).
% Hausdorff distance measures worst-case spatial disagreement between
% two binary masks; it must be in physical units (mm) rather than voxel
% indices so that results are comparable across different scan protocols
% with varying voxel sizes.
rpt_vox_dims = data_vectors_gtvp(j,1,1).vox_dims;
if isempty(rpt_vox_dims) || ~isnumeric(rpt_vox_dims) || numel(rpt_vox_dims) ~= 3
    % Fallback: estimate isotropic voxel dimensions from voxel volume.
    % vox_vol is stored in cm^3, so convert to mm^3 (* 1000) before
    % taking the cube root to get the side length in mm.
    rpt_vox_vol = data_vectors_gtvp(j,1,1).vox_vol;
    if ~isempty(rpt_vox_vol) && ~isnan(rpt_vox_vol) && rpt_vox_vol > 0
        side_mm = (rpt_vox_vol * 1000) ^ (1/3);
        rpt_vox_dims = [side_mm, side_mm, side_mm];
    else
        % Last resort: assume unit voxel spacing (results in voxel-space distances)
        rpt_vox_dims = [1, 1, 1];
    end
end

% Collect valid repeat indices and their parameter vectors.
% Fx1 (fraction 1) scans may have multiple repeat acquisitions (same day)
% to assess test-retest reproducibility of diffusion measurements.
% The 3rd dimension of data_vectors_gtvp indexes these repeat scans.
valid_rpis = [];
rpt_vecs = struct('adc', {{}}, 'd', {{}}, 'f', {{}}, 'dstar', {{}});
for rpi2 = 1:size(data_vectors_gtvp, 3)
    % Extract DWI-type-specific parameter vectors (Standard/DnCNN/IVIMnet)
    [rv_adc, rv_d, rv_f, rv_ds] = select_dwi_vectors(data_vectors_gtvp, j, 1, rpi2, dwi_type);
    if ~isempty(rv_adc) || ~isempty(rv_d)
        valid_rpis(end+1) = rpi2; %#ok<AGROW>
        rpt_vecs.adc{end+1} = rv_adc;
        rpt_vecs.d{end+1} = rv_d;
        rpt_vecs.f{end+1} = rv_f;
        rpt_vecs.dstar{end+1} = rv_ds;
    end
end

% Need at least 2 valid repeats to compute pairwise spatial overlap
if numel(valid_rpis) < 2
    return;
end

% Load 3D GTV masks for each valid repeat.
% The 3D mask is needed to embed 1D voxel parameter vectors back into
% their spatial positions for morphological cleanup and Dice/Hausdorff
% computation. Without 3D masks, spatial repeatability cannot be assessed.
%
% Shared-mask fallback: many sites contour a single GTV per Fx1 session
% (one indexed file, typically GTV1.mat or GTV.mat) and reuse it across
% every back-to-back repeat because the tumour does not move between
% acquisitions.  In that case gtv_locations{j,1,1} is populated but
% {j,1,2..N} are empty, and without a fallback this function would return
% NaN for every patient in the cohort.  Find the first non-empty Fx1 path
% and reuse it for any repeat that is missing its own path.
shared_fx1_path = '';
for ri = 1:size(gtv_locations, 3)
    candidate = gtv_locations{j, 1, ri};
    if ~isempty(candidate); shared_fx1_path = candidate; break; end
end

rpt_masks_3d = cell(numel(valid_rpis), 1);
rpt_has_3d = true;
for ri = 1:numel(valid_rpis)
    rpi_idx = valid_rpis(ri);
    gtv_mat_path = gtv_locations{j, 1, rpi_idx};
    if isempty(gtv_mat_path); gtv_mat_path = shared_fx1_path; end
    if ~isempty(gtv_mat_path)
        % Normalize path separators for cross-platform compatibility
        % (strsplit on both / and \ handles Windows and Unix paths)
        path_parts = strsplit(gtv_mat_path, {'/', '\'});
        gtv_mat_path = fullfile(path_parts{:});
        % Restore leading slash for absolute Unix paths (strsplit produces
        % an empty first element for paths starting with /)
        if isunix && ~startsWith(gtv_mat_path, filesep) && isempty(path_parts{1})
            gtv_mat_path = [filesep gtv_mat_path]; %#ok<AGROW>
        end
        if exist(gtv_mat_path, 'file')
            % Cache the last loaded mask to avoid redundant disk I/O when
            % multiple repeats share the same GTV contour file
            if strcmp(gtv_mat_path, last_rpt_gtv_mat)
                rpt_masks_3d{ri} = last_rpt_gtv_mask;
            else
                % safe_load_mask validates variable type before loading
                % to prevent arbitrary code execution from malicious .mat files
                rpt_masks_3d{ri} = safe_load_mask(gtv_mat_path, 'Stvol3d');
                last_rpt_gtv_mat = gtv_mat_path;
                last_rpt_gtv_mask = rpt_masks_3d{ri};
            end
        end
    end
    if isempty(rpt_masks_3d{ri})
        rpt_has_3d = false;
    end
end

% All repeats must have valid 3D masks for spatial comparison
if ~rpt_has_3d
    return;
end

% Accumulate pairwise Dice and Hausdorff metrics across all repeat pairs.
% With N repeats, there are N*(N-1)/2 unique pairs. Each pair gets its
% own row; columns correspond to the 4 diffusion parameters.
param_names = {'adc', 'd', 'f', 'dstar'};
param_thresholds = [adc_thresh, d_thresh, f_thresh, dstar_thresh];
pair_dice = nan(numel(valid_rpis)*(numel(valid_rpis)-1)/2, 4);
pair_hd_max = nan(size(pair_dice));
pair_hd95 = nan(size(pair_dice));
pi_count = 0;  % running pair counter

for ri1 = 1:numel(valid_rpis)-1
    for ri2 = ri1+1:numel(valid_rpis)
        pi_count = pi_count + 1;
        mask_3d_1 = rpt_masks_3d{ri1};
        mask_3d_2 = rpt_masks_3d{ri2};

        % Verify compatible 3D mask dimensions (mismatched grids indicate
        % different reconstruction matrices and cannot be directly compared)
        if ~isequal(size(mask_3d_1), size(mask_3d_2))
            continue;
        end

        for pi = 1:4
            vec_1 = rpt_vecs.(param_names{pi}){ri1};
            vec_2 = rpt_vecs.(param_names{pi}){ri2};
            if isempty(vec_1) || isempty(vec_2)
                continue;
            end

            % Verify that the 1D parameter vector length matches the number
            % of GTV voxels in the 3D mask (data integrity check)
            n_gtv_1 = sum(mask_3d_1(:) == 1);
            n_gtv_2 = sum(mask_3d_2(:) == 1);
            if numel(vec_1) ~= n_gtv_1 || numel(vec_2) ~= n_gtv_2
                continue;
            end

            % Threshold parameter vector to create a binary sub-volume.
            % For ADC and D: values below threshold indicate restricted
            % diffusion (high cellularity / treatment resistance).
            % For f: values below threshold indicate low perfusion fraction.
            % Embed the binary vector back into the 3D GTV mask positions.
            subvol_3d_1 = false(size(mask_3d_1));
            subvol_3d_1(mask_3d_1 == 1) = vec_1 < param_thresholds(pi);
            % Morphological open-then-close removes small noise speckles
            % (open) and fills small gaps (close) for cleaner sub-volumes.
            % bwareaopen removes connected components smaller than the
            % minimum size threshold to eliminate isolated noise voxels.
            subvol_3d_1 = imclose(imopen(subvol_3d_1, morph_se), morph_se);
            subvol_3d_1 = bwareaopen(subvol_3d_1, morph_min_cc);

            subvol_3d_2 = false(size(mask_3d_2));
            subvol_3d_2(mask_3d_2 == 1) = vec_2 < param_thresholds(pi);
            subvol_3d_2 = imclose(imopen(subvol_3d_2, morph_se), morph_se);
            subvol_3d_2 = bwareaopen(subvol_3d_2, morph_min_cc);

            % Compute spatial overlap (Dice) and distance (Hausdorff) metrics.
            % Dice ranges [0,1] (1 = perfect overlap); Hausdorff is in mm
            % (0 = identical boundaries). HD95 is the 95th-percentile
            % Hausdorff distance, less sensitive to single-voxel outliers.
            [d_val, hm_val, h95_val] = compute_dice_hausdorff( ...
                subvol_3d_1, subvol_3d_2, rpt_vox_dims);
            pair_dice(pi_count, pi) = d_val;
            pair_hd_max(pi_count, pi) = hm_val;
            pair_hd95(pi_count, pi) = h95_val;
        end
    end
end

% Average across all repeat pairs to produce a single repeatability
% estimate per parameter. NaN-safe mean handles pairs where one repeat
% had missing data for a particular parameter.
if exist('OCTAVE_VERSION', 'builtin')
    mean_dice = mean(pair_dice(1:pi_count,:), 1, 'omitnan');
    mean_hd_max = mean(pair_hd_max(1:pi_count,:), 1, 'omitnan');
    mean_hd95 = mean(pair_hd95(1:pi_count,:), 1, 'omitnan');
else
    mean_dice = nanmean(pair_dice(1:pi_count,:), 1);
    mean_hd_max = nanmean(pair_hd_max(1:pi_count,:), 1);
    mean_hd95 = nanmean(pair_hd95(1:pi_count,:), 1);
end
% Unpack columns: 1=ADC, 2=D, 3=f, 4=D*
dice_adc = mean_dice(1);
dice_d   = mean_dice(2);
dice_f   = mean_dice(3);
dice_dstar = mean_dice(4);
hd_max_adc = mean_hd_max(1);
hd_max_d   = mean_hd_max(2);
hd_max_f   = mean_hd_max(3);
hd_max_dstar = mean_hd_max(4);
hd95_adc = mean_hd95(1);
hd95_d   = mean_hd95(2);
hd95_f   = mean_hd95(3);
hd95_dstar = mean_hd95(4);

end

%% --- Local helper function ---

function [adc_vec, d_vec, f_vec, dstar_vec] = select_dwi_vectors(data_vectors_gtvp, j, k, rpi, dwi_type)
% SELECT_DWI_VECTORS  Extracts the appropriate parameter vectors based on DWI type.
%   Standard (1): uses native ADC and IVIM vectors from conventional fitting.
%   DnCNN (2): uses vectors from DnCNN-denoised images (suffix _dncnn).
%   IVIMnet (3): uses native ADC but deep-learning IVIM estimates (suffix _ivimnet).
%   Note: IVIMnet uses native ADC because IVIMnet only re-estimates the IVIM
%   parameters (D, f, D*) via its neural network, not ADC.
switch dwi_type
    case 1
        adc_vec   = data_vectors_gtvp(j,k,rpi).adc_vector;
        d_vec     = data_vectors_gtvp(j,k,rpi).d_vector;
        f_vec     = data_vectors_gtvp(j,k,rpi).f_vector;
        dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector;
    case 2
        adc_vec   = data_vectors_gtvp(j,k,rpi).adc_vector_dncnn;
        d_vec     = data_vectors_gtvp(j,k,rpi).d_vector_dncnn;
        f_vec     = data_vectors_gtvp(j,k,rpi).f_vector_dncnn;
        dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector_dncnn;
    case 3
        adc_vec   = data_vectors_gtvp(j,k,rpi).adc_vector;
        d_vec     = data_vectors_gtvp(j,k,rpi).d_vector_ivimnet;
        f_vec     = data_vectors_gtvp(j,k,rpi).f_vector_ivimnet;
        dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector_ivimnet;
end
end
