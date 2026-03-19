function core_mask = extract_tumor_core(config_struct, adc_vec, d_vec, f_vec, dstar_vec, has_3d, gtv_mask_3d, opts)
% EXTRACT_TUMOR_CORE — Delineates the viable tumor core from parameter maps.
%
% Integrates multiple algorithmic approaches for isolating the restricted,
% densely-cellular sub-volume of a pancreatic tumor.
%
% Inputs:
%   config_struct  - Configuration settings including 'core_method'
%                    Options: 'adc_threshold', 'd_threshold', 'df_intersection',
%                    'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours',
%                    'percentile', 'spectral', 'fdm'
%   adc_vec        - 1D vector of ADC values
%   d_vec          - 1D vector of true diffusion (D) values
%   f_vec          - 1D vector of perfusion fraction (f) values
%   dstar_vec      - 1D vector of pseudo-diffusion (D*) values
%   has_3d         - Logical, true if a 3D GTV mask is available
%   gtv_mask_3d    - 3D logical mask of the total GTV
%   opts           - (Optional) struct with additional options:
%                    .timepoint_index    - Current fraction index (for fDM)
%                    .baseline_adc_vec   - Baseline ADC voxel vector (for fDM)
%                    .baseline_d_vec     - Baseline D voxel vector (for fDM)
%                    .repeatability_cor  - Coefficient of Reproducibility (for fDM)
%                    .voxel_dims         - [x,y,z] voxel dimensions in mm
%                    .tumor_volume_mm3   - Total tumor volume in mm³
%
% Outputs:
%   core_mask      - A 1D logical mask (same size as input vectors) where
%                    true indicates a voxel belongs to the core.

    if nargin < 8, opts = struct(); end

    % Validate inputs and get method-specific parameters
    [valid_params, fallback_mask] = validate_and_prepare(config_struct, adc_vec, d_vec, f_vec, dstar_vec, gtv_mask_3d, opts);
    if valid_params.should_return_fallback
        core_mask = fallback_mask;
        return;
    end

    % Method dispatch table
    method_handlers = containers.Map();
    method_handlers('adc_threshold') = @handle_adc_threshold;
    method_handlers('d_threshold') = @handle_d_threshold;
    method_handlers('df_intersection') = @handle_df_intersection;
    method_handlers('otsu') = @handle_otsu;
    method_handlers('gmm') = @handle_gmm;
    method_handlers('kmeans') = @handle_kmeans;
    method_handlers('region_growing') = @handle_region_growing;
    method_handlers('active_contours') = @handle_active_contours;
    method_handlers('percentile') = @handle_percentile;
    method_handlers('spectral') = @handle_spectral;
    method_handlers('fdm') = @handle_fdm;

    % Execute the selected method
    handler = method_handlers(valid_params.core_method);
    core_mask = handler(valid_params, adc_vec, d_vec, f_vec, dstar_vec, has_3d, gtv_mask_3d, opts);
end

function [valid_params, fallback_mask] = validate_and_prepare(config_struct, adc_vec, d_vec, f_vec, dstar_vec, gtv_mask_3d, opts)
    % Common parameter validation and preparation
    
    % Validate required config field
    if ~isfield(config_struct, 'core_method')
        error('extract_tumor_core:missingField', ...
            'config_struct must contain a ''core_method'' field.');
    end

    VALID_METHODS = {'adc_threshold', 'd_threshold', 'df_intersection', ...
        'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
        'percentile', 'spectral', 'fdm'};
    core_method = lower(config_struct.core_method);
    if ~ismember(core_method, VALID_METHODS)
        error('extract_tumor_core:invalidMethod', ...
            'Unrecognized core_method ''%s''. Valid options: %s', ...
            core_method, strjoin(VALID_METHODS, ', '));
    end

    n_voxels = length(adc_vec);
    fallback_mask = false(n_voxels, 1);

    % Check for early return conditions
    if sum(~isnan(adc_vec)) == 0
        valid_params = struct('should_return_fallback', true);
        return;
    end

    % Prepare common parameters
    valid_params = struct();
    valid_params.core_method = core_method;
    valid_params.n_voxels = n_voxels;
    valid_params.should_return_fallback = false;
    
    % Copy relevant config fields
    if isfield(config_struct, 'adc_thresh'), valid_params.adc_thresh = config_struct.adc_thresh; end
    if isfield(config_struct, 'd_thresh'), valid_params.d_thresh = config_struct.d_thresh; end
    if isfield(config_struct, 'f_thresh'), valid_params.f_thresh = config_struct.f_thresh; end
    if isfield(config_struct, 'min_vox_hist'), valid_params.min_vox_hist = config_struct.min_vox_hist; end
    if isfield(config_struct, 'core_percentile'), valid_params.core_percentile = config_struct.core_percentile; end
    if isfield(config_struct, 'core_n_clusters'), valid_params.core_n_clusters = config_struct.core_n_clusters; end
    if isfield(config_struct, 'spectral_min_voxels'), valid_params.spectral_min_voxels = config_struct.spectral_min_voxels; end
    if isfield(config_struct, 'fdm_parameter'), valid_params.fdm_parameter = config_struct.fdm_parameter; end
    if isfield(config_struct, 'fdm_thresh'), valid_params.fdm_thresh = config_struct.fdm_thresh; end
    
    % Load morphological processing parameters with adaptive sizing
    valid_params.morphology = get_morphology_params(config_struct, gtv_mask_3d, opts);
end

function morphology_params = get_morphology_params(config_struct, gtv_mask_3d, opts)
    % Extract morphological processing parameters with adaptive sizing
    
    % Default configuration parameters
    defaults = struct();
    defaults.erosion_disk_radius = 2;
    defaults.min_region_size = 50;
    defaults.active_contour_iterations = 50;
    defaults.region_growing_std_multiplier = 2.0;
    defaults.enable_adaptive_sizing = true;
    defaults.base_voxel_size_mm = 1.0;  % Reference voxel size for scaling
    defaults.base_tumor_volume_mm3 = 8000;  % Reference tumor volume for scaling
    
    % Override with config values
    morphology_params = defaults;
    if isfield(config_struct, 'erosion_disk_radius')
        morphology_params.erosion_disk_radius = config_struct.erosion_disk_radius;
    end
    if isfield(config_struct, 'min_region_size')
        morphology_params.min_region_size = config_struct.min_region_size;
    end
    if isfield(config_struct, 'active_contour_iterations')
        morphology_params.active_contour_iterations = config_struct.active_contour_iterations;
    end
    if isfield(config_struct, 'region_growing_std_multiplier')
        morphology_params.region_growing_std_multiplier = config_struct.region_growing_std_multiplier;
    end
    if isfield(config_struct, 'enable_adaptive_sizing')
        morphology_params.enable_adaptive_sizing = config_struct.enable_adaptive_sizing;
    end
    if isfield(config_struct, 'base_voxel_size_mm')
        morphology_params.base_voxel_size_mm = config_struct.base_voxel_size_mm;
    end
    if isfield(config_struct, 'base_tumor_volume_mm3')
        morphology_params.base_tumor_volume_mm3 = config_struct.base_tumor_volume_mm3;
    end
    
    % Apply adaptive sizing if enabled
    if morphology_params.enable_adaptive_sizing
        morphology_params = apply_adaptive_sizing(morphology_params, gtv_mask_3d, opts);
    end
end

function morphology_params = apply_adaptive_sizing(morphology_params, gtv_mask_3d, opts)
    % Apply adaptive sizing based on voxel dimensions and tumor volume
    
    % Get voxel dimensions
    voxel_dims = [1.0, 1.0, 1.0];  % Default 1mm isotropic
    if isfield(opts, 'voxel_dims') && ~isempty(opts.voxel_dims)
        voxel_dims = opts.voxel_dims;
    end
    
    % Calculate effective voxel size (geometric mean of dimensions)
    effective_voxel_size = (voxel_dims(1) * voxel_dims(2) * voxel_dims(3))^(1/3);
    
    % Get tumor volume
    tumor_volume_mm3 = morphology_params.base_tumor_volume_mm3;  % Default
    if isfield(opts, 'tumor_volume_mm3') && ~isempty(opts.tumor_volume_mm3)
        tumor_volume_mm3 = opts.tumor_volume_mm3;
    elseif ~isempty(gtv_mask_3d)
        % Estimate volume from mask if not provided
        voxel_volume = prod(voxel_dims);
        tumor_volume_mm3 = sum(gtv_mask_3d(:)) * voxel_volume;
    end
    
    % Calculate scaling factors
    voxel_scale_factor = effective_voxel_size / morphology_params.base_voxel_size_mm;
    volume_scale_factor = (tumor_volume_mm3 / morphology_params.base_tumor_volume_mm3)^(1/3);
    
    % Combined scaling factor (geometric mean)
    combined_scale_factor = sqrt(voxel_scale_factor * volume_scale_factor);
    
    % Apply scaling with reasonable bounds
    min_scale = 0.5;
    max_scale = 3.0;
    combined_scale_factor = max(min_scale, min(max_scale, combined_scale_factor));
    
    % Scale parameters
    morphology_params.erosion_disk_radius = max(1, round(morphology_params.erosion_disk_radius * combined_scale_factor));
    morphology_params.min_region_size = max(10, round(morphology_params.min_region_size * combined_scale_factor^3));
    morphology_params.active_contour_iterations = max(10, round(morphology_params.active_contour_iterations * combined_scale_factor));
    
    % Store scaling info for debugging
    morphology_params.applied_scale_factor = combined_scale_factor;
    morphology_params.voxel_scale_factor = voxel_scale_factor;
    morphology_params.volume_scale_factor = volume_scale_factor;
end

function core_mask = apply_fallback_threshold(valid_params, adc_vec)
    % Common fallback to ADC threshold
    core_mask = adc_vec <= valid_params.adc_thresh;
end

function core_mask = handle_adc_threshold(valid_params, adc_vec, ~, ~, ~, ~, ~, ~)
    core_mask = adc_vec <= valid_params.adc_thresh;
end

function core_mask = handle_d_threshold(valid_params, adc_vec, d_vec, ~, ~, ~, ~, ~)
    if isempty(d_vec) || sum(~isnan(d_vec)) == 0
        warning('extract_tumor_core:noDValues', 'D values missing. Falling back to ADC threshold.');
        core_mask = apply_fallback_threshold(valid_params, adc_vec);
    else
        core_mask = d_vec <= valid_params.d_thresh;
    end
end

function core_mask = handle_df_intersection(valid_params, adc_vec, d_vec, f_vec, ~, ~, ~, ~)
    if isempty(d_vec) || sum(~isnan(d_vec)) == 0 || isempty(f_vec) || sum(~isnan(f_vec)) == 0
        warning('extract_tumor_core:noIVIMValues', 'IVIM D/f values missing. Falling back to ADC threshold.');
        core_mask = apply_fallback_threshold(valid_params, adc_vec);
    else
        core_mask = (d_vec <= valid_params.d_thresh) & (f_vec <= valid_params.f_thresh);
    end
end

function core_mask = handle_otsu(valid_params, adc_vec, ~, ~, ~, ~, ~, ~)
    valid_idx = find(~isnan(adc_vec));
    core_mask = false(valid_params.n_voxels, 1);
    
    if isempty(valid_idx)
        return;
    end
    
    valid_vals = adc_vec(valid_idx);
    if length(valid_vals) > valid_params.min_vox_hist
        try
            level = multithresh(valid_vals, 1);
            core_mask = adc_vec <= level;
        catch
            warning('extract_tumor_core:otsuFailed', 'Otsu method failed. Falling back to ADC threshold.');
            core_mask = apply_fallback_threshold(valid_params, adc_vec);
        end
    else
        core_mask = apply_fallback_threshold(valid_params, adc_vec);
    end
end

function core_mask = handle_gmm(valid_params, adc_vec, ~, ~, ~, ~, ~, ~)
    valid_idx = find(~isnan(adc_vec));
    core_mask = false(valid_params.n_voxels, 1);
    
    if isempty(valid_idx)
        return;
    end
    
    valid_vals = adc_vec(valid_idx);
    if length(valid_vals) > valid_params.min_vox_hist
        try
            if exist('OCTAVE_VERSION', 'builtin')
                % Octave lacks fitgmdist; fall back to k-means
                [idx, C] = kmeans(valid_vals, 2);
                [~, min_cluster] = min(C);
                core_mask_valid = (idx == min_cluster);
            else
                gm = fitgmdist(valid_vals, 2, 'Replicates', 3, 'Options', statset('MaxIter', 500));
                idx = cluster(gm, valid_vals);
                [~, min_c_idx] = min(gm.mu);
                core_mask_valid = (idx == min_c_idx);
            end
            core_mask(valid_idx) = core_mask_valid;
        catch
            warning('extract_tumor_core:gmmFailed', 'GMM failed. Falling back to ADC threshold.');
            core_mask = apply_fallback_threshold(valid_params, adc_vec);
        end
    else
        core_mask = apply_fallback_threshold(valid_params, adc_vec);
    end
end

function core_mask = handle_kmeans(valid_params, adc_vec, ~, ~, ~, ~, ~, ~)
    valid_idx = find(~isnan(adc_vec));
    core_mask = false(valid_params.n_voxels, 1);
    
    if isempty(valid_idx)
        return;
    end
    
    valid_vals = adc_vec(valid_idx);
    if length(valid_vals) > valid_params.min_vox_hist
        try
            [idx, C] = kmeans(valid_vals, 2, 'Replicates', 3);
            [~, min_cluster] = min(C);
            core_mask(valid_idx) = (idx == min_cluster);
        catch
            warning('extract_tumor_core:kmeansFailed', 'K-Means failed. Falling back to ADC threshold.');
            core_mask = apply_fallback_threshold(valid_params, adc_vec);
        end
    else
        core_mask = apply_fallback_threshold(valid_params, adc_vec);
    end
end

function core_mask = handle_region_growing(valid_params, adc_vec, ~, ~, ~, has_3d, gtv_mask_3d, ~)
    if ~has_3d || isempty(gtv_mask_3d)
        warning('extract_tumor_core:no3DForRegionGrowing', '3D mask required for Region Growing. Falling back to ADC threshold.');
        core_mask = apply_fallback_threshold(valid_params, adc_vec);
        return;
    end
    
    valid_idx = find(~isnan(adc_vec));
    if isempty(valid_idx)
        core_mask = false(valid_params.n_voxels, 1);
        return;
    end
    
    % Reconstruct the 3D ADC map inside the GTV
    adc_map_3d = nan(size(gtv_mask_3d));
    adc_map_3d(gtv_mask_3d == 1) = adc_vec;
    
    % Find absolute minimum ADC (the most restricted voxel as seed)
    [min_val, min_idx] = min(adc_map_3d(:));
    
    if isnan(min_val)
        core_mask = apply_fallback_threshold(valid_params, adc_vec);
    else
        [sx, sy, sz] = ind2sub(size(adc_map_3d), min_idx);
        
        sorted_adc = sort(adc_vec(valid_idx));
        bot_10_idx = max(1, round(0.1 * length(sorted_adc)));
        bot_10 = sorted_adc(1:bot_10_idx);
        if length(bot_10) < 2
            tol_upper = valid_params.adc_thresh;
        else
            tol_upper = mean(bot_10) + std(bot_10) * valid_params.morphology.region_growing_std_multiplier;
        end
        if isnan(tol_upper) || tol_upper == 0
            tol_upper = valid_params.adc_thresh;
        end
        
        % BFS-based region growing from the seed voxel
        rg_mask = false(size(adc_map_3d));
        rg_mask(sx, sy, sz) = true;
        n_gtv_vox = sum(gtv_mask_3d(:));
        queue = zeros(n_gtv_vox, 3);
        queue(1,:) = [sx, sy, sz];
        q_head = 1;
        q_tail = 1;

        shifts = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];

        while q_head <= q_tail
            cx = queue(q_head,1); cy = queue(q_head,2); cz = queue(q_head,3);
            q_head = q_head + 1;

            for i=1:6
                nx = cx + shifts(i,1);
                ny = cy + shifts(i,2);
                nz = cz + shifts(i,3);

                if nx >= 1 && nx <= size(adc_map_3d,1) && ...
                   ny >= 1 && ny <= size(adc_map_3d,2) && ...
                   nz >= 1 && nz <= size(adc_map_3d,3)

                    if ~rg_mask(nx, ny, nz) && ~isnan(adc_map_3d(nx,ny,nz)) ...
                       && gtv_mask_3d(nx,ny,nz) == 1

                        if adc_map_3d(nx,ny,nz) <= tol_upper
                            rg_mask(nx,ny,nz) = true;
                            q_tail = q_tail + 1;
                            queue(q_tail,:) = [nx, ny, nz];
                        end
                    end
                end
            end
        end
        
        % Apply morphological post-processing to remove small regions
        if sum(rg_mask(:)) > valid_params.morphology.min_region_size
            rg_mask = apply_morphological_cleanup(rg_mask, valid_params.morphology);
        end
        
        % Map 3D region-grown mask back to 1D voxel vector
        core_mask = rg_mask(gtv_mask_3d == 1);
    end
end

function core_mask = handle_active_contours(valid_params, adc_vec, ~, ~, ~, has_3d, gtv_mask_3d, ~)
    if ~has_3d || isempty(gtv_mask_3d)
        warning('extract_tumor_core:no3DForActiveContours', '3D mask required for Active Contours. Falling back to ADC threshold.');
        core_mask = apply_fallback_threshold(valid_params, adc_vec);
        return;
    end
    
    if exist('OCTAVE_VERSION', 'builtin')
        warning('extract_tumor_core:octaveContours', 'Active contours not supported in Octave natively. Falling back to Otsu.');
        try
            valid_vals = adc_vec(~isnan(adc_vec));
            if length(valid_vals) > valid_params.min_vox_hist
                level = multithresh(valid_vals, 1);
                core_mask = adc_vec <= level;
            else
                core_mask = apply_fallback_threshold(valid_params, adc_vec);
            end
        catch
            core_mask = apply_fallback_threshold(valid_params, adc_vec);
        end
    else
        % Reconstruct 3D ADC map, normalize for activecontour
        adc_map_3d = nan(size(gtv_mask_3d));
        adc_map_3d(gtv_mask_3d == 1) = adc_vec;
        
        init_mask = false(size(gtv_mask_3d));
        init_mask(adc_map_3d <= valid_params.adc_thresh) = true;
        
        max_adc = max(adc_vec);
        if isnan(max_adc) || max_adc == 0
            warning('extract_tumor_core:nanMaxAdc', 'All-NaN or zero ADC values. Falling back to ADC threshold.');
            core_mask = apply_fallback_threshold(valid_params, adc_vec);
            return;
        end
        norm_img = zeros(size(gtv_mask_3d));
        norm_img(gtv_mask_3d == 1) = max_adc - adc_vec;
        norm_max = max(norm_img(:));
        if norm_max == 0 || isnan(norm_max)
            warning('extract_tumor_core:constantAdc', 'Constant ADC within GTV. Falling back to ADC threshold.');
            core_mask = apply_fallback_threshold(valid_params, adc_vec);
            return;
        end
        norm_img = norm_img / norm_max;

        try
            ac_mask = activecontour(norm_img, init_mask, valid_params.morphology.active_contour_iterations, 'Chan-Vese');
            
            % Apply morphological post-processing
            ac_mask = apply_morphological_cleanup(ac_mask, valid_params.morphology);
            
            core_mask = ac_mask(gtv_mask_3d == 1);
        catch
            warning('extract_tumor_core:activeContourFailed', 'Active contour failed. Falling back to ADC threshold.');
            core_mask = apply_fallback_threshold(valid_params, adc_vec);
        end
    end
end

function core_mask = handle_percentile(valid_params, adc_vec, ~, ~, ~, ~, ~, ~)
    pct = valid_params.core_percentile;
    valid_adc = adc_vec(~isnan(adc_vec));
    if isempty(valid_adc)
        core_mask = false(valid_params.n_voxels, 1);
        return;
    end
    pct_thresh = prctile(valid_adc, pct);
    pct_thresh = min(pct_thresh, valid_params.adc_thresh);
    core_mask = adc_vec <= pct_thresh;
end

function core_mask = handle_spectral(valid_params, adc_vec, d_vec, f_vec, dstar_vec, ~, ~, ~)
    % Build feature matrix from available parameters
    features = adc_vec;
    if ~isempty(d_vec) && sum(~isnan(d_vec)) > 0
        features = [features, d_vec];
    end
    if ~isempty(f_vec) && sum(~isnan(f_vec)) > 0
        features = [features, f_vec];
    end
    if ~isempty(dstar_vec) && sum(~isnan(dstar_vec)) > 0
        features = [features, dstar_vec];
    end

    % Only keep voxels where all features are finite
    valid_idx = find(all(~isnan(features), 2));
    spectral_min = valid_params.min_vox_hist;
    if isfield(valid_params, 'spectral_min_voxels')
        spectral_min = valid_params.spectral_min_voxels;
    end
    
    core_mask = false(valid_params.n_voxels, 1);
    
    if numel(valid_idx) < spectral_min
        warning('extract_tumor_core:tooFewForSpectral', ...
            'Too few valid voxels (%d) for spectral clustering. Falling back to ADC threshold.', ...
            numel(valid_idx));
        core_mask = apply_fallback_threshold(valid_params, adc_vec);
        return;
    end

    n_clusters = valid_params.core_n_clusters;
    valid_features = features(valid_idx, :);

    % Z-score normalize
    feat_mean = mean(valid_features, 1);
    feat_std = std(valid_features, 0, 1);
    feat_std(feat_std == 0) = 1;
    valid_features_z = (valid_features - feat_mean) ./ feat_std;

    try
        if exist('spectralcluster', 'file')
            idx = spectralcluster(valid_features_z, n_clusters);
        else
            warning('extract_tumor_core:noSpectralCluster', ...
                'spectralcluster not available (requires R2019b+). Using k-means fallback.');
            [idx, ~] = kmeans(valid_features_z, n_clusters, 'Replicates', 5);
        end

        % Core cluster = the one with lowest mean ADC (first column)
        cluster_means = zeros(n_clusters, 1);
        for ci = 1:n_clusters
            cluster_means(ci) = mean(valid_features(idx == ci, 1));
        end
        [~, core_cluster] = min(cluster_means);
        core_mask(valid_idx) = (idx == core_cluster);
    catch ME
        warning('extract_tumor_core:spectralFailed', ...
            'Spectral clustering failed: %s. Falling back to ADC threshold.', ME.message);
        core_mask = apply_fallback_threshold(valid_params, adc_vec);
    end
end

function core_mask = handle_fdm(valid_params, adc_vec, d_vec, ~, ~, ~, ~, opts)
    % Timepoint index: 1 = baseline (Fx1), 2+ = follow-up fractions
    k_idx = 0;
    if isfield(opts, 'timepoint_index')
        k_idx = opts.timepoint_index;
    end

    if k_idx <= 1
        warning('extract_tumor_core:fdmBaseline', ...
            'fDM requires post-baseline timepoint. Falling back to ADC threshold for baseline.');
        core_mask = apply_fallback_threshold(valid_params, adc_vec);
        return;
    end

    % Select parameter for voxel-level change
    fdm_param = lower(valid_params.fdm_parameter);
    switch fdm_param
        case 'adc'
            current_vec = adc_vec;
            if isfield(opts, 'baseline_adc_vec')
                baseline_vec = opts.baseline_adc_vec;
            else
                warning('extract_tumor_core:fdmNoBaseline', ...
                    'fDM: No baseline ADC vector provided. Falling back to ADC threshold.');
                core_mask = apply_fallback_threshold(valid_params, adc_vec);
                return;
            end
        case 'd'
            current_vec = d_vec;
            if isfield(opts, 'baseline_d_vec')
                baseline_vec = opts.baseline_d_vec;
            else
                warning('extract_tumor_core:fdmNoBaseline', ...
                    'fDM: No baseline D vector provided. Falling back to ADC threshold.');
                core_mask = apply_fallback_threshold(valid_params, adc_vec);
                return;
            end
        otherwise
            error('extract_tumor_core:invalidFdmParam', ...
                'fdm_parameter must be ''adc'' or ''d'', got ''%s''.', fdm_param);
    end

    fdm_significance = valid_params.fdm_thresh;
    if isfield(opts, 'repeatability_cor') && ~isnan(opts.repeatability_cor)
        fdm_significance = opts.repeatability_cor;
    end

    if isempty(baseline_vec) || numel(baseline_vec) ~= numel(current_vec)
        warning('extract_tumor_core:fdmVectorMismatch', ...
            'fDM: Baseline vector length (%d) does not match current (%d). Falling back to ADC threshold.', ...
            numel(baseline_vec), numel(current_vec));
        core_mask = apply_fallback_threshold(valid_params, adc_vec);
        return;
    end

    delta = current_vec - baseline_vec;
    core_mask = delta <= -fdm_significance;
end

function clean_mask = apply_morphological_cleanup(input_mask, morphology_params)
    % Apply morphological operations to clean up the mask
    
    clean_mask = input_mask;
    
    % Skip processing if mask is too small
    if sum(clean_mask(:)) < morphology_params.min_region_size
        return;
    end
    
    % Create structuring element for morphological operations
    se_radius = morphology_params.erosion_disk_radius;
    if ndims(input_mask) == 3
        se = strel('sphere', se_radius);
    else
        se = strel('disk', se_radius);
    end
    
    try
        % Light erosion followed by dilation to remove small connections
        clean_mask = imerode(clean_mask, se);
        clean_mask = imdilate(clean_mask, se);
        
        % Remove small connected components
        if ndims(input_mask) == 3
            cc = bwconncomp(clean_mask, 26);
        else
            cc = bwconncomp(clean_mask, 8);
        end
        
        % Keep only regions above minimum size
        for i = 1:cc.NumObjects