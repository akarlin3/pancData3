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
%
% Outputs:
%   core_mask      - A 1D logical mask (same size as input vectors) where
%                    true indicates a voxel belongs to the core.

    if nargin < 8, opts = struct(); end

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
    % Default: no voxels belong to the core (conservative)
    core_mask = false(n_voxels, 1);

    % Fallback if ADC vector is entirely empty or NaN
    if sum(~isnan(adc_vec)) == 0
        return;
    end

    switch core_method
        case 'adc_threshold'
            % ADC (Apparent Diffusion Coefficient) reflects water mobility.
            % In densely-cellular tumor core, tightly-packed cells restrict
            % water diffusion → low ADC.  Necrotic or edematous regions
            % have high ADC due to unimpeded water motion.  A fixed scalar
            % threshold separates the restricted (core) from unrestricted
            % (margin/necrosis) compartments.
            core_mask = adc_vec <= config_struct.adc_thresh;

        case 'd_threshold'
            % True diffusion (D) from the IVIM biexponential model isolates
            % tissue cellularity more accurately than ADC, which conflates
            % true diffusion with microvascular perfusion (pseudo-diffusion).
            % By removing the perfusion contribution, D provides a purer
            % measure of cell-packing density.
            if isempty(d_vec) || sum(~isnan(d_vec)) == 0
                warning('extract_tumor_core:noDValues', 'D values missing. Falling back to ADC threshold.');
                core_mask = adc_vec <= config_struct.adc_thresh;
            else
                core_mask = d_vec <= config_struct.d_thresh;
            end

        case 'df_intersection'
            % Pancreatic ductal adenocarcinoma (PDAC) is characteristically
            % both hypocellular-stroma-rich (low D from dense desmoplasia)
            % and hypovascular (low f from poor microvasculature).  The
            % intersection of low-D AND low-f identifies the desmoplastic
            % core compartment, distinguishing it from well-perfused normal
            % pancreas (high f) or necrotic regions (high D).
            if isempty(d_vec) || sum(~isnan(d_vec)) == 0 || isempty(f_vec) || sum(~isnan(f_vec)) == 0
                warning('extract_tumor_core:noIVIMValues', 'IVIM D/f values missing. Falling back to ADC threshold.');
                core_mask = adc_vec <= config_struct.adc_thresh;
            else
                core_mask = (d_vec <= config_struct.d_thresh) & (f_vec <= config_struct.f_thresh);
            end

        case 'otsu'
            % Otsu's method finds the threshold that maximizes inter-class
            % variance, assuming the ADC histogram is bimodal (core + non-
            % core).  Data-driven: adapts to each patient's distribution
            % without requiring a fixed threshold.
            valid_idx = find(~isnan(adc_vec));
            if isempty(valid_idx)
                return;
            end
            valid_vals = adc_vec(valid_idx);
            if length(valid_vals) > config_struct.min_vox_hist
                try
                    level = multithresh(valid_vals, 1);
                    % Core is the darker class
                    core_mask = adc_vec <= level;
                catch
                    warning('extract_tumor_core:otsuFailed', 'Otsu method failed. Falling back to ADC threshold.');
                    core_mask = adc_vec <= config_struct.adc_thresh;
                end
            else
                core_mask = adc_vec <= config_struct.adc_thresh;
            end

        case 'gmm'
            % Fits a 2-component Gaussian mixture to the ADC distribution.
            % Unlike Otsu (hard threshold), GMM models each sub-population
            % (core/non-core) as a Gaussian, assigning voxels via posterior
            % probability.  The component with the lower mean (mu) is the
            % restricted-diffusion core.  More robust to skewed or
            % overlapping distributions than a single threshold.
            valid_idx = find(~isnan(adc_vec));
            if isempty(valid_idx)
                return;
            end
            valid_vals = adc_vec(valid_idx);
            if length(valid_vals) > config_struct.min_vox_hist
                try
                    if exist('OCTAVE_VERSION', 'builtin')
                        % Octave lacks fitgmdist; fall back to k-means
                        [idx, C] = kmeans(valid_vals, 2);
                        [~, min_cluster] = min(C);
                        core_mask_valid = (idx == min_cluster);
                    else
                        % 3 replicates with 500 max iterations guards against
                        % local optima in the EM algorithm
                        gm = fitgmdist(valid_vals, 2, 'Replicates', 3, 'Options', statset('MaxIter', 500));
                        idx = cluster(gm, valid_vals);
                        % The component with the smaller mean (mu) corresponds
                        % to the restricted-diffusion (core) population
                        [~, min_c_idx] = min(gm.mu);
                        core_mask_valid = (idx == min_c_idx);
                    end
                    core_mask(valid_idx) = core_mask_valid;
                catch
                    warning('extract_tumor_core:gmmFailed', 'GMM failed. Falling back to ADC threshold.');
                    core_mask = adc_vec <= config_struct.adc_thresh;
                end
            else
                core_mask = adc_vec <= config_struct.adc_thresh;
            end

        case 'kmeans'
            % K-means partitions ADC values into 2 clusters by minimizing
            % within-cluster variance.  The cluster with the lower centroid
            % represents restricted diffusion (core).  Simpler than GMM
            % (assumes spherical clusters), but faster and deterministic
            % with fixed random seed.
            valid_idx = find(~isnan(adc_vec));
            if isempty(valid_idx)
                return;
            end
            valid_vals = adc_vec(valid_idx);
            if length(valid_vals) > config_struct.min_vox_hist
                try
                    [idx, C] = kmeans(valid_vals, 2, 'Replicates', 3);
                    [~, min_cluster] = min(C);
                    core_mask(valid_idx) = (idx == min_cluster);
                catch
                    warning('extract_tumor_core:kmeansFailed', 'K-Means failed. Falling back to ADC threshold.');
                    core_mask = adc_vec <= config_struct.adc_thresh;
                end
            else
                core_mask = adc_vec <= config_struct.adc_thresh;
            end

        case 'region_growing'
            % Spatially-contiguous core: starts from the most restricted
            % voxel (lowest ADC = densest cellularity) and grows outward
            % through 6-connected neighbors, accepting voxels whose ADC
            % falls within a data-driven tolerance.  Enforces spatial
            % coherence — isolated low-ADC voxels (noise) are excluded
            % because they aren't reachable from the seed.
            if ~has_3d || isempty(gtv_mask_3d)
                warning('extract_tumor_core:no3DForRegionGrowing', '3D mask required for Region Growing. Falling back to ADC threshold.');
                core_mask = adc_vec <= config_struct.adc_thresh;
            else
                valid_idx = find(~isnan(adc_vec));
                if isempty(valid_idx), return; end
                
                % Reconstruct the 3D ADC map inside the GTV
                adc_map_3d = nan(size(gtv_mask_3d));
                adc_map_3d(gtv_mask_3d == 1) = adc_vec;
                
                % Find absolute minimum ADC (the most restricted voxel as seed)
                [min_val, min_idx] = min(adc_map_3d(:));
                
                if isnan(min_val)
                    core_mask = adc_vec <= config_struct.adc_thresh;
                else
                    [sx, sy, sz] = ind2sub(size(adc_map_3d), min_idx);
                    
                    % Adaptive tolerance derived from the most restricted
                    % 10% of voxels (the tail of the core distribution).
                    % mean + 2*SD of that tail captures ~95% of the core
                    % sub-population's variability while rejecting margin
                    % tissue whose ADC lies well above that bound.
                    sorted_adc = sort(adc_vec(valid_idx));
                    bot_10_idx = max(1, round(0.1 * length(sorted_adc)));
                    bot_10 = sorted_adc(1:bot_10_idx);
                    tol_upper = mean(bot_10) + std(bot_10) * 2;
                    if isnan(tol_upper) || tol_upper == 0
                        tol_upper = config_struct.adc_thresh;
                    end
                    
                    % BFS-based region growing from the seed voxel.
                    % Pre-allocate queue to the number of GTV voxels (upper
                    % bound on reachable voxels) to avoid dynamic growth
                    % inside the loop, which causes O(n^2) memory copies.
                    rg_mask = false(size(adc_map_3d));
                    rg_mask(sx, sy, sz) = true;
                    n_gtv_vox = sum(gtv_mask_3d(:));
                    queue = zeros(n_gtv_vox, 3);
                    queue(1,:) = [sx, sy, sz];
                    q_head = 1;  % next index to dequeue
                    q_tail = 1;  % last enqueued index

                    % 6-connected neighborhood shifts (face-adjacent only;
                    % excludes diagonal neighbors to prevent "leaking"
                    % through thin tissue boundaries)
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
                    
                    % Map 3D region-grown mask back to 1D voxel vector
                    % (same ordering as the GTV mask linearization)
                    core_mask = rg_mask(gtv_mask_3d == 1);
                end
            end

        case 'active_contours'
            % Level-set segmentation (Chan-Vese model) finds the contour
            % that best separates the image into two intensity-homogeneous
            % regions.  Unlike threshold methods, it produces spatially
            % smooth boundaries that respect tumor geometry.
            if ~has_3d || isempty(gtv_mask_3d)
                warning('extract_tumor_core:no3DForActiveContours', '3D mask required for Active Contours. Falling back to ADC threshold.');
                core_mask = adc_vec <= config_struct.adc_thresh;
            else
                if exist('OCTAVE_VERSION', 'builtin')
                    warning('extract_tumor_core:octaveContours', 'Active contours not supported in Octave natively. Falling back to Otsu.');
                    try
                        valid_vals = adc_vec(~isnan(adc_vec));
                        if length(valid_vals) > config_struct.min_vox_hist
                            level = multithresh(valid_vals, 1);
                            core_mask = adc_vec <= level;
                        else
                            core_mask = adc_vec <= config_struct.adc_thresh;
                        end
                    catch
                        core_mask = adc_vec <= config_struct.adc_thresh;
                    end
                else
                    % Reconstruct 3D ADC map, normalize for activecontour
                    adc_map_3d = nan(size(gtv_mask_3d));
                    adc_map_3d(gtv_mask_3d == 1) = adc_vec;
                    
                    % Initialize the level-set with a threshold-based seed mask.
                    % Chan-Vese evolves this initial contour toward the
                    % optimal two-region partition; the seed need not be
                    % perfect but should overlap with the core.
                    init_mask = false(size(gtv_mask_3d));
                    init_mask(adc_map_3d <= config_struct.adc_thresh) = true;
                    
                    % Invert contrast: Chan-Vese evolves toward bright
                    % regions.  Low ADC (core) is dark, so we invert so
                    % the core becomes the brightest region and the contour
                    % converges around it.
                    max_adc = max(adc_vec);

                    % Build a normalized inverted ADC image within the GTV
                    norm_img = zeros(size(gtv_mask_3d));
                    norm_img(gtv_mask_3d == 1) = max_adc - adc_vec;
                    norm_img = norm_img / max(norm_img(:));

                    try
                        % 50 iterations balances convergence against
                        % over-smoothing; Chan-Vese is chosen over edge-
                        % based models because ADC maps lack sharp edges.
                        ac_mask = activecontour(norm_img, init_mask, 50, 'Chan-Vese');
                        core_mask = ac_mask(gtv_mask_3d == 1);
                    catch
                        warning('extract_tumor_core:activeContourFailed', 'Active contour failed. Falling back to ADC threshold.');
                        core_mask = adc_vec <= config_struct.adc_thresh;
                    end
                end
            end

        case 'percentile'
            % Patient-adaptive threshold: Nth percentile of the patient's
            % own ADC distribution within the GTV.  Adapts to inter-patient
            % variability in baseline tissue properties and scanner
            % calibration, unlike a fixed threshold.
            pct = config_struct.core_percentile;
            valid_adc = adc_vec(~isnan(adc_vec));
            if isempty(valid_adc)
                return;
            end
            pct_thresh = prctile(valid_adc, pct);
            % Safety floor: cap at adc_thresh to prevent labelling normal
            % tissue as core in uniformly necrotic/high-ADC tumors.
            pct_thresh = min(pct_thresh, config_struct.adc_thresh);
            core_mask = adc_vec <= pct_thresh;

        case 'spectral'
            % Multi-parameter spectral clustering using available IVIM
            % parameters.  Discovers tissue sub-populations via the graph
            % Laplacian eigenspace, capturing non-linear relationships
            % between diffusion and perfusion parameters.

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
            spectral_min = config_struct.min_vox_hist;
            if isfield(config_struct, 'spectral_min_voxels')
                spectral_min = config_struct.spectral_min_voxels;
            end
            if numel(valid_idx) < spectral_min
                warning('extract_tumor_core:tooFewForSpectral', ...
                    'Too few valid voxels (%d) for spectral clustering. Falling back to ADC threshold.', ...
                    numel(valid_idx));
                core_mask = adc_vec <= config_struct.adc_thresh;
                return;
            end

            n_clusters = config_struct.core_n_clusters;  % typically 2 (core vs non-core)
            valid_features = features(valid_idx, :);

            % Z-score normalize so parameters on different scales
            % (ADC ~1e-3, f ~0.1) contribute equally
            feat_mean = mean(valid_features, 1);
            feat_std = std(valid_features, 0, 1);
            feat_std(feat_std == 0) = 1;  % prevent division by zero
            valid_features_z = (valid_features - feat_mean) ./ feat_std;

            try
                if exist('spectralcluster', 'file')
                    idx = spectralcluster(valid_features_z, n_clusters);
                else
                    % Fallback: k-means when spectralcluster unavailable
                    % (MATLAB < R2019b or Octave)
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
                core_mask = adc_vec <= config_struct.adc_thresh;
            end

        case 'fdm'
            % Functional Diffusion Maps: voxel-level change classification
            % between baseline and current timepoint.  Identifies the
            % treatment-resistant core as voxels where diffusivity
            % decreased (progressing = increased cellularity).

            % Timepoint index: 1 = baseline (Fx1), 2+ = follow-up fractions
            k_idx = 0;
            if isfield(opts, 'timepoint_index')
                k_idx = opts.timepoint_index;
            end

            if k_idx <= 1
                % fDM is undefined at baseline — requires delta from Fx1
                warning('extract_tumor_core:fdmBaseline', ...
                    'fDM requires post-baseline timepoint. Falling back to ADC threshold for baseline.');
                core_mask = adc_vec <= config_struct.adc_thresh;
                return;
            end

            % Select parameter for voxel-level change
            fdm_param = lower(config_struct.fdm_parameter);
            switch fdm_param
                case 'adc'
                    current_vec = adc_vec;
                    if isfield(opts, 'baseline_adc_vec')
                        baseline_vec = opts.baseline_adc_vec;
                    else
                        warning('extract_tumor_core:fdmNoBaseline', ...
                            'fDM: No baseline ADC vector provided. Falling back to ADC threshold.');
                        core_mask = adc_vec <= config_struct.adc_thresh;
                        return;
                    end
                case 'd'
                    current_vec = d_vec;
                    if isfield(opts, 'baseline_d_vec')
                        baseline_vec = opts.baseline_d_vec;
                    else
                        warning('extract_tumor_core:fdmNoBaseline', ...
                            'fDM: No baseline D vector provided. Falling back to ADC threshold.');
                        core_mask = adc_vec <= config_struct.adc_thresh;
                        return;
                    end
                otherwise
                    error('extract_tumor_core:invalidFdmParam', ...
                        'fdm_parameter must be ''adc'' or ''d'', got ''%s''.', fdm_param);
            end

            % Determine significance threshold: prefer repeatability-
            % derived CoR when available, else use config fallback
            fdm_significance = config_struct.fdm_thresh;
            if isfield(opts, 'repeatability_cor') && ~isnan(opts.repeatability_cor)
                fdm_significance = opts.repeatability_cor;
            end

            % Guard against vector length mismatch (different GTV contour
            % sizes across fractions after deformable registration)
            if isempty(baseline_vec) || numel(baseline_vec) ~= numel(current_vec)
                warning('extract_tumor_core:fdmVectorMismatch', ...
                    'fDM: Baseline vector length (%d) does not match current (%d). Falling back to ADC threshold.', ...
                    numel(baseline_vec), numel(current_vec));
                core_mask = adc_vec <= config_struct.adc_thresh;
                return;
            end

            % Voxel-wise change: positive delta = increased diffusivity (responding),
            % negative delta = decreased diffusivity (progressing/resistant)
            delta = current_vec - baseline_vec;

            % Classification:
            %   delta <= -threshold  →  "progressing" (ADC decreased = core)
            %   delta >=  threshold  →  "responding"  (ADC increased = non-core)
            %   |delta| < threshold  →  "stable"      (within noise floor)
            % Core = progressing voxels (treatment-resistant)
            core_mask = delta <= -fdm_significance;

        otherwise
            error('extract_tumor_core:invalidMethod', 'Unrecognized core_method: %s', config_struct.core_method);
    end
end
