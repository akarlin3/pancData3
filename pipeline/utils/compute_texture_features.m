function features = compute_texture_features(param_map, mask, n_levels, voxel_spacing, quantization_method, texture_3d)
% COMPUTE_TEXTURE_FEATURES  Extract first-order, GLCM, GLRLM, and shape features.
%
%   Computes first-order statistics, gray-level co-occurrence matrix (GLCM),
%   gray-level run-length matrix (GLRLM), and shape features from a 2D or
%   3D parameter map within a binary mask.
%
% IBSI Compliance Notes (Image Biomarker Standardisation Initiative, v11):
% -----------------------------------------------------------------------
%   This function implements a subset of IBSI-defined radiomics features.
%   The following describes compliance status per feature family:
%
%   IBSI-COMPLIANT features:
%     First-order: energy (IBSI id: N2K1), kurtosis (3RF5), skewness (88K1),
%       10th/90th percentile (QG58/GBY2), interquartile range (WL9B),
%       mean absolute deviation (D2ZX), robust MAD (1128).
%     GLCM: contrast (ACUI), correlation (NI2N), energy/ASM (8ZQL),
%       homogeneity/inverse difference moment (IB1Z).
%       Multi-offset averaging (0/45/90/135 degrees) for 2D rotational
%       invariance per IBSI recommendation (Section 3.4.2).
%     GLRLM: short run emphasis (22OV), long run emphasis (W4KF),
%       grey level non-uniformity (FP8K), run length non-uniformity (IC23),
%       run percentage (9ZK5). Multi-direction averaging at 4 angles (2D)
%       or 13 directions (3D, when texture_3d=true).
%     Shape: volume (RNU0), surface area (C0JK), sphericity (QCFX),
%       compactness (XPCE).
%
%   APPROXIMATIONS (deviations from strict IBSI definitions):
%     - Entropy: uses histogram-based Shannon entropy with equal-probability
%       bins (histcounts) rather than the IBSI intensity histogram entropy
%       definition which uses a fixed bin width scheme. Result varies with
%       quantization_method choice.
%     - Uniformity: histogram-based sum-of-squared-probabilities; equivalent
%       to IBSI only under fixed_bin_number quantization.
%     - Elongation: ratio of smallest to largest eigenvalue of the inertia
%       tensor. IBSI defines this as sqrt(lambda_minor/lambda_major); our
%       implementation matches this definition but uses the covariance
%       matrix rather than the inertia tensor, which is equivalent for
%       binary masks.
%     - GLCM is computed on a single 2D slice (largest mask cross-section)
%       rather than full 3D. This is standard for pancreatic DWI where
%       slice thickness (5-7mm) >> in-plane resolution (1-2mm), making
%       inter-slice co-occurrence physically less meaningful.
%     - GLRLM uses 2D (4 directions) or 3D (13 directions) depending on
%       the texture_3d parameter.
%
%   QUANTIZATION:
%     The quantization_method parameter controls how continuous parameter
%     values are discretised into grey levels for GLCM/GLRLM computation.
%     IBSI specifies both approaches and notes that the choice affects
%     feature values (Section 3.4.1):
%       'fixed_bin_number' (default): Rescales to [1, n_levels] using the
%         min/max of the ROI. Number of bins = n_levels. Suitable when
%         comparing texture across patients with different value ranges.
%       'fixed_bin_width': Bins are spaced at fixed intervals of
%         (range / n_levels). Number of bins varies with value range.
%         Suitable when absolute parameter values are meaningful (e.g.,
%         ADC in mm^2/s). Bin width = (max - min) / n_levels.
%
% Inputs:
%   param_map            - 2D or 3D numeric array (parameter map, e.g., ADC)
%   mask                 - Binary mask (same size as param_map)
%   n_levels             - Number of quantization levels (default: 32)
%   voxel_spacing        - [dx dy dz] voxel spacing in mm for physical-unit
%                          shape features (default: [1 1 1])
%   quantization_method  - 'fixed_bin_number' (default) or 'fixed_bin_width'
%   texture_3d           - Use 3D GLRLM (13 directions) for 3D volumes
%                          (default: true). When false, uses 2D GLRLM on
%                          the largest mask slice.
%
% Outputs:
%   features       - Struct with named texture feature fields:
%                    First-order: energy, uniformity, entropy, kurtosis, skewness,
%                                 p10, p90, iqr, mad, rmad
%                    GLCM: glcm_contrast, glcm_correlation, glcm_energy, glcm_homogeneity
%                    GLRLM: glrlm_sre, glrlm_lre, glrlm_gln, glrlm_rln, glrlm_rp
%                    Shape: shape_volume, shape_surface_area, shape_sphericity,
%                           shape_elongation, shape_compactness

    if nargin < 3 || isempty(n_levels), n_levels = 32; end
    if nargin < 4 || isempty(voxel_spacing), voxel_spacing = [1 1 1]; end
    % Backward compatibility: 5th arg may be a logical (legacy texture_3d usage)
    if nargin >= 5 && islogical(quantization_method)
        texture_3d = quantization_method;
        quantization_method = 'fixed_bin_number';
    end
    if nargin < 5 || isempty(quantization_method), quantization_method = 'fixed_bin_number'; end
    if nargin < 6 || isempty(texture_3d), texture_3d = true; end

    % Extract masked voxels and validate input data
    vals = double(param_map(mask > 0));
    
    % Remove NaN/Inf values and handle negative values
    vals = vals(isfinite(vals));
    if isempty(vals)
        % Handle case where all values are NaN/Inf
        features = initialize_nan_features();
        return;
    end
    
    % Handle negative values by taking absolute value for GLRLM quantization
    has_negative = any(vals < 0);
    vals_abs = abs(vals);

    % Initialize output struct
    features = struct();

    % NaN defaults for all features
    nan_fields = {'energy', 'uniformity', 'entropy', 'kurtosis', 'skewness', ...
        'p10', 'p90', 'iqr', 'mad', 'rmad', ...
        'glcm_contrast', 'glcm_correlation', 'glcm_energy', 'glcm_homogeneity', ...
        'glrlm_sre', 'glrlm_lre', 'glrlm_gln', 'glrlm_rln', 'glrlm_rp', ...
        'shape_volume', 'shape_surface_area', 'shape_sphericity', ...
        'shape_elongation', 'shape_compactness'};

    if length(vals) < 2
        for fi = 1:numel(nan_fields)
            features.(nan_fields{fi}) = NaN;
        end
        return;
    end

    % --- First-order features (use original values including negative) ---
    features.energy = sum(vals.^2);
    features.kurtosis = kurtosis(vals);
    features.skewness = skewness(vals);
    features.p10 = prctile(vals, 10);
    features.p90 = prctile(vals, 90);
    features.iqr = prctile(vals, 75) - prctile(vals, 25);
    features.mad = mean(abs(vals - mean(vals)));

    % Robust mean absolute deviation (10th-90th percentile)
    p10_val = features.p10;
    p90_val = features.p90;
    robust_vals = vals(vals >= p10_val & vals <= p90_val);
    if ~isempty(robust_vals)
        features.rmad = mean(abs(robust_vals - mean(robust_vals)));
    else
        features.rmad = features.mad;
    end

    % Entropy (histogram-based)
    [counts, ~] = histcounts(vals, n_levels);
    prob = counts / sum(counts);
    prob_nz = prob(prob > 0);
    features.entropy = -sum(prob_nz .* log2(prob_nz));

    % Uniformity (sum of squared histogram bin probabilities)
    features.uniformity = sum(prob.^2);

    % --- Determine whether to use 3D GLRLM ---
    is_3d_volume = ndims(param_map) == 3 && size(mask, 3) > 1;
    use_3d_glrlm = texture_3d && is_3d_volume;

    % --- GLCM features (2D in-plane, always uses best slice for 3D) ---
    if ndims(param_map) == 3
        n_slices = size(mask, 3);
        slice_counts = zeros(n_slices, 1);
        for s = 1:n_slices
            slice_counts(s) = sum(sum(mask(:,:,s)));
        end
        [~, best_slice] = max(slice_counts);
        img_2d = double(param_map(:,:,best_slice));
        mask_2d = mask(:,:,best_slice);
    else
        img_2d = double(param_map);
        mask_2d = mask;
    end

    % Validate and clean the 2D image data
    img_2d_clean = validate_and_clean_image(img_2d, mask_2d);
    
    % Quantize to n_levels within mask
    masked_img = img_2d_clean;
    masked_img(~mask_2d) = NaN;
    valid_pixels = masked_img(isfinite(masked_img));

    has_quantized = false;
    if length(valid_pixels) >= 4
        min_val = min(valid_pixels);
        max_val = max(valid_pixels);
        if min_val == max_val
            features.glcm_contrast = 0;
            features.glcm_correlation = NaN;
            features.glcm_energy = 1;
            features.glcm_homogeneity = 1;
            % GLRLM: all runs of length=image_width for constant image
            features.glrlm_sre = NaN;
            features.glrlm_lre = NaN;
            features.glrlm_gln = NaN;
            features.glrlm_rln = NaN;
            features.glrlm_rp = NaN;
        else
            has_quantized = true;

            % Quantize using the selected method (IBSI Section 3.4.1)
            if strcmpi(quantization_method, 'fixed_bin_width')
                % Fixed bin width: bin_width = (max - min) / n_levels
                % Number of occupied bins may be <= n_levels
                bin_width = (max_val - min_val) / n_levels;
                quantized = floor((img_2d_clean - min_val) / bin_width) + 1;
                quantized(~mask_2d) = 0;
                actual_n_levels = max(quantized(mask_2d));
                quantized = max(1, min(actual_n_levels, quantized));
                n_levels_glcm = actual_n_levels;
            else
                % Fixed bin number (default): rescale to [1, n_levels]
                quantized = round((img_2d_clean - min_val) / (max_val - min_val) * (n_levels - 1)) + 1;
                quantized(~mask_2d) = 0;
                quantized = max(1, min(n_levels, quantized));
                n_levels_glcm = n_levels;
            end

            % Compute GLCM at 4 angles and average using built-in functions
            [glcm_contrast_val, glcm_correlation_val, glcm_energy_val, glcm_homogeneity_val] = ...
                compute_glcm_features(quantized, n_levels_glcm);
            
            features.glcm_contrast = glcm_contrast_val;
            features.glcm_correlation = glcm_correlation_val;
            features.glcm_energy = glcm_energy_val;
            features.glcm_homogeneity = glcm_homogeneity_val;

            % --- GLRLM features ---
            if use_3d_glrlm
                % 3D GLRLM: validate and quantize the full 3D volume
                param_map_clean = validate_and_clean_image_3d(param_map, mask);
                
                if strcmpi(quantization_method, 'fixed_bin_width')
                    all_vals = double(param_map_clean(mask > 0));
                    all_vals = all_vals(isfinite(all_vals));
                    min_val_3d = min(all_vals);
                    max_val_3d = max(all_vals);
                    bw = (max_val_3d - min_val_3d) / n_levels;
                    quantized_3d = floor((double(param_map_clean) - min_val_3d) / bw) + 1;
                    quantized_3d(~mask) = 0;
                    actual_n_3d = max(quantized_3d(mask));
                    quantized_3d = max(1, min(actual_n_3d, quantized_3d));
                    n_levels_3d = actual_n_3d;
                else
                    all_vals = double(param_map_clean(mask > 0));
                    all_vals = all_vals(isfinite(all_vals));
                    min_val_3d = min(all_vals);
                    max_val_3d = max(all_vals);
                    quantized_3d = round((double(param_map_clean) - min_val_3d) / (max_val_3d - min_val_3d) * (n_levels - 1)) + 1;
                    quantized_3d(~mask) = 0;
                    quantized_3d = max(1, min(n_levels, quantized_3d));
                    n_levels_3d = n_levels;
                end
                [features.glrlm_sre, features.glrlm_lre, features.glrlm_gln, ...
                    features.glrlm_rln, features.glrlm_rp] = compute_glrlm_3d(quantized_3d, mask, n_levels_3d);
            else
                [features.glrlm_sre, features.glrlm_lre, features.glrlm_gln, ...
                    features.glrlm_rln, features.glrlm_rp] = compute_glrlm(quantized, mask_2d, n_levels_glcm);
            end
        end
    else
        features.glcm_contrast = NaN;
        features.glcm_correlation = NaN;
        features.glcm_energy = NaN;
        features.glcm_homogeneity = NaN;
        features.glrlm_sre = NaN;
        features.glrlm_lre = NaN;
        features.glrlm_gln = NaN;
        features.glrlm_rln = NaN;
        features.glrlm_rp = NaN;
    end

    % --- Shape features from binary mask ---
    [features.shape_volume, features.shape_surface_area, features.shape_sphericity, ...
        features.shape_elongation, features.shape_compactness] = compute_shape_features(mask, voxel_spacing);
end


%% ===== Input validation and cleaning functions =====

function features = initialize_nan_features()
%INITIALIZE_NAN_FEATURES  Initialize feature struct with NaN values.
    nan_fields = {'energy', 'uniformity', 'entropy', 'kurtosis', 'skewness', ...
        'p10', 'p90', 'iqr', 'mad', 'rmad', ...
        'glcm_contrast', 'glcm_correlation', 'glcm_energy', 'glcm_homogeneity', ...
        'glrlm_sre', 'glrlm_lre', 'glrlm_gln', 'glrlm_rln', 'glrlm_rp', ...
        'shape_volume', 'shape_surface_area', 'shape_sphericity', ...
        'shape_elongation', 'shape_compactness'};
    
    features = struct();
    for fi = 1:numel(nan_fields)
        features.(nan_fields{fi}) = NaN;
    end
end


function img_clean = validate_and_clean_image(img, mask)
%VALIDATE_AND_CLEAN_IMAGE  Remove NaN/Inf values and handle negative values for 2D image.
    img_clean = double(img);
    
    % Replace NaN/Inf values outside mask with 0
    invalid_mask = ~isfinite(img_clean);
    img_clean(invalid_mask & ~mask) = 0;
    
    % For values inside mask, replace NaN/Inf with median of valid values
    if any(invalid_mask & mask)
        valid_vals = img_clean(mask & isfinite(img_clean));
        if ~isempty(valid_vals)
            replacement_val = median(valid_vals);
        else
            replacement_val = 0;
        end
        img_clean(invalid_mask & mask) = replacement_val;
    end
    
    % Handle negative values by taking absolute value for texture calculations
    % (preserves original values for first-order features computed earlier)
    img_clean = abs(img_clean);
end


function img_clean = validate_and_clean_image_3d(img, mask)
%VALIDATE_AND_CLEAN_IMAGE_3D  Remove NaN/Inf values and handle negative values for 3D image.
    img_clean = double(img);
    
    % Replace NaN/Inf values outside mask with 0
    invalid_mask = ~isfinite(img_clean);
    img_clean(invalid_mask & ~mask) = 0;
    
    % For values inside mask, replace NaN/Inf with median of valid values
    if any(invalid_mask(:) & mask(:))
        valid_vals = img_clean(mask & isfinite(img_clean));
        if ~isempty(valid_vals)
            replacement_val = median(valid_vals);
        else
            replacement_val = 0;
        end
        img_clean(invalid_mask & mask) = replacement_val;
    end
    
    % Handle negative values by taking absolute value for texture calculations
    img_clean = abs(img_clean);
end


%% ===== GLCM computation with built-in functions and fallback =====

function [contrast_val, correlation_val, energy_val, homogeneity_val] = compute_glcm_features(quantized, n_levels)
%COMPUTE_GLCM_FEATURES  Compute GLCM features using built-in graycomatrix/graycoprops with fallback.

    % Check if Image Processing Toolbox is available
    has_ipt = license('test', 'image_toolbox') && exist('graycomatrix', 'file') == 2 && exist('graycoprops', 'file') == 2;
    
    if has_ipt
        % Use vectorized built-in functions
        try
            offsets = [0 1; -1 1; -1 0; -1 -1];
            contrast_vals = zeros(4, 1);
            correlation_vals = zeros(4, 1);
            energy_vals = zeros(4, 1);
            homogeneity_vals = zeros(4, 1);
            valid_angles = false(4, 1);

            for a = 1:size(offsets, 1)
                try
                    glcm = graycomatrix(quantized, 'Offset', offsets(a,:), ...
                        'NumLevels', n_levels, 'GrayLimits', [1 n_levels], ...
                        'Symmetric', true);
                    props = graycoprops(glcm, {'Contrast', 'Correlation', 'Energy', 'Homogeneity'});
                    contrast_vals(a) = props.Contrast;
                    correlation_vals(a) = props.Correlation;
                    energy_vals(a) = props.Energy;
                    homogeneity_vals(a) = props.Homogeneity;
                    valid_angles(a) = true;
                catch ME_glcm
                    warning('compute_texture_features:glcmFailed', ...
                        'GLCM computation failed for angle %d: %s', a, ME_glcm.message);
                end
            end

            if any(valid_angles)
                contrast_val = mean(contrast_vals(valid_angles));
                correlation_val = mean(correlation_vals(valid_angles));
                energy_val = mean(energy_vals(valid_angles));
                homogeneity_val = mean(homogeneity_vals(valid_angles));
                return;
            end
        catch ME_ipt
            warning('compute_texture_features:iptFallback', ...
                'Image Processing Toolbox functions failed (%s), using manual implementation', ME_ipt.message);
        end
    end
    
    % Fallback: manual GLCM computation
    [contrast_val, correlation_val, energy_val, homogeneity_val] = compute_glcm_manual(quantized, n_levels);
end


%% ===== Manual GLCM computation (fallback) =====

function [contrast_val, correlation_val, energy_val, homogeneity_val] = compute_glcm_manual(quantized, n_levels)
%COMPUTE_GLCM_MANUAL  Manual GLCM computation for systems without Image Processing Toolbox.

    [nrows, ncols] = size(quantized);
    offsets = [0 1; -1 1; -1 0; -1 -1];
    
    contrast_vals = zeros(4, 1);
    correlation_vals = zeros(4, 1);
    energy_vals = zeros(4, 1);
    homogeneity_vals = zeros(4, 1);
    valid_angles = false(4, 1);

    for a = 1:size(offsets, 1)
        dr = offsets(a, 1);
        dc = offsets(a, 2);
        
        % Initialize GLCM matrix
        glcm = zeros(n_levels, n_levels);
        
        % Compute co-occurrences
        for r = 1:nrows
            for c = 1:ncols
                if quantized(r, c) == 0, continue; end
                
                r2 = r + dr;
                c2 = c + dc;
                
                if r2 >= 1 && r2 <= nrows && c2 >= 1 && c2 <= ncols && quantized(r2, c2) > 0
                    i = quantized(r, c);
                    j = quantized(r2, c2);
                    glcm(i, j) = glcm(i, j) + 1;
                    glcm(j, i) = glcm(j, i) + 1; % Symmetric
                end
            end
        end
        
        % Normalize GLCM
        total = sum(glcm(:));
        if total == 0
            continue;
        end
        valid_angles(a) = true;
        
        glcm_norm = glcm / total;
        
        % Compute features
        [i_idx, j_idx] = meshgrid(1:n_levels, 1:n_levels);
        
        % Contrast: sum of |i-j|^2 * P(i,j)
        contrast_vals(a) = sum(sum((i_idx - j_idx).^2 .* glcm_norm));
        
        % Energy: sum of P(i,j)^2
        energy_vals(a) = sum(sum(glcm_norm.^2));
        
        % Homogeneity: sum of P(i,j) / (1 + |i-j|)
        homogeneity_vals(a) = sum(sum(glcm_norm ./ (1 + abs(i_idx - j_idx))));
        
        % Correlation: more complex computation
        mu_i = sum(sum(i_idx .* glcm_norm));
        mu_j = sum(sum(j_idx .* glcm_norm));
        sigma_i = sqrt(sum(sum((i_idx - mu_i).^2 .* glcm_norm)));
        sigma_j = sqrt(sum(sum((j_idx - mu_j).^2 .* glcm_norm)));
        
        if sigma_i > 0 && sigma_j > 0
            correlation_vals(a) = sum(sum((i_idx - mu_i) .* (j_idx - mu_j) .* glcm_norm)) / (sigma_i * sigma_j);
        else
            correlation_vals(a) = NaN;
        end
    end

    if any(valid_angles)
        contrast_val = mean(contrast_vals(valid_angles));
        correlation_val = mean(correlation_vals(valid_angles));
        energy_val = mean(energy_vals(valid_angles));
        homogeneity_val = mean(homogeneity_vals(valid_angles));
    else
        contrast_val = NaN;
        correlation_val = NaN;
        energy_val = NaN;
        homogeneity_val = NaN;
    end
end


%% ===== 3D GLRLM computation (13 directions) =====

function [sre, lre, gln, rln, rp] = compute_glrlm_3d(quantized, mask, n_levels)
%COMPUTE_GLRLM_3D  Gray-Level Run-Length Matrix features in 13 3D directions, averaged.
%
%   Scans lines through the 3D quantized volume for each of 13 unique
%   directions: 3 axis-aligned, 6 face-diagonal, 4 body-diagonal.
%   For each direction, records run lengths per gray level and computes
%   SRE, LRE, GLN, RLN, RP.  Features are averaged across valid directions.

    sz = size(quantized);
    nrows = sz(1); ncols = sz(2); nslices = sz(3);
    max_run = max([nrows, ncols, nslices]);

    % 13 unique 3D directions: [dr, dc, ds]
    directions = [
        0  0  1;   % z-axis
        0  1  0;   % y-axis
        1  0  0;   % x-axis
        0  1  1;   % yz face-diagonal
        0  1 -1;   % yz face-diagonal
        1  0  1;   % xz face-diagonal
        1  0 -1;   % xz face-diagonal
        1  1  0;   % xy face-diagonal
        1 -1  0;   % xy face-diagonal
        1  1  1;   % body diagonal
        1  1 -1;   % body diagonal
        1 -1  1;   % body diagonal
        1 -1 -1;   % body diagonal
    ];

    % Pre-allocate feature arrays
    sre_vals = zeros(13, 1);
    lre_vals = zeros(13, 1);
    gln_vals = zeros(13, 1);
    rln_vals = zeros(13, 1);
    rp_vals = zeros(13, 1);
    valid_dirs = false(13, 1);
    
    n_pixels = sum(mask(:));

    for d = 1:size(directions, 1)
        dr = directions(d, 1);
        dc = directions(d, 2);
        ds = directions(d, 3);

        % Pre-allocate GLRLM matrix
        rlm = zeros(n_levels, max_run);
        total_runs = 0;

        % Find starting positions using vectorized approach
        starts = find_3d_starts_vectorized(nrows, ncols, nslices, dr, dc, ds);

        for si = 1:size(starts, 1)
            r = starts(si, 1);
            c = starts(si, 2);
            s = starts(si, 3);
            run_len = 0;
            run_val = 0;

            while r >= 1 && r <= nrows && c >= 1 && c <= ncols && s >= 1 && s <= nslices
                pix = quantized(r, c, s);
                in_mask = mask(r, c, s);

                if in_mask && pix > 0
                    if pix == run_val
                        run_len = run_len + 1;
                    else
                        if run_len > 0 && run_val > 0
                            rl = min(run_len, max_run);
                            rlm(run_val, rl) = rlm(run_val, rl) + 1;
                            total_runs = total_runs + 1;
                        end
                        run_val = pix;
                        run_len = 1;
                    end
                else
                    if run_len > 0 && run_val > 0
                        rl = min(run_len, max_run);
                        rlm(run_val, rl) = rlm(run_val, rl) + 1;
                        total_runs = total_runs + 1;
                    end
                    run_val = 0;
                    run_len = 0;
                end
                r = r + dr;
                c = c + dc;
                s = s + ds;
            end
            % End of line: record final run
            if run_len > 0 && run_val > 0
                rl = min(run_len, max_run);
                rlm(run_val, rl) = rlm(run_val, rl) + 1;
                total_runs = total_runs + 1;
            end
        end

        if total_runs == 0
            continue;
        end
        valid_dirs(d) = true;

        % Vectorized GLRLM feature computation
        n_r = total_runs;
        run_lengths = 1:max_run;
        
        % Short Run Emphasis - vectorized
        run_weights = 1 ./ (run_lengths.^2);
        col_sums = sum(rlm, 1);
        sre_vals(d) = sum(col_sums .* run_weights) / n_r;

        % Long Run Emphasis - vectorized
        run_weights_lre = run_lengths.^2;
        lre_vals(d) = sum(col_sums .* run_weights_lre) / n_r;

        % Gray-Level Non-Uniformity - vectorized
        gi_sums = sum(rlm, 2