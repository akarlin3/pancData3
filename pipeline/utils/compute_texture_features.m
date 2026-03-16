function features = compute_texture_features(param_map, mask, n_levels, voxel_spacing, quantization_method)
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
%       run percentage (9ZK5). Multi-direction averaging at 4 angles.
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
%     - GLRLM is computed per 2D slice (4 in-plane directions). Full 3D
%       GLRLM (13 directions) is not implemented.
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
    if nargin < 5 || isempty(quantization_method), quantization_method = 'fixed_bin_number'; end

    % Extract masked voxels
    vals = double(param_map(mask > 0));
    vals = vals(isfinite(vals));

    % Initialize output struct
    features = struct();

    % NaN defaults for all features
    nan_fields = {'energy', 'uniformity', 'entropy', 'kurtosis', 'skewness', ...
        'p10', 'p90', 'iqr', 'mad', 'rmad', ...
        'glcm_contrast', 'glcm_correlation', 'glcm_energy', 'glcm_homogeneity', ...
        'glrlm_sre', 'glrlm_lre', 'glrlm_gln', 'glrlm_rln', 'glrlm_rp', ...
        'shape_volume', 'shape_surface_area', 'shape_sphericity', ...
        'shape_elongation', 'shape_compactness'};

    if isempty(vals) || length(vals) < 2
        for fi = 1:numel(nan_fields)
            features.(nan_fields{fi}) = NaN;
        end
        return;
    end

    % --- First-order features ---
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

    % --- GLCM features (2D in-plane) ---
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

    % Quantize to n_levels within mask
    masked_img = img_2d;
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
                quantized = floor((img_2d - min_val) / bin_width) + 1;
                quantized(~mask_2d) = 0;
                actual_n_levels = max(quantized(mask_2d));
                quantized = max(1, min(actual_n_levels, quantized));
                n_levels_glcm = actual_n_levels;
            else
                % Fixed bin number (default): rescale to [1, n_levels]
                quantized = round((img_2d - min_val) / (max_val - min_val) * (n_levels - 1)) + 1;
                quantized(~mask_2d) = 0;
                quantized = max(1, min(n_levels, quantized));
                n_levels_glcm = n_levels;
            end

            % Compute GLCM at 4 angles and average
            offsets = [0 1; -1 1; -1 0; -1 -1];
            contrast_sum = 0; correlation_sum = 0;
            energy_sum = 0; homogeneity_sum = 0;
            n_valid_angles = 0;

            for a = 1:size(offsets, 1)
                try
                    glcm = graycomatrix(quantized, 'Offset', offsets(a,:), ...
                        'NumLevels', n_levels_glcm, 'GrayLimits', [1 n_levels_glcm], ...
                        'Symmetric', true);
                    props = graycoprops(glcm, {'Contrast', 'Correlation', 'Energy', 'Homogeneity'});
                    contrast_sum = contrast_sum + props.Contrast;
                    correlation_sum = correlation_sum + props.Correlation;
                    energy_sum = energy_sum + props.Energy;
                    homogeneity_sum = homogeneity_sum + props.Homogeneity;
                    n_valid_angles = n_valid_angles + 1;
                catch
                end
            end

            if n_valid_angles > 0
                features.glcm_contrast = contrast_sum / n_valid_angles;
                features.glcm_correlation = correlation_sum / n_valid_angles;
                features.glcm_energy = energy_sum / n_valid_angles;
                features.glcm_homogeneity = homogeneity_sum / n_valid_angles;
            else
                features.glcm_contrast = NaN;
                features.glcm_correlation = NaN;
                features.glcm_energy = NaN;
                features.glcm_homogeneity = NaN;
            end

            % --- GLRLM features ---
            [features.glrlm_sre, features.glrlm_lre, features.glrlm_gln, ...
                features.glrlm_rln, features.glrlm_rp] = compute_glrlm(quantized, mask_2d, n_levels_glcm);
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


%% ===== GLRLM computation (no toolbox dependency) =====

function [sre, lre, gln, rln, rp] = compute_glrlm(quantized, mask_2d, n_levels)
%COMPUTE_GLRLM  Gray-Level Run-Length Matrix features at 4 angles, averaged.
%
%   Scans rows of the quantized image for each of 4 directions:
%   0° (horizontal), 90° (vertical), 45° (diagonal), 135° (anti-diagonal).
%   For each direction, records run lengths per gray level.

    [nrows, ncols] = size(quantized);
    max_run = max(nrows, ncols);

    % Directions: [row_step, col_step] for scanning
    directions = [0 1; 1 0; 1 1; 1 -1];  % 0, 90, 45, 135 degrees

    sre_sum = 0; lre_sum = 0; gln_sum = 0; rln_sum = 0; rp_sum = 0;
    n_valid = 0;

    for d = 1:size(directions, 1)
        dr = directions(d, 1);
        dc = directions(d, 2);

        % Build GLRLM for this direction
        rlm = zeros(n_levels, max_run);
        total_runs = 0;

        % Determine starting positions
        if dr == 0 && dc == 1
            % Horizontal: start from each row, col 1
            starts_r = 1:nrows; starts_c = ones(1, nrows);
        elseif dr == 1 && dc == 0
            % Vertical: start from row 1, each col
            starts_r = ones(1, ncols); starts_c = 1:ncols;
        elseif dr == 1 && dc == 1
            % 45 degrees: start from left column and top row
            starts_r = [1:nrows, ones(1, ncols-1)];
            starts_c = [ones(1, nrows), 2:ncols];
        elseif dr == 1 && dc == -1
            % 135 degrees: start from right column and top row
            starts_r = [1:nrows, ones(1, ncols-1)];
            starts_c = [ncols*ones(1, nrows), (ncols-1):-1:1];
        end

        for si = 1:numel(starts_r)
            r = starts_r(si);
            c = starts_c(si);
            run_len = 0;
            run_val = 0;

            while r >= 1 && r <= nrows && c >= 1 && c <= ncols
                pix = quantized(r, c);
                in_mask = mask_2d(r, c);

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
        n_valid = n_valid + 1;

        % Compute GLRLM features from this direction's matrix
        n_r = total_runs;
        run_lengths = 1:max_run;

        % Short Run Emphasis: sum(rlm(i,j)/j^2) / n_r
        sre_d = 0;
        for j = 1:max_run
            sre_d = sre_d + sum(rlm(:, j)) / (j^2);
        end
        sre_sum = sre_sum + sre_d / n_r;

        % Long Run Emphasis: sum(rlm(i,j)*j^2) / n_r
        lre_d = 0;
        for j = 1:max_run
            lre_d = lre_d + sum(rlm(:, j)) * (j^2);
        end
        lre_sum = lre_sum + lre_d / n_r;

        % Gray-Level Non-Uniformity: sum(sum_j(rlm(i,j))^2) / n_r
        gi_sums = sum(rlm, 2);  % sum over run lengths per gray level
        gln_sum = gln_sum + sum(gi_sums.^2) / n_r;

        % Run-Length Non-Uniformity: sum(sum_i(rlm(i,j))^2) / n_r
        rj_sums = sum(rlm, 1);  % sum over gray levels per run length
        rln_sum = rln_sum + sum(rj_sums.^2) / n_r;

        % Run Percentage: n_r / n_pixels
        n_pixels = sum(mask_2d(:));
        rp_sum = rp_sum + n_r / n_pixels;
    end

    if n_valid > 0
        sre = sre_sum / n_valid;
        lre = lre_sum / n_valid;
        gln = gln_sum / n_valid;
        rln = rln_sum / n_valid;
        rp = rp_sum / n_valid;
    else
        sre = NaN; lre = NaN; gln = NaN; rln = NaN; rp = NaN;
    end
end


%% ===== Shape features =====

function [volume, sa, sphericity, elongation, compactness] = compute_shape_features(mask, voxel_spacing)
%COMPUTE_SHAPE_FEATURES  Shape features from a binary 3D mask.

    mask = logical(mask);
    n_voxels = sum(mask(:));

    if n_voxels < 1
        volume = NaN; sa = NaN; sphericity = NaN;
        elongation = NaN; compactness = NaN;
        return;
    end

    dx = voxel_spacing(1);
    dy = voxel_spacing(2);
    dz = 1;
    if numel(voxel_spacing) >= 3
        dz = voxel_spacing(3);
    end

    % Volume (physical units)
    voxel_vol = dx * dy * dz;
    volume = n_voxels * voxel_vol;

    % Surface area via marching cubes (isosurface)
    sa = NaN;
    try
        if ndims(mask) == 3 && all(size(mask) >= 2)
            % Pad to avoid edge artifacts
            padded = false(size(mask) + 2);
            padded(2:end-1, 2:end-1, 2:end-1) = mask;
            fv = isosurface(padded, 0.5);
            if ~isempty(fv.vertices) && ~isempty(fv.faces)
                % Scale vertices to physical coordinates
                v = fv.vertices;
                v(:,1) = v(:,1) * dx;
                v(:,2) = v(:,2) * dy;
                v(:,3) = v(:,3) * dz;
                % Triangle areas
                v1 = v(fv.faces(:,1), :);
                v2 = v(fv.faces(:,2), :);
                v3 = v(fv.faces(:,3), :);
                cross_prod = cross(v2 - v1, v3 - v1, 2);
                tri_areas = 0.5 * sqrt(sum(cross_prod.^2, 2));
                sa = sum(tri_areas);
            end
        elseif ndims(mask) == 2
            % 2D: perimeter approximation
            perim = bwperim(mask);
            sa = sum(perim(:)) * dx;
        end
    catch
        sa = NaN;
    end

    % Sphericity: (36*pi*V^2 / SA^3)^(1/3)
    if isfinite(sa) && sa > 0 && isfinite(volume) && volume > 0
        sphericity = (36 * pi * volume^2 / sa^3)^(1/3);
    else
        sphericity = NaN;
    end

    % Elongation: ratio of eigenvalues of inertia tensor
    elongation = NaN;
    try
        [r, c, s] = ind2sub(size(mask), find(mask));
        coords = [r * dx, c * dy, s * dz];
        if size(coords, 1) >= 3
            coords_centered = coords - mean(coords, 1);
            cov_mat = (coords_centered' * coords_centered) / size(coords_centered, 1);
            eigenvalues = sort(eig(cov_mat), 'descend');
            eigenvalues = eigenvalues(eigenvalues > 0);
            if numel(eigenvalues) >= 2
                elongation = sqrt(eigenvalues(end) / eigenvalues(1));
            end
        end
    catch
        elongation = NaN;
    end

    % Compactness: V / (sqrt(pi) * SA^1.5)
    if isfinite(sa) && sa > 0 && isfinite(volume) && volume > 0
        compactness = volume / (sqrt(pi) * sa^1.5);
    else
        compactness = NaN;
    end
end
