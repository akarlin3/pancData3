function features = compute_texture_features(param_map, mask, n_levels)
% COMPUTE_TEXTURE_FEATURES  Extract first-order and GLCM texture features.
%
%   Computes first-order statistics and gray-level co-occurrence matrix
%   (GLCM) features from a 2D or 3D parameter map within a binary mask.
%
% Inputs:
%   param_map  - 2D or 3D numeric array (parameter map, e.g., ADC)
%   mask       - Binary mask (same size as param_map)
%   n_levels   - Number of quantization levels for GLCM (default: 32)
%
% Outputs:
%   features   - Struct with named texture feature fields

    if nargin < 3 || isempty(n_levels), n_levels = 32; end

    % Extract masked voxels
    vals = double(param_map(mask > 0));
    vals = vals(isfinite(vals));

    % Initialize output struct
    features = struct();

    if isempty(vals) || length(vals) < 2
        % Return NaN for all features
        features.energy = NaN;
        features.entropy = NaN;
        features.kurtosis = NaN;
        features.skewness = NaN;
        features.p10 = NaN;
        features.p90 = NaN;
        features.iqr = NaN;
        features.mad = NaN;
        features.rmad = NaN;
        features.glcm_contrast = NaN;
        features.glcm_correlation = NaN;
        features.glcm_energy = NaN;
        features.glcm_homogeneity = NaN;
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
    prob = prob(prob > 0);
    features.entropy = -sum(prob .* log2(prob));

    % --- GLCM features (2D in-plane) ---
    % For 3D volumes, compute GLCM on the central slice with most mask voxels
    if ndims(param_map) == 3
        % Find slice with most mask voxels
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
    if length(valid_pixels) < 4
        features.glcm_contrast = NaN;
        features.glcm_correlation = NaN;
        features.glcm_energy = NaN;
        features.glcm_homogeneity = NaN;
        return;
    end

    min_val = min(valid_pixels);
    max_val = max(valid_pixels);
    if min_val == max_val
        features.glcm_contrast = 0;
        features.glcm_correlation = NaN;
        features.glcm_energy = 1;
        features.glcm_homogeneity = 1;
        return;
    end

    % Quantize
    quantized = round((img_2d - min_val) / (max_val - min_val) * (n_levels - 1)) + 1;
    quantized(~mask_2d) = 0;
    quantized = max(1, min(n_levels, quantized));

    % Compute GLCM at 4 angles and average
    offsets = [0 1; -1 1; -1 0; -1 -1];  % 0, 45, 90, 135 degrees
    contrast_sum = 0;
    correlation_sum = 0;
    energy_sum = 0;
    homogeneity_sum = 0;
    n_valid_angles = 0;

    for a = 1:size(offsets, 1)
        try
            glcm = graycomatrix(quantized, 'Offset', offsets(a,:), ...
                'NumLevels', n_levels, 'GrayLimits', [1 n_levels], ...
                'Symmetric', true);
            props = graycoprops(glcm, {'Contrast', 'Correlation', 'Energy', 'Homogeneity'});
            contrast_sum = contrast_sum + props.Contrast;
            correlation_sum = correlation_sum + props.Correlation;
            energy_sum = energy_sum + props.Energy;
            homogeneity_sum = homogeneity_sum + props.Homogeneity;
            n_valid_angles = n_valid_angles + 1;
        catch
            % Skip this angle if GLCM computation fails
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
end
