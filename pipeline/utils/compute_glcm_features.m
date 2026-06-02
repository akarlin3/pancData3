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
