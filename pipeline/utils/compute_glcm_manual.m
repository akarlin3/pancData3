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
