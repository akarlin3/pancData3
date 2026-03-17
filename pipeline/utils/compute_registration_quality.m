function quality = compute_registration_quality(reference_volume, warped_volume, deformation_field, voxel_spacing)
% COMPUTE_REGISTRATION_QUALITY  Quantify deformable image registration fidelity.
%
%   Computes Jacobian determinant statistics, normalized cross-correlation,
%   and mutual information between reference and warped volumes to assess
%   registration quality.
%
% Inputs:
%   reference_volume  - 3D reference image volume
%   warped_volume     - 3D warped (registered) image volume
%   deformation_field - 4D deformation field [x, y, z, 3] or empty
%   voxel_spacing     - [dx, dy, dz] voxel dimensions in mm (default: [1 1 1])
%                       Used for physically correct Jacobian gradient computation.
%
% Outputs:
%   quality           - Struct with registration quality metrics

    if nargin < 4 || isempty(voxel_spacing)
        voxel_spacing = [1 1 1];
    end

    quality = struct();

    ref = double(reference_volume);
    wrp = double(warped_volume);

    % --- Jacobian determinant from deformation field ---
    if ~isempty(deformation_field) && ndims(deformation_field) >= 4
        % Compute Jacobian determinant at each voxel
        dx = deformation_field(:,:,:,1);
        dy = deformation_field(:,:,:,2);
        dz = deformation_field(:,:,:,3);

        % Numerical gradient of deformation field
        % MATLAB gradient() uses (dy, dx, dz) ordering for 3D arrays
        [dxx, dxy, dxz] = gradient(dx, voxel_spacing(2), voxel_spacing(1), voxel_spacing(3));
        [dyx, dyy, dyz] = gradient(dy, voxel_spacing(2), voxel_spacing(1), voxel_spacing(3));
        [dzx, dzy, dzz] = gradient(dz, voxel_spacing(2), voxel_spacing(1), voxel_spacing(3));

        % Add identity (deformation = displacement + identity)
        dxx = dxx + 1;
        dyy = dyy + 1;
        dzz = dzz + 1;

        % Jacobian determinant
        jac_det = dxx .* (dyy .* dzz - dyz .* dzy) ...
                - dxy .* (dyx .* dzz - dyz .* dzx) ...
                + dxz .* (dyx .* dzy - dyy .* dzx);

        jac_vals = jac_det(:);
        jac_vals = jac_vals(isfinite(jac_vals));

        quality.jacobian_mean = mean(jac_vals);
        quality.jacobian_std = std(jac_vals);
        quality.jacobian_min = min(jac_vals);
        quality.jacobian_max = max(jac_vals);
        quality.jacobian_folding_pct = 100 * sum(jac_vals < 0) / numel(jac_vals);
    else
        quality.jacobian_mean = NaN;
        quality.jacobian_std = NaN;
        quality.jacobian_min = NaN;
        quality.jacobian_max = NaN;
        quality.jacobian_folding_pct = NaN;
    end

    % --- Normalized Cross-Correlation (NCC) ---
    % Compute within overlapping non-zero regions
    valid = isfinite(ref) & isfinite(wrp) & (ref ~= 0 | wrp ~= 0);
    if sum(valid(:)) > 1
        r = ref(valid) - mean(ref(valid));
        w = wrp(valid) - mean(wrp(valid));
        denom = sqrt(sum(r.^2) * sum(w.^2));
        if denom > 0
            quality.ncc = sum(r .* w) / denom;
        else
            quality.ncc = NaN;
        end
    else
        quality.ncc = NaN;
    end

    % --- Mutual Information ---
    quality.mutual_information = compute_mi(ref(valid), wrp(valid));
end


function mi = compute_mi(x, y)
% Compute mutual information between two vectors using histogram method
    n_bins = 32;
    if isempty(x) || isempty(y) || length(x) < 10
        mi = NaN;
        return;
    end

    % Joint histogram
    x_edges = linspace(min(x), max(x) + eps, n_bins + 1);
    y_edges = linspace(min(y), max(y) + eps, n_bins + 1);

    joint_hist = zeros(n_bins, n_bins);
    for i = 1:length(x)
        xi = min(n_bins, max(1, sum(x(i) >= x_edges(1:end-1))));
        yi = min(n_bins, max(1, sum(y(i) >= y_edges(1:end-1))));
        joint_hist(xi, yi) = joint_hist(xi, yi) + 1;
    end

    % Normalize
    joint_prob = joint_hist / sum(joint_hist(:));
    px = sum(joint_prob, 2);
    py = sum(joint_prob, 1);

    % MI = sum p(x,y) * log(p(x,y) / (p(x)*p(y)))
    mi = 0;
    for i = 1:n_bins
        for j = 1:n_bins
            if joint_prob(i,j) > 0 && px(i) > 0 && py(j) > 0
                mi = mi + joint_prob(i,j) * log2(joint_prob(i,j) / (px(i) * py(j)));
            end
        end
    end
end
