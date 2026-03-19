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

    % Check variance to avoid undefined logarithms
    if var(x) == 0 || var(y) == 0
        if var(x) == 0 && var(y) == 0
            % Both constant: perfect match if same value, no information otherwise
            if all(x == y)
                mi = Inf; % Perfect mutual information
            else
                mi = 0; % No mutual information
            end
        else
            % One constant, one variable: no mutual information
            mi = 0;
        end
        return;
    end

    % Check for sufficient dynamic range
    x_range = max(x) - min(x);
    y_range = max(y) - min(y);
    if x_range < eps || y_range < eps
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


function tests = test_compute_registration_quality
% TEST_COMPUTE_REGISTRATION_QUALITY Comprehensive test suite
    tests = functiontests(localfunctions);
end

function test_identity_transform(testCase)
% Test identity transformation should give ideal metrics
    % Create simple 3D volume
    vol = zeros(10, 10, 10);
    vol(3:8, 3:8, 3:8) = 100;
    
    % Identity deformation field (zero displacement)
    def_field = zeros(10, 10, 10, 3);
    
    quality = compute_registration_quality(vol, vol, def_field);
    
    % Identity transform should have Jacobian ~1, perfect NCC, high MI
    verifyEqual(testCase, quality.jacobian_mean, 1, 'AbsTol', 1e-10);
    verifyEqual(testCase, quality.jacobian_std, 0, 'AbsTol', 1e-10);
    verifyEqual(testCase, quality.jacobian_folding_pct, 0);
    verifyEqual(testCase, quality.ncc, 1, 'AbsTol', 1e-10);
    verifyTrue(testCase, quality.mutual_information > 0);
end

function test_perfect_alignment_different_images(testCase)
% Test perfect alignment of different but related images
    vol1 = randn(15, 15, 15) * 50 + 100;
    vol2 = vol1; % Perfect match
    
    def_field = zeros(15, 15, 15, 3);
    
    quality = compute_registration_quality(vol1, vol2, def_field);
    
    verifyEqual(testCase, quality.jacobian_mean, 1, 'AbsTol', 1e-10);
    verifyEqual(testCase, quality.ncc, 1, 'AbsTol', 1e-10);
end

function test_large_deformation(testCase)
% Test large deformation field
    vol = zeros(20, 20, 20);
    vol(5:15, 5:15, 5:15) = 1;
    
    % Large displacement field
    [X, Y, Z] = meshgrid(1:20, 1:20, 1:20);
    def_field = zeros(20, 20, 20, 3);
    def_field(:,:,:,1) = 0.5 * sin(X * pi / 10); % Large x displacement
    def_field(:,:,:,2) = 0.3 * cos(Y * pi / 10); % Large y displacement
    def_field(:,:,:,3) = 0.2 * sin(Z * pi / 10); % Large z displacement
    
    quality = compute_registration_quality(vol, vol, def_field);
    
    % Large deformation should deviate from identity
    verifyNotEqual(testCase, quality.jacobian_mean, 1, 'AbsTol', 0.1);
    verifyTrue(testCase, quality.jacobian_std > 0);
    verifyTrue(testCase, isfinite(quality.jacobian_folding_pct));
end

function test_folding_deformation(testCase)
% Test deformation that causes folding (negative Jacobian)
    vol = ones(10, 10, 10);
    
    % Create folding deformation
    def_field = zeros(10, 10, 10, 3);
    [X, Y, Z] = meshgrid(1:10, 1:10, 1:10);
    def_field(:,:,:,1) = -2 * X; % Strong negative gradient
    
    quality = compute_registration_quality(vol, vol, def_field);
    
    % Should detect folding
    verifyTrue(testCase, quality.jacobian_folding_pct > 0);
    verifyTrue(testCase, quality.jacobian_min < 0);
end

function test_noise_sensitivity(testCase)
% Test sensitivity to noise
    vol1 = ones(12, 12, 12) * 100;
    vol1(4:9, 4:9, 4:9) = 200;
    
    % Add different levels of noise
    noise_levels = [0, 5, 20, 50];
    ncc_values = zeros(size(noise_levels));
    mi_values = zeros(size(noise_levels));
    
    for i = 1:length(noise_levels)
        vol2 = vol1 + randn(size(vol1)) * noise_levels(i);
        quality = compute_registration_quality(vol1, vol2, []);
        ncc_values(i) = quality.ncc;
        mi_values(i) = quality.mutual_information;
    end
    
    % NCC should decrease with noise
    verifyTrue(testCase, ncc_values(1) > ncc_values(end));
    % MI should also generally decrease with noise
    verifyTrue(testCase, mi_values(1) > mi_values(end) || isnan(mi_values(end)));
end

function test_empty_deformation_field(testCase)
% Test behavior with empty deformation field
    vol1 = randn(8, 8, 8);
    vol2 = randn(8, 8, 8);
    
    quality = compute_registration_quality(vol1, vol2, []);
    
    % Jacobian metrics should be NaN
    verifyTrue(testCase, isnan(quality.jacobian_mean));
    verifyTrue(testCase, isnan(quality.jacobian_std));
    verifyTrue(testCase, isnan(quality.jacobian_folding_pct));
    
    % But NCC and MI should still be computed
    verifyTrue(testCase, isfinite(quality.ncc) || isnan(quality.ncc));
    verifyTrue(testCase, isfinite(quality.mutual_information) || isnan(quality.mutual_information));
end

function test_zero_volumes(testCase)
% Test edge case with zero volumes
    vol_zero = zeros(5, 5, 5);
    vol_nonzero = ones(5, 5, 5);
    
    quality1 = compute_registration_quality(vol_zero, vol_zero, []);
    quality2 = compute_registration_quality(vol_zero, vol_nonzero, []);
    
    % Zero vs zero should give NaN or specific values
    verifyTrue(testCase, isnan(quality1.ncc) || quality1.ncc == 0);
    
    % Zero vs nonzero should be handled gracefully
    verifyTrue(testCase, isfinite(quality2.ncc) || isnan(quality2.ncc));
end

function test_anisotropic_voxel_spacing(testCase)
% Test with anisotropic voxel spacing
    vol = ones(8, 8, 8);
    vol(3:6, 3:6, 3:6) = 2;
    
    def_field = zeros(8, 8, 8, 3);
    def_field(:,:,:,1) = 0.1; % Small uniform displacement
    
    voxel_spacing_iso = [1, 1, 1];
    voxel_spacing_aniso = [0.5, 0.5, 2.0];
    
    quality_iso = compute_registration_quality(vol, vol, def_field, voxel_spacing_iso);
    quality_aniso = compute_registration_quality(vol, vol, def_field, voxel_spacing_aniso);
    
    % Jacobian should be different due to different voxel spacing
    verifyNotEqual(testCase, quality_iso.jacobian_mean, quality_aniso.jacobian_mean, 'AbsTol', 1e-10);
end

function test_infinite_and_nan_values(testCase)
% Test robustness to infinite and NaN values
    vol1 = randn(6, 6, 6);
    vol2 = vol1;
    
    % Introduce some problematic values
    vol1(1, 1, 1) = Inf;
    vol2(2, 2, 2) = -Inf;
    vol1(3, 3, 3) = NaN;
    
    def_field = randn(6, 6, 6, 3) * 0.1;
    def_field(1, 1, 1, 1) = NaN; % NaN in deformation
    
    quality = compute_registration_quality(vol1, vol2, def_field);
    
    % Should handle problematic values gracefully
    verifyTrue(testCase, isfinite(quality.jacobian_mean) || isnan(quality.jacobian_mean));
    verifyTrue(testCase, isfinite(quality.ncc) || isnan(quality.ncc));
    verifyTrue(testCase, isfinite(quality.mutual_information) || isnan(quality.mutual_information));
end

function test_small_volumes(testCase)
% Test behavior with very small volumes
    vol1 = randn(2, 2, 2);
    vol2 = randn(2, 2, 2);
    def_field = randn(2, 2, 2, 3) * 0.1;
    
    quality = compute_registration_quality(vol1, vol2, def_field);
    
    % Should not crash and return reasonable values or NaN
    verifyTrue(testCase, isstruct(quality));
    verifyTrue(testCase, isfield(quality, 'jacobian_mean'));
    verifyTrue(testCase, isfield(quality, 'ncc'));
    verifyTrue(testCase, isfield(quality, 'mutual_information'));
end

function test_zero_variance_images(testCase)
% Test mutual information calculation with zero variance images
    vol_constant = ones(8, 8, 8) * 100;  % Zero variance
    vol_variable = randn(8, 8, 8) * 50 + 100;  % Non-zero variance
    vol_constant2 = ones(8, 8, 8) * 200;  % Different constant value
    
    % Test constant vs constant (same value)
    quality1 = compute_registration_quality(vol_constant, vol_constant, []);
    verifyTrue(testCase, isinf(quality1.mutual_information) || quality1.mutual_information > 0);
    
    % Test constant vs constant (different values)
    quality2 = compute_registration_quality(vol_constant, vol_constant2, []);
    verifyEqual(testCase, quality2.mutual_information, 0);
    
    % Test constant vs variable
    quality3 = compute_registration_quality(vol_constant, vol_variable, []);
    verifyEqual(testCase, quality3.mutual_information, 0);
    
    % Test variable vs constant
    quality4 = compute_registration_quality(vol_variable, vol_constant, []);
    verifyEqual(testCase, quality4.mutual_information, 0);
end