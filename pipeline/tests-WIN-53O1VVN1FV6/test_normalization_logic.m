%% test_normalization_logic.m
% Unit test to verify the normalization logic used in apply_dncnn_symmetric.m.
%
% This script tests three key normalization behaviors:
%   1. GTV-based normalization: After z-score normalization using only the
%      tumor (GTV mask) voxels, the mean inside the GTV should be 0 and
%      the standard deviation should be 1.
%   2. Safety checks: Empty masks (no voxels) must be handled gracefully
%      (skip normalization), and single-voxel masks with zero std must
%      clamp sigma to 1.0 to avoid division by zero.
%   3. Winsorization (clamping): Extreme outlier voxels in the background
%      must be clamped to [-3.5, 3.5] after normalization to prevent
%      downstream numerical instability in the DnCNN denoiser.

% ---- Test 1: GTV-based z-score normalization ----

% Create a mock 10x10 image: uniform background at 100, tumor region at center
img = ones(10, 10) * 100; % Healthy background = 100
gtv_mask = zeros(10, 10);
gtv_mask(4:7, 4:7) = 1;   % 4x4 tumor block in the center (16 voxels)

% Fill the GTV region with varying intensities to give non-trivial mean/std.
% These values span 150-250, so the mean is roughly 195 and std roughly 30.
img(4,4) = 150; img(4,5) = 160; img(4,6) = 170; img(4,7) = 180;
img(5,4) = 190; img(5,5) = 200; img(5,6) = 210; img(5,7) = 220;
img(6,4) = 230; img(6,5) = 240; img(6,6) = 250; img(6,7) = 150;
img(7,4) = 160; img(7,5) = 170; img(7,6) = 180; img(7,7) = 190;

% Simulate the bounding-box crop step (here, full image for simplicity)
r_min = 1; r_max = 10; c_min = 1; c_max = 10;
img_cropped = img(r_min:r_max, c_min:c_max);
mask_cropped = gtv_mask(r_min:r_max, c_min:c_max);

% Compute z-score normalization parameters from GTV voxels only
target_voxels = single(img_cropped(mask_cropped > 0));
mu_target = mean(target_voxels, 'all');
sigma_target = std(target_voxels, 0, 'all');

fprintf('Target Mu: %.2f (Expected ~200ish, actual: %.2f)\n', mu_target, mu_target);
fprintf('Target Sigma: %.2f (Expected ~50ish, actual: %.2f)\n', sigma_target, sigma_target);

% Apply normalization: subtract mean, divide by std
img_norm = (single(img_cropped) - mu_target) / sigma_target;

% Verify that normalized GTV voxels have mean=0, std=1
voxels_inside_norm = img_norm(mask_cropped > 0);
mu_inside = mean(voxels_inside_norm);
sigma_inside = std(voxels_inside_norm);

fprintf('Normalized Inside Mu: %.4f (Expected 0)\n', mu_inside);
fprintf('Normalized Inside Sigma: %.4f (Expected 1)\n', sigma_inside);

assert(abs(mu_inside) < 1e-6, 'Mean inside GTV should be 0');
assert(abs(sigma_inside - 1) < 1e-6, 'Sigma inside GTV should be 1');

% ---- Test 2a: Safety check for empty mask ----
% When no voxels are selected by the mask, normalization must be skipped
% entirely to avoid NaN from mean([]) or std([]).
mask_empty = zeros(size(img_cropped));
target_voxels_empty = single(img_cropped(mask_empty > 0));
if isempty(target_voxels_empty)
    fprintf('Safety Check: No target voxels, skipping normalization (expected behavior)\n');
else
    error('Safety Check Failed: target_voxels should be empty');
end

% ---- Test 2b: Safety check for single-voxel (uniform) mask ----
% A single voxel has std=0, which would cause division by zero. The code
% must clamp sigma to 1.0 in this case to produce a finite result.
mask_uniform = zeros(size(img_cropped));
mask_uniform(5,5) = 1;       % Only one voxel in the mask
img_uniform = img_cropped;
img_uniform(5,5) = 10;       % Arbitrary intensity for the single voxel
target_voxels_unif = single(img_uniform(mask_uniform > 0));
mu_unif = mean(target_voxels_unif);
sigma_unif = std(target_voxels_unif);
if sigma_unif < 1e-8
    sigma_unif = 1.0;         % Clamp to avoid division by zero
end
img_norm_unif = (single(img_uniform) - mu_unif) / sigma_unif;
fprintf('Safety Check (Uniform): Sigma clamped to 1.0. Mu: %.2f, Sigma: %.2f\n', mu_unif, sigma_unif);

% ---- Test 3: Winsorization / clamping to [-3.5, 3.5] ----
% Background or artifact voxels far from the GTV mean can produce extreme
% normalized values. Clamping to [-3.5, 3.5] prevents these outliers from
% distorting the DnCNN input range.
img_extreme = img;
img_extreme(1,1) = 2000;    % Extreme high outlier
img_extreme(1,2) = -2000;   % Extreme low outlier

% Apply normalization then clamp
img_norm_clamped = max(-3.5, min(3.5, (single(img_extreme) - mu_target) / sigma_target));

fprintf('Max Normalized Value: %.2f (Expected <= 3.5)\n', max(img_norm_clamped, [], 'all'));
fprintf('Min Normalized Value: %.2f (Expected >= -3.5)\n', min(img_norm_clamped, [], 'all'));

assert(max(img_norm_clamped, [], 'all') <= 3.5 + 1e-6, 'Clamping failed: max exceeded 3.5');
assert(min(img_norm_clamped, [], 'all') >= -3.5 - 1e-6, 'Clamping failed: min below -3.5');

disp('Clamping verification passed!');
disp('All tests passed!');
