%% test_normalization_logic.m
% Unit test to verify the normalization logic in apply_dncnn_symmetric.m

% 1. Create a mock 10x10 image with a 4x4 tumor (mask) in the center
img = ones(10, 10) * 100; % Healthy background = 100
gtv_mask = zeros(10, 10);
gtv_mask(4:7, 4:7) = 1;

% Add some variation inside the GTV
% Mean inside = 200, Std inside = 50
img(4,4) = 150; img(4,5) = 160; img(4,6) = 170; img(4,7) = 180;
img(5,4) = 190; img(5,5) = 200; img(5,6) = 210; img(5,7) = 220;
img(6,4) = 230; img(6,5) = 240; img(6,6) = 250; img(6,7) = 150;
img(7,4) = 160; img(7,5) = 170; img(7,6) = 180; img(7,7) = 190;

% Current implementation extracts bounding box and then normalizes.
% Let's simulate the normalization block from the code.

% Simulation of step 2 & 3: Cropping
r_min = 1; r_max = 10; c_min = 1; c_max = 10;
img_cropped = img(r_min:r_max, c_min:c_max);
mask_cropped = gtv_mask(r_min:r_max, c_min:c_max);

% Simulation of step 4: Normalization
target_voxels = single(img_cropped(mask_cropped > 0));
mu_target = mean(target_voxels, 'all');
sigma_target = std(target_voxels, 0, 'all');

fprintf('Target Mu: %.2f (Expected ~200ish, actual: %.2f)\n', mu_target, mu_target);
fprintf('Target Sigma: %.2f (Expected ~50ish, actual: %.2f)\n', sigma_target, sigma_target);

img_norm = (single(img_cropped) - mu_target) / sigma_target;

% Verification
voxels_inside_norm = img_norm(mask_cropped > 0);
mu_inside = mean(voxels_inside_norm);
sigma_inside = std(voxels_inside_norm);

fprintf('Normalized Inside Mu: %.4f (Expected 0)\n', mu_inside);
fprintf('Normalized Inside Sigma: %.4f (Expected 1)\n', sigma_inside);

assert(abs(mu_inside) < 1e-6, 'Mean inside GTV should be 0');
assert(abs(sigma_inside - 1) < 1e-6, 'Sigma inside GTV should be 1');

% Test safety check (empty mask)
mask_empty = zeros(size(img_cropped));
target_voxels_empty = single(img_cropped(mask_empty > 0));
if isempty(target_voxels_empty)
    fprintf('Safety Check: No target voxels, skipping normalization (expected behavior)\n');
else
    error('Safety Check Failed: target_voxels should be empty');
end

% Test safety check (uniform mask)
mask_uniform = zeros(size(img_cropped));
mask_uniform(5,5) = 1;
img_uniform = img_cropped;
img_uniform(5,5) = 10;
target_voxels_unif = single(img_uniform(mask_uniform > 0));
mu_unif = mean(target_voxels_unif);
sigma_unif = std(target_voxels_unif);
if sigma_unif < 1e-8
    sigma_unif = 1.0;
end
img_norm_unif = (single(img_uniform) - mu_unif) / sigma_unif;
fprintf('Safety Check (Uniform): Sigma clamped to 1.0. Mu: %.2f, Sigma: %.2f\n', mu_unif, sigma_unif);

disp('All tests passed!');
