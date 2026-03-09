%% test_core_methods.m - Unit testing for tumor core extraction methods
%
% Exercises all 11 tumor core delineation algorithms in extract_tumor_core()
% using a synthetic bimodal voxel population. The data consists of 30 "core"
% voxels (low ADC/D, low f) and 70 "margin" voxels (high ADC/D, high f),
% mimicking the expected pattern in treatment-resistant pancreatic tumours.
%
% For each method, the test verifies:
%   - Output is logical (boolean mask)
%   - Output length matches input (100 voxels)
%   - Mask is neither all-true nor all-false (algorithm found some structure)
%
% Additionally tests:
%   - fDM method: baseline fallback (timepoint 1) and longitudinal mode (timepoint 2)
%   - Edge case: all-NaN input should yield an all-false mask
%
% Usage:
%   test_core_methods   % from MATLAB command window (run as function)

function test_core_methods()
    disp('==== Running test_core_methods ====');
    [dir_path, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(dir_path, '..', 'utils'));

    % Setup mock configuration with clinically typical thresholds.
    % adc_thresh and d_thresh at 0.001 mm^2/s capture restricted diffusion.
    % f_thresh at 0.1 (10%) separates perfusion-dominant from diffusion-dominant.
    config_struct = struct();
    config_struct.adc_thresh = 0.001;
    config_struct.d_thresh = 0.001;
    config_struct.f_thresh = 0.1;
    config_struct.dstar_thresh = 0.01;
    config_struct.min_vox_hist = 50; % Lower than default (100) to work with 100-voxel test data
    config_struct.core_percentile = 25;
    config_struct.core_n_clusters = 2; % Binary clustering (core vs. margin)
    config_struct.fdm_parameter = 'adc';
    config_struct.fdm_thresh = 0.0002;

    % Setup mock parameter maps (100 voxels total).
    % The bimodal distribution is designed so that threshold-based, clustering-based,
    % and data-driven methods should all identify roughly 30 core voxels.
    rng(42); % Fix seed for reproducibility across test runs
    % 30 "core" voxels: low ADC (~0.0005), low D (~0.0004), low f (~0.05)
    % 70 "margin" voxels: high ADC (~0.0015), high D (~0.0014), high f (~0.20)
    adc_vec = [0.0005 + 0.0001*randn(30,1); 0.0015 + 0.0002*randn(70,1)];
    d_vec   = [0.0004 + 0.0001*randn(30,1); 0.0014 + 0.0002*randn(70,1)];
    f_vec   = [0.05 + 0.01*randn(30,1);     0.2 + 0.05*randn(70,1)];
    dstar_vec = nan(100,1); % D* is intentionally NaN to test NaN tolerance

    % Setup mock 3D geometry: 10x10x1 slab (single-slice).
    % Required by methods that use spatial connectivity (region_growing,
    % active_contours, spectral).
    has_3d = true;
    gtv_mask_3d = true(10, 10, 1); % Entire slab is within GTV

    % All 10 non-longitudinal core methods to test. The fDM method requires
    % additional opts (baseline vectors) and is tested separately below.
    methods = {'adc_threshold', 'd_threshold', 'df_intersection', 'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', 'percentile', 'spectral'};

    for i=1:length(methods)
        config_struct.core_method = methods{i};
        try
            core_mask = extract_tumor_core(config_struct, adc_vec, d_vec, f_vec, dstar_vec, has_3d, gtv_mask_3d);

            % Basic contract checks: output must be logical, correct length
            assert(islogical(core_mask), sprintf('Method %s: Output is not logical', methods{i}));
            assert(length(core_mask) == 100, sprintf('Method %s: Output length is incorrect', methods{i}));

            % Sanity: with a clearly bimodal distribution, every method should
            % find at least some core voxels and leave some as margin.
            % Note: GMM/kmeans may occasionally flip cluster labels, and
            % active_contours can behave unpredictably on small 10x10 grids,
            % but the bimodal gap is large enough to prevent degenerate results.
            n_core = sum(core_mask);
            assert(n_core > 0 && n_core < 100, sprintf('Method %s: Mask is all or nothing (%d voxels)', methods{i}, n_core));

            fprintf('[PASS] %s: Found %d core voxels\n', methods{i}, n_core);

        catch ME
            fprintf('[FAIL] %s: %s\n', methods{i}, ME.message);
            rethrow(ME);
        end
    end

    % Test fDM method separately (requires opts with baseline vectors)
    config_struct.core_method = 'fdm';
    try
        % fDM at baseline (k=1) should fall back to threshold
        opts_baseline = struct('timepoint_index', 1);
        core_mask = extract_tumor_core(config_struct, adc_vec, d_vec, f_vec, dstar_vec, has_3d, gtv_mask_3d, opts_baseline);
        assert(islogical(core_mask), 'fDM baseline: Output is not logical');
        assert(length(core_mask) == 100, 'fDM baseline: Output length is incorrect');
        fprintf('[PASS] fdm (baseline fallback): Found %d core voxels\n', sum(core_mask));

        % fDM at k=2 with simulated baseline
        opts_fx2 = struct('timepoint_index', 2);
        % Simulate treatment response: some voxels increase ADC, some decrease
        opts_fx2.baseline_adc_vec = adc_vec + 0.0003 * randn(100, 1);
        core_mask = extract_tumor_core(config_struct, adc_vec, d_vec, f_vec, dstar_vec, has_3d, gtv_mask_3d, opts_fx2);
        assert(islogical(core_mask), 'fDM fx2: Output is not logical');
        assert(length(core_mask) == 100, 'fDM fx2: Output length is incorrect');
        fprintf('[PASS] fdm (fx2 with baseline): Found %d core voxels\n', sum(core_mask));

    catch ME
        fprintf('[FAIL] fdm: %s\n', ME.message);
        rethrow(ME);
    end

    % Edge case test: empty inputs
    config_struct.core_method = 'otsu';
    empty_vec = nan(10,1);
    % Should return gracefully with all-false mask
    empty_mask = extract_tumor_core(config_struct, empty_vec, empty_vec, empty_vec, empty_vec, false, []);
    assert(sum(empty_mask) == 0, 'Empty input should yield empty mask');
    fprintf('[PASS] edge_case: Handled NaN/Empty input successfully\n');

    disp('==== All core method tests passed! ====');
end
