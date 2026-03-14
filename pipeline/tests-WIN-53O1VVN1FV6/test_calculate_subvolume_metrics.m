classdef test_calculate_subvolume_metrics < matlab.unittest.TestCase
    % TEST_CALCULATE_SUBVOLUME_METRICS Unit tests for calculate_subvolume_metrics.
    %
    % calculate_subvolume_metrics computes dose-coverage metrics (D95, V50)
    % for diffusion-defined resistant sub-volumes within a GTV. The sub-volume
    % consists of voxels whose diffusion parameter (e.g., ADC) is strictly
    % below a threshold, indicating treatment-resistant tissue.
    %
    % Two code paths are tested:
    %   1. 1D path (has_3d = false): sub-volume defined by simple vector thresholding
    %   2. 3D path (has_3d = true): sub-volume undergoes morphological processing
    %      (imopen, imclose, bwareaopen) to remove small disconnected clusters
    %
    % Key metrics:
    %   D95: dose covering 95% of the sub-volume (5th percentile of dose)
    %   V50: percentage of sub-volume voxels receiving >= 50 Gy
    %
    % A minimum of 100 voxels (min_subvol_voxels) is required for valid output;
    % smaller sub-volumes return NaN.

    properties
        OriginalPath  % Saved MATLAB path for restoration in teardown
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            % Save the current path and add utils for calculate_subvolume_metrics.
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            % Restore the original MATLAB path.
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testBasicD95AndV50(testCase)
            % All 200 voxels are below the threshold (0 < 0.5), so the
            % entire GTV is the sub-volume. D95 and V50 are computed from
            % a dose ramp spanning 10-80 Gy.
            n = 200;
            vector   = zeros(n, 1);          % all below threshold = 0.5
            threshold = 0.5;
            dose_vec = linspace(10, 80, n)'; % linear dose ramp 10-80 Gy

            [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, false, []);

            % D95 = prctile(dose_vec, 5): the 5th percentile of dose
            expected_d95 = prctile(dose_vec, 5);
            % V50 = fraction of voxels receiving >= 50 Gy, as a percentage
            expected_v50 = sum(dose_vec >= 50) / n * 100;

            testCase.verifyEqual(d95, expected_d95, 'AbsTol', 1e-6, ...
                'D95 should equal the 5th percentile of dose in the sub-volume.');
            testCase.verifyEqual(v50, expected_v50, 'AbsTol', 1e-6, ...
                'V50 should equal the percentage of sub-volume voxels receiving >= 50 Gy.');
        end

        function testEmptySubvolumeReturnsNaN(testCase)
            % All voxels are above the threshold (1.0 > 0.5), so the
            % sub-volume is empty. Both D95 and V50 should be NaN.
            n = 200;
            vector    = ones(n, 1);    % all ABOVE threshold
            threshold = 0.5;
            dose_vec  = 60 * ones(n, 1);

            [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, false, []);

            testCase.verifyTrue(isnan(d95), 'D95 should be NaN when sub-volume is empty.');
            testCase.verifyTrue(isnan(v50), 'V50 should be NaN when sub-volume is empty.');
        end

        function testInsufficientVoxelsBelowMinReturnsNaN(testCase)
            % Only 50 voxels are below threshold (< min_subvol_voxels = 100).
            % Too few voxels for reliable dose metrics, so NaN is returned.
            n = 200;
            vector    = ones(n, 1);
            threshold = 0.5;
            vector(1:50) = 0;  % only 50 below threshold
            dose_vec = 60 * ones(n, 1);

            [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, false, []);

            testCase.verifyTrue(isnan(d95), 'D95 should be NaN when sub-volume has < 100 voxels.');
            testCase.verifyTrue(isnan(v50), 'V50 should be NaN when sub-volume has < 100 voxels.');
        end

        function testV50AllAbove50Gy(testCase)
            % All 150 sub-volume voxels receive 60 Gy (>= 50), so V50 = 100%.
            n = 150;
            vector    = zeros(n, 1);
            threshold = 1.0;
            dose_vec  = 60 * ones(n, 1);

            [~, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, false, []);

            testCase.verifyEqual(v50, 100.0, 'AbsTol', 1e-9, ...
                'V50 should be 100 % when all sub-volume voxels receive >= 50 Gy.');
        end

        function testV50NoneAbove50Gy(testCase)
            % All 150 sub-volume voxels receive 30 Gy (< 50), so V50 = 0%.
            n = 150;
            vector    = zeros(n, 1);
            threshold = 1.0;
            dose_vec  = 30 * ones(n, 1);

            [~, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, false, []);

            testCase.verifyEqual(v50, 0.0, 'AbsTol', 1e-9, ...
                'V50 should be 0 % when no sub-volume voxels receive >= 50 Gy.');
        end

        function testExactThresholdBoundary(testCase)
            % Tests the strict-less-than comparison: voxels exactly equal to
            % the threshold (1.0) are NOT included in the sub-volume.
            % Only the 150 voxels at 0.5 (< 1.0) qualify.
            n = 200;
            threshold = 1.0;
            vector = [0.5 * ones(150, 1); threshold * ones(50, 1)];
            dose_vec = 55 * ones(n, 1);

            [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, false, []);

            % 150 voxels in sub-volume (>= 100 min) -> valid computation
            testCase.verifyFalse(isnan(d95), 'D95 should be computed for 150-voxel sub-volume.');
            % All 150 voxels receive 55 Gy (>= 50), so V50 = 100%
            testCase.verifyEqual(v50, 100.0, 'AbsTol', 1e-9, ...
                'V50 should be 100 % since all 150 sub-volume voxels receive 55 Gy.');
        end

        function testWith3DMaskRunsWithoutError(testCase)
            % Exercises the 3D morphological processing path (has_3d = true).
            % A 7x7x7 = 343-voxel cube is placed in the center of a 12^3
            % volume. This is large enough to survive morphological opening,
            % closing, and bwareaopen(10), yielding valid (non-NaN) metrics.
            vol_size = [12, 12, 12];
            gtv_mask_3d = false(vol_size);
            gtv_mask_3d(3:9, 3:9, 3:9) = true;  % 343-voxel contiguous block

            n_gtv = sum(gtv_mask_3d(:));   % 343 voxels
            vector    = zeros(n_gtv, 1);   % all below threshold
            threshold = 0.5;
            dose_vec  = 55 * ones(n_gtv, 1);

            [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, true, gtv_mask_3d);

            % Block may shrink slightly after morphological ops but should
            % remain above the 100-voxel minimum for valid output
            testCase.verifyFalse(isnan(d95), ...
                'D95 should not be NaN for a large 3D contiguous sub-volume.');
            testCase.verifyFalse(isnan(v50), ...
                'V50 should not be NaN for a large 3D contiguous sub-volume.');
            testCase.verifyGreaterThanOrEqual(v50, 0.0);
            testCase.verifyLessThanOrEqual(v50, 100.0);
        end

        function testWith3DMaskSmallBlockReturnsNaN(testCase)
            % A 3x3x3 = 27-voxel cube is too small: even if it survives
            % bwareaopen(10), the 27 voxels are below the min_subvol_voxels
            % threshold of 100, so NaN is returned.
            vol_size = [10, 10, 10];
            gtv_mask_3d = false(vol_size);
            gtv_mask_3d(5:7, 5:7, 5:7) = true;  % 27-voxel block

            n_gtv = sum(gtv_mask_3d(:));   % 27 voxels
            vector    = zeros(n_gtv, 1);
            threshold = 0.5;
            dose_vec  = 55 * ones(n_gtv, 1);

            [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, true, gtv_mask_3d);

            testCase.verifyTrue(isnan(d95), ...
                'D95 should be NaN when 3D sub-volume has fewer than 100 voxels.');
            testCase.verifyTrue(isnan(v50), ...
                'V50 should be NaN when 3D sub-volume has fewer than 100 voxels.');
        end

    end
end
