classdef test_calculate_subvolume_metrics < matlab.unittest.TestCase
    % TEST_CALCULATE_SUBVOLUME_METRICS Unit tests for calculate_subvolume_metrics.
    %
    % Tests D95 and V50 dose-metric computation for diffusion-defined
    % resistant sub-volumes, covering the no-3D-mask (1D) path and the
    % 3D morphological-processing path.

    properties
        OriginalPath
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testBasicD95AndV50(testCase)
            % 200 voxels all below threshold; verifies D95 and V50 are computed.
            n = 200;
            vector   = zeros(n, 1);          % all below threshold = 0.5
            threshold = 0.5;
            dose_vec = linspace(10, 80, n)'; % spans 10–80 Gy

            [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, false, []);

            % D95 = prctile(dose_vec, 5) — dose covering 95 % of sub-volume
            expected_d95 = prctile(dose_vec, 5);
            % V50 = fraction receiving >= 50 Gy
            expected_v50 = sum(dose_vec >= 50) / n * 100;

            testCase.verifyEqual(d95, expected_d95, 'AbsTol', 1e-6, ...
                'D95 should equal the 5th percentile of dose in the sub-volume.');
            testCase.verifyEqual(v50, expected_v50, 'AbsTol', 1e-6, ...
                'V50 should equal the percentage of sub-volume voxels receiving ≥ 50 Gy.');
        end

        function testEmptySubvolumeReturnsNaN(testCase)
            % No voxels below threshold → sub-volume is empty → NaN outputs.
            n = 200;
            vector    = ones(n, 1);    % all ABOVE threshold = 0.5
            threshold = 0.5;
            dose_vec  = 60 * ones(n, 1);

            [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, false, []);

            testCase.verifyTrue(isnan(d95), 'D95 should be NaN when sub-volume is empty.');
            testCase.verifyTrue(isnan(v50), 'V50 should be NaN when sub-volume is empty.');
        end

        function testInsufficientVoxelsBelowMinReturnsNaN(testCase)
            % Fewer than 100 voxels in sub-volume → NaN outputs (min_subvol_voxels = 100).
            n = 200;
            vector    = ones(n, 1);
            threshold = 0.5;
            % Only 50 voxels are below threshold
            vector(1:50) = 0;
            dose_vec = 60 * ones(n, 1);

            [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, false, []);

            testCase.verifyTrue(isnan(d95), 'D95 should be NaN when sub-volume has < 100 voxels.');
            testCase.verifyTrue(isnan(v50), 'V50 should be NaN when sub-volume has < 100 voxels.');
        end

        function testV50AllAbove50Gy(testCase)
            % All sub-volume voxels receive ≥ 50 Gy → V50 = 100 %.
            n = 150;
            vector    = zeros(n, 1);
            threshold = 1.0;
            dose_vec  = 60 * ones(n, 1);  % all 60 Gy

            [~, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, false, []);

            testCase.verifyEqual(v50, 100.0, 'AbsTol', 1e-9, ...
                'V50 should be 100 % when all sub-volume voxels receive ≥ 50 Gy.');
        end

        function testV50NoneAbove50Gy(testCase)
            % No sub-volume voxels receive ≥ 50 Gy → V50 = 0 %.
            n = 150;
            vector    = zeros(n, 1);
            threshold = 1.0;
            dose_vec  = 30 * ones(n, 1);  % all 30 Gy

            [~, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, false, []);

            testCase.verifyEqual(v50, 0.0, 'AbsTol', 1e-9, ...
                'V50 should be 0 % when no sub-volume voxels receive ≥ 50 Gy.');
        end

        function testExactThresholdBoundary(testCase)
            % Voxels equal to threshold are NOT in sub-volume (strict < comparison).
            n = 200;
            threshold = 1.0;
            % First 150 strictly below, last 50 equal to threshold
            vector = [0.5 * ones(150, 1); threshold * ones(50, 1)];
            dose_vec = 55 * ones(n, 1);

            [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, false, []);

            % 150 voxels in sub-volume (>= 100 min) → should compute
            testCase.verifyFalse(isnan(d95), 'D95 should be computed for 150-voxel sub-volume.');
            testCase.verifyEqual(v50, 100.0, 'AbsTol', 1e-9, ...
                'V50 should be 100 % since all 150 sub-volume voxels receive 55 Gy.');
        end

        function testWith3DMaskRunsWithoutError(testCase)
            % Exercises the 3D morphological processing path (has_3d = true).
            % Creates a volume with a large contiguous GTV block to survive
            % imopen/imclose and bwareaopen(10).
            vol_size = [12, 12, 12];
            gtv_mask_3d = false(vol_size);
            % 7×7×7 = 343-voxel cube in the centre of the volume
            gtv_mask_3d(3:9, 3:9, 3:9) = true;

            n_gtv = sum(gtv_mask_3d(:));   % 343 voxels
            vector    = zeros(n_gtv, 1);   % all below threshold
            threshold = 0.5;
            dose_vec  = 55 * ones(n_gtv, 1);

            [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, true, gtv_mask_3d);

            % After morphological ops the block may shrink slightly, but
            % should remain >= 100 voxels → expect valid (non-NaN) results.
            testCase.verifyFalse(isnan(d95), ...
                'D95 should not be NaN for a large 3D contiguous sub-volume.');
            testCase.verifyFalse(isnan(v50), ...
                'V50 should not be NaN for a large 3D contiguous sub-volume.');
            testCase.verifyGreaterThanOrEqual(v50, 0.0);
            testCase.verifyLessThanOrEqual(v50, 100.0);
        end

        function testWith3DMaskSmallBlockReturnsNaN(testCase)
            % If the 3D block is too small to survive bwareaopen(10), sub-volume
            % drops to zero voxels → NaN outputs.
            vol_size = [10, 10, 10];
            gtv_mask_3d = false(vol_size);
            % 3×3×3 = 27-voxel cube — will survive bwareaopen(10) but may not
            % survive the minimum sub-volume check (100 voxels).
            gtv_mask_3d(5:7, 5:7, 5:7) = true;

            n_gtv = sum(gtv_mask_3d(:));   % 27 voxels
            vector    = zeros(n_gtv, 1);
            threshold = 0.5;
            dose_vec  = 55 * ones(n_gtv, 1);

            [d95, v50] = calculate_subvolume_metrics(vector, threshold, dose_vec, true, gtv_mask_3d);

            % 27 voxels < min_subvol_voxels (100) → expect NaN
            testCase.verifyTrue(isnan(d95), ...
                'D95 should be NaN when 3D sub-volume has fewer than 100 voxels.');
            testCase.verifyTrue(isnan(v50), ...
                'V50 should be NaN when 3D sub-volume has fewer than 100 voxels.');
        end

    end
end
