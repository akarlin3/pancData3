classdef test_compute_texture_features < matlab.unittest.TestCase
    % TEST_COMPUTE_TEXTURE_FEATURES  Tests for texture feature extraction.
    %
    % Covers:
    %   - Checkerboard pattern: high contrast, low homogeneity
    %   - Uniform image: zero contrast
    %   - Output struct field count
    %   - 3D input handling
    %   - 3D GLRLM (13 directions) vs 2D GLRLM (4 directions)
    %   - Empty mask edge case

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

        function testCheckerboardHighContrast(testCase)
            % Checkerboard pattern should have high GLCM contrast and low homogeneity.
            img = zeros(64, 64);
            for i = 1:64
                for j = 1:64
                    if mod(i + j, 2) == 0
                        img(i, j) = 1;
                    else
                        img(i, j) = 0;
                    end
                end
            end
            mask = true(64, 64);

            features = compute_texture_features(img, mask, 32);

            testCase.verifyGreaterThan(features.glcm_contrast, 0, ...
                'Checkerboard should have non-zero contrast.');
            testCase.verifyTrue(isfinite(features.glcm_homogeneity));
        end

        function testUniformImageZeroContrast(testCase)
            % Uniform image should have zero GLCM contrast.
            img = ones(32, 32) * 0.5;
            mask = true(32, 32);

            features = compute_texture_features(img, mask, 32);

            testCase.verifyEqual(features.glcm_contrast, 0, 'AbsTol', 1e-10, ...
                'Uniform image should have zero contrast.');
            testCase.verifyEqual(features.glcm_energy, 1, 'AbsTol', 1e-6, ...
                'Uniform image GLCM energy should be 1.');
        end

        function testOutputFieldCount(testCase)
            % Verify the output struct has all expected fields (now 24).
            img = rand(20, 20);
            mask = true(20, 20);

            features = compute_texture_features(img, mask);

            expected_fields = {'energy', 'uniformity', 'entropy', 'kurtosis', 'skewness', ...
                'p10', 'p90', 'iqr', 'mad', 'rmad', ...
                'glcm_contrast', 'glcm_correlation', 'glcm_energy', 'glcm_homogeneity', ...
                'glrlm_sre', 'glrlm_lre', 'glrlm_gln', 'glrlm_rln', 'glrlm_rp', ...
                'shape_volume', 'shape_surface_area', 'shape_sphericity', ...
                'shape_elongation', 'shape_compactness'};
            for i = 1:length(expected_fields)
                testCase.verifyTrue(isfield(features, expected_fields{i}), ...
                    sprintf('Missing field: %s', expected_fields{i}));
            end
            testCase.verifyEqual(length(fieldnames(features)), 24, ...
                'Should have exactly 24 texture features.');
        end

        function test3DInput(testCase)
            % 3D volume input should work (uses central slice for GLCM).
            img_3d = rand(20, 20, 10);
            mask_3d = true(20, 20, 10);

            features = compute_texture_features(img_3d, mask_3d, 16);

            testCase.verifyTrue(isfinite(features.glcm_contrast), ...
                'Should compute GLCM from 3D input.');
            testCase.verifyTrue(isfinite(features.energy));
        end

        function testEmptyMask(testCase)
            % Empty mask should return all NaN features.
            img = rand(20, 20);
            mask = false(20, 20);

            features = compute_texture_features(img, mask);

            testCase.verifyTrue(isnan(features.energy));
            testCase.verifyTrue(isnan(features.glcm_contrast));
        end

        function testFirstOrderStatistics(testCase)
            % Verify first-order statistics are computed correctly.
            rng(42);
            img = randn(50, 50) * 0.5 + 2;  % mean ~2, std ~0.5
            mask = true(50, 50);

            features = compute_texture_features(img, mask);

            testCase.verifyGreaterThan(features.energy, 0, 'Energy should be positive.');
            testCase.verifyGreaterThan(features.entropy, 0, 'Entropy should be positive.');
            testCase.verifyTrue(isfinite(features.iqr), 'IQR should be finite.');
            testCase.verifyTrue(features.p90 > features.p10, 'P90 should exceed P10.');
            testCase.verifyTrue(isfinite(features.uniformity), 'Uniformity should be finite.');
        end

        function testCheckerboardGLRLM(testCase)
            % Checkerboard pattern: all runs should be length 1, so
            % SRE should be high (near 1) and LRE should be low (near 1).
            img = zeros(64, 64);
            for i = 1:64
                for j = 1:64
                    if mod(i + j, 2) == 0
                        img(i, j) = 1;
                    end
                end
            end
            mask = true(64, 64);
            features = compute_texture_features(img, mask, 32);

            testCase.verifyGreaterThan(features.glrlm_sre, 0.8, ...
                'Checkerboard should have high SRE (all runs length 1).');
            testCase.verifyLessThan(features.glrlm_lre, 2.0, ...
                'Checkerboard should have low LRE (all runs length 1).');
        end

        function testSphereMaskSphericity(testCase)
            % A sphere mask should have sphericity near 1.0.
            sz = 51;
            center = (sz + 1) / 2;
            radius = 20;
            mask = false(sz, sz, sz);
            for x = 1:sz
                for y = 1:sz
                    for z = 1:sz
                        if (x-center)^2 + (y-center)^2 + (z-center)^2 <= radius^2
                            mask(x,y,z) = true;
                        end
                    end
                end
            end
            img = rand(sz, sz, sz);
            features = compute_texture_features(img, mask, 16, [1 1 1]);

            testCase.verifyGreaterThan(features.shape_sphericity, 0.8, ...
                'Sphere mask should have sphericity near 1.0.');
            testCase.verifyLessThanOrEqual(features.shape_sphericity, 1.0 + 0.05, ...
                'Sphericity should not exceed 1.0 significantly.');
        end

        function testElongatedEllipsoidLowSphericity(testCase)
            % An elongated ellipsoid should have low sphericity.
            sz_x = 10; sz_y = 10; sz_z = 60;
            mask = false(sz_x, sz_y, sz_z);
            cx = 5; cy = 5; cz = 30;
            rx = 3; ry = 3; rz = 25;
            for x = 1:sz_x
                for y = 1:sz_y
                    for z = 1:sz_z
                        if ((x-cx)/rx)^2 + ((y-cy)/ry)^2 + ((z-cz)/rz)^2 <= 1
                            mask(x,y,z) = true;
                        end
                    end
                end
            end
            img = rand(sz_x, sz_y, sz_z);
            features = compute_texture_features(img, mask, 16, [1 1 1]);

            testCase.verifyLessThan(features.shape_sphericity, 0.8, ...
                'Elongated ellipsoid should have low sphericity.');
        end

        function test3DCheckerboardGLRLM(testCase)
            % 3D checkerboard pattern: 3D GLRLM (13 directions) should
            % produce different features than 2D-only (4 directions)
            % because the inter-slice alternation creates additional
            % short runs visible only in the z and diagonal directions.
            sz = 16;
            img_3d = zeros(sz, sz, sz);
            for i = 1:sz
                for j = 1:sz
                    for k = 1:sz
                        if mod(i + j + k, 2) == 0
                            img_3d(i, j, k) = 1;
                        end
                    end
                end
            end
            mask_3d = true(sz, sz, sz);

            % 3D GLRLM (13 directions, default)
            f3d = compute_texture_features(img_3d, mask_3d, 32, [1 1 1], true);

            % 2D GLRLM (4 directions, forced)
            f2d = compute_texture_features(img_3d, mask_3d, 32, [1 1 1], false);

            % Both should produce finite GLRLM features
            testCase.verifyTrue(isfinite(f3d.glrlm_sre), ...
                '3D GLRLM SRE should be finite.');
            testCase.verifyTrue(isfinite(f2d.glrlm_sre), ...
                '2D GLRLM SRE should be finite.');

            % 3D and 2D should differ because 3D captures additional
            % directions (z-axis, face-diagonals, body-diagonals) that
            % contribute different run-length statistics.
            testCase.verifyNotEqual(f3d.glrlm_sre, f2d.glrlm_sre, ...
                '3D GLRLM SRE should differ from 2D for 3D checkerboard.');
            testCase.verifyNotEqual(f3d.glrlm_lre, f2d.glrlm_lre, ...
                '3D GLRLM LRE should differ from 2D for 3D checkerboard.');
            testCase.verifyNotEqual(f3d.glrlm_rp, f2d.glrlm_rp, ...
                '3D GLRLM RP should differ from 2D for 3D checkerboard.');
        end

        function test3DGLRLMSingleSliceFallback(testCase)
            % Single-slice 3D input should use 2D GLRLM even with texture_3d=true.
            rng(99);
            img_3d = rand(20, 20, 1);
            mask_3d = true(20, 20, 1);

            f_3d_flag = compute_texture_features(img_3d, mask_3d, 16, [1 1 1], true);
            f_2d_flag = compute_texture_features(img_3d, mask_3d, 16, [1 1 1], false);

            % With only 1 slice, both should produce identical GLRLM results
            testCase.verifyEqual(f_3d_flag.glrlm_sre, f_2d_flag.glrlm_sre, 'AbsTol', 1e-10, ...
                'Single-slice should fall back to 2D GLRLM regardless of texture_3d flag.');
            testCase.verifyEqual(f_3d_flag.glrlm_lre, f_2d_flag.glrlm_lre, 'AbsTol', 1e-10, ...
                'Single-slice LRE should match between 3D and 2D flags.');
        end

        function test3DGLRLMDisabledFlag(testCase)
            % When texture_3d=false, a multi-slice 3D input should use 2D GLRLM.
            rng(77);
            img_3d = rand(16, 16, 8);
            mask_3d = true(16, 16, 8);

            f_off = compute_texture_features(img_3d, mask_3d, 16, [1 1 1], false);

            % Should still produce finite features
            testCase.verifyTrue(isfinite(f_off.glrlm_sre), ...
                '2D fallback GLRLM should produce finite SRE for 3D input.');
            testCase.verifyTrue(isfinite(f_off.glrlm_gln), ...
                '2D fallback GLRLM should produce finite GLN for 3D input.');
        end

        function testShapeVolumePhysicalUnits(testCase)
            % Volume should scale with voxel spacing.
            mask = true(10, 10, 10);  % 1000 voxels
            img = rand(10, 10, 10);

            % Unit spacing
            f1 = compute_texture_features(img, mask, 16, [1 1 1]);
            testCase.verifyEqual(f1.shape_volume, 1000, 'AbsTol', 1e-6, ...
                'Volume should be 1000 mm^3 at unit spacing.');

            % 2mm spacing
            f2 = compute_texture_features(img, mask, 16, [2 2 2]);
            testCase.verifyEqual(f2.shape_volume, 8000, 'AbsTol', 1e-6, ...
                'Volume should be 8000 mm^3 at 2mm spacing.');
        end

    end
end
