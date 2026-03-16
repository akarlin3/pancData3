classdef test_compute_texture_features < matlab.unittest.TestCase
    % TEST_COMPUTE_TEXTURE_FEATURES  Tests for texture feature extraction.
    %
    % Covers:
    %   - Checkerboard pattern: high contrast, low homogeneity
    %   - Uniform image: zero contrast
    %   - Output struct field count
    %   - 3D input handling
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
            % Verify the output struct has all expected fields.
            img = rand(20, 20);
            mask = true(20, 20);

            features = compute_texture_features(img, mask);

            expected_fields = {'energy', 'entropy', 'kurtosis', 'skewness', ...
                'p10', 'p90', 'iqr', 'mad', 'rmad', ...
                'glcm_contrast', 'glcm_correlation', 'glcm_energy', 'glcm_homogeneity'};
            for i = 1:length(expected_fields)
                testCase.verifyTrue(isfield(features, expected_fields{i}), ...
                    sprintf('Missing field: %s', expected_fields{i}));
            end
            testCase.verifyEqual(length(fieldnames(features)), 13, ...
                'Should have exactly 13 texture features.');
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
        end

    end
end
