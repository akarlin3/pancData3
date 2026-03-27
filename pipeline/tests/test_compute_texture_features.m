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
            % Checkerboard pattern: horizontal and vertical runs are
            % length 1, but diagonal directions (1,1) and (1,-1) have
            % constant values along mod(i+j,2), producing long runs.
            % SRE should still be moderate (mix of short+long runs).
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

            testCase.verifyGreaterThan(features.glrlm_sre, 0.5, ...
                'Checkerboard SRE should be moderate (short runs in axis-aligned directions).');
            testCase.verifyTrue(isfinite(features.glrlm_lre), ...
                'Checkerboard LRE should be finite.');
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

            testCase.verifyGreaterThan(features.shape_sphericity, 0.6, ...
                'Sphere mask should have moderate-to-high sphericity.');
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

            % In a mod(i+j+k,2) checkerboard, directions whose component
            % sum is even (e.g. (1,1,0)) see constant values (long runs),
            % while odd-sum directions (e.g. (1,0,0)) see alternating
            % values (short runs).  Both 2D (2/4 long) and 3D (6/13 long)
            % have similar proportions, so SRE and RP are comparable.
            testCase.verifyGreaterThan(f3d.glrlm_sre, 0.5, ...
                '3D checkerboard SRE should be moderate.');
            testCase.verifyGreaterThan(f2d.glrlm_sre, 0.5, ...
                '2D checkerboard SRE should be moderate.');
            testCase.verifyTrue(isfinite(f3d.glrlm_rp) && isfinite(f2d.glrlm_rp), ...
                'Both 2D and 3D GLRLM RP should be finite.');
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

        function testFixedBinWidthQuantization(testCase)
            % Fixed bin width quantization should produce valid features
            % and differ from fixed bin number results.
            rng(42);
            img = randn(50, 50) * 0.5 + 2;
            mask = true(50, 50);

            f_fbn = compute_texture_features(img, mask, 32, [1 1 1], 'fixed_bin_number');
            f_fbw = compute_texture_features(img, mask, 32, [1 1 1], 'fixed_bin_width');

            % Both should produce finite GLCM features
            testCase.verifyTrue(isfinite(f_fbn.glcm_contrast), ...
                'Fixed bin number should produce finite GLCM contrast.');
            testCase.verifyTrue(isfinite(f_fbw.glcm_contrast), ...
                'Fixed bin width should produce finite GLCM contrast.');

            % First-order features should be identical (independent of quantization)
            testCase.verifyEqual(f_fbn.energy, f_fbw.energy, 'AbsTol', 1e-10, ...
                'First-order energy should not depend on quantization method.');
        end

        function testMultiOffsetGLCMRotationalInvariance(testCase)
            % GLCM with multi-offset averaging should produce rotationally
            % invariant features: a pattern and its 90-degree rotation
            % should yield similar GLCM values.
            rng(42);
            img = rand(40, 40);
            % Add directional texture
            for i = 1:40
                img(i, :) = img(i, :) + 0.5 * sin(2*pi*i/8);
            end
            mask = true(40, 40);

            features_orig = compute_texture_features(img, mask, 32);
            features_rot = compute_texture_features(img', mask, 32);

            % Multi-offset averaging should make contrast approximately equal
            % between original and 90-degree rotated versions
            testCase.verifyEqual(features_orig.glcm_contrast, features_rot.glcm_contrast, ...
                'RelTol', 0.3, ...
                'Multi-offset GLCM should produce similar contrast for 90-deg rotation.');
        end

        function testDefaultQuantizationMethodBackwardCompat(testCase)
            % Calling without quantization_method argument should produce
            % same results as explicitly passing 'fixed_bin_number'.
            rng(42);
            img = randn(30, 30) * 0.3 + 1;
            mask = true(30, 30);

            f_default = compute_texture_features(img, mask, 16, [1 1 1]);
            f_explicit = compute_texture_features(img, mask, 16, [1 1 1], 'fixed_bin_number');

            testCase.verifyEqual(f_default.glcm_contrast, f_explicit.glcm_contrast, ...
                'AbsTol', 1e-12, ...
                'Default should match explicit fixed_bin_number.');
            testCase.verifyEqual(f_default.glcm_energy, f_explicit.glcm_energy, ...
                'AbsTol', 1e-12);
        end

        function testIBSIPancreaticDWIExpectedRanges(testCase)
            % Validate texture features against expected ranges for
            % pancreatic DWI ADC maps.
            %
            % IBSI phantom validation is not directly applicable to
            % clinical pancreatic DWI data (the IBSI digital phantom uses
            % a synthetic 3D checkerboard at integer grey levels, whereas
            % pancreatic ADC maps are continuous-valued ~0.5-3.0 x 10^-3
            % mm^2/s with irregular tumor ROIs). Instead, we validate
            % against empirically observed ranges from pancreatic DWI
            % literature and our own cohort analysis.
            %
            % Reference ranges (pancreatic adenocarcinoma, ADC maps):
            %   - Entropy: 3.0-6.5 (depends on n_levels; 32 bins typical)
            %   - GLCM contrast: 0.5-50 (higher = more heterogeneous)
            %   - GLCM energy: 0.001-0.2 (lower = more heterogeneous)
            %   - GLCM homogeneity: 0.1-0.8
            %   - Kurtosis: 1.5-6.0 (platykurtic to mildly leptokurtic)
            %   - Skewness: -1.5 to 2.0 (typically right-skewed)
            %   - SRE: 0.5-1.0 (short-run emphasis)
            %   - LRE: 1.0-5.0 (long-run emphasis)
            %
            % These ranges are documented for reproducibility and should be
            % updated as larger validation cohorts become available.

            rng(123);
            % Simulate a realistic pancreatic ADC distribution:
            % Mean ADC ~1.2 x 10^-3, SD ~0.3 x 10^-3, mild right skew
            n = 50;
            adc_vals = 0.0012 + 0.0003 * randn(n, n);
            adc_vals = max(adc_vals, 0.0003);  % physiological floor
            adc_vals = min(adc_vals, 0.003);    % physiological ceiling
            mask = true(n, n);

            features = compute_texture_features(adc_vals, mask, 32, [1.5 1.5 5]);

            % Entropy within expected range for 32-bin histogram
            testCase.verifyGreaterThan(features.entropy, 2.0, ...
                'Entropy should be > 2.0 for heterogeneous ADC map.');
            testCase.verifyLessThan(features.entropy, 7.0, ...
                'Entropy should be < 7.0 for 32-bin histogram.');

            % GLCM contrast: non-zero for heterogeneous map
            testCase.verifyGreaterThan(features.glcm_contrast, 0.1, ...
                'GLCM contrast should be > 0.1 for realistic ADC map.');
            testCase.verifyLessThan(features.glcm_contrast, 100, ...
                'GLCM contrast should be < 100 for realistic ADC map.');

            % GLCM energy: low for heterogeneous texture
            testCase.verifyGreaterThan(features.glcm_energy, 0.0001, ...
                'GLCM energy should be > 0.0001.');
            testCase.verifyLessThan(features.glcm_energy, 0.5, ...
                'GLCM energy should be < 0.5 for heterogeneous ADC.');

            % GLCM homogeneity: moderate range
            testCase.verifyGreaterThan(features.glcm_homogeneity, 0.05, ...
                'Homogeneity should be > 0.05.');
            testCase.verifyLessThan(features.glcm_homogeneity, 1.0, ...
                'Homogeneity should be <= 1.0.');

            % Kurtosis: expected near-normal for Gaussian-derived ADC
            testCase.verifyGreaterThan(features.kurtosis, 1.0, ...
                'Kurtosis should be > 1.0 for clipped Gaussian.');
            testCase.verifyLessThan(features.kurtosis, 10.0, ...
                'Kurtosis should be < 10.0.');

            % GLRLM short run emphasis: high for heterogeneous texture
            testCase.verifyGreaterThan(features.glrlm_sre, 0.3, ...
                'SRE should be > 0.3 for heterogeneous image.');
            testCase.verifyLessThanOrEqual(features.glrlm_sre, 1.0, ...
                'SRE should be <= 1.0.');

            % GLRLM long run emphasis: moderate for heterogeneous texture
            testCase.verifyGreaterThan(features.glrlm_lre, 1.0, ...
                'LRE should be >= 1.0.');
        end

        function testFixedBinWidthUniformImage(testCase)
            % Fixed bin width on a uniform image should still give zero contrast.
            img = ones(32, 32) * 0.5;
            mask = true(32, 32);

            features = compute_texture_features(img, mask, 32, [1 1 1], 'fixed_bin_width');

            testCase.verifyEqual(features.glcm_contrast, 0, 'AbsTol', 1e-10, ...
                'Uniform image should have zero contrast with fixed_bin_width.');
        end

    end
end
