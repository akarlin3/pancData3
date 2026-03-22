classdef test_compute_registration_quality < matlab.unittest.TestCase
    % TEST_COMPUTE_REGISTRATION_QUALITY  Tests for registration quality metrics.
    %
    % Covers:
    %   - Identity transform (perfect NCC, J=1)
    %   - Known affine shift
    %   - Degenerate (zero) deformation field
    %   - Missing deformation field

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

        function testIdentityTransform(testCase)
            % Identity transform: NCC should be 1, Jacobian should be 1 everywhere.
            rng(42);
            vol = rand(20, 20, 10);

            % Identity deformation field (zero displacement)
            def_field = zeros(20, 20, 10, 3);

            quality = compute_registration_quality(vol, vol, def_field);

            testCase.verifyEqual(quality.ncc, 1.0, 'AbsTol', 1e-10, ...
                'NCC should be 1.0 for identical volumes.');
            testCase.verifyEqual(quality.jacobian_mean, 1.0, 'AbsTol', 0.2, ...
                'Mean Jacobian should be ~1 for identity transform.');
            testCase.verifyEqual(quality.jacobian_folding_pct, 0, 'AbsTol', 1, ...
                'No Jacobian folding expected for identity transform.');
        end

        function testKnownShift(testCase)
            % Known uniform shift: NCC should still be high for shifted versions.
            rng(42);
            vol = rand(20, 20, 10);
            shifted = circshift(vol, [1, 0, 0]);

            quality = compute_registration_quality(vol, shifted, []);

            testCase.verifyGreaterThan(quality.ncc, 0.5, ...
                'NCC should be reasonably high for small shift.');
            testCase.verifyTrue(isfinite(quality.mutual_information), ...
                'MI should be finite.');
        end

        function testDegenerateDeformationField(testCase)
            % Zero deformation field: all displacements are zero.
            vol = rand(10, 10, 5);
            def_field = zeros(10, 10, 5, 3);

            quality = compute_registration_quality(vol, vol, def_field);

            testCase.verifyTrue(isfinite(quality.jacobian_mean), ...
                'Should compute Jacobian from zero deformation field.');
            testCase.verifyEqual(quality.ncc, 1.0, 'AbsTol', 1e-10);
        end

        function testNoDeformationField(testCase)
            % Missing deformation field: Jacobian metrics should be NaN.
            vol = rand(10, 10, 5);

            quality = compute_registration_quality(vol, vol, []);

            testCase.verifyTrue(isnan(quality.jacobian_mean), ...
                'Jacobian mean should be NaN without deformation field.');
            testCase.verifyTrue(isnan(quality.jacobian_folding_pct));
            testCase.verifyEqual(quality.ncc, 1.0, 'AbsTol', 1e-10, ...
                'NCC should still be computed.');
        end

        function testMutualInformation(testCase)
            % MI between identical volumes should be higher than between random ones.
            rng(42);
            vol = rand(15, 15, 8);
            random_vol = rand(15, 15, 8);

            q_same = compute_registration_quality(vol, vol, []);
            q_diff = compute_registration_quality(vol, random_vol, []);

            testCase.verifyGreaterThan(q_same.mutual_information, q_diff.mutual_information, ...
                'MI for identical volumes should exceed MI for random volumes.');
        end

        function testIdentityAnisotropicSpacing(testCase)
            % Identity displacement field with anisotropic spacing:
            % Jacobian determinant should be exactly 1.0 everywhere.
            sz = [20, 15, 10];
            rng(42);
            vol = rand(sz);

            % Zero displacement field (identity transform)
            def_field = zeros([sz, 3]);

            % Anisotropic voxel spacing: 0.5mm x 1.0mm x 2.5mm
            voxel_spacing = [0.5, 1.0, 2.5];

            quality = compute_registration_quality(vol, vol, def_field, voxel_spacing);

            testCase.verifyEqual(quality.jacobian_mean, 1.0, 'AbsTol', 1e-10, ...
                'Mean Jacobian should be exactly 1.0 for identity with anisotropic spacing.');
            testCase.verifyEqual(quality.jacobian_std, 0.0, 'AbsTol', 1e-10, ...
                'Jacobian std should be 0 for identity with anisotropic spacing.');
            testCase.verifyEqual(quality.jacobian_min, 1.0, 'AbsTol', 1e-10, ...
                'Min Jacobian should be 1.0 for identity with anisotropic spacing.');
            testCase.verifyEqual(quality.jacobian_max, 1.0, 'AbsTol', 1e-10, ...
                'Max Jacobian should be 1.0 for identity with anisotropic spacing.');
            testCase.verifyEqual(quality.jacobian_folding_pct, 0, ...
                'No folding for identity transform.');
        end

    end
end
