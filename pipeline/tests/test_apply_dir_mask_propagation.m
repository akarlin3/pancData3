classdef test_apply_dir_mask_propagation < matlab.unittest.TestCase
    % TEST_APPLY_DIR_MASK_PROPAGATION Unit tests for Deformable Image Registration
    %
    % apply_dir_mask_propagation performs deformable image registration (DIR)
    % between two b0 DWI volumes and warps a GTV mask from the fixed frame
    % to the moving frame. It uses MATLAB's imregdemons with percentile-
    % based intensity normalization for robustness.
    %
    % Validates:
    %   - Identity registration (identical images -> near-zero displacement)
    %   - Mask warp produces a non-empty logical output
    %   - Empty input handling (each of the three inputs)
    %   - Dimension mismatch rejection (fixed vs moving, fixed vs mask)
    %   - Correct output types (logical mask, numeric displacement, imref3d)
    %   - Robust normalization handles bright outlier voxels

    methods(TestMethodSetup)
        function addPaths(testCase)
            % Add the utils directory containing apply_dir_mask_propagation.
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'utils'));
        end
    end

    methods(Test)

        function testIdentityRegistrationSmallDisplacement(testCase)
            % When fixed and moving images are identical, the deformable
            % registration displacement field should be approximately zero,
            % and the warped mask should closely overlap the original.
            rng(1);
            img = rand(32, 32, 8) * 100;  % 32x32x8 random volume
            mask = false(32, 32, 8);
            mask(10:20, 10:20, 3:6) = true;  % 11x11x4 rectangular mask

            [warped, D_forward, ref3d] = apply_dir_mask_propagation(img, img, mask);

            testCase.verifyNotEmpty(warped, 'Warped mask should not be empty.');
            testCase.verifyClass(warped, 'logical', 'Warped mask should be logical.');

            % Maximum displacement should be small (< 5 voxels) for identical images
            max_disp = max(abs(D_forward(:)));
            testCase.verifyLessThan(max_disp, 5, ...
                'Displacement field should be near zero for identical images.');

            % Dice coefficient between original and warped mask should be high
            intersection = sum(warped(:) & mask(:));
            union = sum(warped(:)) + sum(mask(:));
            dice = 2 * intersection / union;
            testCase.verifyGreaterThan(dice, 0.8, ...
                'Dice coefficient should be high for identity registration.');
        end

        function testEmptyFixedReturnsEmpty(testCase)
            % Empty b0_fixed input should trigger a warning and return
            % empty outputs rather than crashing.
            mask = true(5, 5, 5);
            moving = rand(5, 5, 5);

            testCase.verifyWarning(@() apply_dir_mask_propagation([], moving, mask), ...
                'apply_dir_mask_propagation:emptyInput');
        end

        function testEmptyMovingReturnsEmpty(testCase)
            % Empty b0_moving input should trigger a warning and return
            % empty outputs.
            fixed = rand(5, 5, 5);
            mask = true(5, 5, 5);

            testCase.verifyWarning(@() apply_dir_mask_propagation(fixed, [], mask), ...
                'apply_dir_mask_propagation:emptyInput');
        end

        function testEmptyMaskReturnsEmpty(testCase)
            % Empty gtv_mask input should trigger a warning and return
            % empty outputs (nothing to warp).
            fixed = rand(5, 5, 5);
            moving = rand(5, 5, 5);

            testCase.verifyWarning(@() apply_dir_mask_propagation(fixed, moving, []), ...
                'apply_dir_mask_propagation:emptyInput');
        end

        function testDimensionMismatchReturnsEmpty(testCase)
            % When fixed (5x5x5) and moving (6x6x6) have different
            % dimensions, registration cannot proceed. All outputs should
            % be empty. The sizeMismatch warning is suppressed for cleaner
            % test output and restored via onCleanup.
            fixed = rand(5, 5, 5);
            moving = rand(6, 6, 6);
            mask = true(5, 5, 5);

            warning('off', 'apply_dir_mask_propagation:sizeMismatch');
            cleanup = onCleanup(@() warning('on', 'apply_dir_mask_propagation:sizeMismatch'));

            [warped, D_forward, ref3d] = apply_dir_mask_propagation(fixed, moving, mask);

            testCase.verifyEmpty(warped, 'Warped mask should be empty for mismatched sizes.');
            testCase.verifyEmpty(D_forward, 'D_forward should be empty for mismatched sizes.');
            testCase.verifyEmpty(ref3d, 'ref3d should be empty for mismatched sizes.');
        end

        function testMaskFixedSizeMismatchReturnsEmpty(testCase)
            % When fixed and moving match (5x5x5) but the mask is a different
            % size (4x4x4), the function should return empty outputs because
            % the mask cannot be applied to the registration grid.
            fixed = rand(5, 5, 5);
            moving = rand(5, 5, 5);
            mask = true(4, 4, 4);

            warning('off', 'apply_dir_mask_propagation:sizeMismatch');
            cleanup = onCleanup(@() warning('on', 'apply_dir_mask_propagation:sizeMismatch'));

            [warped, ~, ~] = apply_dir_mask_propagation(fixed, moving, mask);
            testCase.verifyEmpty(warped, ...
                'Warped mask should be empty when mask size differs from fixed image.');
        end

        function testOutputTypes(testCase)
            % Verify that the three outputs have the correct MATLAB types:
            % warped mask -> logical, displacement field -> numeric array,
            % spatial reference -> imref3d object.
            rng(2);
            img = rand(16, 16, 8) * 100;
            mask = false(16, 16, 8);
            mask(5:10, 5:10, 2:6) = true;

            [warped, D_forward, ref3d] = apply_dir_mask_propagation(img, img, mask);

            testCase.verifyClass(warped, 'logical');
            testCase.verifyTrue(isnumeric(D_forward), 'D_forward should be numeric.');
            testCase.verifyClass(ref3d, 'imref3d');
        end

        function testRobustNormalizationWithOutliers(testCase)
            % A single extreme outlier voxel (1e6 vs ~100 range) should not
            % break registration. The percentile-based normalization in
            % apply_dir_mask_propagation clips such outliers before running
            % imregdemons, so registration should still succeed.
            rng(3);
            img = rand(20, 20, 8) * 100;
            img(1, 1, 1) = 1e6;  % extreme bright outlier
            mask = false(20, 20, 8);
            mask(5:15, 5:15, 2:7) = true;

            [warped, ~, ~] = apply_dir_mask_propagation(img, img, mask);

            testCase.verifyNotEmpty(warped, ...
                'Registration should succeed despite bright outlier voxels.');
        end

    end
end
