classdef test_apply_dir_mask_propagation < matlab.unittest.TestCase
    % TEST_APPLY_DIR_MASK_PROPAGATION Unit tests for Deformable Image Registration
    %
    % Validates:
    %   - Identity registration (identical images -> near-zero displacement)
    %   - Mask warp produces a non-empty logical output
    %   - Empty input handling
    %   - Dimension mismatch rejection
    %   - Robust normalization handles bright outlier voxels

    methods(TestMethodSetup)
        function addPaths(testCase)
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'utils'));
        end
    end

    methods(Test)

        function testIdentityRegistrationSmallDisplacement(testCase)
            % When fixed and moving images are identical, the displacement
            % field should be approximately zero.
            rng(1);
            img = rand(32, 32, 8) * 100;
            mask = false(32, 32, 8);
            mask(10:20, 10:20, 3:6) = true;

            [warped, D_forward, ref3d] = apply_dir_mask_propagation(img, img, mask);

            testCase.verifyNotEmpty(warped, 'Warped mask should not be empty.');
            testCase.verifyClass(warped, 'logical', 'Warped mask should be logical.');

            % Displacement field should be near zero for identical images
            max_disp = max(abs(D_forward(:)));
            testCase.verifyLessThan(max_disp, 5, ...
                'Displacement field should be near zero for identical images.');

            % Warped mask should have good Dice overlap with original
            intersection = sum(warped(:) & mask(:));
            union = sum(warped(:)) + sum(mask(:));
            dice = 2 * intersection / union;
            testCase.verifyGreaterThan(dice, 0.8, ...
                'Dice coefficient should be high for identity registration.');
        end

        function testEmptyFixedReturnsEmpty(testCase)
            % Empty b0_fixed should return empty outputs with a warning
            mask = true(5, 5, 5);
            moving = rand(5, 5, 5);

            testCase.verifyWarning(@() apply_dir_mask_propagation([], moving, mask), ...
                'apply_dir_mask_propagation:emptyInput');
        end

        function testEmptyMovingReturnsEmpty(testCase)
            % Empty b0_moving should return empty outputs with a warning
            fixed = rand(5, 5, 5);
            mask = true(5, 5, 5);

            testCase.verifyWarning(@() apply_dir_mask_propagation(fixed, [], mask), ...
                'apply_dir_mask_propagation:emptyInput');
        end

        function testEmptyMaskReturnsEmpty(testCase)
            % Empty gtv_mask should return empty outputs with a warning
            fixed = rand(5, 5, 5);
            moving = rand(5, 5, 5);

            testCase.verifyWarning(@() apply_dir_mask_propagation(fixed, moving, []), ...
                'apply_dir_mask_propagation:emptyInput');
        end

        function testDimensionMismatchReturnsEmpty(testCase)
            % Mismatched dimensions between fixed and moving should return empty
            fixed = rand(5, 5, 5);
            moving = rand(6, 6, 6);
            mask = true(5, 5, 5);

            [warped, D_forward, ref3d] = apply_dir_mask_propagation(fixed, moving, mask);

            testCase.verifyEmpty(warped, 'Warped mask should be empty for mismatched sizes.');
            testCase.verifyEmpty(D_forward, 'D_forward should be empty for mismatched sizes.');
            testCase.verifyEmpty(ref3d, 'ref3d should be empty for mismatched sizes.');
        end

        function testMaskFixedSizeMismatchReturnsEmpty(testCase)
            % Mismatched dimensions between fixed and mask should return empty
            fixed = rand(5, 5, 5);
            moving = rand(5, 5, 5);
            mask = true(4, 4, 4);

            [warped, ~, ~] = apply_dir_mask_propagation(fixed, moving, mask);
            testCase.verifyEmpty(warped, ...
                'Warped mask should be empty when mask size differs from fixed image.');
        end

        function testOutputTypes(testCase)
            % Verify correct types for all outputs
            rng(2);
            img = rand(16, 16, 4) * 100;
            mask = false(16, 16, 4);
            mask(5:10, 5:10, 2:3) = true;

            [warped, D_forward, ref3d] = apply_dir_mask_propagation(img, img, mask);

            testCase.verifyClass(warped, 'logical');
            testCase.verifyTrue(isnumeric(D_forward), 'D_forward should be numeric.');
            testCase.verifyClass(ref3d, 'imref3d');
        end

        function testRobustNormalizationWithOutliers(testCase)
            % An image with extreme outlier voxels should still register
            % successfully (percentile normalization clips outliers).
            rng(3);
            img = rand(20, 20, 6) * 100;
            % Add extreme outlier
            img(1, 1, 1) = 1e6;
            mask = false(20, 20, 6);
            mask(5:15, 5:15, 2:5) = true;

            [warped, ~, ~] = apply_dir_mask_propagation(img, img, mask);

            testCase.verifyNotEmpty(warped, ...
                'Registration should succeed despite bright outlier voxels.');
        end

    end
end
