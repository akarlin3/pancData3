classdef test_compute_dice_hausdorff < matlab.unittest.TestCase
    % TEST_COMPUTE_DICE_HAUSDORFF Unit tests for Dice and Hausdorff utility.

    methods(TestMethodSetup)
        function addPaths(testCase)
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'utils'));
        end
    end

    methods(Test)

        function testIdenticalMasks(testCase)
            mask = false(10, 10, 10);
            mask(3:7, 3:7, 3:7) = true;
            [dice, hd_max, hd95] = compute_dice_hausdorff(mask, mask, [1 1 1]);
            testCase.verifyEqual(dice, 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(hd_max, 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(hd95, 0, 'AbsTol', 1e-12);
        end

        function testBothEmpty(testCase)
            mask = false(5, 5, 5);
            [dice, hd_max, hd95] = compute_dice_hausdorff(mask, mask, [1 1 1]);
            testCase.verifyTrue(isnan(dice), 'Both empty should give NaN Dice.');
            testCase.verifyTrue(isnan(hd_max), 'Both empty should give NaN HD max.');
            testCase.verifyTrue(isnan(hd95), 'Both empty should give NaN HD95.');
        end

        function testOneEmptyOneNonEmpty(testCase)
            mask_a = false(5, 5, 5);
            mask_b = false(5, 5, 5);
            mask_b(2:4, 2:4, 2:4) = true;
            [dice, hd_max, hd95] = compute_dice_hausdorff(mask_a, mask_b, [1 1 1]);
            testCase.verifyEqual(dice, 0);
            testCase.verifyEqual(hd_max, Inf);
            testCase.verifyEqual(hd95, Inf);
        end

        function testCompletelyDisjoint(testCase)
            mask_a = false(10, 10, 10);
            mask_b = false(10, 10, 10);
            mask_a(1:3, 1:3, 1:3) = true;
            mask_b(8:10, 8:10, 8:10) = true;
            [dice, hd_max, hd95] = compute_dice_hausdorff(mask_a, mask_b, [1 1 1]);
            testCase.verifyEqual(dice, 0);
            testCase.verifyGreaterThan(hd_max, 0, 'Disjoint masks should have positive HD.');
            testCase.verifyGreaterThan(hd95, 0);
        end

        function testPartialOverlap(testCase)
            % Two blocks sharing a 2-voxel overlap in one dimension
            mask_a = false(10, 10, 10);
            mask_b = false(10, 10, 10);
            mask_a(3:6, 3:6, 3:6) = true;  % 4x4x4 = 64 voxels
            mask_b(5:8, 3:6, 3:6) = true;  % 4x4x4 = 64 voxels
            % Overlap: rows 5:6, cols 3:6, slices 3:6 = 2x4x4 = 32 voxels
            [dice, ~, ~] = compute_dice_hausdorff(mask_a, mask_b, [1 1 1]);
            expected_dice = 2 * 32 / (64 + 64);
            testCase.verifyEqual(dice, expected_dice, 'AbsTol', 1e-12);
        end

        function testAnisotropicVoxels(testCase)
            % Two single-voxel masks separated by 1 voxel in the z direction
            mask_a = false(5, 5, 5);
            mask_b = false(5, 5, 5);
            mask_a(3, 3, 2) = true;
            mask_b(3, 3, 4) = true;
            % Separation = 2 voxels in z. With dz=5mm, distance = 10mm
            [~, hd_max, ~] = compute_dice_hausdorff(mask_a, mask_b, [2 2 5]);
            testCase.verifyEqual(hd_max, 10, 'AbsTol', 1e-10, ...
                'HD should be 2 voxels * 5mm/voxel = 10mm in z direction.');
        end

        function testAnisotropicVsIsotropic(testCase)
            % Same voxel-space separation, different physical distances
            mask_a = false(5, 5, 5);
            mask_b = false(5, 5, 5);
            mask_a(2, 3, 3) = true;
            mask_b(4, 3, 3) = true;
            % Separation = 2 voxels in x
            [~, hd_iso, ~] = compute_dice_hausdorff(mask_a, mask_b, [1 1 1]);
            [~, hd_aniso, ~] = compute_dice_hausdorff(mask_a, mask_b, [3 1 1]);
            testCase.verifyEqual(hd_iso, 2, 'AbsTol', 1e-10);
            testCase.verifyEqual(hd_aniso, 6, 'AbsTol', 1e-10, ...
                'HD should scale with voxel dimension.');
        end

        function testSingleVoxelSameLocation(testCase)
            mask = false(5, 5, 5);
            mask(3, 3, 3) = true;
            [dice, hd_max, hd95] = compute_dice_hausdorff(mask, mask, [2 2 5]);
            testCase.verifyEqual(dice, 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(hd_max, 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(hd95, 0, 'AbsTol', 1e-12);
        end

        function testDefaultVoxDims(testCase)
            % Verify default vox_dims = [1 1 1] when omitted
            mask_a = false(5, 5, 5);
            mask_b = false(5, 5, 5);
            mask_a(2, 3, 3) = true;
            mask_b(4, 3, 3) = true;
            [~, hd_default, ~] = compute_dice_hausdorff(mask_a, mask_b);
            [~, hd_explicit, ~] = compute_dice_hausdorff(mask_a, mask_b, [1 1 1]);
            testCase.verifyEqual(hd_default, hd_explicit, 'AbsTol', 1e-12);
        end

        function testSymmetry(testCase)
            mask_a = false(8, 8, 8);
            mask_b = false(8, 8, 8);
            mask_a(2:4, 2:4, 2:4) = true;
            mask_b(3:6, 3:6, 3:6) = true;
            [dice_ab, hd_ab, hd95_ab] = compute_dice_hausdorff(mask_a, mask_b, [1 1 1]);
            [dice_ba, hd_ba, hd95_ba] = compute_dice_hausdorff(mask_b, mask_a, [1 1 1]);
            testCase.verifyEqual(dice_ab, dice_ba, 'AbsTol', 1e-12, ...
                'Dice should be symmetric.');
            testCase.verifyEqual(hd_ab, hd_ba, 'AbsTol', 1e-12, ...
                'Hausdorff should be symmetric.');
            testCase.verifyEqual(hd95_ab, hd95_ba, 'AbsTol', 1e-12, ...
                'HD95 should be symmetric.');
        end

    end
end
