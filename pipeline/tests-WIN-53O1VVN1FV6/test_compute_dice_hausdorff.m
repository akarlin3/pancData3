classdef test_compute_dice_hausdorff < matlab.unittest.TestCase
    % TEST_COMPUTE_DICE_HAUSDORFF Unit tests for Dice and Hausdorff utility.
    %
    % compute_dice_hausdorff computes three spatial overlap/distance metrics
    % between two 3D binary masks:
    %   - Dice coefficient: 2*|A & B| / (|A| + |B|), range [0, 1]
    %   - Hausdorff distance (max): worst-case surface-to-surface distance
    %   - HD95: 95th percentile of surface distances (robust to outliers)
    %
    % These metrics are used in compare_core_methods.m to quantify agreement
    % between different tumor core delineation algorithms. The vox_dims
    % parameter converts voxel-space distances to physical (mm) distances,
    % which is critical for anisotropic MRI acquisitions.
    %
    % Tests cover:
    %   - Perfect overlap (identical masks)
    %   - Both masks empty (NaN results)
    %   - One empty, one non-empty (Dice=0, HD=Inf)
    %   - Disjoint masks (Dice=0, positive HD)
    %   - Partial overlap with known Dice value
    %   - Anisotropic voxel dimensions
    %   - Single-voxel masks
    %   - Default vox_dims when argument omitted
    %   - Symmetry: swap(A,B) must give identical results
    %
    % Run tests with:
    %   results = runtests('tests/test_compute_dice_hausdorff.m');

    methods(TestMethodSetup)
        function addPaths(testCase)
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'utils'));
        end
    end

    methods(Test)

        function testIdenticalMasks(testCase)
            % Perfect overlap: Dice must be 1, both HD metrics must be 0.
            mask = false(10, 10, 10);
            mask(3:7, 3:7, 3:7) = true; % 5x5x5 = 125-voxel cube
            [dice, hd_max, hd95] = compute_dice_hausdorff(mask, mask, [1 1 1]);
            testCase.verifyEqual(dice, 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(hd_max, 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(hd95, 0, 'AbsTol', 1e-12);
        end

        function testBothEmpty(testCase)
            % When both masks are empty, Dice is undefined (0/0 = NaN).
            % HD is also undefined since there are no surface points.
            mask = false(5, 5, 5);
            [dice, hd_max, hd95] = compute_dice_hausdorff(mask, mask, [1 1 1]);
            testCase.verifyTrue(isnan(dice), 'Both empty should give NaN Dice.');
            testCase.verifyTrue(isnan(hd_max), 'Both empty should give NaN HD max.');
            testCase.verifyTrue(isnan(hd95), 'Both empty should give NaN HD95.');
        end

        function testOneEmptyOneNonEmpty(testCase)
            % One empty mask and one non-empty: zero overlap, infinite
            % distance (every point in B is infinitely far from A's empty surface).
            mask_a = false(5, 5, 5);
            mask_b = false(5, 5, 5);
            mask_b(2:4, 2:4, 2:4) = true;
            [dice, hd_max, hd95] = compute_dice_hausdorff(mask_a, mask_b, [1 1 1]);
            testCase.verifyEqual(dice, 0);
            testCase.verifyEqual(hd_max, Inf);
            testCase.verifyEqual(hd95, Inf);
        end

        function testCompletelyDisjoint(testCase)
            % Two non-overlapping cubes at opposite corners of the volume.
            % Dice must be 0 (no intersection). HD must be positive (the
            % distance between the closest surfaces of the two cubes).
            mask_a = false(10, 10, 10);
            mask_b = false(10, 10, 10);
            mask_a(1:3, 1:3, 1:3) = true; % Corner near origin
            mask_b(8:10, 8:10, 8:10) = true; % Opposite corner
            [dice, hd_max, hd95] = compute_dice_hausdorff(mask_a, mask_b, [1 1 1]);
            testCase.verifyEqual(dice, 0);
            testCase.verifyGreaterThan(hd_max, 0, 'Disjoint masks should have positive HD.');
            testCase.verifyGreaterThan(hd95, 0);
        end

        function testPartialOverlap(testCase)
            % Two 4x4x4 cubes shifted by 2 voxels in the row direction.
            % Overlap region: rows 5:6, cols 3:6, slices 3:6 = 2x4x4 = 32 voxels.
            % Expected Dice = 2 * 32 / (64 + 64) = 0.5
            mask_a = false(10, 10, 10);
            mask_b = false(10, 10, 10);
            mask_a(3:6, 3:6, 3:6) = true;  % 4x4x4 = 64 voxels
            mask_b(5:8, 3:6, 3:6) = true;  % 4x4x4 = 64 voxels
            [dice, ~, ~] = compute_dice_hausdorff(mask_a, mask_b, [1 1 1]);
            expected_dice = 2 * 32 / (64 + 64); % = 0.5
            testCase.verifyEqual(dice, expected_dice, 'AbsTol', 1e-12);
        end

        function testAnisotropicVoxels(testCase)
            % Two single-voxel masks separated by 2 voxels in the z direction.
            % With anisotropic voxel size [2, 2, 5] mm, the z-separation of
            % 2 voxels corresponds to 2 * 5 = 10 mm physical distance.
            mask_a = false(5, 5, 5);
            mask_b = false(5, 5, 5);
            mask_a(3, 3, 2) = true;
            mask_b(3, 3, 4) = true;
            [~, hd_max, ~] = compute_dice_hausdorff(mask_a, mask_b, [2 2 5]);
            testCase.verifyEqual(hd_max, 10, 'AbsTol', 1e-10, ...
                'HD should be 2 voxels * 5mm/voxel = 10mm in z direction.');
        end

        function testAnisotropicVsIsotropic(testCase)
            % Demonstrates that voxel dimensions scale the HD correctly.
            % Same 2-voxel separation in x, but [1 1 1] gives HD=2 while
            % [3 1 1] gives HD=6 (2 voxels * 3 mm/voxel).
            mask_a = false(5, 5, 5);
            mask_b = false(5, 5, 5);
            mask_a(2, 3, 3) = true;
            mask_b(4, 3, 3) = true;
            [~, hd_iso, ~] = compute_dice_hausdorff(mask_a, mask_b, [1 1 1]);
            [~, hd_aniso, ~] = compute_dice_hausdorff(mask_a, mask_b, [3 1 1]);
            testCase.verifyEqual(hd_iso, 2, 'AbsTol', 1e-10);
            testCase.verifyEqual(hd_aniso, 6, 'AbsTol', 1e-10, ...
                'HD should scale with voxel dimension.');
        end

        function testSingleVoxelSameLocation(testCase)
            % A single voxel compared to itself: perfect Dice, zero distance.
            % Also verifies correct behaviour with anisotropic voxel dims.
            mask = false(5, 5, 5);
            mask(3, 3, 3) = true;
            [dice, hd_max, hd95] = compute_dice_hausdorff(mask, mask, [2 2 5]);
            testCase.verifyEqual(dice, 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(hd_max, 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(hd95, 0, 'AbsTol', 1e-12);
        end

        function testDefaultVoxDims(testCase)
            % When vox_dims is omitted, the function should default to
            % [1 1 1] (isotropic unit voxels). Results must match an
            % explicit [1 1 1] call.
            mask_a = false(5, 5, 5);
            mask_b = false(5, 5, 5);
            mask_a(2, 3, 3) = true;
            mask_b(4, 3, 3) = true;
            [~, hd_default, ~] = compute_dice_hausdorff(mask_a, mask_b);
            [~, hd_explicit, ~] = compute_dice_hausdorff(mask_a, mask_b, [1 1 1]);
            testCase.verifyEqual(hd_default, hd_explicit, 'AbsTol', 1e-12);
        end

        function testSymmetry(testCase)
            % Dice, HD, and HD95 are all symmetric metrics: swapping
            % the two mask arguments must produce identical results.
            mask_a = false(8, 8, 8);
            mask_b = false(8, 8, 8);
            mask_a(2:4, 2:4, 2:4) = true; % 3x3x3 = 27 voxels
            mask_b(3:6, 3:6, 3:6) = true; % 4x4x4 = 64 voxels (asymmetric sizes)
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
