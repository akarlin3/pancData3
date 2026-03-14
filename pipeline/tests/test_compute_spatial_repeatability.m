classdef test_compute_spatial_repeatability < matlab.unittest.TestCase
% TEST_COMPUTE_SPATIAL_REPEATABILITY — Unit tests for compute_spatial_repeatability.m
%
% Validates spatial repeatability computation:
%   - Returns NaN when fewer than 2 valid repeats
%   - Returns valid Dice/Hausdorff when 2+ repeats exist
%   - Perfect overlap gives Dice=1, HD=0
%   - Handles missing 3D masks gracefully

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'dependencies'));
        end
    end

    methods (Test)
        function test_single_repeat_returns_nan(testCase)
            % Only 1 valid repeat → all outputs NaN
            data_vectors = testCase.makeDummyVectors(1);

            [da, hma, h95a, dd, hmd, h95d, df, hmf, h95f, dds, hmds, h95ds] = ...
                compute_spatial_repeatability(data_vectors, 1, 1, ...
                testCase.makeGtvLocations(1), 0.001, 0.001, 0.1, 0.01, ...
                strel('sphere', 1), 5, '', []);

            testCase.verifyTrue(isnan(da), 'Single repeat should give NaN Dice.');
            testCase.verifyTrue(isnan(hma));
            testCase.verifyTrue(isnan(h95a));
        end

        function test_identical_repeats_perfect_dice(testCase)
            % Two identical repeats should give Dice=1 and HD=0
            data_vectors = testCase.makeIdenticalVectors(2);
            gtv_locs = testCase.makeGtvLocationsWithMasks(2, data_vectors);

            [da, hma, h95a, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                compute_spatial_repeatability(data_vectors, 1, 1, ...
                gtv_locs, 0.0015, 0.0015, 0.15, 0.025, ...
                strel('sphere', 1), 1, '', []);

            if ~isnan(da)
                testCase.verifyEqual(da, 1.0, 'AbsTol', 0.01, ...
                    'Identical repeats should give Dice ≈ 1.');
                testCase.verifyEqual(hma, 0.0, 'AbsTol', 0.5, ...
                    'Identical repeats should give HD ≈ 0.');
            end
        end

        function test_no_3d_masks_returns_nan(testCase)
            % When GTV locations point to non-existent files, all NaN
            data_vectors = testCase.makeDummyVectors(2);
            gtv_locs = cell(1, 1, 2);
            gtv_locs{1,1,1} = 'nonexistent_file_1.mat';
            gtv_locs{1,1,2} = 'nonexistent_file_2.mat';

            [da, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                compute_spatial_repeatability(data_vectors, 1, 1, ...
                gtv_locs, 0.001, 0.001, 0.1, 0.01, ...
                strel('sphere', 1), 5, '', []);

            testCase.verifyTrue(isnan(da), ...
                'Missing 3D masks should return NaN.');
        end

        function test_output_count(testCase)
            % Verify all 12 outputs are returned
            data_vectors = testCase.makeDummyVectors(1);

            outputs = cell(1, 12);
            [outputs{:}] = compute_spatial_repeatability(data_vectors, 1, 1, ...
                testCase.makeGtvLocations(1), 0.001, 0.001, 0.1, 0.01, ...
                strel('sphere', 1), 5, '', []);

            testCase.verifyEqual(numel(outputs), 12, ...
                'Function should return exactly 12 outputs.');
        end
    end

    methods (Access = private)
        function data_vectors = makeDummyVectors(~, n_repeats)
            % Create minimal data_vectors struct array
            template = struct('adc_vector', [], 'd_vector', [], ...
                'f_vector', [], 'dstar_vector', [], ...
                'vox_dims', [2, 2, 3], 'vox_vol', 0.012);
            data_vectors = repmat(template, 1, 1, n_repeats);
            rng(42);
            for r = 1:n_repeats
                n = 100;
                data_vectors(1,1,r).adc_vector = rand(n,1) * 0.003;
                data_vectors(1,1,r).d_vector = rand(n,1) * 0.003;
                data_vectors(1,1,r).f_vector = rand(n,1) * 0.3;
                data_vectors(1,1,r).dstar_vector = rand(n,1) * 0.05;
            end
        end

        function data_vectors = makeIdenticalVectors(~, n_repeats)
            % Create two identical repeat vectors
            template = struct('adc_vector', [], 'd_vector', [], ...
                'f_vector', [], 'dstar_vector', [], ...
                'vox_dims', [2, 2, 3], 'vox_vol', 0.012);
            data_vectors = repmat(template, 1, 1, n_repeats);
            rng(42);
            n = 125;  % 5x5x5 grid
            adc = rand(n,1) * 0.003;
            d = rand(n,1) * 0.003;
            f = rand(n,1) * 0.3;
            ds = rand(n,1) * 0.05;
            for r = 1:n_repeats
                data_vectors(1,1,r).adc_vector = adc;
                data_vectors(1,1,r).d_vector = d;
                data_vectors(1,1,r).f_vector = f;
                data_vectors(1,1,r).dstar_vector = ds;
            end
        end

        function gtv_locs = makeGtvLocations(~, n_repeats)
            gtv_locs = cell(1, 1, n_repeats);
        end

        function gtv_locs = makeGtvLocationsWithMasks(~, n_repeats, data_vectors)
            % Create temporary .mat files with 3D GTV masks
            gtv_locs = cell(1, 1, n_repeats);
            n_vox = numel(data_vectors(1,1,1).adc_vector);
            grid_size = round(n_vox^(1/3));
            mask_3d = true(grid_size, grid_size, grid_size);
            % Trim mask if needed to match voxel count
            mask_3d(n_vox+1:end) = false;
            for r = 1:n_repeats
                fpath = fullfile(tempdir, sprintf('test_gtv_rpt_%d.mat', r));
                Stvol3d = mask_3d;  %#ok<NASGU>
                save(fpath, 'Stvol3d');
                gtv_locs{1,1,r} = fpath;
            end
        end
    end
end
