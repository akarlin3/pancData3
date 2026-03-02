classdef test_parsave_dir_cache < matlab.unittest.TestCase
    % TEST_PARSAVE_DIR_CACHE Unit tests for parsave_dir_cache.
    %
    % parsave_dir_cache is a thin wrapper around save() that allows
    % variable-passing inside parfor loops.  Tests verify:
    %   - The output file is created on disk
    %   - All three variables (gtv_mask_warped, D_forward, ref3d) are
    %     faithfully written to the .mat file
    %   - A second call to the same path overwrites the previous file
    %   - The function accepts arrays of varied sizes

    properties
        TempDir
        OriginalPath
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testFileIsCreated(testCase)
            % The .mat file should exist after the call.
            fname = fullfile(testCase.TempDir, 'cache_test.mat');
            parsave_dir_cache(fname, rand(5,5,5) > 0.5, rand(5,5,5,3), rand(5,5,5));

            testCase.verifyTrue(exist(fname, 'file') > 0, ...
                'parsave_dir_cache should create the output .mat file.');
        end

        function testVariablesSavedWithCorrectValues(testCase)
            % Load the saved file and verify all three arrays round-trip correctly.
            fname = fullfile(testCase.TempDir, 'cache_vars.mat');
            rng(7);
            gtv_mask_warped = rand(4,4,4) > 0.5;
            D_forward       = rand(4,4,4,3);
            ref3d           = rand(4,4,4);

            parsave_dir_cache(fname, gtv_mask_warped, D_forward, ref3d);
            loaded = load(fname);

            testCase.verifyEqual(loaded.gtv_mask_warped, gtv_mask_warped, ...
                'gtv_mask_warped should match the value passed to parsave_dir_cache.');
            testCase.verifyEqual(loaded.D_forward, D_forward, ...
                'D_forward should match the value passed to parsave_dir_cache.');
            testCase.verifyEqual(loaded.ref3d, ref3d, ...
                'ref3d should match the value passed to parsave_dir_cache.');
        end

        function testVariableNamesInFile(testCase)
            % The saved file must expose variables named exactly
            % 'gtv_mask_warped', 'D_forward', and 'ref3d'.
            fname = fullfile(testCase.TempDir, 'cache_names.mat');
            parsave_dir_cache(fname, true(2,2,2), zeros(2,2,2,3), ones(2,2,2));

            info = whos('-file', fname);
            names = {info.name};

            testCase.verifyTrue(ismember('gtv_mask_warped', names), ...
                'Saved file must contain variable ''gtv_mask_warped''.');
            testCase.verifyTrue(ismember('D_forward', names), ...
                'Saved file must contain variable ''D_forward''.');
            testCase.verifyTrue(ismember('ref3d', names), ...
                'Saved file must contain variable ''ref3d''.');
        end

        function testOverwritesPreviousFile(testCase)
            % A second call with different data should overwrite the first.
            fname = fullfile(testCase.TempDir, 'cache_overwrite.mat');

            % First write
            parsave_dir_cache(fname, true(3,3,3), ones(3,3,3,3), zeros(3,3,3));
            % Second write with different data
            gtv2  = false(3,3,3);
            Dfwd2 = 2 * ones(3,3,3,3);
            ref2  = 5 * ones(3,3,3);
            parsave_dir_cache(fname, gtv2, Dfwd2, ref2);

            loaded = load(fname);
            testCase.verifyEqual(loaded.gtv_mask_warped, gtv2, ...
                'Second call should overwrite gtv_mask_warped.');
            testCase.verifyEqual(loaded.D_forward, Dfwd2, ...
                'Second call should overwrite D_forward.');
            testCase.verifyEqual(loaded.ref3d, ref2, ...
                'Second call should overwrite ref3d.');
        end

        function testAcceptsVariedArraySizes(testCase)
            % Arrays with different spatial dimensions should all save correctly.
            fname = fullfile(testCase.TempDir, 'cache_sizes.mat');
            gtv   = rand(8, 10, 6) > 0.5;
            Dfwd  = rand(8, 10, 6, 3);
            ref   = rand(8, 10, 6);

            parsave_dir_cache(fname, gtv, Dfwd, ref);
            loaded = load(fname);

            testCase.verifyEqual(size(loaded.gtv_mask_warped), size(gtv), ...
                'Saved gtv_mask_warped size should match input.');
            testCase.verifyEqual(size(loaded.D_forward), size(Dfwd), ...
                'Saved D_forward size should match input.');
            testCase.verifyEqual(size(loaded.ref3d), size(ref), ...
                'Saved ref3d size should match input.');
        end

        function testEmptyArraysSavedCorrectly(testCase)
            % Edge case: all-empty arrays should be saved and loaded without error.
            fname = fullfile(testCase.TempDir, 'cache_empty.mat');
            gtv   = false(0,0,0);
            Dfwd  = zeros(0,0,0,3);
            ref   = zeros(0,0,0);

            parsave_dir_cache(fname, gtv, Dfwd, ref);
            loaded = load(fname);

            testCase.verifyEqual(size(loaded.gtv_mask_warped), [0 0 0]);
            testCase.verifyEqual(size(loaded.D_forward),       [0 0 0 3]);
            testCase.verifyEqual(size(loaded.ref3d),           [0 0 0]);
        end

    end
end
