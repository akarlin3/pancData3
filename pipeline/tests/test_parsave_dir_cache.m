classdef test_parsave_dir_cache < matlab.unittest.TestCase
    % TEST_PARSAVE_DIR_CACHE Unit tests for the parallel-safe save wrapper.
    %
    % parsave_dir_cache.m (in utils/) is a thin wrapper around MATLAB's save()
    % that enables saving variables inside parfor loops. MATLAB's save() cannot
    % be called directly in parfor because it requires variable names to be
    % known at compile time. This wrapper accepts the data as function arguments
    % and saves them with fixed variable names: gtv_mask_warped, D_forward, ref3d.
    %
    % Tests verify:
    %   - The output .mat file is created on disk
    %   - All three variables are faithfully round-tripped (saved then loaded)
    %   - Variable names in the .mat file match the expected identifiers
    %   - A second call to the same path overwrites the previous file
    %   - Arrays of varied spatial dimensions are handled correctly
    %   - Edge case: empty (0x0x0) arrays save and load without error
    %   - Concurrent parfor access writes all files with correct content

    properties
        TempDir
        OriginalPath
    end

    methods(TestMethodSetup)
        function setup(testCase)
            % Create an isolated temp directory for each test and add utils/ to path
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            % Remove the temp directory and restore the original MATLAB path
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testFileIsCreated(testCase)
            % Basic smoke test: verify the .mat file is created on disk after calling
            % parsave_dir_cache with random 5x5x5 arrays (typical small volume size).
            fname = fullfile(testCase.TempDir, 'cache_test.mat');
            parsave_dir_cache(fname, rand(5,5,5) > 0.5, rand(5,5,5,3), rand(5,5,5));

            testCase.verifyTrue(exist(fname, 'file') > 0, ...
                'parsave_dir_cache should create the output .mat file.');
        end

        function testVariablesSavedWithCorrectValues(testCase)
            % Verify data integrity: save known arrays, reload, and compare.
            % Uses a fixed random seed so the test is deterministic.
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
            % Verify idempotent overwrite: a second call with different data
            % should completely replace the previous file contents. This is
            % critical for checkpoint recovery where stale data must not persist.
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
            % Verify that non-cubic arrays (8x10x6) save correctly. Real patient
            % volumes are rarely cubic, so this tests dimensional generality.
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
            % Edge case: empty (0x0x0) arrays must save and reload without error.
            % This can occur when a patient has no valid GTV mask or no
            % deformation field data.
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

        function testConcurrentParforAccess(testCase)
            % Primary use-case validation: invoke parsave_dir_cache from
            % within a parfor loop with multiple iterations writing to
            % different file paths simultaneously. Verify every output file
            % is created and contains the correct content.
            %
            % This guards against race conditions, workspace-scoping issues,
            % and any other problems that arise only under parallel execution.

            nIter = 8;  % enough iterations to exercise concurrent workers

            % Pre-generate deterministic test data for each iteration so we
            % can verify round-trip fidelity after the parfor completes.
            gtv_cells = cell(nIter, 1);
            Dfwd_cells = cell(nIter, 1);
            ref_cells  = cell(nIter, 1);
            fnames     = cell(nIter, 1);

            for k = 1:nIter
                rng(100 + k);
                sz = [3 + k, 4 + k, 5 + k];  % unique dimensions per iteration
                gtv_cells{k}  = rand(sz) > 0.5;
                Dfwd_cells{k} = rand([sz, 3]);
                ref_cells{k}  = rand(sz);
                fnames{k}     = fullfile(testCase.TempDir, sprintf('parfor_iter_%02d.mat', k));
            end

            % --- parfor block (falls back to for-loop when no pool) ---
            parfor k = 1:nIter
                parsave_dir_cache(fnames{k}, gtv_cells{k}, Dfwd_cells{k}, ref_cells{k}); %#ok<PFBNS>
            end

            % Verify every file was created and contains the correct data
            for k = 1:nIter
                testCase.verifyTrue(exist(fnames{k}, 'file') > 0, ...
                    sprintf('parfor iteration %d: output file must exist.', k));

                loaded = load(fnames{k});

                testCase.verifyEqual(loaded.gtv_mask_warped, gtv_cells{k}, ...
                    sprintf('parfor iteration %d: gtv_mask_warped mismatch.', k));
                testCase.verifyEqual(loaded.D_forward, Dfwd_cells{k}, ...
                    sprintf('parfor iteration %d: D_forward mismatch.', k));
                testCase.verifyEqual(loaded.ref3d, ref_cells{k}, ...
                    sprintf('parfor iteration %d: ref3d mismatch.', k));
            end
        end

    end
end