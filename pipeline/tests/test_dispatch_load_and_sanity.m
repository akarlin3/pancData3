classdef test_dispatch_load_and_sanity < matlab.unittest.TestCase
% TEST_DISPATCH_LOAD_AND_SANITY — Unit tests for dispatch_load_and_sanity.m
%
% Validates the load/sanity dispatch logic:
%   - Skip-both path loads data from disk correctly
%   - Missing .mat file causes abort=true
%   - Skipping sanity passes data through unchanged
%   - Auxiliary biomarker loading is skipped when not configured
%   - abort=false on successful skip paths
%   - All 4 outputs are returned

    properties
        OriginalPath
        TempDir
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path;
            utilsDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils');
            coreDir  = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'core');
            addpath(utilsDir);
            addpath(coreDir);
        end

        function createTempDir(testCase)
            testCase.TempDir = fullfile(tempdir, ['test_dispatch_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempDir);
        end
    end

    methods (TestMethodTeardown)
        function teardownEnv(testCase)
            diary off;
            path(testCase.OriginalPath);
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods (Test)
        function testSkipBothStepsWithData(testCase)
            % Skip both load and sanity steps, provide .mat files on disk.
            % Verify data loads correctly and abort=false.
            session = testCase.makeSession({});

            % Create voxel-cache .mat
            data_vectors_gtvp = struct('adc', rand(5,1), 'D', rand(5,1));
            data_vectors_gtvn = struct('adc', rand(3,1), 'D', rand(3,1));
            save(session.voxel_cache_file, 'data_vectors_gtvp', 'data_vectors_gtvn');

            % Create summary_metrics .mat
            summary_metrics = struct('mean_adc', 0.001, 'n_patients', 10);
            save(session.summary_metrics_file, 'summary_metrics');

            evalc('[gtvp, gtvn, sm, ab] = dispatch_load_and_sanity(session)');

            testCase.verifyFalse(ab, 'abort should be false when data loads successfully.');
            testCase.verifyEqual(gtvp.adc, data_vectors_gtvp.adc, ...
                'GTVp data should match saved data.');
            testCase.verifyEqual(gtvn.adc, data_vectors_gtvn.adc, ...
                'GTVn data should match saved data.');
            testCase.verifyEqual(sm.mean_adc, 0.001, ...
                'Summary metrics should match saved data.');
        end

        function testSkipLoadMissingFile(testCase)
            % Skip load step but point to non-existent .mat file.
            % Verify abort=true.
            session = testCase.makeSession({});
            session.voxel_cache_file = fullfile(testCase.TempDir, 'nonexistent.mat');
            session.voxel_cache_fallback_file = fullfile(testCase.TempDir, 'also_nonexistent.mat');

            evalc('[~, ~, ~, ab] = dispatch_load_and_sanity(session)');

            testCase.verifyTrue(ab, 'abort should be true when data files are missing.');
        end

        function testSkipSanityPassesDataThrough(testCase)
            % Skip sanity step, verify data from load step passes through unchanged.
            % Here we skip both load and sanity; data loaded from disk should
            % appear in validated_data outputs unchanged.
            session = testCase.makeSession({});

            data_vectors_gtvp = struct('adc', [1.1; 2.2; 3.3]);
            data_vectors_gtvn = struct('adc', [4.4; 5.5]);
            save(session.voxel_cache_file, 'data_vectors_gtvp', 'data_vectors_gtvn');

            summary_metrics = struct('ok', true);
            save(session.summary_metrics_file, 'summary_metrics');

            evalc('[gtvp, gtvn, ~, ab] = dispatch_load_and_sanity(session)');

            testCase.verifyFalse(ab);
            testCase.verifyEqual(gtvp.adc, [1.1; 2.2; 3.3], ...
                'Validated GTVp should be identical to loaded data when sanity is skipped.');
            testCase.verifyEqual(gtvn.adc, [4.4; 5.5], ...
                'Validated GTVn should be identical to loaded data when sanity is skipped.');
        end

        function testAuxiliaryBiomarkerSkippedWhenNotConfigured(testCase)
            % Verify no error when auxiliary_biomarker_csv is not set.
            session = testCase.makeSession({});
            session.config_struct.auxiliary_biomarker_csv = '';

            data_vectors_gtvp = struct('adc', rand(3,1));
            data_vectors_gtvn = struct('adc', rand(3,1));
            save(session.voxel_cache_file, 'data_vectors_gtvp', 'data_vectors_gtvn');

            summary_metrics = struct('mean_adc', 0.002);
            save(session.summary_metrics_file, 'summary_metrics');

            evalc('[~, ~, sm, ab] = dispatch_load_and_sanity(session)');

            testCase.verifyFalse(ab, 'abort should be false.');
            testCase.verifyFalse(isfield(sm, 'auxiliary_biomarkers'), ...
                'auxiliary_biomarkers should not be added when CSV is empty.');
        end

        function testAbortReturnsFalseOnSuccess(testCase)
            % Verify abort=false when all skip steps succeed with valid data.
            session = testCase.makeSession({});

            data_vectors_gtvp = struct('adc', rand(8,1));
            data_vectors_gtvn = struct('adc', rand(4,1));
            save(session.voxel_cache_file, 'data_vectors_gtvp', 'data_vectors_gtvn');

            summary_metrics = struct('count', 8);
            save(session.summary_metrics_file, 'summary_metrics');

            evalc('[~, ~, ~, ab] = dispatch_load_and_sanity(session)');

            testCase.verifyFalse(ab, 'abort should be false on successful data load.');
        end

        function testOutputDimensions(testCase)
            % Verify all 4 outputs are returned.
            session = testCase.makeSession({});

            data_vectors_gtvp = struct('adc', rand(5,1));
            data_vectors_gtvn = struct('adc', rand(5,1));
            save(session.voxel_cache_file, 'data_vectors_gtvp', 'data_vectors_gtvn');

            summary_metrics = struct('val', 42);
            save(session.summary_metrics_file, 'summary_metrics');

            evalc('[gtvp, gtvn, sm, ab] = dispatch_load_and_sanity(session)');

            testCase.verifyNotEmpty(gtvp, 'validated_data_gtvp should not be empty.');
            testCase.verifyNotEmpty(gtvn, 'validated_data_gtvn should not be empty.');
            testCase.verifyNotEmpty(sm, 'summary_metrics should not be empty.');
            testCase.verifyClass(ab, 'logical', 'abort should be a logical.');
        end
    end

    methods (Access = private)
        function session = makeSession(testCase, steps)
            % Build a minimal session struct for testing skip paths.
            %   steps - cell array of step names to include in steps_to_run
            session = struct();
            session.steps_to_run = steps;
            session.config_struct = struct( ...
                'output_folder', testCase.TempDir, ...
                'dwi_type', 'Standard', ...
                'auxiliary_biomarker_csv', '' ...
            );
            session.pipeGUI = [];
            session.log_fid = -1;
            session.master_diary_file = fullfile(testCase.TempDir, 'test_diary.log');
            session.current_name = 'Standard';
            session.current_dtype = 1;
            session.voxel_cache_file = fullfile(testCase.TempDir, 'pipeline_voxels_Standard.mat');
            session.voxel_cache_fallback_file = fullfile(testCase.TempDir, 'pipeline_voxels.mat');
            session.summary_metrics_file = fullfile(testCase.TempDir, 'summary_metrics_Standard.mat');
        end
    end
end
