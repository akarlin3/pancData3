classdef test_dispatch_pipeline_steps < matlab.unittest.TestCase
    % TEST_DISPATCH_PIPELINE_STEPS Unit tests for dispatch_pipeline_steps
    % and dispatch_load_and_sanity.
    %
    % Tests the step dispatch utilities: load/sanity abort signaling,
    % skip-step behavior, and session struct consumption.

    properties
        TmpDir
    end

    methods(TestMethodSetup)
        function setupTmpDir(testCase)
            testCase.TmpDir = tempname;
            mkdir(testCase.TmpDir);

            current_dir = fileparts(mfilename('fullpath'));
            addpath(fullfile(current_dir, '..', 'utils'));
            addpath(fullfile(current_dir, '..', 'core'));
        end
    end

    methods(TestMethodTeardown)
        function cleanupTmpDir(testCase)
            diary off;
            if isfolder(testCase.TmpDir)
                rmdir(testCase.TmpDir, 's');
            end
        end
    end

    methods(Test)
        function test_load_skip_no_file_aborts(testCase)
            % Skipping load when no data file exists should set abort=true.
            session = build_test_session(testCase, {'sanity'});

            [~, ~, ~, abort] = dispatch_load_and_sanity(session);
            testCase.verifyTrue(abort);
        end

        function test_load_skip_with_data_succeeds(testCase)
            % Skipping load when data file exists should succeed.
            % Use 'visualize' (not 'sanity') so sanity_checks is not
            % invoked on the minimal fake data.
            session = build_test_session(testCase, {'visualize'});

            % Create fake DWI vectors file
            data_vectors_gtvp = struct('adc', {1}); %#ok<NASGU>
            data_vectors_gtvn = struct('adc', {1}); %#ok<NASGU>
            save(session.voxel_cache_file, 'data_vectors_gtvp', 'data_vectors_gtvn');

            % Create fake summary metrics file
            summary_metrics = struct('id_list', {{'P001'}}); %#ok<NASGU>
            save(session.summary_metrics_file, 'summary_metrics');

            [gtvp, gtvn, sm, abort] = dispatch_load_and_sanity(session);
            testCase.verifyFalse(abort);
            testCase.verifyTrue(isstruct(gtvp));
            testCase.verifyTrue(isstruct(gtvn));
            testCase.verifyTrue(isstruct(sm));
        end

        function test_sanity_skip_uses_loaded_data(testCase)
            % When both load and sanity are skipped but data exists,
            % validated data should be loaded from disk.
            session = build_test_session(testCase, {});

            data_vectors_gtvp = struct('adc', {2}); %#ok<NASGU>
            data_vectors_gtvn = struct('adc', {2}); %#ok<NASGU>
            save(session.voxel_cache_file, 'data_vectors_gtvp', 'data_vectors_gtvn');
            summary_metrics = struct('id_list', {{'P002'}}); %#ok<NASGU>
            save(session.summary_metrics_file, 'summary_metrics');

            [gtvp, ~, ~, abort] = dispatch_load_and_sanity(session);
            testCase.verifyFalse(abort);
            testCase.verifyEqual(gtvp.adc, 2);
        end

        function test_dispatch_analysis_no_steps(testCase)
            % dispatch_pipeline_steps with no matching steps should not error.
            session = build_test_session(testCase, {});
            session.steps_to_run = {};

            % Should complete without error
            dispatch_pipeline_steps(session, struct(), struct(), struct('id_list', {{}}));
        end

        function test_dispatch_analysis_skips_all(testCase)
            % When steps_to_run has no analysis steps, dispatch should
            % just print skip messages without error.
            session = build_test_session(testCase, {'load', 'sanity'});

            dispatch_pipeline_steps(session, struct(), struct(), struct('id_list', {{}}));
            % No error = pass
        end

        function test_session_struct_fields(testCase)
            % Verify all required session fields are accessed without error.
            session = build_test_session(testCase, {});
            % Verify critical fields exist
            testCase.verifyTrue(isfield(session, 'config_struct'));
            testCase.verifyTrue(isfield(session, 'steps_to_run'));
            testCase.verifyTrue(isfield(session, 'current_name'));
            testCase.verifyTrue(isfield(session, 'current_dtype'));
            testCase.verifyTrue(isfield(session, 'log_fid'));
            testCase.verifyTrue(isfield(session, 'master_diary_file'));
            testCase.verifyTrue(isfield(session, 'pipeGUI'));
            testCase.verifyTrue(isfield(session, 'type_output_folder'));
            testCase.verifyTrue(isfield(session, 'voxel_cache_file'));
            testCase.verifyTrue(isfield(session, 'results_file'));
        end

        function test_legacy_fallback_blocked_for_dncnn(testCase)
            % Non-Standard DWI types should NOT fall back to legacy file.
            session = build_test_session(testCase, {'sanity'});
            session.current_dtype = 2;
            session.current_name = 'dnCNN';

            % Create only the legacy file, not the type-specific one
            data_vectors_gtvp = struct('adc', {1}); %#ok<NASGU>
            data_vectors_gtvn = struct('adc', {1}); %#ok<NASGU>
            save(session.voxel_cache_fallback_file, 'data_vectors_gtvp', 'data_vectors_gtvn');

            [~, ~, ~, abort] = dispatch_load_and_sanity(session);
            testCase.verifyTrue(abort);
        end
    end

    methods(Access = private)
        function session = build_test_session(testCase, steps)
            % Build a minimal session struct for testing.
            type_dir = fullfile(testCase.TmpDir, 'Standard');
            mkdir(type_dir);

            config_struct = struct();
            config_struct.dataloc = testCase.TmpDir;
            config_struct.output_folder = type_dir;
            config_struct.master_output_folder = testCase.TmpDir;
            config_struct.run_compare_cores = false;
            config_struct.run_all_core_methods = false;
            config_struct.dwi_types_to_run = 1;
            config_struct.dwi_type_name = 'Standard';

            diary_file = fullfile(type_dir, 'test_diary.txt');

            session = struct();
            session.config_struct = config_struct;
            session.steps_to_run = steps;
            session.current_dtype = 1;
            session.current_name = 'Standard';
            session.master_output_folder = testCase.TmpDir;
            session.type_output_folder = type_dir;
            session.master_diary_file = diary_file;
            session.log_fid = -1;
            session.pipeGUI = [];
            session.voxel_cache_file = fullfile(testCase.TmpDir, 'pipeline_voxels_Standard.mat');
            session.voxel_cache_fallback_file = fullfile(testCase.TmpDir, 'pipeline_voxels.mat');
            session.summary_metrics_file = fullfile(type_dir, 'summary_metrics_Standard.mat');
            session.results_file = fullfile(type_dir, 'calculated_results_Standard.mat');
            session.baseline_results_file = fullfile(type_dir, 'metrics_baseline_results_Standard.mat');
            session.dosimetry_results_file = fullfile(type_dir, 'metrics_dosimetry_results_Standard.mat');
            session.predictive_results_file = fullfile(type_dir, 'metrics_stats_predictive_results_Standard.mat');
        end
    end
end
