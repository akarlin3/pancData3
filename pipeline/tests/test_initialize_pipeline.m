classdef test_initialize_pipeline < matlab.unittest.TestCase
% TEST_INITIALIZE_PIPELINE — Unit tests for initialize_pipeline.m
%
% Validates pipeline initialization including:
%   - Config path resolution (absolute and relative)
%   - Path setup (core, utils, dependencies added to MATLAB path)
%   - Toolbox license verification
%   - Pre-flight test skip logic (env var, config flag, already-passed)
%   - Test staleness detection

    properties
        TempDir
        OrigPath
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
        end

        function createTempDir(testCase)
            testCase.TempDir = fullfile(tempdir, ['test_init_pipe_' char(java.util.UUID.randomUUID())]);
            mkdir(testCase.TempDir);
            testCase.OrigPath = path();
        end
    end

    methods (TestMethodTeardown)
        function cleanupTempDir(testCase)
            diary off;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
            path(testCase.OrigPath);
        end
    end

    methods (Test)
        function test_resolves_relative_config_path(testCase)
            % When config_path is relative, it should be resolved against pipeline_dir
            pipeline_dir = testCase.TempDir;
            % Create a minimal config.json
            cfg_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(cfg_path, 'w');
            fprintf(fid, '{"skip_tests": true}');
            fclose(fid);

            % Create required subdirs
            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));

            [resolved, ~, ~] = initialize_pipeline( ...
                pipeline_dir, 'config.json', {}, '', false, []);

            testCase.verifyEqual(resolved, cfg_path, ...
                'Relative config path should resolve to full path.');
        end

        function test_absolute_config_path_unchanged(testCase)
            pipeline_dir = testCase.TempDir;
            cfg_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(cfg_path, 'w');
            fprintf(fid, '{"skip_tests": true}');
            fclose(fid);

            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));

            [resolved, ~, ~] = initialize_pipeline( ...
                pipeline_dir, cfg_path, {}, '', false, []);

            testCase.verifyEqual(resolved, cfg_path, ...
                'Absolute config path should be returned unchanged.');
        end

        function test_adds_core_utils_deps_to_path(testCase)
            pipeline_dir = testCase.TempDir;
            cfg_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(cfg_path, 'w');
            fprintf(fid, '{"skip_tests": true}');
            fclose(fid);

            core_dir = fullfile(pipeline_dir, 'core');
            utils_dir = fullfile(pipeline_dir, 'utils');
            deps_dir = fullfile(pipeline_dir, 'dependencies');
            mkdir(core_dir);
            mkdir(utils_dir);
            mkdir(deps_dir);

            initialize_pipeline(pipeline_dir, cfg_path, {}, '', false, []);

            p = path();
            testCase.verifyTrue(contains(p, core_dir), ...
                'core/ should be on the MATLAB path.');
            testCase.verifyTrue(contains(p, utils_dir), ...
                'utils/ should be on the MATLAB path.');
            testCase.verifyTrue(contains(p, deps_dir), ...
                'dependencies/ should be on the MATLAB path.');
        end

        function test_skip_preflight_via_env_var(testCase)
            pipeline_dir = testCase.TempDir;
            cfg_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(cfg_path, 'w');
            fprintf(fid, '{}');
            fclose(fid);

            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));

            % Set env var to skip pre-flight
            setenv('SKIP_PIPELINE_PREFLIGHT', '1');
            testCase.addTeardown(@() setenv('SKIP_PIPELINE_PREFLIGHT', ''));

            % With 'test' in steps but env skip, should not run tests
            [~, tp, ~] = initialize_pipeline( ...
                pipeline_dir, cfg_path, {'test'}, '', false, []);

            % tests_passed should remain false (not run)
            testCase.verifyFalse(tp, ...
                'Pre-flight should be skipped via env var.');
        end

        function test_skip_preflight_via_config_flag(testCase)
            pipeline_dir = testCase.TempDir;
            cfg_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(cfg_path, 'w');
            fprintf(fid, '{"skip_tests": true}');
            fclose(fid);

            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));

            [~, tp, ~] = initialize_pipeline( ...
                pipeline_dir, cfg_path, {'test'}, '', false, []);

            testCase.verifyFalse(tp, ...
                'Pre-flight should be skipped via config skip_tests flag.');
        end

        function test_already_passed_skips_rerun(testCase)
            % When tests_passed_in is true and timestamp is recent, skip tests
            pipeline_dir = testCase.TempDir;
            cfg_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(cfg_path, 'w');
            fprintf(fid, '{}');
            fclose(fid);

            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));
            mkdir(fullfile(pipeline_dir, 'tests'));

            % Pass in true with a very recent timestamp (future)
            [~, tp, ts] = initialize_pipeline( ...
                pipeline_dir, cfg_path, {'test'}, '', true, now + 1);

            testCase.verifyTrue(tp, ...
                'Already-passed tests should not be re-run.');
        end

        function test_no_test_step_skips_preflight(testCase)
            % When 'test' is not in steps_to_run, pre-flight is skipped
            pipeline_dir = testCase.TempDir;
            cfg_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(cfg_path, 'w');
            fprintf(fid, '{}');
            fclose(fid);

            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));

            [~, tp, ~] = initialize_pipeline( ...
                pipeline_dir, cfg_path, {'load', 'sanity'}, '', false, []);

            testCase.verifyFalse(tp, ...
                'Tests should not run when test step is not requested.');
        end
    end
end
