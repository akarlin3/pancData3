classdef test_initialize_pipeline < matlab.unittest.TestCase
% TEST_INITIALIZE_PIPELINE — Unit tests for initialize_pipeline.m
%
% Validates pipeline initialization including:
%   - Config path resolution (absolute and relative)
%   - Path setup (core, utils, dependencies added)
%   - Toolbox license verification
%   - Pre-flight test skip via config flag
%   - Pre-flight test skip via environment variable
%   - Test staleness detection

    properties
        TempDir
        OrigPath
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            addpath(fullfile(testDir, '..', 'utils'));
        end

        function setupTempDir(testCase)
            testCase.TempDir = fullfile(tempdir, ['test_init_' char(datetime('now','Format','yyyyMMddHHmmss'))]);
            mkdir(testCase.TempDir);
            testCase.OrigPath = path();
        end
    end

    methods (TestMethodTeardown)
        function cleanupTempDir(testCase)
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods (Test)
        function test_config_path_resolution_absolute(testCase)
            % Absolute path to config should be returned as-is if file exists
            config_path = fullfile(testCase.TempDir, 'config.json');
            fid = fopen(config_path, 'w');
            fprintf(fid, '{"skip_tests": true}\n');
            fclose(fid);

            pipeline_dir = testCase.TempDir;
            % Create subdirs to avoid errors
            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));

            [resolved, ~, ~] = initialize_pipeline( ...
                pipeline_dir, config_path, {}, '', false, []);

            testCase.verifyEqual(resolved, config_path);
        end

        function test_config_path_resolution_relative(testCase)
            % Relative path should resolve against pipeline_dir
            pipeline_dir = testCase.TempDir;
            config_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(config_path, 'w');
            fprintf(fid, '{"skip_tests": true}\n');
            fclose(fid);

            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));

            [resolved, ~, ~] = initialize_pipeline( ...
                pipeline_dir, 'config.json', {}, '', false, []);

            testCase.verifyTrue(isfile(resolved), ...
                'Resolved config path should point to an existing file.');
        end

        function test_skip_tests_env_variable(testCase)
            % SKIP_PIPELINE_PREFLIGHT=1 should skip pre-flight tests
            pipeline_dir = testCase.TempDir;
            config_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(config_path, 'w');
            fprintf(fid, '{"skip_tests": false}\n');
            fclose(fid);
            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));

            setenv('SKIP_PIPELINE_PREFLIGHT', '1');
            testCase.addTeardown(@() setenv('SKIP_PIPELINE_PREFLIGHT', ''));

            % With 'test' in steps but env skip set, should not run tests
            [~, passed, ~] = initialize_pipeline( ...
                pipeline_dir, config_path, {'test'}, '', false, []);

            % Should return false (tests not run, not verified)
            testCase.verifyFalse(passed);
        end

        function test_skip_tests_config_flag(testCase)
            % skip_tests: true in config should skip pre-flight tests
            pipeline_dir = testCase.TempDir;
            config_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(config_path, 'w');
            fprintf(fid, '{"skip_tests": true}\n');
            fclose(fid);
            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));

            setenv('SKIP_PIPELINE_PREFLIGHT', '');

            [~, passed, ~] = initialize_pipeline( ...
                pipeline_dir, config_path, {'test'}, '', false, []);

            testCase.verifyFalse(passed);
        end

        function test_no_test_step_skips_tests(testCase)
            % When 'test' is not in steps_to_run, pre-flight should be skipped
            pipeline_dir = testCase.TempDir;
            config_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(config_path, 'w');
            fprintf(fid, '{"skip_tests": false}\n');
            fclose(fid);
            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));

            setenv('SKIP_PIPELINE_PREFLIGHT', '');

            [~, passed, ~] = initialize_pipeline( ...
                pipeline_dir, config_path, {'load'}, '', false, []);

            % Tests not run, flag should remain false
            testCase.verifyFalse(passed);
        end

        function test_already_passed_tests_cached(testCase)
            % When tests_passed is already true, should skip re-running
            pipeline_dir = testCase.TempDir;
            config_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(config_path, 'w');
            fprintf(fid, '{"skip_tests": false}\n');
            fclose(fid);
            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));
            mkdir(fullfile(pipeline_dir, 'tests'));

            setenv('SKIP_PIPELINE_PREFLIGHT', '');

            % Pass in true + future timestamp so staleness check passes
            [~, passed, ts] = initialize_pipeline( ...
                pipeline_dir, config_path, {'test'}, '', true, now + 1);

            testCase.verifyTrue(passed, ...
                'Should preserve cached test pass.');
        end

        function test_toolbox_check_no_error(testCase)
            % On a properly licensed MATLAB installation, toolbox checks should pass
            if exist('OCTAVE_VERSION', 'builtin'), return; end

            pipeline_dir = testCase.TempDir;
            config_path = fullfile(pipeline_dir, 'config.json');
            fid = fopen(config_path, 'w');
            fprintf(fid, '{"skip_tests": true}\n');
            fclose(fid);
            mkdir(fullfile(pipeline_dir, 'core'));
            mkdir(fullfile(pipeline_dir, 'utils'));
            mkdir(fullfile(pipeline_dir, 'dependencies'));

            % Should not throw on a licensed installation
            [~, ~, ~] = initialize_pipeline( ...
                pipeline_dir, config_path, {}, '', false, []);
        end
    end
end
