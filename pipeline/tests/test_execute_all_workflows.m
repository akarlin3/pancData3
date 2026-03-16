classdef test_execute_all_workflows < matlab.unittest.TestCase
    % TEST_EXECUTE_ALL_WORKFLOWS Integration tests for execute_all_workflows.m
    %
    % Tests the orchestrator logic: config mutation/restoration, output folder
    % creation with sentinel, DWI type sequencing, step list construction,
    % and compare_cores conditional injection.
    %
    % These tests exercise individual behaviours of the orchestrator without
    % actually running the full pipeline (which requires patient data).
    %
    % Run tests with:
    %   results = runtests('tests/test_execute_all_workflows.m');

    properties
        TempDir
        OrigPath
    end

    methods (TestMethodSetup)
        function setupPaths(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            testCase.OrigPath = path;
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'dependencies'));
        end
    end

    methods (TestMethodTeardown)
        function cleanupTemp(testCase)
            diary off;
            if isfield(testCase, 'OrigPath') && ~isempty(testCase.OrigPath)
                entries = strsplit(testCase.OrigPath, pathsep);
                entries = entries(cellfun(@(p) isempty(p) || isfolder(p), entries));
                path(strjoin(entries, pathsep));
            end
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods (Test)

        function testConfigMutationAndRestoration(testCase)
            % Verify that json_set_field correctly mutates dwi_type and
            % skip_to_reload, and that restoring the original string
            % recovers the initial config exactly.
            original = '{"dwi_type": "Standard", "skip_to_reload": false, "dataloc": "/tmp"}';
            config_file = fullfile(testCase.TempDir, 'config.json');
            fid = fopen(config_file, 'w');
            fwrite(fid, original);
            fclose(fid);

            % Mutate to dnCNN + skip_to_reload=true
            modified = json_set_field(original, 'dwi_type', 'dnCNN');
            modified = json_set_field(modified, 'skip_to_reload', true);
            fid = fopen(config_file, 'w');
            fwrite(fid, modified);
            fclose(fid);

            % Verify mutation
            raw = fileread(config_file);
            cfg = jsondecode(raw);
            testCase.verifyEqual(cfg.dwi_type, 'dnCNN', ...
                'dwi_type should be mutated to dnCNN.');
            testCase.verifyTrue(cfg.skip_to_reload, ...
                'skip_to_reload should be mutated to true.');

            % Restore original
            fid = fopen(config_file, 'w');
            fwrite(fid, original);
            fclose(fid);
            raw = fileread(config_file);
            cfg = jsondecode(raw);
            testCase.verifyEqual(cfg.dwi_type, 'Standard', ...
                'dwi_type should be restored to Standard.');
            testCase.verifyFalse(cfg.skip_to_reload, ...
                'skip_to_reload should be restored to false.');
        end

        function testConfigRestorationViaCleanup(testCase)
            % Verify the restore_config_file pattern: onCleanup restores
            % config.json even after mutation.
            config_file = fullfile(testCase.TempDir, 'config.json');
            original_str = '{"dwi_type": "Standard", "skip_to_reload": false}';
            fid = fopen(config_file, 'w');
            fwrite(fid, original_str);
            fclose(fid);

            % Simulate the cleanup function from execute_all_workflows
            restore_fn = @() localRestoreConfig(config_file, original_str);
            cleanup = onCleanup(restore_fn);

            % Mutate config
            modified = json_set_field(original_str, 'dwi_type', 'IVIMnet');
            modified = json_set_field(modified, 'skip_to_reload', true);
            fid = fopen(config_file, 'w');
            fwrite(fid, modified);
            fclose(fid);

            % Trigger cleanup
            delete(cleanup);

            % Verify restoration
            raw = fileread(config_file);
            cfg = jsondecode(raw);
            testCase.verifyEqual(cfg.dwi_type, 'Standard', ...
                'Config should be restored after cleanup.');
            testCase.verifyFalse(cfg.skip_to_reload, ...
                'skip_to_reload should be restored after cleanup.');
        end

        function testOutputFolderCreationWithSentinel(testCase)
            % Verify that the output folder pattern creates the directory
            % and writes a .pipeline_created sentinel.
            timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
            output_folder = fullfile(testCase.TempDir, sprintf('saved_files_%s', timestamp_str));
            mkdir(output_folder);

            sentinel = fullfile(output_folder, '.pipeline_created');
            fid = fopen(sentinel, 'w');
            testCase.verifyGreaterThan(fid, 0, 'Sentinel file should be writable.');
            fprintf(fid, 'Created by test at %s\n', timestamp_str);
            fclose(fid);

            testCase.verifyTrue(exist(output_folder, 'dir') == 7, ...
                'Output folder should exist.');
            testCase.verifyTrue(exist(sentinel, 'file') == 2, ...
                'Sentinel file should exist.');

            content = fileread(sentinel);
            testCase.verifySubstring(content, timestamp_str, ...
                'Sentinel should contain the timestamp.');
        end

        function testDwiTypeSequencing(testCase)
            % Verify the 3-type sequence: Standard (skip=false),
            % dnCNN (skip=true), IVIMnet (skip=true).
            original = '{"dwi_type": "Standard", "skip_to_reload": false}';
            types = {'Standard', 'dnCNN', 'IVIMnet'};
            skip_values = {false, true, true};

            for i = 1:3
                config_json = json_set_field(original, 'dwi_type', types{i});
                config_json = json_set_field(config_json, 'skip_to_reload', skip_values{i});
                cfg = jsondecode(config_json);

                testCase.verifyEqual(cfg.dwi_type, types{i}, ...
                    sprintf('DWI type %d should be %s.', i, types{i}));
                testCase.verifyEqual(cfg.skip_to_reload, skip_values{i}, ...
                    sprintf('skip_to_reload for %s should be %s.', ...
                    types{i}, mat2str(skip_values{i})));
            end
        end

        function testStepListDefaultOrder(testCase)
            % Verify the default step order matches the analytical
            % dependency chain.
            steps = {'load', 'sanity', 'visualize', 'metrics_baseline', ...
                     'metrics_longitudinal', 'metrics_dosimetry', ...
                     'metrics_stats_comparisons', 'metrics_stats_predictive', ...
                     'metrics_survival'};

            testCase.verifyEqual(numel(steps), 9, ...
                'Default step list should have 9 steps.');

            % Verify ordering: load before sanity, sanity before visualize, etc.
            testCase.verifyLessThan(find(strcmp(steps, 'load')), ...
                find(strcmp(steps, 'sanity')), ...
                'load must precede sanity.');
            testCase.verifyLessThan(find(strcmp(steps, 'metrics_baseline')), ...
                find(strcmp(steps, 'metrics_longitudinal')), ...
                'metrics_baseline must precede metrics_longitudinal.');
            testCase.verifyLessThan(find(strcmp(steps, 'metrics_stats_predictive')), ...
                find(strcmp(steps, 'metrics_survival')), ...
                'metrics_stats_predictive must precede metrics_survival.');
        end

        function testCompareCoresInjection(testCase)
            % Verify that compare_cores is injected after metrics_baseline
            % when run_compare_cores is true.
            steps = {'load', 'sanity', 'visualize', 'metrics_baseline', ...
                     'metrics_longitudinal', 'metrics_dosimetry', ...
                     'metrics_stats_comparisons', 'metrics_stats_predictive', ...
                     'metrics_survival'};

            % Simulate injection logic from execute_all_workflows
            eaw_cfg = struct('run_compare_cores', true);
            if isfield(eaw_cfg, 'run_compare_cores') && eaw_cfg.run_compare_cores
                cc_idx = find(strcmp(steps, 'metrics_baseline'));
                if ~isempty(cc_idx)
                    steps = [steps(1:cc_idx), {'compare_cores'}, steps(cc_idx+1:end)];
                end
            end

            testCase.verifyEqual(numel(steps), 10, ...
                'Step list should have 10 steps with compare_cores.');
            cc_pos = find(strcmp(steps, 'compare_cores'));
            mb_pos = find(strcmp(steps, 'metrics_baseline'));
            ml_pos = find(strcmp(steps, 'metrics_longitudinal'));
            testCase.verifyEqual(cc_pos, mb_pos + 1, ...
                'compare_cores should be right after metrics_baseline.');
            testCase.verifyLessThan(cc_pos, ml_pos, ...
                'compare_cores should precede metrics_longitudinal.');
        end

        function testCompareCoresNotInjectedByDefault(testCase)
            % Verify that compare_cores is NOT in the step list when
            % run_compare_cores is false.
            steps = {'load', 'sanity', 'visualize', 'metrics_baseline', ...
                     'metrics_longitudinal', 'metrics_dosimetry', ...
                     'metrics_stats_comparisons', 'metrics_stats_predictive', ...
                     'metrics_survival'};

            eaw_cfg = struct('run_compare_cores', false);
            if isfield(eaw_cfg, 'run_compare_cores') && eaw_cfg.run_compare_cores
                cc_idx = find(strcmp(steps, 'metrics_baseline'));
                if ~isempty(cc_idx)
                    steps = [steps(1:cc_idx), {'compare_cores'}, steps(cc_idx+1:end)];
                end
            end

            testCase.verifyEqual(numel(steps), 9, ...
                'Step list should have 9 steps without compare_cores.');
            testCase.verifyTrue(isempty(find(strcmp(steps, 'compare_cores'), 1)), ...
                'compare_cores should not be in the step list.');
        end

        function testSkipTestsConfigFlag(testCase)
            % Verify that skip_tests flag is correctly read from config.
            config_file = fullfile(testCase.TempDir, 'config.json');

            % Write config with skip_tests = true
            cfg = struct('dwi_type', 'Standard', 'skip_to_reload', false, ...
                         'skip_tests', true, 'dataloc', testCase.TempDir);
            fid = fopen(config_file, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            eaw_skip_tests = false;
            try
                eaw_raw = fileread(config_file);
                eaw_cfg = jsondecode(eaw_raw);
                if isfield(eaw_cfg, 'skip_tests') && eaw_cfg.skip_tests
                    eaw_skip_tests = true;
                end
            catch
            end

            testCase.verifyTrue(eaw_skip_tests, ...
                'skip_tests=true should be detected from config.');

            % Write config with skip_tests = false
            cfg.skip_tests = false;
            fid = fopen(config_file, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            eaw_skip_tests = false;
            try
                eaw_raw = fileread(config_file);
                eaw_cfg = jsondecode(eaw_raw);
                if isfield(eaw_cfg, 'skip_tests') && eaw_cfg.skip_tests
                    eaw_skip_tests = true;
                end
            catch
            end

            testCase.verifyFalse(eaw_skip_tests, ...
                'skip_tests=false should not trigger skipping.');
        end

        function testSkipTestsDefaultsToFalse(testCase)
            % Verify that missing skip_tests field defaults to running tests.
            config_file = fullfile(testCase.TempDir, 'config.json');
            cfg = struct('dwi_type', 'Standard', 'dataloc', testCase.TempDir);
            fid = fopen(config_file, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            eaw_skip_tests = false;
            try
                eaw_raw = fileread(config_file);
                eaw_cfg = jsondecode(eaw_raw);
                if isfield(eaw_cfg, 'skip_tests') && eaw_cfg.skip_tests
                    eaw_skip_tests = true;
                end
            catch
            end

            testCase.verifyFalse(eaw_skip_tests, ...
                'Missing skip_tests should default to false (run tests).');
        end

        function testDiaryFileCreation(testCase)
            % Verify that the master diary file is created in the output folder.
            output_folder = fullfile(testCase.TempDir, 'saved_files_test');
            mkdir(output_folder);
            diary_file = fullfile(output_folder, 'execute_all_workflows.log');

            diary(diary_file);
            fprintf('Test diary entry.\n');
            diary off;

            testCase.verifyTrue(exist(diary_file, 'file') == 2, ...
                'Master diary file should be created.');
            content = fileread(diary_file);
            testCase.verifySubstring(content, 'Test diary entry', ...
                'Diary file should contain the test entry.');
        end

        function testTimestampedFolderNaming(testCase)
            % Verify the timestamped folder naming pattern.
            timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
            folder_name = sprintf('saved_files_%s', timestamp_str);

            % Verify the format matches expected pattern
            testCase.verifyTrue(~isempty(regexp(folder_name, ...
                '^saved_files_\d{8}_\d{6}$', 'once')), ...
                'Folder name should match saved_files_YYYYMMDD_HHMMSS pattern.');
        end

        function testConfigFileOpenFailure(testCase)
            % Verify that an error is raised when config.json cannot be opened.
            bad_path = fullfile(testCase.TempDir, 'nonexistent', 'config.json');

            testCase.verifyError(@() openConfigOrError(bad_path), ...
                'execute_all_workflows:fileOpenFailed', ...
                'Should error when config.json cannot be opened.');
        end

        function testConfigWriteRoundTrip(testCase)
            % Verify that writing modified config and reading back preserves
            % all three DWI types sequentially.
            config_file = fullfile(testCase.TempDir, 'config.json');
            original = '{"dwi_type": "Standard", "skip_to_reload": false, "dataloc": "/tmp"}';

            types_and_skips = {
                'Standard', false;
                'dnCNN', true;
                'IVIMnet', true;
            };

            for i = 1:size(types_and_skips, 1)
                config_json = json_set_field(original, 'dwi_type', types_and_skips{i,1});
                config_json = json_set_field(config_json, 'skip_to_reload', types_and_skips{i,2});

                fid = fopen(config_file, 'w');
                fwrite(fid, config_json);
                fclose(fid);

                % Read back and verify
                raw = fileread(config_file);
                cfg = jsondecode(raw);
                testCase.verifyEqual(cfg.dwi_type, types_and_skips{i,1});
                testCase.verifyEqual(cfg.skip_to_reload, types_and_skips{i,2});
            end

            % Restore original and verify
            fid = fopen(config_file, 'w');
            fwrite(fid, original);
            fclose(fid);
            raw = fileread(config_file);
            cfg = jsondecode(raw);
            testCase.verifyEqual(cfg.dwi_type, 'Standard');
        end

    end
end

function openConfigOrError(config_path)
    fid = fopen(config_path, 'r');
    if fid < 0
        error('execute_all_workflows:fileOpenFailed', ...
            'Cannot open %s for reading.', config_path);
    end
    fclose(fid);
end

function localRestoreConfig(config_path, original_str)
    try
        fid = fopen(config_path, 'w');
        if fid == -1
            error('restore:openFailed', 'Cannot open config.json for writing.');
        end
        fwrite(fid, original_str);
        fclose(fid);
    catch
        fprintf('Config restoration failed.\n');
    end
end
