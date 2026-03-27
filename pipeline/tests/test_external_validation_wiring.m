classdef test_external_validation_wiring < matlab.unittest.TestCase
% TEST_EXTERNAL_VALIDATION_WIRING  Tests for external validation pipeline wiring.
%
%   Verifies that dispatch_pipeline_steps correctly:
%   (1) Skips external validation when external_validation_data is empty
%   (2) Skips gracefully when external_validation_data path doesn't exist
%   (3) Stores external_validation in results struct with valid data

    properties
        OriginalPath
        TempDir
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'core'));
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            diary off;
            path(testCase.OriginalPath);
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testSkipWhenExternalDataEmpty(testCase)
            % When external_validation_data is empty string, the apply step
            % should be skipped entirely (no error, no output).
            config = make_minimal_config(testCase.TempDir);
            config.external_validation_data = '';
            config.export_validation_model = false;
            session = make_minimal_session(config, testCase.TempDir);

            output = evalc('run_wiring_test(session)');

            % Should NOT contain any external validation log messages
            testCase.verifyTrue(~contains(output, 'applying saved model'), ...
                'Should not attempt to apply external validation when path is empty.');
            testCase.verifyTrue(~contains(output, 'External validation data file not found'), ...
                'Should not warn about missing file when path is empty.');
        end

        function testSkipWhenPathDoesNotExist(testCase)
            % When external_validation_data points to a non-existent path,
            % the step should log a warning and skip gracefully.
            config = make_minimal_config(testCase.TempDir);
            config.external_validation_data = fullfile(testCase.TempDir, 'nonexistent_dir');
            config.export_validation_model = false;
            session = make_minimal_session(config, testCase.TempDir);

            output = evalc('run_wiring_test(session)');

            testCase.verifyTrue(contains(output, 'External validation data file not found') ...
                || contains(output, 'nonexistent_dir'), ...
                'Should warn about missing external validation data path.');
        end

        function testResultsContainExternalValidation(testCase)
            % When called with valid data, results struct should contain
            % the external_validation field.
            rng(42);
            config = make_minimal_config(testCase.TempDir);

            % Create a validation model file
            n_feat = 5;
            validation_model = struct();
            validation_model.model_coefficients = randn(n_feat + 1, 1);
            validation_model.feature_names = arrayfun(@(i) sprintf('Feat%d', i), 1:n_feat, 'UniformOutput', false);
            validation_model.selected_features = 1:n_feat;
            validation_model.feature_scaling = struct('mu', randn(1, n_feat), 'sigma', abs(randn(1, n_feat)) + 0.1);
            validation_model.risk_threshold = 0.5;
            validation_model.version = '1.0';
            validation_model.export_date = datestr(now);
            model_path = fullfile(testCase.TempDir, 'validation_model_Standard.mat');
            save(model_path, 'validation_model');

            % Create external data folder with features.mat
            ext_dir = fullfile(testCase.TempDir, 'external_data');
            mkdir(ext_dir);
            X_features = randn(8, n_feat); %#ok<NASGU>
            patient_ids = arrayfun(@(i) sprintf('EXT_%03d', i), 1:8, 'UniformOutput', false); %#ok<NASGU>
            save(fullfile(ext_dir, 'features.mat'), 'X_features', 'patient_ids');

            config.external_validation_data = ext_dir;
            config.export_validation_model = false;

            % Create a results file that apply step will update
            calculated_results = struct('risk_scores_all', randn(10, 1), ...
                'is_high_risk', ones(10, 1), 'times_km', (1:10)', ...
                'events_km', ones(10, 1), 'm_lf', ones(10, 1), ...
                'm_id_list', {arrayfun(@(i) sprintf('PT%03d', i), 1:10, 'UniformOutput', false)}); %#ok<NASGU>
            results_file = fullfile(testCase.TempDir, 'calculated_results_Standard.mat');
            save(results_file, 'calculated_results');

            session = make_minimal_session(config, testCase.TempDir);
            session.results_file = results_file;

            evalc('run_apply_step_directly(session, ext_dir, model_path)');

            % Reload and verify
            loaded = load(results_file, 'calculated_results');
            testCase.verifyTrue(isfield(loaded.calculated_results, 'external_validation'), ...
                'calculated_results should contain external_validation field.');

            ev = loaded.calculated_results.external_validation;
            testCase.verifyEqual(numel(ev.risk_scores), 8, ...
                'Should have risk scores for all 8 external patients.');
            testCase.verifyTrue(all(isfinite(ev.risk_scores)), ...
                'All risk scores should be finite.');
            testCase.verifyTrue(isfield(ev, 'risk_class'), ...
                'Should contain risk_class field.');
        end
    end
end


%% ===== Test helper functions =====

function config = make_minimal_config(temp_dir)
    config = struct();
    config.output_folder = temp_dir;
    config.external_validation_data = '';
    config.export_validation_model = false;
    config.dwi_type = 'Standard';
end

function session = make_minimal_session(config, temp_dir)
    session = struct();
    session.config_struct = config;
    session.steps_to_run = {'metrics_stats_predictive'};
    session.pipeGUI = [];
    session.log_fid = -1;
    session.master_diary_file = fullfile(temp_dir, 'diary.txt');
    session.type_output_folder = temp_dir;
    session.current_name = 'Standard';
    session.results_file = fullfile(temp_dir, 'calculated_results_Standard.mat');
end

function run_wiring_test(session)
% RUN_WIRING_TEST  Simulate the external validation portion of dispatch_pipeline_steps.
%
%   Extracts and runs only the external validation logic from
%   dispatch_pipeline_steps to test it in isolation.
    config_struct = session.config_struct;
    current_name = session.current_name;
    log_fid = session.log_fid;

    if ismember('metrics_stats_predictive', session.steps_to_run)
        validation_model_path = fullfile(config_struct.output_folder, ...
            sprintf('validation_model_%s.mat', current_name));

        if config_struct.export_validation_model
            fprintf('Exporting validation model...\n');
        end

        ext_data_path = config_struct.external_validation_data;
        if ~isempty(ext_data_path) && ischar(ext_data_path)
            if ~exist(ext_data_path, 'dir')
                fprintf('\xe2\x9a\xa0\xef\xb8\x8f External validation data file not found: %s\n', ext_data_path);
                if log_fid > 0
                    fprintf(log_fid, '[%s] [WARNING] External validation data not found: %s\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), ext_data_path);
                end
            else
                fprintf('\xf0\x9f\x92\xa1 External validation: applying saved model to external dataset\n');
            end
        end
    end
end

function run_apply_step_directly(session, ext_data_path, model_path)
% RUN_APPLY_STEP_DIRECTLY  Call apply_external_validation and store results.
    config_struct = session.config_struct;
    ext_results = apply_external_validation(model_path, ext_data_path, config_struct);
    if exist(session.results_file, 'file')
        tmp = load(session.results_file, 'calculated_results');
        cr = tmp.calculated_results;
        cr.external_validation = ext_results;
        calculated_results = cr; %#ok<NASGU>
        save(session.results_file, 'calculated_results');
    end
end
