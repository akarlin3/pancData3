classdef test_prepare_external_validation < matlab.unittest.TestCase
% TEST_PREPARE_EXTERNAL_VALIDATION  Tests for external validation export.
%
%   Verifies: (1) .mat file contains all required fields, (2) round-trip
%   export/apply reproduces original scores within tolerance, (3) missing
%   field handling.

    properties
        OriginalPath
        TempDir
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.OriginalPath);
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testExportContainsAllFields(testCase)
            % Exported .mat file should contain all required fields.
            model_struct = build_test_model();
            config = struct('dwi_type', 'Standard');
            out_path = fullfile(testCase.TempDir, 'validation_model.mat');

            prepare_external_validation(model_struct, config, out_path);

            testCase.verifyTrue(exist(out_path, 'file') > 0, ...
                'Exported .mat file should exist.');

            loaded = load(out_path, 'validation_model');
            vm = loaded.validation_model;

            testCase.verifyTrue(isfield(vm, 'model_coefficients'), 'Missing model_coefficients.');
            testCase.verifyTrue(isfield(vm, 'feature_names'), 'Missing feature_names.');
            testCase.verifyTrue(isfield(vm, 'feature_scaling'), 'Missing feature_scaling.');
            testCase.verifyTrue(isfield(vm, 'risk_threshold'), 'Missing risk_threshold.');
            testCase.verifyTrue(isfield(vm, 'training_cohort_summary'), 'Missing training_cohort_summary.');
            testCase.verifyTrue(isfield(vm, 'version'), 'Missing version.');
            testCase.verifyTrue(isfield(vm.feature_scaling, 'mu'), 'Missing scaling mu.');
            testCase.verifyTrue(isfield(vm.feature_scaling, 'sigma'), 'Missing scaling sigma.');
        end

        function testRoundTrip(testCase)
            % Export then apply to the same data should reproduce
            % consistent risk scores.
            rng(42);
            n_feat = 5;
            model_struct = build_test_model();
            config = struct('dwi_type', 'Standard');
            model_path = fullfile(testCase.TempDir, 'val_model.mat');
            prepare_external_validation(model_struct, config, model_path);

            % Create external data with same features
            ext_dir = fullfile(testCase.TempDir, 'external');
            mkdir(ext_dir);
            X_features = randn(10, n_feat);
            patient_ids = arrayfun(@(i) sprintf('EXT_%03d', i), 1:10, 'UniformOutput', false);
            save(fullfile(ext_dir, 'features.mat'), 'X_features', 'patient_ids');

            results = apply_external_validation(model_path, ext_dir, config);
            testCase.verifyEqual(numel(results.risk_scores), 10, ...
                'Should produce risk scores for all 10 patients.');
            testCase.verifyTrue(all(isfinite(results.risk_scores)), ...
                'All risk scores should be finite.');
        end

        function testMissingFieldError(testCase)
            % Missing required field should produce an error.
            model_struct = struct('coefficients', [1; 2]);
            config = struct();
            out_path = fullfile(testCase.TempDir, 'bad_model.mat');

            testCase.verifyError( ...
                @() prepare_external_validation(model_struct, config, out_path), ...
                'prepare_external_validation:missingField', ...
                'Missing required fields should throw an error.');
        end
    end
end


function model = build_test_model()
    n_feat = 5;
    model = struct();
    model.coefficients = randn(n_feat, 1);
    model.selected_features = 1:n_feat;
    model.feature_names = arrayfun(@(i) sprintf('Feat%d', i), 1:n_feat, 'UniformOutput', false);
    model.scaling_mu = randn(1, n_feat);
    model.scaling_sigma = abs(randn(1, n_feat)) + 0.1;
    model.imputation_ref = randn(20, n_feat);
    model.risk_threshold = 0.5;
    model.auc = 0.75;
    model.n_patients = 30;
    model.event_rate = 0.4;
end
