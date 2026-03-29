classdef test_apply_external_validation < matlab.unittest.TestCase
% TEST_APPLY_EXTERNAL_VALIDATION — Unit tests for apply_external_validation.m
%
% Validates external validation model application:
%   - File/directory not found errors
%   - Missing features.mat and X_features variable errors
%   - Basic risk scoring with known coefficients
%   - Feature mismatch warning
%   - Zero-variance feature handling (no NaN)
%   - Feature padding with zeros for missing features
%   - Risk classification using model threshold
%   - Auto-generated patient IDs when not provided

    properties
        OriginalPath
        TempDir
        ModelPath
        ExtDir
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path;
            utilsDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils');
            addpath(utilsDir);
        end

        function createTempDir(testCase)
            testCase.TempDir = fullfile(tempdir, ['test_extval_' char(java.util.UUID.randomUUID)]);
            mkdir(testCase.TempDir);
            testCase.ModelPath = fullfile(testCase.TempDir, 'validation_model.mat');
            testCase.ExtDir = fullfile(testCase.TempDir, 'external');
            mkdir(testCase.ExtDir);
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
        function testFileNotFoundError(testCase)
            % Verify error when model_path doesn't exist.
            bad_path = fullfile(testCase.TempDir, 'no_such_model.mat');
            cfg = struct();
            ext_dir = testCase.ExtDir;

            testCase.verifyError( ...
                @() apply_external_validation(bad_path, ext_dir, cfg), ...
                'apply_external_validation:fileNotFound');
        end

        function testDirNotFoundError(testCase)
            % Verify error when external_data_path doesn't exist.
            testCase.saveModel(3);
            bad_dir = fullfile(testCase.TempDir, 'nonexistent_dir');
            cfg = struct();
            model_path = testCase.ModelPath;

            testCase.verifyError( ...
                @() apply_external_validation(model_path, bad_dir, cfg), ...
                'apply_external_validation:dirNotFound');
        end

        function testMissingFeaturesMatError(testCase)
            % Verify error when features.mat doesn't exist in external dir.
            testCase.saveModel(3);
            % ExtDir exists but has no features.mat
            cfg = struct();
            model_path = testCase.ModelPath;
            ext_dir = testCase.ExtDir;

            testCase.verifyError( ...
                @() apply_external_validation(model_path, ext_dir, cfg), ...
                'apply_external_validation:noFeatures');
        end

        function testMissingXFeaturesError(testCase)
            % Verify error when features.mat exists but lacks X_features.
            testCase.saveModel(3);
            dummy_var = 42; %#ok<NASGU>
            save(fullfile(testCase.ExtDir, 'features.mat'), 'dummy_var');
            cfg = struct();
            model_path = testCase.ModelPath;
            ext_dir = testCase.ExtDir;

            testCase.verifyError( ...
                @() apply_external_validation(model_path, ext_dir, cfg), ...
                'apply_external_validation:missingFeatures');
        end

        function testBasicScoring(testCase)
            % Create a mock model and features, verify risk_scores are computed.
            n_feat = 3;
            n_pat = 5;
            testCase.saveModel(n_feat);

            X_features = ones(n_pat, n_feat); %#ok<NASGU>
            save(fullfile(testCase.ExtDir, 'features.mat'), 'X_features');

            cfg = struct();
            evalc('results = apply_external_validation(testCase.ModelPath, testCase.ExtDir, cfg)');

            testCase.verifyEqual(numel(results.risk_scores), n_pat, ...
                'Should have one risk score per patient.');
            testCase.verifyTrue(isnumeric(results.risk_scores), ...
                'Risk scores should be numeric.');
            testCase.verifyTrue(all(results.risk_probs >= 0 & results.risk_probs <= 1), ...
                'Risk probabilities should be in [0,1].');
        end

        function testFeatureMismatchWarning(testCase)
            % When X_features has different number of columns than expected,
            % verify warning is issued.
            n_model_feat = 4;
            testCase.saveModel(n_model_feat);

            % Provide only 2 features instead of 4
            X_features = rand(3, 2); %#ok<NASGU>
            save(fullfile(testCase.ExtDir, 'features.mat'), 'X_features');

            cfg = struct();
            model_path = testCase.ModelPath;
            ext_dir = testCase.ExtDir;
            testCase.verifyWarning( ...
                @() apply_external_validation(model_path, ext_dir, cfg), ...
                'apply_external_validation:featureMismatch');
        end

        function testZeroVarianceFeatureHandled(testCase)
            % Verify zero-variance features are handled (center only, no NaN).
            n_feat = 3;
            n_pat = 4;

            % Model with sigma=0 for feature 2
            model = testCase.makeModel(n_feat);
            model.feature_scaling.sigma(2) = 0;
            validation_model = model; %#ok<NASGU>
            save(testCase.ModelPath, 'validation_model');

            X_features = rand(n_pat, n_feat); %#ok<NASGU>
            save(fullfile(testCase.ExtDir, 'features.mat'), 'X_features');

            cfg = struct();
            evalc('results = apply_external_validation(testCase.ModelPath, testCase.ExtDir, cfg)');

            testCase.verifyTrue(all(isfinite(results.risk_scores)), ...
                'Risk scores should be finite even with zero-variance features.');
        end

        function testFeaturePaddingWithZeros(testCase)
            % When external has fewer features than model expects, verify
            % padding with zeros and warning.
            n_model_feat = 5;
            testCase.saveModel(n_model_feat);

            % Provide only 3 features
            X_features = rand(4, 3); %#ok<NASGU>
            save(fullfile(testCase.ExtDir, 'features.mat'), 'X_features');

            cfg = struct();
            evalc('results = apply_external_validation(testCase.ModelPath, testCase.ExtDir, cfg)');

            testCase.verifyEqual(numel(results.risk_scores), 4, ...
                'Should produce scores for all patients despite fewer features.');
            testCase.verifyTrue(all(isfinite(results.risk_scores)), ...
                'Padded features should not introduce NaN.');
        end

        function testRiskClassification(testCase)
            % Verify risk_class uses model.risk_threshold correctly.
            n_feat = 2;
            n_pat = 6;

            model = testCase.makeModel(n_feat);
            % Set coefficients so that risk_score = sum of features
            % (mu=0, sigma=1, coefs=[1;1])
            model.feature_scaling.mu = zeros(1, n_feat);
            model.feature_scaling.sigma = ones(1, n_feat);
            model.model_coefficients = ones(n_feat, 1);
            model.risk_threshold = 0.5;  % logistic(0) = 0.5
            validation_model = model; %#ok<NASGU>
            save(testCase.ModelPath, 'validation_model');

            % Patients with large positive scores -> high risk
            % Patients with large negative scores -> low risk
            X_features = [5 5; 3 3; 1 1; -1 -1; -3 -3; -5 -5]; %#ok<NASGU>
            save(fullfile(testCase.ExtDir, 'features.mat'), 'X_features');

            cfg = struct();
            evalc('results = apply_external_validation(testCase.ModelPath, testCase.ExtDir, cfg)');

            testCase.verifyEqual(results.risk_class(1), 1, ...
                'Patient with large positive score should be high risk.');
            testCase.verifyEqual(results.risk_class(end), 0, ...
                'Patient with large negative score should be low risk.');
            testCase.verifyTrue(all(results.risk_class == 0 | results.risk_class == 1), ...
                'Risk class should be binary (0 or 1).');
        end

        function testPatientIdsGenerated(testCase)
            % When features.mat lacks patient_ids, verify auto-generated IDs.
            n_feat = 2;
            n_pat = 3;
            testCase.saveModel(n_feat);

            X_features = rand(n_pat, n_feat); %#ok<NASGU>
            % Deliberately do not save patient_ids
            save(fullfile(testCase.ExtDir, 'features.mat'), 'X_features');

            cfg = struct();
            evalc('results = apply_external_validation(testCase.ModelPath, testCase.ExtDir, cfg)');

            testCase.verifyEqual(numel(results.patient_ids), n_pat, ...
                'Should generate one ID per patient.');
            testCase.verifySubstring(results.patient_ids{1}, 'EXT_', ...
                'Auto-generated IDs should start with EXT_.');
        end
    end

    methods (Access = private)
        function model = makeModel(~, n_feat)
            % Build a minimal validation model struct.
            model = struct();
            model.version = '1.0';
            model.export_date = datestr(now);
            model.feature_names = arrayfun(@(i) sprintf('feat_%d', i), ...
                1:n_feat, 'UniformOutput', false);
            model.feature_scaling.mu = rand(1, n_feat);
            model.feature_scaling.sigma = ones(1, n_feat);
            model.model_coefficients = randn(n_feat, 1);
            model.risk_threshold = 0.5;
        end

        function saveModel(testCase, n_feat)
            % Save a default model to testCase.ModelPath.
            validation_model = testCase.makeModel(n_feat); %#ok<NASGU>
            save(testCase.ModelPath, 'validation_model');
        end
    end
end
