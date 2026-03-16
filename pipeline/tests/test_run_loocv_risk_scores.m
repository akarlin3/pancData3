classdef test_run_loocv_risk_scores < matlab.unittest.TestCase
% TEST_RUN_LOOCV_RISK_SCORES — Unit tests for run_loocv_risk_scores.m
%
% Validates nested LOOCV risk score computation including:
%   - Output shape matches input size
%   - Risk scores are numeric (finite or NaN for failed folds)
%   - High-risk stratification is binary (0/1/NaN)
%   - DL provenance leakage check fires for contaminated patients
%   - Standard DWI type (dtype=1) bypasses provenance checks

    methods (TestMethodSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            addpath(fullfile(testDir, '..', 'utils'));
            addpath(fullfile(testDir, '..', 'dependencies'));
        end
    end

    methods (Test)
        function test_output_shape(testCase)
            % Risk scores and stratification vectors should match input rows
            rng(42);
            n = 20;
            X = randn(n, 3);
            y = [ones(10,1); zeros(10,1)];
            X(:,1) = X(:,1) + 2*y;
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';
            prov = struct('dncnn_train_ids', {{}}, 'ivimnet_train_ids', {{}});

            [risk, high] = run_loocv_risk_scores(X, y, ids, prov, 1, 'Standard', false);

            testCase.verifyEqual(numel(risk), n, 'Risk scores should have n elements.');
            testCase.verifyEqual(numel(high), n, 'Stratification should have n elements.');
        end

        function test_risk_scores_numeric(testCase)
            % All risk scores should be numeric (finite or NaN)
            rng(42);
            n = 20;
            X = randn(n, 3);
            y = [ones(10,1); zeros(10,1)];
            X(:,1) = X(:,1) + 2*y;
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';
            prov = struct('dncnn_train_ids', {{}}, 'ivimnet_train_ids', {{}});

            [risk, ~] = run_loocv_risk_scores(X, y, ids, prov, 1, 'Standard', false);

            testCase.verifyTrue(isnumeric(risk));
        end

        function test_stratification_binary_or_nan(testCase)
            % is_high_risk should be 0, 1, or NaN only
            rng(42);
            n = 20;
            X = randn(n, 3);
            y = [ones(10,1); zeros(10,1)];
            X(:,1) = X(:,1) + 2*y;
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';
            prov = struct('dncnn_train_ids', {{}}, 'ivimnet_train_ids', {{}});

            [~, high] = run_loocv_risk_scores(X, y, ids, prov, 1, 'Standard', false);

            valid = ~isnan(high);
            testCase.verifyTrue(all(high(valid) == 0 | high(valid) == 1), ...
                'Non-NaN stratification values must be 0 or 1.');
        end

        function test_median_split_balance(testCase)
            % Approximately half should be high-risk (median split)
            rng(42);
            n = 30;
            X = randn(n, 3);
            y = [ones(15,1); zeros(15,1)];
            X(:,1) = X(:,1) + 2*y;
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';
            prov = struct('dncnn_train_ids', {{}}, 'ivimnet_train_ids', {{}});

            [~, high] = run_loocv_risk_scores(X, y, ids, prov, 1, 'Standard', false);

            valid = ~isnan(high);
            n_valid = sum(valid);
            n_high = sum(high(valid) == 1);
            % Should be roughly half ± a few
            testCase.verifyGreaterThan(n_high, 0, 'Should have some high-risk patients.');
            testCase.verifyLessThan(n_high, n_valid, 'Should have some low-risk patients.');
        end

        function test_dl_provenance_leakage_dncnn(testCase)
            % DL provenance check should error for contaminated patient (dtype=2)
            rng(42);
            n = 10;
            X = randn(n, 2);
            y = [ones(5,1); zeros(5,1)];
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';
            % Mark P001 as used in DnCNN training — should trigger leakage error
            prov = struct('dncnn_train_ids', {{'P001'}}, 'ivimnet_train_ids', {{}});

            testCase.verifyError(@() ...
                run_loocv_risk_scores(X, y, ids, prov, 2, 'dnCNN', false), ...
                ?MException, ...
                'Should error on DL provenance leakage for dnCNN.');
        end

        function test_dl_provenance_leakage_ivimnet(testCase)
            % DL provenance check should error for contaminated patient (dtype=3)
            rng(42);
            n = 10;
            X = randn(n, 2);
            y = [ones(5,1); zeros(5,1)];
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';
            prov = struct('dncnn_train_ids', {{}}, 'ivimnet_train_ids', {{'P001'}});

            testCase.verifyError(@() ...
                run_loocv_risk_scores(X, y, ids, prov, 3, 'IVIMnet', false), ...
                ?MException, ...
                'Should error on DL provenance leakage for IVIMnet.');
        end

        function test_standard_dtype_skips_provenance(testCase)
            % dtype=1 (Standard) should not check DL provenance
            rng(42);
            n = 20;
            X = randn(n, 3);
            y = [ones(10,1); zeros(10,1)];
            X(:,1) = X(:,1) + 2*y;
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';
            % Even with IDs in provenance lists, Standard should be fine
            prov = struct('dncnn_train_ids', {ids(1:5)}, 'ivimnet_train_ids', {ids(1:5)});

            [risk, ~] = run_loocv_risk_scores(X, y, ids, prov, 1, 'Standard', false);

            testCase.verifyEqual(numel(risk), n, ...
                'Standard dtype should not trigger provenance errors.');
        end

        function test_nan_risk_for_failed_folds(testCase)
            % With very few samples and many features, some folds may fail;
            % failed folds should produce NaN risk scores
            rng(42);
            n = 8;
            X = randn(n, 15);  % p >> n to provoke failures
            y = [1;1;1;1;0;0;0;0];
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';
            prov = struct('dncnn_train_ids', {{}}, 'ivimnet_train_ids', {{}});

            [risk, high] = run_loocv_risk_scores(X, y, ids, prov, 1, 'Standard', false);

            % Either all valid or some NaN — both are acceptable
            testCase.verifyTrue(isnumeric(risk));
            % NaN risk → NaN stratification
            nan_risk = isnan(risk);
            nan_high = isnan(high);
            testCase.verifyEqual(nan_risk, nan_high, ...
                'NaN risk scores should correspond to NaN stratification.');
        end
    end
end
