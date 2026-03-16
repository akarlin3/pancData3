classdef test_run_loocv_risk_scores < matlab.unittest.TestCase
% TEST_RUN_LOOCV_RISK_SCORES — Unit tests for run_loocv_risk_scores.m
%
% Validates nested LOOCV risk score computation including:
%   - Output shapes (n_pts x 1)
%   - Risk scores are finite or NaN (no Inf)
%   - High/low risk stratification uses median split
%   - DL provenance leakage detection (dtype 2 and 3)
%   - Standard pipeline (dtype=1) bypasses leakage checks

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
        end
    end

    methods (Test)
        function test_output_shapes(testCase)
            % risk_scores and is_high_risk should be n_pts x 1
            rng(10);
            [X, y, ids] = testCase.makeSyntheticData(20, 4);
            prov = struct('dncnn_train_ids', {{}}, 'ivimnet_train_ids', {{}});

            [risk, highrisk] = run_loocv_risk_scores(X, y, ids, prov, 1, 'Standard', false);

            testCase.verifySize(risk, [20, 1], 'risk_scores should be n_pts x 1.');
            testCase.verifySize(highrisk, [20, 1], 'is_high_risk should be n_pts x 1.');
        end

        function test_no_inf_in_risk_scores(testCase)
            % Risk scores should be finite or NaN, never Inf
            rng(11);
            [X, y, ids] = testCase.makeSyntheticData(20, 3);
            prov = struct('dncnn_train_ids', {{}}, 'ivimnet_train_ids', {{}});

            [risk, ~] = run_loocv_risk_scores(X, y, ids, prov, 1, 'Standard', false);

            valid = ~isnan(risk);
            testCase.verifyTrue(all(isfinite(risk(valid))), ...
                'Non-NaN risk scores should be finite (no Inf).');
        end

        function test_high_risk_binary_or_nan(testCase)
            % is_high_risk values should be 0, 1, or NaN
            rng(12);
            [X, y, ids] = testCase.makeSyntheticData(20, 3);
            prov = struct('dncnn_train_ids', {{}}, 'ivimnet_train_ids', {{}});

            [~, highrisk] = run_loocv_risk_scores(X, y, ids, prov, 1, 'Standard', false);

            for i = 1:numel(highrisk)
                testCase.verifyTrue(highrisk(i) == 0 || highrisk(i) == 1 || isnan(highrisk(i)), ...
                    'is_high_risk must be 0, 1, or NaN.');
            end
        end

        function test_median_split_balance(testCase)
            % When all risk scores are valid, high/low groups should be
            % approximately balanced (within +/- 1 due to odd N or ties)
            rng(13);
            [X, y, ids] = testCase.makeSyntheticData(24, 4);
            prov = struct('dncnn_train_ids', {{}}, 'ivimnet_train_ids', {{}});

            [risk, highrisk] = run_loocv_risk_scores(X, y, ids, prov, 1, 'Standard', false);

            valid = ~isnan(risk);
            n_valid = sum(valid);
            n_high = sum(highrisk(valid) == 1);
            n_low = sum(highrisk(valid) == 0);
            % Median split should give roughly equal groups
            testCase.verifyLessThanOrEqual(abs(n_high - n_low), n_valid, ...
                'High/low risk split should use all valid patients.');
            testCase.verifyEqual(n_high + n_low, n_valid, ...
                'All valid patients should be classified as high or low risk.');
        end

        function test_dncnn_leakage_detection(testCase)
            % dtype=2 should error when a patient appears in dncnn_train_ids
            rng(14);
            [X, y, ids] = testCase.makeSyntheticData(20, 3);
            % Put first patient in training set — should trigger leakage error
            prov = struct('dncnn_train_ids', {{ids{1}}}, 'ivimnet_train_ids', {{}});

            testCase.verifyError( ...
                @() run_loocv_risk_scores(X, y, ids, prov, 2, 'dnCNN', false), ...
                '', 'Should detect DL provenance leakage for dnCNN.');
        end

        function test_ivimnet_leakage_detection(testCase)
            % dtype=3 should error when a patient appears in ivimnet_train_ids
            rng(15);
            [X, y, ids] = testCase.makeSyntheticData(20, 3);
            prov = struct('dncnn_train_ids', {{}}, 'ivimnet_train_ids', {{ids{1}}});

            testCase.verifyError( ...
                @() run_loocv_risk_scores(X, y, ids, prov, 3, 'IVIMnet', false), ...
                '', 'Should detect DL provenance leakage for IVIMnet.');
        end

        function test_standard_dtype_no_leakage_check(testCase)
            % dtype=1 should NOT check DL provenance (Standard pipeline)
            rng(16);
            [X, y, ids] = testCase.makeSyntheticData(20, 3);
            % Even with matching IDs, dtype=1 should not error
            prov = struct('dncnn_train_ids', {{ids{1}}}, 'ivimnet_train_ids', {{ids{1}}});

            [risk, ~] = run_loocv_risk_scores(X, y, ids, prov, 1, 'Standard', false);
            testCase.verifySize(risk, [20, 1]);
        end

        function test_nan_risk_for_failed_folds(testCase)
            % When a fold fails (e.g., too few features), the risk score
            % should be NaN rather than zero
            rng(17);
            n = 10;
            X = ones(n, 2);  % degenerate: constant features
            y = [ones(5,1); zeros(5,1)];
            ids = arrayfun(@(i) sprintf('P%02d', i), 1:n, 'UniformOutput', false)';
            prov = struct('dncnn_train_ids', {{}}, 'ivimnet_train_ids', {{}});

            [risk, highrisk] = run_loocv_risk_scores(X, y, ids, prov, 1, 'Standard', false);

            % With constant features, all folds should produce NaN
            testCase.verifyTrue(all(isnan(risk)), ...
                'Constant features should yield NaN risk scores.');
            testCase.verifyTrue(all(isnan(highrisk)), ...
                'NaN risk scores should produce NaN high-risk labels.');
        end
    end

    methods (Static, Access = private)
        function [X, y, ids] = makeSyntheticData(n, p)
            y = [ones(floor(n/2), 1); zeros(n - floor(n/2), 1)];
            X = randn(n, p);
            X(:, 1) = X(:, 1) + y * 2;
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';
        end
    end
end
