classdef test_run_elastic_net_cv < matlab.unittest.TestCase
% TEST_RUN_ELASTIC_NET_CV — Unit tests for run_elastic_net_cv.m
%
% Validates 5-fold elastic net CV including:
%   - Output shapes and types
%   - Lambda grid generation and optimal lambda selection
%   - Graceful failure when data is degenerate
%   - Consensus collinearity filtering via keep_fold_counts
%   - Reproducibility via fixed random seed

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
        end
    end

    methods (Test)
        function test_output_shapes(testCase)
            % Verify that all outputs have expected types and dimensions
            rng(1);
            [X, y, ids] = testCase.makeSyntheticData(40, 6);
            feat_names = arrayfun(@(i) sprintf('F%d', i), 1:6, 'UniformOutput', false);
            orig_idx = 1:6;

            [sel, opt_lam, cLam, cv_fail, kfc, coefs, fin_idx] = ...
                run_elastic_net_cv(X, y, ids, 5, false, orig_idx, feat_names, 'Fx3');

            testCase.verifyFalse(cv_fail, 'CV should not fail on well-conditioned data.');
            testCase.verifyClass(opt_lam, 'double');
            testCase.verifyGreaterThan(opt_lam, 0, 'Optimal lambda should be positive.');
            testCase.verifyFalse(isempty(cLam), 'Lambda grid should not be empty.');
            testCase.verifyTrue(isnumeric(sel), 'selected_indices should be numeric.');
            testCase.verifyTrue(isnumeric(kfc), 'keep_fold_counts should be numeric.');
            testCase.verifyEqual(numel(kfc), 6, 'keep_fold_counts should have one entry per feature.');
        end

        function test_lambda_grid_length(testCase)
            % Lambda grid should have exactly n_lambdas=10 entries (hardcoded in function)
            rng(2);
            [X, y, ids] = testCase.makeSyntheticData(40, 4);
            feat_names = arrayfun(@(i) sprintf('F%d', i), 1:4, 'UniformOutput', false);

            [~, ~, cLam, cv_fail, ~, ~, ~] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1:4, feat_names, 'Fx3');

            testCase.verifyFalse(cv_fail);
            testCase.verifyLessThanOrEqual(numel(cLam), 10, ...
                'Lambda grid should have at most 10 entries.');
            testCase.verifyGreaterThan(numel(cLam), 0);
        end

        function test_selected_indices_subset_of_original(testCase)
            % Selected feature indices must be a subset of the original indices
            rng(3);
            [X, y, ids] = testCase.makeSyntheticData(50, 8);
            feat_names = arrayfun(@(i) sprintf('F%d', i), 1:8, 'UniformOutput', false);
            orig_idx = [2, 5, 7, 10, 12, 15, 18, 20];

            [sel, ~, ~, cv_fail, ~, ~, ~] = ...
                run_elastic_net_cv(X, y, ids, 5, false, orig_idx, feat_names, 'Fx3');

            if ~cv_fail && ~isempty(sel)
                testCase.verifyTrue(all(ismember(sel, orig_idx)), ...
                    'Selected indices must be a subset of original_feature_indices.');
            end
        end

        function test_keep_fold_counts_bounded(testCase)
            % keep_fold_counts entries should be between 0 and n_folds
            rng(4);
            [X, y, ids] = testCase.makeSyntheticData(40, 5);
            feat_names = arrayfun(@(i) sprintf('F%d', i), 1:5, 'UniformOutput', false);

            [~, ~, ~, cv_fail, kfc, ~, ~] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1:5, feat_names, 'Fx3');

            if ~cv_fail
                testCase.verifyGreaterThanOrEqual(kfc, 0);
                testCase.verifyLessThanOrEqual(kfc, 5, ...
                    'Each fold count should be at most n_folds.');
            end
        end

        function test_reproducibility(testCase)
            % Two runs with same data should produce identical results
            % (function uses rng(42) internally)
            [X, y, ids] = testCase.makeSyntheticData(40, 4);
            feat_names = arrayfun(@(i) sprintf('F%d', i), 1:4, 'UniformOutput', false);

            [sel1, lam1, ~, fail1, ~, ~, ~] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1:4, feat_names, 'Fx3');
            [sel2, lam2, ~, fail2, ~, ~, ~] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1:4, feat_names, 'Fx3');

            testCase.verifyEqual(fail1, fail2);
            if ~fail1
                testCase.verifyEqual(lam1, lam2, 'AbsTol', 1e-12);
                testCase.verifyEqual(sel1, sel2);
            end
        end

        function test_cv_failed_flag_on_constant_data(testCase)
            % When all features are constant, elastic net should fail gracefully
            n = 30;
            X = ones(n, 3);
            y = [ones(15,1); zeros(15,1)];
            ids = arrayfun(@(i) sprintf('P%d', i), 1:n, 'UniformOutput', false)';
            feat_names = {'F1', 'F2', 'F3'};

            [sel, ~, ~, cv_fail, ~, ~, ~] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1:3, feat_names, 'Fx3');

            % Either cv_failed is true or no features are selected
            testCase.verifyTrue(cv_fail || isempty(sel), ...
                'Constant data should yield failure or empty selection.');
        end

        function test_final_feature_indices_consensus(testCase)
            % final_feature_indices should only include features kept in
            % majority of folds (keep_fold_counts > n_folds/2)
            rng(5);
            [X, y, ids] = testCase.makeSyntheticData(50, 6);
            feat_names = arrayfun(@(i) sprintf('F%d', i), 1:6, 'UniformOutput', false);
            orig_idx = 1:6;

            [~, ~, ~, cv_fail, kfc, ~, fin_idx] = ...
                run_elastic_net_cv(X, y, ids, 5, false, orig_idx, feat_names, 'Fx3');

            if ~cv_fail && ~isempty(fin_idx)
                % final_feature_indices correspond to features where
                % keep_fold_counts > n_folds_en/2
                consensus_mask = (kfc > 2.5);  % > 5/2
                expected_idx = orig_idx(consensus_mask);
                testCase.verifyTrue(all(ismember(fin_idx, expected_idx)), ...
                    'final_feature_indices should only include consensus features.');
            end
        end

        function test_coefs_match_selected_count(testCase)
            % When cv succeeds, coefs_en length should match the number of
            % consensus-filtered features
            rng(6);
            [X, y, ids] = testCase.makeSyntheticData(50, 5);
            feat_names = arrayfun(@(i) sprintf('F%d', i), 1:5, 'UniformOutput', false);

            [~, ~, ~, cv_fail, kfc, coefs, fin_idx] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1:5, feat_names, 'Fx3');

            if ~cv_fail && ~isempty(coefs)
                testCase.verifyEqual(numel(coefs), numel(fin_idx), ...
                    'coefs_en length should equal final_feature_indices length.');
            end
        end
    end

    methods (Static, Access = private)
        function [X, y, ids] = makeSyntheticData(n, p)
            % Generate synthetic classification data with signal
            y = [ones(floor(n/2), 1); zeros(n - floor(n/2), 1)];
            X = randn(n, p);
            % Add signal to first feature
            X(:, 1) = X(:, 1) + y * 2;
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';
        end
    end
end
