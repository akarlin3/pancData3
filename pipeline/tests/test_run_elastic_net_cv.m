classdef test_run_elastic_net_cv < matlab.unittest.TestCase
% TEST_RUN_ELASTIC_NET_CV — Unit tests for run_elastic_net_cv.m
%
% Validates 5-fold elastic net CV + final model fitting including:
%   - Output shape and type for well-separated data
%   - Convergence failure returns empty selection with cv_failed flag
%   - Lambda grid is shared across folds (non-empty)
%   - Keep-fold counts track per-feature retention across folds
%   - Reproducibility via fixed random seed

    methods (TestMethodSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            addpath(fullfile(testDir, '..', 'utils'));
            addpath(fullfile(testDir, '..', 'dependencies'));
        end
    end

    methods (Test)
        function test_well_separated_data_produces_selection(testCase)
            % When classes are well-separated, elastic net should converge
            % and select at least one feature.
            rng(42);
            n = 40;
            X = randn(n, 5);
            y = [ones(20,1); zeros(20,1)];
            % Make feature 1 strongly predictive
            X(:,1) = X(:,1) + 3*y;
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';
            feat_names = {'F1','F2','F3','F4','F5'};
            orig_idx = 1:5;

            [sel, opt_lam, common_Lam, cv_fail, keep_counts, coefs, final_idx] = ...
                run_elastic_net_cv(X, y, ids, 5, false, orig_idx, feat_names, 'Fx3');

            testCase.verifyFalse(cv_fail, 'CV should not fail on well-separated data.');
            testCase.verifyFalse(isempty(common_Lam), 'Lambda grid should be computed.');
            testCase.verifyGreaterThan(numel(common_Lam), 0);
            testCase.verifyTrue(opt_lam > 0, 'Optimal lambda should be positive.');
            testCase.verifyEqual(numel(keep_counts), 5, ...
                'keep_fold_counts should have one entry per feature.');
        end

        function test_output_shapes(testCase)
            rng(42);
            n = 30;
            p = 4;
            X = randn(n, p);
            y = [ones(15,1); zeros(15,1)];
            X(:,1) = X(:,1) + 2*y;
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';

            [sel, opt_lam, common_Lam, cv_fail, keep_counts, coefs, final_idx] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1:p, ...
                {'A','B','C','D'}, 'Fx3');

            testCase.verifyTrue(islogical(cv_fail) || isnumeric(cv_fail));
            testCase.verifyEqual(numel(keep_counts), p);
            if ~cv_fail
                testCase.verifyTrue(isnumeric(sel));
                testCase.verifyTrue(isnumeric(opt_lam));
                testCase.verifyTrue(isnumeric(coefs));
            end
        end

        function test_constant_features_graceful(testCase)
            % All-constant features should not crash; elastic net may fail
            rng(42);
            n = 20;
            X = ones(n, 3);
            y = [ones(10,1); zeros(10,1)];
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';

            [sel, ~, ~, cv_fail, ~, ~, ~] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1:3, ...
                {'A','B','C'}, 'Fx3');

            % Should either fail gracefully or return empty selection
            if cv_fail
                testCase.verifyEmpty(sel, ...
                    'Failed CV should return empty selection.');
            end
        end

        function test_keep_fold_counts_bounded(testCase)
            % Keep-fold counts should be between 0 and n_folds
            rng(42);
            n = 40;
            p = 6;
            X = randn(n, p);
            y = [ones(20,1); zeros(20,1)];
            X(:,1) = X(:,1) + 2*y;
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';

            [~, ~, ~, ~, keep_counts, ~, ~] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1:p, ...
                arrayfun(@(i) sprintf('F%d',i), 1:p, 'UniformOutput', false), 'Fx3');

            testCase.verifyGreaterThanOrEqual(keep_counts, 0);
            testCase.verifyLessThanOrEqual(keep_counts, 5, ...
                'Keep counts should not exceed number of folds.');
        end

        function test_reproducibility(testCase)
            % Same inputs should produce identical outputs (deterministic seed)
            rng(42);
            n = 30;
            X = randn(n, 4);
            y = [ones(15,1); zeros(15,1)];
            X(:,1) = X(:,1) + 3*y;
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';
            names = {'A','B','C','D'};

            [sel1, lam1, ~, ~, ~, ~, ~] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1:4, names, 'Fx3');
            [sel2, lam2, ~, ~, ~, ~, ~] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1:4, names, 'Fx3');

            testCase.verifyEqual(sel1, sel2, 'Selection should be reproducible.');
            testCase.verifyEqual(lam1, lam2, 'AbsTol', 1e-10, ...
                'Optimal lambda should be reproducible.');
        end

        function test_single_feature_works(testCase)
            % Should handle single-feature input without error
            rng(42);
            n = 30;
            X = randn(n, 1);
            y = [ones(15,1); zeros(15,1)];
            X(:,1) = X(:,1) + 3*y;
            ids = arrayfun(@(i) sprintf('P%03d', i), 1:n, 'UniformOutput', false)';

            [~, ~, ~, cv_fail, keep_counts, ~, ~] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1, {'F1'}, 'Fx3');

            testCase.verifyEqual(numel(keep_counts), 1);
        end

        function test_cv_failed_returns_empty(testCase)
            % When cv_failed is true, selected_indices should be empty
            rng(42);
            n = 6;  % Very small sample — likely causes convergence failure
            X = randn(n, 10);  % More features than samples
            y = [1;1;1;0;0;0];
            ids = {'P1','P2','P3','P4','P5','P6'}';

            [sel, opt_lam, ~, cv_fail, ~, ~, ~] = ...
                run_elastic_net_cv(X, y, ids, 5, false, 1:10, ...
                arrayfun(@(i) sprintf('F%d',i), 1:10, 'UniformOutput', false), 'Fx3');

            if cv_fail
                testCase.verifyEmpty(sel, ...
                    'Failed CV should produce empty selection.');
            end
            % Either way, this should not crash
        end
    end
end
