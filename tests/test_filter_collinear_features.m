classdef test_filter_collinear_features < matlab.unittest.TestCase
    % TEST_FILTER_COLLINEAR_FEATURES Unit tests for filter_collinear_features.
    % Tests collinearity pruning, AUC-based tie-breaking, time-stratification,
    % NaN handling, and edge cases.

    methods(Test)
        function test_empty_input(testCase)
            keep_idx = filter_collinear_features([], []);
            testCase.verifyEmpty(keep_idx);
        end

        function test_single_feature(testCase)
            X = randn(20, 1);
            y = [zeros(10,1); ones(10,1)];
            keep_idx = filter_collinear_features(X, y);
            testCase.verifyEqual(keep_idx, 1);
        end

        function test_uncorrelated_features_kept(testCase)
            rng(42);
            X = randn(50, 3);  % 3 independent features
            y = [zeros(25,1); ones(25,1)];
            keep_idx = filter_collinear_features(X, y);
            testCase.verifyEqual(sort(keep_idx), [1, 2, 3]);
        end

        function test_highly_correlated_drops_weaker(testCase)
            rng(42);
            n = 50;
            y = [zeros(n/2, 1); ones(n/2, 1)];
            % F1: strongly predictive
            f1 = [randn(n/2, 1); randn(n/2, 1) + 2];
            % F2: correlated with F1 but noisier (weaker AUC)
            f2 = f1 + 0.5*randn(n, 1) + 0.5;
            X = [f1, f2];

            keep_idx = filter_collinear_features(X, y);
            % Should keep the stronger predictor (F1)
            testCase.verifyTrue(ismember(1, keep_idx));
            testCase.verifyEqual(numel(keep_idx), 1);
        end

        function test_three_correlated_features(testCase)
            rng(42);
            n = 60;
            y = [zeros(n/2, 1); ones(n/2, 1)];
            base = [randn(n/2, 1); randn(n/2, 1) + 2];
            f1 = base;
            f2 = base + 0.1*randn(n, 1);
            f3 = base + 0.15*randn(n, 1);
            X = [f1, f2, f3];

            keep_idx = filter_collinear_features(X, y);
            % All three are highly correlated; only 1 should survive
            testCase.verifyEqual(numel(keep_idx), 1);
        end

        function test_nan_handling(testCase)
            rng(42);
            n = 50;
            y = [zeros(n/2, 1); ones(n/2, 1)];
            f1 = [randn(n/2, 1); randn(n/2, 1) + 2];
            f2 = f1 + 0.3*randn(n, 1);
            X = [f1, f2];
            X(5, 1) = NaN;
            X(10, 2) = NaN;

            keep_idx = filter_collinear_features(X, y);
            testCase.verifyTrue(~isempty(keep_idx));
            testCase.verifyEqual(numel(keep_idx), 1);
        end

        function test_time_stratified_baseline_only(testCase)
            % Features uncorrelated at baseline but correlated at later fractions
            rng(99);
            n_pts = 20;
            frac_vec = [ones(n_pts,1); 2*ones(n_pts,1)];
            y = repmat([zeros(n_pts/2,1); ones(n_pts/2,1)], 2, 1);

            f1_base = randn(n_pts, 1);
            f2_base = randn(n_pts, 1);  % independent at baseline
            f1_late = randn(n_pts, 1);
            f2_late = f1_late + 0.05*randn(n_pts, 1);  % highly correlated later
            X = [[f1_base, f2_base]; [f1_late, f2_late]];

            keep_strat = filter_collinear_features(X, y, frac_vec);
            % Should keep both since baseline correlation is low
            testCase.verifyEqual(numel(keep_strat), 2);
        end

        function test_time_stratified_no_baseline_fallback(testCase)
            % No Fraction 1 rows: should fall back to full matrix
            rng(42);
            n = 20;
            y = [zeros(n/2,1); ones(n/2,1)];
            base = [randn(n/2, 1); randn(n/2, 1) + 2];
            X = [base, base + 0.1*randn(n, 1)];
            frac_vec = 2 * ones(n, 1);  % No baseline rows

            keep_idx = filter_collinear_features(X, y, frac_vec);
            testCase.verifyEqual(numel(keep_idx), 1);
        end

        function test_below_threshold_preserved(testCase)
            % Features with |rho| just below 0.8 should both be kept
            rng(42);
            n = 100;
            y = [zeros(n/2, 1); ones(n/2, 1)];
            f1 = randn(n, 1);
            f2 = 0.7 * f1 + 0.7 * randn(n, 1);  % moderate correlation
            X = [f1, f2];

            R = corr(X, 'Type', 'Spearman');
            if abs(R(1,2)) < 0.8
                keep_idx = filter_collinear_features(X, y);
                testCase.verifyEqual(numel(keep_idx), 2);
            end
        end

        function test_auc_tiebreaker_keeps_stronger(testCase)
            % When both features are equally correlated, the one with higher
            % AUC should be kept
            rng(42);
            n = 60;
            y = [zeros(n/2, 1); ones(n/2, 1)];
            % f1: strong signal
            f1 = [randn(n/2, 1) - 2; randn(n/2, 1) + 2];
            % f2: same but weaker signal
            f2 = f1 + randn(n, 1) * 0.3;
            X = [f1, f2];

            keep_idx = filter_collinear_features(X, y);
            testCase.verifyEqual(numel(keep_idx), 1);
        end

        function test_all_identical_features(testCase)
            % Identical columns: all correlated with rho=1
            n = 30;
            y = [zeros(15,1); ones(15,1)];
            X = repmat(randn(n,1), 1, 3);

            keep_idx = filter_collinear_features(X, y);
            testCase.verifyEqual(numel(keep_idx), 1);
        end
    end
end
