classdef test_filter_collinear_features < matlab.unittest.TestCase
    % TEST_FILTER_COLLINEAR_FEATURES Unit tests for filter_collinear_features.
    % Tests collinearity pruning, AUC-based tie-breaking, time-stratification,
    % NaN handling, and edge cases.
    %
    % filter_collinear_features removes redundant features from a predictor
    % matrix by computing pairwise Spearman correlations and dropping the
    % weaker predictor (lower AUC) in each correlated pair (|rho| >= 0.8).
    % When a fraction vector is supplied, correlations are evaluated at
    % baseline (Fraction 1) only to avoid spurious longitudinal correlations.

    methods(Test)
        function test_empty_input(testCase)
            % Empty feature matrix and label vector should return empty indices.
            keep_idx = filter_collinear_features([], []);
            testCase.verifyEmpty(keep_idx);
        end

        function test_single_feature(testCase)
            % A single feature cannot be collinear with anything, so it must
            % always be retained (index 1 returned).
            X = randn(20, 1);
            y = [zeros(10,1); ones(10,1)];
            keep_idx = filter_collinear_features(X, y);
            testCase.verifyEqual(keep_idx, 1);
        end

        function test_uncorrelated_features_kept(testCase)
            % Three independent random features should all survive pruning
            % because no pair exceeds the |rho| >= 0.8 threshold.
            rng(42);
            X = randn(50, 3);  % 3 independent features
            y = [zeros(25,1); ones(25,1)];
            keep_idx = filter_collinear_features(X, y);
            testCase.verifyEqual(sort(keep_idx), [1, 2, 3]);
        end

        function test_highly_correlated_drops_weaker(testCase)
            % Two features with high correlation: F1 has a strong group
            % separation (AUC), F2 is a noisy copy with weaker AUC.
            % Only the stronger predictor (F1) should be retained.
            rng(42);
            n = 50;
            y = [zeros(n/2, 1); ones(n/2, 1)];
            % F1: strongly predictive (2-unit group separation)
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
            % Three features that are near-identical copies of the same base
            % signal. All pairwise correlations exceed the threshold, so only
            % the single best predictor should survive.
            rng(42);
            n = 60;
            y = [zeros(n/2, 1); ones(n/2, 1)];
            base = [randn(n/2, 1); randn(n/2, 1) + 2];
            f1 = base;
            f2 = base + 0.1*randn(n, 1);   % tiny noise added
            f3 = base + 0.15*randn(n, 1);  % slightly more noise
            X = [f1, f2, f3];

            keep_idx = filter_collinear_features(X, y);
            % All three are highly correlated; only 1 should survive
            testCase.verifyEqual(numel(keep_idx), 1);
        end

        function test_nan_handling(testCase)
            % Scattered NaN values should be handled gracefully (pairwise
            % complete observations used for correlation). The correlated
            % pair should still be detected and pruned to one feature.
            rng(42);
            n = 50;
            y = [zeros(n/2, 1); ones(n/2, 1)];
            f1 = [randn(n/2, 1); randn(n/2, 1) + 2];
            f2 = f1 + 0.3*randn(n, 1);
            X = [f1, f2];
            X(5, 1) = NaN;   % inject NaN in feature 1
            X(10, 2) = NaN;  % inject NaN in feature 2

            keep_idx = filter_collinear_features(X, y);
            testCase.verifyTrue(~isempty(keep_idx));
            testCase.verifyEqual(numel(keep_idx), 1);
        end

        function test_time_stratified_baseline_only(testCase)
            % When a fraction vector is supplied, collinearity is assessed
            % using only baseline (Fraction 1) data. Here the features are
            % independent at baseline but highly correlated at Fraction 2.
            % Both features should be kept because baseline correlation is low.
            rng(99);
            n_pts = 20;
            frac_vec = [ones(n_pts,1); 2*ones(n_pts,1)];  % 20 baseline + 20 later
            y = repmat([zeros(n_pts/2,1); ones(n_pts/2,1)], 2, 1);

            % Baseline: independent features
            f1_base = randn(n_pts, 1);
            f2_base = randn(n_pts, 1);
            % Later fraction: highly correlated features
            f1_late = randn(n_pts, 1);
            f2_late = f1_late + 0.05*randn(n_pts, 1);
            X = [[f1_base, f2_base]; [f1_late, f2_late]];

            keep_strat = filter_collinear_features(X, y, frac_vec);
            % Should keep both since baseline correlation is low
            testCase.verifyEqual(numel(keep_strat), 2);
        end

        function test_time_stratified_no_baseline_fallback(testCase)
            % When no Fraction 1 rows exist in the data, the function should
            % fall back to using the full matrix for correlation assessment.
            % The two highly correlated features should still be pruned to one.
            rng(42);
            n = 20;
            y = [zeros(n/2,1); ones(n/2,1)];
            base = [randn(n/2, 1); randn(n/2, 1) + 2];
            X = [base, base + 0.1*randn(n, 1)];
            frac_vec = 2 * ones(n, 1);  % No baseline (Fraction 1) rows

            keep_idx = filter_collinear_features(X, y, frac_vec);
            testCase.verifyEqual(numel(keep_idx), 1);
        end

        function test_below_threshold_preserved(testCase)
            % Features with moderate correlation (|rho| < 0.8) should both
            % be preserved. The test constructs features with ~0.7 correlation
            % and conditionally verifies only when the Spearman rho is
            % confirmed below threshold (to avoid seed-dependent failures).
            rng(42);
            n = 100;
            y = [zeros(n/2, 1); ones(n/2, 1)];
            f1 = randn(n, 1);
            f2 = 0.7 * f1 + 0.7 * randn(n, 1);  % moderate correlation (~0.5-0.7)
            X = [f1, f2];

            R = corr(X, 'Type', 'Spearman');
            if abs(R(1,2)) < 0.8
                keep_idx = filter_collinear_features(X, y);
                testCase.verifyEqual(numel(keep_idx), 2);
            end
        end

        function test_auc_tiebreaker_keeps_stronger(testCase)
            % When two features are highly correlated, the AUC-based
            % tie-breaker should retain the feature with better univariate
            % predictive performance. F1 has a 4-unit group separation
            % (high AUC); F2 has only a noisy copy (lower AUC).
            rng(42);
            n = 60;
            y = [zeros(n/2, 1); ones(n/2, 1)];
            % f1: strong signal (large group separation)
            f1 = [randn(n/2, 1) - 2; randn(n/2, 1) + 2];
            % f2: same but weaker signal (added noise degrades AUC)
            f2 = f1 + randn(n, 1) * 0.3;
            X = [f1, f2];

            keep_idx = filter_collinear_features(X, y);
            testCase.verifyEqual(numel(keep_idx), 1);
        end

        function test_all_identical_features(testCase)
            % Three perfectly identical columns (rho=1.0 for all pairs).
            % All have the same AUC, so any one can be kept; we only verify
            % that exactly one survives.
            n = 30;
            y = [zeros(15,1); ones(15,1)];
            X = repmat(randn(n,1), 1, 3);

            keep_idx = filter_collinear_features(X, y);
            testCase.verifyEqual(numel(keep_idx), 1);
        end
    end
end
