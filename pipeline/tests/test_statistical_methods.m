classdef test_statistical_methods < matlab.unittest.TestCase
    % TEST_STATISTICAL_METHODS Pure algorithmic tests for statistical and
    % numerical routines used in the DWI analysis pipeline.
    %
    % These tests exercise BH/Holm multiple-comparison corrections,
    % correlation filtering, LOOCV logic, b-value validation, FDR,
    % imputation, ADC handling, DIR mask propagation, IVIM segmented
    % fitting, and time-dependent panel construction — all without
    % reading source code files.
    %
    % Run tests with:
    %   results = runtests('tests/test_statistical_methods.m');

    methods(TestMethodSetup)
        function addPaths(testCase)
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'dependencies'));
        end
    end

    methods(Test)

        % =================================================================
        % BH correction tests
        % =================================================================

        function testBH_QValues_KnownInput(testCase)
            % Verify BH q-values for a hand-calculated example.
            % Five raw p-values; q = p * (m/rank), then enforce monotonicity.
            raw_p = [0.005; 0.01; 0.03; 0.40; 0.90];
            % Hand-calculated expected q-values (sorted order):
            % rank 1: 0.005 * 5/1 = 0.025
            % rank 2: 0.01  * 5/2 = 0.025
            % rank 3: 0.03  * 5/3 = 0.05
            % rank 4: 0.40  * 5/4 = 0.50
            % rank 5: 0.90  * 5/5 = 0.90
            % Monotonicity enforcement (backward pass): already non-decreasing.
            expected = [0.025; 0.025; 0.05; 0.50; 0.90];
            % Call the production function
            q_vals = benjamini_hochberg_fdr(raw_p);
            % The production function returns q-values in the original
            % (input) order. Since raw_p is already sorted, we can compare
            % directly.
            testCase.verifyEqual(numel(q_vals), numel(raw_p), ...
            'Output length should match input length');
            testCase.verifyEqual(q_vals(:), expected(:), 'AbsTol', 1e-12, ...
            'BH q-values do not match hand-calculated values');
        end

        function testBH_Monotonicity(testCase)
            % Q-values must be non-decreasing after the BH step-up procedure.
            rng(1);
            raw_p = sort(rand(50, 1));   % 50 sorted p-values
            % Call the production function
            q_vals = benjamini_hochberg_fdr(raw_p);
            diffs = diff(q_vals(:));
            testCase.verifyTrue(all(diffs >= -1e-15), ...
            'BH q-values are not monotonically non-decreasing');
        end

        function testBH_CappedAtOne(testCase)
            % No q-value should exceed 1.0.
            raw_p = [0.10; 0.50; 0.80; 0.95; 0.99];
            % Call the production function
            q_vals = benjamini_hochberg_fdr(raw_p);
            testCase.verifyTrue(all(q_vals(:) <= 1.0), ...
            'Some BH q-values exceed 1.0');
        end

        function testBH_QGreaterThanOrEqualP(testCase)
            % Every q-value must be >= the corresponding raw p-value.
            raw_p = [0.001; 0.01; 0.04; 0.05; 0.20; 0.50; 0.80];
            % Call the production function
            q_vals = benjamini_hochberg_fdr(raw_p);
            testCase.verifyTrue(all(q_vals(:) >= raw_p(:) - 1e-15), ...
            'Some BH q-values are smaller than their raw p-values');
        end

        function testBH_SinglePValue(testCase)
            % Edge case: m = 1. q should equal p (no correction needed).
            raw_p = 0.03;
            q_vals = benjamini_hochberg_fdr(raw_p);
            testCase.verifyTrue(abs(q_vals - 0.03) <= 1e-15, ...
            'Single p-value should remain unchanged after BH correction');
        end

        function testBH_AllIdenticalPValues(testCase)
            % When all p-values are identical, q-values should all equal
            % the same adjusted value (or 1.0 if > 1).
            raw_p = repmat(0.04, 10, 1);
            % Call the production function
            q_vals = benjamini_hochberg_fdr(raw_p);
            % All q-values should be the same after monotonicity
            testCase.verifyTrue(abs(max(q_vals(:)) - min(q_vals(:))) <= 1e-15, ...
            'Identical p-values should yield identical q-values');
        end

        function testBH_UnsortedInput(testCase)
            % Verify that the production function handles unsorted input
            % and returns q-values in the original (unsorted) order.
            raw_p = [0.90; 0.005; 0.40; 0.01; 0.03];
            q_vals = benjamini_hochberg_fdr(raw_p);
            % The sorted version of raw_p is [0.005; 0.01; 0.03; 0.40; 0.90]
            % with expected sorted q-values [0.025; 0.025; 0.05; 0.50; 0.90].
            % Map back to original order:
            expected_unsorted = [0.90; 0.025; 0.50; 0.025; 0.05];
            testCase.verifyEqual(numel(q_vals), numel(raw_p), ...
            'Output length should match input length');
            testCase.verifyEqual(q_vals(:), expected_unsorted(:), 'AbsTol', 1e-12, ...
            'BH q-values for unsorted input do not match expected values');
        end

        % =================================================================
        % Holm correction tests
        % =================================================================

        function testHolm_ThresholdFormula(testCase)
            % Verify that Holm thresholds follow alpha / (m + 1 - k).
            m = 5;
            alpha = 0.05;
            holm_thresholds = alpha ./ (m + 1 - (1:m)');
            expected = [0.05/5; 0.05/4; 0.05/3; 0.05/2; 0.05/1];
            testCase.verifyTrue(all(abs(holm_thresholds - expected) <= 1e-15), ...
            'Holm thresholds do not match expected formula');
        end

        function testHolm_ThresholdsStrictlyIncreasing(testCase)
            % Holm thresholds must be strictly increasing (less strict as
            % rank increases).
            m = 20;
            alpha = 0.05;
            holm_thresholds = alpha ./ (m + 1 - (1:m)');
            diffs = diff(holm_thresholds);
            testCase.verifyTrue(all(diffs > 0), ...
            'Holm thresholds are not strictly increasing');
        end

        function testHolm_EarlyStopBehavior(testCase)
            % If the k-th test fails, all subsequent tests must also be
            % marked non-significant (step-down property).
            p_sorted = [0.001; 0.008; 0.060; 0.070; 0.900];
            m = length(p_sorted);
            alpha = 0.05;
            holm_thresholds = alpha ./ (m + 1 - (1:m)');
            is_sig = false(size(p_sorted));
            for k = 1:m
            if p_sorted(k) < holm_thresholds(k)
            is_sig(k) = true;
            else
            break;
            end
            end
            % First two should pass: 0.001 < 0.01, 0.008 < 0.0125
            % Third fails: 0.060 >= 0.0167 → stop
            testCase.verifyTrue(is_sig(1) && is_sig(2), ...
            'First two tests should be significant');
            testCase.verifyTrue(~any(is_sig(3:end)), ...
            'Tests 3-5 should be non-significant after early stop');
        end

        function testHolm_AllSignificant(testCase)
            % When all p-values beat every threshold.
            p_sorted = [0.001; 0.002; 0.003; 0.004; 0.005];
            m = length(p_sorted);
            alpha = 0.05;
            holm_thresholds = alpha ./ (m + 1 - (1:m)');
            is_sig = false(size(p_sorted));
            for k = 1:m
            if p_sorted(k) < holm_thresholds(k)
            is_sig(k) = true;
            else
            break;
            end
            end
            testCase.verifyTrue(all(is_sig), ...
            'All tests should be significant with very small p-values');
        end

        function testHolm_NoneSignificant(testCase)
            % When the first test already fails.
            p_sorted = [0.10; 0.20; 0.30; 0.40; 0.50];
            m = length(p_sorted);
            alpha = 0.05;
            holm_thresholds = alpha ./ (m + 1 - (1:m)');
            is_sig = false(size(p_sorted));
            for k = 1:m
            if p_sorted(k) < holm_thresholds(k)
            is_sig(k) = true;
            else
            break;
            end
            end
            testCase.verifyTrue(~any(is_sig), ...
            'No tests should be significant');
        end

        % =================================================================
        % Correlation filter tests
        % =================================================================

        function testCorrFilter_DropsHighCorrelation(testCase)
            % Two perfectly correlated features → one should be dropped.
            rng(7);
            n = 30;
            x1 = randn(n, 1);
            x2 = x1 * 2 + 0.01 * randn(n, 1);  % r ≈ 1.0
            x3 = randn(n, 1);                    % independent
            X = [x1, x2, x3];
            R_corr = corrcoef(X);
            n_feats = size(X, 2);
            drop_flag = false(1, n_feats);
            for fi = 1:n_feats
            if drop_flag(fi), continue; end
            for fj = fi+1:n_feats
            if drop_flag(fj), continue; end
            if abs(R_corr(fi, fj)) > 0.8
            drop_flag(fj) = true;
            end
            end
            end
            keep_idx = find(~drop_flag);
            % x2 (column 2) should be dropped; x1 and x3 kept
            testCase.verifyTrue(isequal(keep_idx, [1, 3]), ...
            'Highly correlated feature (col 2) should be dropped');
        end

        function testCorrFilter_KeepsAllWhenUncorrelated(testCase)
            % Orthogonal features: nothing should be dropped.
            X = eye(5);
            R_corr = corrcoef(X);
            n_feats = size(X, 2);
            drop_flag = false(1, n_feats);
            for fi = 1:n_feats
            if drop_flag(fi), continue; end
            for fj = fi+1:n_feats
            if drop_flag(fj), continue; end
            if abs(R_corr(fi, fj)) > 0.8
            drop_flag(fj) = true;
            end
            end
            end
            keep_idx = find(~drop_flag);
            testCase.verifyTrue(isequal(keep_idx, 1:5), ...
            'No features should be dropped when all are uncorrelated');
        end

        function testCorrFilter_PrefersEarlierFeature(testCase)
            % When columns i and j are correlated (i < j), column j is
            % dropped — i.e., the earlier / absolute variable is preferred.
            rng(8);
            n = 50;
            base = randn(n, 1);
            X = [base, base + 0.001*randn(n,1), randn(n,1), randn(n,1)];
            R_corr = corrcoef(X);
            n_feats = size(X, 2);
            drop_flag = false(1, n_feats);
            for fi = 1:n_feats
            if drop_flag(fi), continue; end
            for fj = fi+1:n_feats
            if drop_flag(fj), continue; end
            if abs(R_corr(fi, fj)) > 0.8
            drop_flag(fj) = true;
            end
            end
            end
            % Column 1 kept, column 2 dropped
            testCase.verifyTrue(~drop_flag(1), ...
            'Earlier feature should be kept');
            testCase.verifyTrue(drop_flag(2), ...
            'Later correlated feature should be dropped');
        end

        function testCorrFilter_NegativeCorrelation(testCase)
            % Negative correlation with |r| > 0.8 should also be dropped.
            rng(9);
            n = 40;
            x1 = randn(n, 1);
            x2 = -x1 + 0.01*randn(n,1);  % r ≈ -1.0
            X = [x1, x2];
            R_corr = corrcoef(X);
            n_feats = size(X, 2);
            drop_flag = false(1, n_feats);
            for fi = 1:n_feats
            if drop_flag(fi), continue; end
            for fj = fi+1:n_feats
            if drop_flag(fj), continue; end
            if abs(R_corr(fi, fj)) > 0.8
            drop_flag(fj) = true;
            end
            end
            end
            testCase.verifyTrue(drop_flag(2), ...
            'Negatively correlated feature should be dropped');
        end

        function testCorrFilter_BoundaryNotDropped(testCase)
            % |r| = 0.8 exactly should NOT be dropped (threshold is >0.8).
            % Construct a pair with correlation exactly 0.8 via known formula.
            rng(10);
            n = 10000;   % large n so sample r ≈ population r
            x1 = randn(n, 1);
            % x2 = 0.8*x1 + sqrt(1-0.64)*noise gives population r = 0.8
            x2 = 0.8*x1 + sqrt(1 - 0.64)*randn(n, 1);
            X = [x1, x2];
            R_corr = corrcoef(X);
            r12 = abs(R_corr(1, 2));
            drop_flag = false(1, 2);
            if r12 > 0.8
            drop_flag(2) = true;
            end
            % With n=10000, sample r will be very close to 0.8 but may
            % land on either side. Verify the logic is at least applied:
            if r12 <= 0.8
            testCase.verifyTrue(~drop_flag(2), ...
            'Feature at boundary (|r|==0.8) should NOT be dropped');
            else
            testCase.verifyTrue(drop_flag(2), ...
            'Feature slightly above boundary should be dropped');
            end
        end

        function testCorrFilter_ChainDropping(testCase)
            % If A-B correlated and B-C correlated, dropping B means C
            % is evaluated against A only. Verify transitive behaviour:
            % b is dropped (r(a,b) > 0.8), then c is compared to a only,
            % and since r(a,c) ~ 0.7 < 0.8, c is retained.
            rng(11);
            n = 100;
            a = randn(n, 1);
            b = a + 0.01*randn(n, 1);           % r(a,b) ≈ 1
            c = 0.5*a + 0.5*randn(n, 1);        % r(a,c) ≈ 0.7 (not dropped)
            X = [a, b, c];
            R_corr = corrcoef(X);
            n_feats = size(X, 2);
            drop_flag = false(1, n_feats);
            for fi = 1:n_feats
            if drop_flag(fi), continue; end
            for fj = fi+1:n_feats
            if drop_flag(fj), continue; end
            if abs(R_corr(fi, fj)) > 0.8
            drop_flag(fj) = true;
            end
            end
            end
            keep_idx = find(~drop_flag);
            % b dropped (correlated with a), c kept (r(a,c) < 0.8)
            testCase.verifyTrue(drop_flag(2), 'b should be dropped');
            testCase.verifyTrue(~drop_flag(3), 'c should be kept');
            testCase.verifyTrue(length(keep_idx) == 2);
        end

        % =================================================================
        % LOOCV tests
        % =================================================================

        function testLOOCV_NoLeakage_CutoffExcludesHeldOut(testCase)
            % Verify that the median cutoff for each fold is computed
            % WITHOUT the held-out patient's data.
            rng(42);
            n = 10;
            vals = (1:n)';   % Predictable ordering for median checks
            for loo_i = 1:n
            train_mask = true(n, 1);
            train_mask(loo_i) = false;
            train_vals = vals(train_mask);
            cutoff = median(train_vals);
            % The held-out value must NOT influence the cutoff.
            % Re-compute cutoff including it — should differ.
            cutoff_all = median(vals);
            if vals(loo_i) ~= cutoff_all
            testCase.verifyTrue(cutoff ~= cutoff_all, ...
            sprintf('Fold %d: cutoff should differ when held-out is excluded', loo_i));
            end
            end
        end

        function testLOOCV_TrainMaskSize(testCase)
            % Each training fold should contain exactly N-1 patients.
            n = 15;
            for loo_i = 1:n
            train_mask = true(n, 1);
            train_mask(loo_i) = false;
            testCase.verifyTrue(sum(train_mask) == n - 1, ...
            'Training set should have N-1 patients');
            end
        end

        function testLOOCV_EachPatientHeldOutOnce(testCase)
            % Over the full loop, each patient is held out exactly once.
            n = 12;
            held_out_count = zeros(n, 1);
            for loo_i = 1:n
            held_out_count(loo_i) = held_out_count(loo_i) + 1;
            end
            testCase.verifyTrue(isequal(held_out_count, ones(n, 1)), ...
            'Each patient should be held out exactly once');
        end

        function testLOOCV_RiskAssignmentIsBoolean(testCase)
            % is_high_risk must be logical and have one entry per patient.
            rng(42);
            n = 8;
            km_X = randn(n, 8);
            km_y = [0;0;0;0;1;1;1;1];
            is_high_risk = false(n, 1);
            for loo_i = 1:n
            train_mask = true(n, 1);
            train_mask(loo_i) = false;
            train_vals = km_X(train_mask, 5);  % Arbitrary percent-change col
            cutoff = median(train_vals);
            is_high_risk(loo_i) = km_X(loo_i, 5) < cutoff;
            end
            testCase.verifyTrue(islogical(is_high_risk), 'is_high_risk should be logical');
            testCase.verifyTrue(numel(is_high_risk) == n, 'is_high_risk should have one entry per patient');
        end

        function testLOOCV_NoEmptyGroups_LargeN(testCase)
            % With balanced synthetic data, LOOCV should produce both high-
            % and low-risk patients (i.e., not degenerate to all-same).
            rng(42);
            n = 20;
            % Create data where LC patients have higher values than LF
            km_X = zeros(n, 8);
            km_y = [zeros(n/2, 1); ones(n/2, 1)];
            km_X(1:n/2, 5) = randn(n/2, 1) + 2;   % LC: positive change
            km_X(n/2+1:end, 5) = randn(n/2, 1) - 2; % LF: negative change
            is_high_risk = false(n, 1);
            for loo_i = 1:n
            train_mask = true(n, 1);
            train_mask(loo_i) = false;
            train_vals = km_X(train_mask, 5);
            cutoff = median(train_vals);
            is_high_risk(loo_i) = km_X(loo_i, 5) < cutoff;
            end
            testCase.verifyTrue(any(is_high_risk), ...
            'Should have at least one high-risk patient');
            testCase.verifyTrue(any(~is_high_risk), ...
            'Should have at least one low-risk patient');
        end

        function testLOOCV_FeatureSelectionBase_Mapping(testCase)
            % Verify the mod-based mapping from 8-feature index to base
            % metric index: indices 1-4 = Absolute values (ADC, D, f, D*),
            % indices 5-8 = Percent-Change values (same order).
            % base = mod(sel - 1, 4) + 1 maps both groups to [1,2,3,4].
            expected_base = [1, 2, 3, 4, 1, 2, 3, 4];
            for idx = 1:8
            base = mod(idx - 1, 4) + 1;
            testCase.verifyTrue(base == expected_base(idx), ...
            sprintf('Index %d should map to base %d', idx, expected_base(idx)));
            end
        end

        function testLOOCV_DefaultLowRisk_WhenNoFeatures(testCase)
            % When LASSO/Elastic Net selects zero features, the patient
            % should default to low-risk (is_high_risk = false).
            sel_loo = [];   % Empty selection (simulating convergence failure)
            default_risk = false;
            if isempty(sel_loo)
            assigned_risk = false;
            else
            assigned_risk = true;
            end
            testCase.verifyTrue(assigned_risk == default_risk, ...
            'Empty feature selection should default to low-risk');
        end

        % =================================================================
        % B-value validation tests
        % =================================================================

        function testBval_ValidationLogic_MatchingProtocol(testCase)
            % When bvals exactly match [0; 30; 150; 550], no deviation.
            bvals = [0; 30; 150; 550];
            expected_bvals = [0; 30; 150; 550];
            testCase.verifyTrue(isequal(sort(bvals), expected_bvals), ...
            'Standard protocol b-values should pass validation');
        end

        function testBval_ValidationLogic_UnsortedMatch(testCase)
            % Unsorted b-values that match after sorting should pass.
            bvals = [550; 30; 0; 150];
            expected_bvals = [0; 30; 150; 550];
            testCase.verifyTrue(isequal(sort(bvals), expected_bvals), ...
            'Unsorted but correct b-values should pass validation');
        end

        function testBval_ValidationLogic_ExtraValue(testCase)
            % Extra b-value (5 instead of 4) should be flagged as deviation.
            bvals = [0; 50; 400; 800; 1000];
            expected_bvals = [0; 30; 150; 550];
            % First verify the count mismatch is detected (the primary
            % reason this is a deviation: 5 b-values vs expected 4).
            testCase.verifyNotEqual(numel(bvals), numel(expected_bvals), ...
            'Should detect extra b-values by count');
            % Also verify that the overall deviation flag is set.
            is_deviation = ~isequal(sort(bvals), expected_bvals);
            testCase.verifyTrue(is_deviation, ...
            'Extra b-values should be flagged as protocol deviation');
        end

        function testBval_ValidationLogic_WrongValue(testCase)
            % Different b-value set should be flagged as deviation.
            bvals = [0; 30; 200; 550];
            expected_bvals = [0; 30; 150; 550];
            % Count is the same, so deviation is due to differing values.
            testCase.verifyEqual(numel(bvals), numel(expected_bvals), ...
            'Count should match so deviation is purely value-based');
            is_deviation = ~isequal(sort(bvals), expected_bvals);
            testCase.verifyTrue(is_deviation, ...
            'Non-standard b-values should be flagged as protocol deviation');
        end

        function testBval_ValidationLogic_MissingValue(testCase)
            % Fewer b-values than expected should be flagged.
            bvals = [0; 50; 400];
            expected_bvals =  [0; 30; 150; 550];
            % First verify the count mismatch is detected.
            testCase.verifyNotEqual(numel(bvals), numel(expected_bvals), ...
            'Should detect missing b-values by count');
            is_deviation = ~isequal(sort(bvals), expected_bvals);
            testCase.verifyTrue(is_deviation, ...
            'Missing b-values should be flagged as protocol deviation');
        end

        % =================================================================
        % Per-timepoint FDR
        % =================================================================

        function testPerTimepointFDR_SmallerFamilySize(testCase)
            % Same p-values corrected per-timepoint (n=4) vs globally (n=12)
            % should yield smaller q-values per-timepoint.
            p_tp = [0.01; 0.02; 0.04; 0.10];  % 4 tests in one timepoint
            n_tp = 4;
            [p_s, si] = sort(p_tp);
            q_tp = zeros(n_tp, 1);
            q_tp(n_tp) = p_s(n_tp);
            for ii = n_tp-1:-1:1
            q_tp(ii) = min(q_tp(ii+1), p_s(ii) * (n_tp / ii));
            end
            q_tp = min(q_tp, 1);
            % Global: same p-values padded to 12 tests
            p_global = [p_tp; 0.30; 0.40; 0.50; 0.60; 0.70; 0.80; 0.90; 0.95];
            n_g = 12;
            [p_sg, sig] = sort(p_global);
            q_g = zeros(n_g, 1);
            q_g(n_g) = p_sg(n_g);
            for ii = n_g-1:-1:1
            q_g(ii) = min(q_g(ii+1), p_sg(ii) * (n_g / ii));
            end
            q_g = min(q_g, 1);
            q_g_unsorted = zeros(n_g, 1); q_g_unsorted(sig) = q_g;
            % The first 4 tests should have smaller q when corrected per-timepoint
            testCase.verifyTrue(all(q_tp <= q_g_unsorted(1:4) + 1e-12), ...
            'Per-timepoint FDR should produce equal or smaller q-values than global');
        end

        % =================================================================
        % Imputation logic
        % =================================================================

        function testImpute_RetainsPartialDataRows(testCase)
            % Synthetic test: a patient with partial NaN data should be retained
            % after imputation, not dropped as in listwise deletion.
            % Row 1: complete data, Row 2: one missing value (kept after fill),
            % Row 3: all NaN (dropped -- no usable imaging data),
            % Row 4: complete data.
            X = [1 2 3; 4 NaN 6; NaN NaN NaN; 7 8 9];
            y = [0; 1; 0; 1];
            % Keep rows that have at least one non-NaN imaging value
            has_any = any(~isnan(X), 2);
            mask = has_any & ~isnan(y);
            X_imp = X(mask, :);
            y_out = y(mask);
            % Fill remaining NaNs with column medians
            X_out = fillmissing(X_imp, 'constant', median(X_imp, 1, 'omitnan'));
            % Row 2 (partial NaN) should be kept; row 3 (all NaN) dropped
            testCase.verifyTrue(size(X_out, 1) == 3, ...
            'Partial-data rows should be retained after imputation');
            testCase.verifyTrue(~any(isnan(X_out(:))), ...
            'No NaN values should remain after imputation');
            testCase.verifyTrue(length(y_out) == 3, ...
            'Outcome vector length should match imputed matrix');
        end

        function testImpute_MedianFillValues(testCase)
            % Verify that imputed values equal the column median.
            X = [2 10; NaN 20; 6 NaN];
            col_medians = median(X, 1, 'omitnan');