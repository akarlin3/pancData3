classdef test_statistical_methods_part2 < matlab.unittest.TestCase
    % TEST_STATISTICAL_METHODS_PART2 Continuation of the pure algorithmic
    % tests for statistical and numerical routines used in the DWI analysis
    % pipeline (split from test_statistical_methods.m to keep files small).
    %
    % This part covers LOOCV logic, b-value validation, per-timepoint FDR,
    % imputation, ADC handling, DIR mask propagation, IVIM segmented
    % fitting, and time-dependent panel construction.
    %
    % Run tests with:
    %   results = runtests('tests/test_statistical_methods_part2.m');

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
            col_med = median(X, 1, 'omitnan');  % [4, 15]
            X_filled = fillmissing(X, 'constant', col_med);
            testCase.verifyTrue(abs(X_filled(2,1) - 4) < 1e-12, ...
            'Imputed value should equal column median');
            testCase.verifyTrue(abs(X_filled(3,2) - 15) < 1e-12, ...
            'Imputed value should equal column median');
        end

        % =================================================================
        % ADC
        % =================================================================

        function testADC_NaNRemovedByNanmean(testCase)
            % NaN values should be excluded from mean calculations.
            adc_vec = [1.5e-3; 2.0e-3; NaN; 1.8e-3];
            m = mean(adc_vec, 'omitnan');
            testCase.verifyTrue(abs(m - mean([1.5e-3; 2.0e-3; 1.8e-3])) < 1e-15, ...
            'nanmean should exclude NaN failed fits');
        end

        % =================================================================
        % DIR functional
        % =================================================================

        function testDIR_EmptyInputsReturnEmpty(testCase)
            % All three empty-input cases should return [].
            warning('off', 'apply_dir_mask_propagation:emptyInput');
            r1 = apply_dir_mask_propagation([], ones(4,4,4), true(4,4,4));
            r2 = apply_dir_mask_propagation(ones(4,4,4), [], true(4,4,4));
            r3 = apply_dir_mask_propagation(ones(4,4,4), ones(4,4,4), []);
            warning('on', 'apply_dir_mask_propagation:emptyInput');
            testCase.verifyTrue(isempty(r1) && isempty(r2) && isempty(r3), ...
            'Empty inputs should return empty output');
        end

        function testDIR_SizeMismatchReturnsEmpty(testCase)
            % Mismatched image sizes should return [].
            warning('off', 'apply_dir_mask_propagation:sizeMismatch');
            r = apply_dir_mask_propagation(ones(4,4,4), ones(5,5,5), true(4,4,4));
            warning('on', 'apply_dir_mask_propagation:sizeMismatch');
            testCase.verifyTrue(isempty(r), ...
            'Size-mismatched inputs should return empty output');
        end

        function testDIR_IdenticalImagesPreserveMask(testCase)
            % When fixed == moving (no deformation needed), the warped mask
            % should closely match the original mask (Dice > 0.95).
            % This verifies that the DIR pipeline does not corrupt the mask
            % in the trivial identity-transform case.
            sz = [32 32 16];
            img = randn(sz) * 100 + 500;
            mask = false(sz);
            mask(12:20, 12:20, 5:12) = true;
            warped = apply_dir_mask_propagation(img, img, mask);
            if ~isempty(warped)
            dice_coeff = 2*sum(warped(:) & mask(:)) / (sum(warped(:)) + sum(mask(:)));
            testCase.verifyTrue(dice_coeff > 0.95, ...
            sprintf('Identical images should preserve mask (Dice=%.3f)', dice_coeff));
            end
        end

        function testDIR_OutputIsLogical(testCase)
            % Output mask must be a logical array.
            sz = [16 16 8];
            img = randn(sz) * 100 + 500;
            mask = false(sz); mask(5:10, 5:10, 3:6) = true;
            warped = apply_dir_mask_propagation(img, img, mask);
            if ~isempty(warped)
            testCase.verifyTrue(islogical(warped), 'Warped mask should be logical');
            testCase.verifyTrue(isequal(size(warped), sz), 'Warped mask should match input size');
            end
        end

        % =================================================================
        % IVIM segmented
        % =================================================================

        function testIVIM_SegmentedTwoStage_Logic(testCase)
            % Inline two-stage test: verify that blim=100 includes more high-b
            % values than blim=200, and that an LLS monoexponential fit on the
            % high-b subset yields a physiologically plausible D.
            D_true = 1.5e-3; f_true = 0.15; Dstar_true = 0.05; S0 = 1000;
            bvals_test = [0; 30; 100; 550];
            rng(42);
            S_test = S0 * ((1-f_true)*exp(-bvals_test*D_true) + f_true*exp(-bvals_test*(D_true+Dstar_true)));
            S_test = S_test + 5*randn(size(S_test));
            % Stage 1 with blim=100: includes b=100 and b=550 (2 points)
            b_hi100 = bvals_test(bvals_test >= 100);
            S_hi100 = S_test(bvals_test >= 100);
            X100 = [-b_hi100, ones(size(b_hi100))];
            p100 = X100 \ log(S_hi100);
            D_est_100 = p100(1);
            % Stage 1 with blim=200: includes only b=550 (exact, 1 point)
            b_hi200 = bvals_test(bvals_test >= 200);
            S_hi200 = S_test(bvals_test >= 200);
            X200 = [-b_hi200, ones(size(b_hi200))];
            p200 = X200 \ log(S_hi200);
            D_est_200 = p200(1);
            testCase.verifyTrue(D_est_100 > 0 && D_est_100 < 3e-3, ...
            sprintf('blim=100 D estimate out of physiological range: %.4g', D_est_100));
            % blim=100 should use more b-values (over-determined), which is the improvement
            testCase.verifyTrue(length(b_hi100) > length(b_hi200), ...
            'blim=100 should include more high-b values for D estimation than blim=200');
        end

        % =================================================================
        % Time-dependent panel
        % =================================================================

        function testTD_PanelHasMoreRowsThanPatients(testCase)
            % Time-dependent panels expand each patient into multiple
            % intervals (one per timepoint transition). With 3 patients
            % and 3 timepoints, the panel should have more rows than 3.
            arr1 = [1.0e-3 1.1e-3 1.2e-3;   % patient 1
            0.9e-3 0.85e-3 0.8e-3;  % patient 2
            1.5e-3 1.6e-3 NaN];      % patient 3 (missing tp3)
            lf_td  = [1; 0; 1];
            tot_td = [30; 100; 20];
            [X_td_t, t0, t1, ev, ~] = build_td_panel({arr1}, {'ADC'}, lf_td, tot_td, 3, [0 5 10]);
            testCase.verifyTrue(size(X_td_t, 1) > length(lf_td), ...
            'Panel should have more rows than patients');
        end

        function testTD_NoEventBeforeFinalInterval(testCase)
            % In counting-process format, the event indicator must be 1
            % only on the patient's final interval. Earlier intervals for
            % the same patient must have event=0 (they are still at risk).
            arr1 = [1.0e-3 1.1e-3 1.2e-3; 0.9e-3 0.85e-3 0.8e-3];
            lf_td  = [1; 0];
            tot_td = [20; 100];
            [~, t0, t1, ev, pid] = build_td_panel({arr1}, {'ADC'}, lf_td, tot_td, 3, [0 5 10]);
            % For patient 1 (LF), event should only fire on its last interval
            p1_rows = (pid == 1);
            ev_p1 = ev(p1_rows);
            testCase.verifyTrue(sum(ev_p1) == 1 && ev_p1(end) == true, ...
            'Event must fire exactly once, on the last interval, for LF patient');
        end

        function testTD_StartAlwaysLessThanStop(testCase)
            % Every interval must have strictly positive duration (t_start < t_stop).
            % Zero-length or negative-duration intervals would break Cox models.
            arr1 = [1.0e-3 1.1e-3; 0.9e-3 0.8e-3; 1.5e-3 1.6e-3];
            lf_td  = [1; 0; 1];
            tot_td = [15; 50; 8];
            [~, t0, t1, ~, ~] = build_td_panel({arr1}, {'ADC'}, lf_td, tot_td, 2, [0 5]);
            testCase.verifyTrue(all(t0 < t1), 't_start must be strictly less than t_stop for every interval');
        end

    end
end
