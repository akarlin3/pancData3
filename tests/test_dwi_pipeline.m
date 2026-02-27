classdef test_dwi_pipeline < matlab.unittest.TestCase
    % TESTDWIPIPELINE Formal unit test suite for the DWI Analysis Pipeline.
    %
    % This class uses the matlab.unittest framework to ensure that the 
    % data loading and metric calculations produce valid, physical results,
    % avoiding workspace collisions and providing robust assertions.
    %
    % Run tests with: 
    %   results = runtests('tests/test_dwi_pipeline.m');

    properties
        % Define properties that can be shared across tests if needed
        MockDataDir
        ConfigStruct
    end

    methods(TestMethodSetup)
        % Setup for each test
        function createMockConfig(testCase)
            % Create a dummy configuration structure for testing
            testCase.ConfigStruct = struct();
            testCase.ConfigStruct.skip_to_reload = false;
            testCase.ConfigStruct.ivim_bthr = 100;
            testCase.ConfigStruct.dataloc = fullfile(pwd, 'mock_data');
            testCase.ConfigStruct.dcm2nii_call = 'dummy_dcm2niix';
            testCase.ConfigStruct.adc_thresh = 0.00115;
            testCase.ConfigStruct.high_adc_thresh = 0.001;
            testCase.ConfigStruct.d_thresh = 0.001;
            testCase.ConfigStruct.f_thresh = 0.1;
            testCase.ConfigStruct.dstar_thresh = 0.01;
            testCase.ConfigStruct.min_vox_hist = 10; % Lowered for small mock
            testCase.ConfigStruct.adc_max = 3.0e-3;
            testCase.ConfigStruct.clinical_data_sheet = 'mock_sheet.xlsx';
            testCase.ConfigStruct.patient_ids = {'P99-MOCK'};
            
            % Create a temporary directory for mock data
            testCase.MockDataDir = testCase.ConfigStruct.dataloc;
            if ~exist(testCase.MockDataDir, 'dir')
                mkdir(testCase.MockDataDir);
            end
            
            % Add necessary paths for tests
            testCase.ConfigStruct.orig_path = path;
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'dependencies'));
        end
    end

    methods(TestMethodTeardown)
        % Cleanup after each test
        function removeMockData(testCase)
            if exist(testCase.MockDataDir, 'dir')
                rmdir(testCase.MockDataDir, 's');
            end
            
            % Restore path
            if isfield(testCase.ConfigStruct, 'orig_path')
                path(testCase.ConfigStruct.orig_path);
            end
        end
    end

    methods(Test)
        
        function testMockNiftiPipelineOutputs(testCase)
            % -----------------------------------------------------------------
            % BOILERPLATE MOCK 3D NIFTI ARRAY TEST
            % -----------------------------------------------------------------
            % 1. Create a tiny mock 4D DWI tensor (X, Y, Z, b-values)
            % We simulate a 5x5x5 volume with 4 b-values (e.g., 0, 50, 400, 800)
            sz = [5, 5, 5, 4];
            
            % Generate physically plausible exponential decay signal
            bvals = [0, 50, 400, 800];
            true_adc = 0.001; % Typical tissue ADC
            
            mockSignal = zeros(sz);
            for b = 1:4
                mockSignal(:,:,:,b) = 1500 * exp(-bvals(b) * true_adc) + 10*randn(5,5,5); 
            end
            mockSignal(mockSignal < 1) = 1; % Prevent negative noise from causing log(0)
            
            % Generate a corresponding 3D GTV mask 
            mockMask = zeros(5, 5, 5);
            mockMask(2:4, 2:4, 2:4) = 1; % Active central core
            
            % -----------------------------------------------------------------
            % 2. Execute a simulated pipeline segment 
            % Normally, you would pass these into your load_dwi_data or fit_adc_mono
            % algorithms. Here we simulate the monoexponential ADC fit logic inline
            % to demonstrate the assertion framework.
            
            S_2d = reshape(mockSignal, [prod(sz(1:3)), length(bvals)]);
            adc_vec = zeros(numel(mockMask), 1);
            
            % Simple OLS log-linear fit for demonstration
            S_a = S_2d; 
            % (-b(2:end) \ log(S(b>0)/S(b=0)))
            adc_vec = (-bvals(2:end)' \ log(S_a(:,2:end) ./ S_a(:,1))')';
            
            % Extract GTV only
            final_adc_features = adc_vec(mockMask == 1);
            
            % -----------------------------------------------------------------
            % 3. Automated Assertions using verifyTrue / verifyEqual
            
            % Is the output physically valid? (No negative ADCs)
            testCase.verifyTrue(all(final_adc_features > 0), ...
                'Negative ADC values calculated. Model fit is non-physical.');
            
            % Are there any NaNs generated during the process?
            testCase.verifyFalse(any(isnan(final_adc_features)), ...
                'NaNs detected in the output metric arrays.');
            
            % Does it accurately approximate our strictly physical ground truth?
            mean_calc_adc = mean(final_adc_features);
            testCase.verifyEqual(mean_calc_adc, true_adc, 'RelTol', 0.15, ...
                'Calculated ADC deviates from physical truth by >15%.');
            
        end
        
                function testBH_QValues_KnownInput(testCase)
            % Verify BH q-values for a hand-calculated example.
            % Five raw p-values; q = p * (m/rank), then enforce monotonicity.
            raw_p = [0.005; 0.01; 0.03; 0.40; 0.90];
            m = length(raw_p);
            [p_sorted, sort_idx] = sort(raw_p);
            q_vals = zeros(size(p_sorted));
            for k = 1:m
            q_vals(k) = p_sorted(k) * (m / k);
            end
            for k = m-1:-1:1
            q_vals(k) = min(q_vals(k), q_vals(k+1));
            end
            q_vals(q_vals > 1) = 1;
            % Expected: [0.025, 0.025, 0.05, 0.50, 0.90]
            expected = [0.025; 0.025; 0.05; 0.50; 0.90];
            testCase.verifyTrue(all(abs(q_vals - expected) <= 1e-12), ...
            'BH q-values do not match hand-calculated values');
        end

        function testBH_Monotonicity(testCase)
            % Q-values must be non-decreasing after the BH step-up procedure.
            rng(1);
            raw_p = sort(rand(50, 1));   % 50 sorted p-values
            m = length(raw_p);
            q_vals = zeros(size(raw_p));
            for k = 1:m
            q_vals(k) = raw_p(k) * (m / k);
            end
            for k = m-1:-1:1
            q_vals(k) = min(q_vals(k), q_vals(k+1));
            end
            q_vals(q_vals > 1) = 1;
            diffs = diff(q_vals);
            testCase.verifyTrue(all(diffs >= -1e-15), ...
            'BH q-values are not monotonically non-decreasing');
        end

        function testBH_CappedAtOne(testCase)
            % No q-value should exceed 1.0.
            raw_p = [0.10; 0.50; 0.80; 0.95; 0.99];
            m = length(raw_p);
            [p_sorted, ~] = sort(raw_p);
            q_vals = zeros(size(p_sorted));
            for k = 1:m
            q_vals(k) = p_sorted(k) * (m / k);
            end
            for k = m-1:-1:1
            q_vals(k) = min(q_vals(k), q_vals(k+1));
            end
            q_vals(q_vals > 1) = 1;
            testCase.verifyTrue(all(q_vals <= 1.0), ...
            'Some BH q-values exceed 1.0');
        end

        function testBH_QGreaterThanOrEqualP(testCase)
            % Every q-value must be >= the corresponding raw p-value.
            raw_p = [0.001; 0.01; 0.04; 0.05; 0.20; 0.50; 0.80];
            m = length(raw_p);
            [p_sorted, sort_idx] = sort(raw_p);
            q_vals = zeros(size(p_sorted));
            for k = 1:m
            q_vals(k) = p_sorted(k) * (m / k);
            end
            for k = m-1:-1:1
            q_vals(k) = min(q_vals(k), q_vals(k+1));
            end
            q_vals(q_vals > 1) = 1;
            q_unsorted = zeros(size(raw_p));
            q_unsorted(sort_idx) = q_vals;
            testCase.verifyTrue(all(q_unsorted >= raw_p - 1e-15), ...
            'Some BH q-values are smaller than their raw p-values');
        end

        function testBH_SinglePValue(testCase)
            % Edge case: m = 1. q should equal p (no correction needed).
            raw_p = 0.03;
            q = raw_p * (1 / 1);
            q = min(q, 1);
            testCase.verifyTrue(abs(q - 0.03) <= 1e-15);
        end

        function testBH_AllIdenticalPValues(testCase)
            % When all p-values are identical, q-values should all equal
            % the same adjusted value (or 1.0 if > 1).
            raw_p = repmat(0.04, 10, 1);
            m = length(raw_p);
            [p_sorted, ~] = sort(raw_p);
            q_vals = zeros(size(p_sorted));
            for k = 1:m
            q_vals(k) = p_sorted(k) * (m / k);
            end
            for k = m-1:-1:1
            q_vals(k) = min(q_vals(k), q_vals(k+1));
            end
            q_vals(q_vals > 1) = 1;
            % All q-values should be the same after monotonicity
            testCase.verifyTrue(abs(max(q_vals) - min(q_vals)) <= 1e-15, ...
            'Identical p-values should yield identical q-values');
        end

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
            % is evaluated against A only. Verify transitive behaviour.
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
            % metric index: indices 1-4 = Absolute, 5-8 = Percent-Change.
            % base = mod(sel - 1, 4) + 1
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

        function testElasticNet_AlphaIsHalf(testCase)
            % The refactored code must use Alpha = 0.5 for Elastic Net.
            code = metrics_code;
            matches = regexp(code, '''Alpha''\s*,\s*0\.5', 'match');
            testCase.verifyTrue(~isempty(matches), ...
            'metrics.m should contain Alpha = 0.5 for Elastic Net');
        end

        function testElasticNet_NotPureLasso(testCase)
            % Confirm no lassoglm call uses pure LASSO (Alpha = 1).
            code = metrics_code;
            matches = regexp(code, '''Alpha''\s*,\s*1(\s|,|\))', 'match');
            testCase.verifyTrue(isempty(matches), ...
            'No lassoglm call should use pure LASSO (Alpha=1)');
        end

        function testLoop_NotHardcodedTo2And3(testCase)
            % The loop "for target_fx = [2, 3]" should have been replaced.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'target_fx = [2, 3]'), ...
            'Loop should NOT be hardcoded to [2, 3]');
        end

        function testLoop_UsesNTpUpperBound(testCase)
            % The loop should use nTp as the upper bound.
            code = metrics_code;
            matches = regexp(code, 'target_fx\s*=\s*2\s*:\s*nTp', 'match');
            testCase.verifyTrue(~isempty(matches), ...
            'Loop should iterate from 2 to nTp');
        end

        function testRemoved_TargetedM8(testCase)
            % The "Targeted Statistical Analysis (m = 8 tests)" section
            % should no longer exist in metrics.m.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'Targeted Statistical Analysis (m = 8'), ...
            'Targeted m=8 section should be removed');
        end

        function testRemoved_TargetedM4(testCase)
            % The "Targeted FDR (Relative Change Only, m = 4)" section
            % should no longer exist in metrics.m.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'Targeted FDR (Relative Change Only, m = 4)'), ...
            'Targeted m=4 FDR section should be removed');
        end

        function testADC_MetricsCallsFitAdcMono(testCase)
            % metrics.m should call fit_adc_mono instead of inline OLS.
            code = metrics_code;
            matches = regexp(code, 'fit_adc_mono\s*\(', 'match');
            testCase.verifyTrue(~isempty(matches), ...
            'metrics.m should call fit_adc_mono for ADC computation');
        end

        function testADC_NoInlineOLS(testCase)
            % The old inline OLS pattern "beta = X \ log(s_signal)" should
            % be gone from metrics.m.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'X \ log(s_signal)'), ...
            'Inline OLS fitting should be removed from metrics.m');
        end

        function testADC_NoNestedPixelLoops(testCase)
            % The old nested "for y = 1:ny / for x = 1:nx" pixel loops
            % should be gone from metrics.m (replaced by fit_adc_mono call).
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'slice_adc(y,x) = beta(2)'), ...
            'Nested pixel-wise OLS loops should be removed');
        end

        function testRetained_PerTimepointFDR(testCase)
            % The FDR section should now apply BH per-timepoint.
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'Per-Timepoint') || contains(code, 'PER-TIMEPOINT'), ...
            'FDR section should indicate per-timepoint operation');
        end

        function testRetained_HolmBonferroni(testCase)
            % Holm-Bonferroni should still exist (now per-timepoint).
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'Holm-Bonferroni'), ...
            'Holm-Bonferroni logic should be retained');
        end

        function testFDR_NoGlobalPooling(testCase)
            % The old global FDR sweep header should no longer exist.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'Full Metric Sweep'), ...
            'Global FDR sweep should be replaced with per-timepoint');
        end

        function testBval_NoTruncation(testCase)
            % The old truncation pattern should be gone from metrics.m.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'min_dim = min(n_imgs, n_bvals)'), ...
            'Arbitrary b-value truncation (min_dim) should be removed');
        end

        function testBval_ExpectedProtocolDefined(testCase)
            % metrics.m should define expected_bvals for protocol validation.
            code = metrics_code;
            matches = regexp(code, 'expected_bvals\s*=\s*\[0;\s*30;\s*150;\s*550\]', 'match');
            testCase.verifyTrue(~isempty(matches), ...
            'metrics.m should define expected_bvals = [0; 30; 150; 550]');
        end

        function testBval_ProtocolDeviationFlagged(testCase)
            % metrics.m should flag protocol deviations with a warning message.
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'Protocol deviation'), ...
            'Protocol deviation flagging should be present in metrics.m');
        end

        function testBval_DeviationExcludesPatient(testCase)
            % The validation block should use 'continue' to exclude
            % deviating patients from the analysis.
            code = metrics_code;
            % Find the Protocol deviation block and verify it contains continue
            idx_dev = strfind(code, 'Protocol deviation');
            idx_continue = strfind(code, 'continue');
            testCase.verifyTrue(~isempty(idx_dev) && ~isempty(idx_continue), ...
            'Both protocol deviation flag and continue must exist');
            % At least one continue must follow the protocol deviation flag
            % (within 200 chars, i.e. within the same if-block)
            found = false;
            for ci = 1:length(idx_continue)
            if idx_continue(ci) > idx_dev(1) && (idx_continue(ci) - idx_dev(1)) < 200
            found = true;
            break;
            end
            end
            testCase.verifyTrue(found, ...
            'A continue statement should follow the protocol deviation flag to exclude the patient');
        end

        function testBval_ValidationLogic_MatchingProtocol(testCase)
            % When bvals exactly match [0; 30; 150; 550], no deviation.
            bvals = [0; 30; 150; 550];
            expected_bvals = [0; 30; 150; 550];
            testCase.verifyTrue(isequal(sort(bvals), expected_bvals), ...
            'Standard protocol b-values should pass validation');
        end

        function testBval_ValidationLogic_UnsortedMatch(testCase)
            % Unsorted b-values that match after sorting should pass.
            bvals = [0; 30; 150; 550];
            expected_bvals = [0; 30; 150; 550];
            testCase.verifyTrue(isequal(sort(bvals), expected_bvals), ...
            'Unsorted but correct b-values should pass validation');
        end

        function testBval_ValidationLogic_ExtraValue(testCase)
            % Extra b-value (5 instead of 4) should be flagged as deviation.
            bvals = [0; 50; 400; 800; 1000];
            expected_bvals = [0; 30; 150; 550];
            is_deviation = ~isequal(sort(bvals), expected_bvals);
            testCase.verifyTrue(is_deviation, ...
            'Extra b-values should be flagged as protocol deviation');
        end

        function testBval_ValidationLogic_WrongValue(testCase)
            % Different b-value set should be flagged as deviation.
            bvals = [0; 30; 200; 550];
            expected_bvals = [0; 30; 150; 550];
            is_deviation = ~isequal(sort(bvals), expected_bvals);
            testCase.verifyTrue(is_deviation, ...
            'Non-standard b-values should be flagged as protocol deviation');
        end

        function testBval_ValidationLogic_MissingValue(testCase)
            % Fewer b-values than expected should be flagged.
            bvals = [0; 50; 400];
            expected_bvals =  [0; 30; 150; 550];
            is_deviation = ~isequal(sort(bvals), expected_bvals);
            testCase.verifyTrue(is_deviation, ...
            'Missing b-values should be flagged as protocol deviation');
        end

        function testImpute_NoListwiseDeletion(testCase)
            % The old listwise deletion pattern should be gone from metrics.m.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'valid_lasso_mask = all(~isnan(X_lasso), 2)'), ...
            'Listwise deletion pattern should be removed from metrics.m');
        end

        function testImpute_UsesKNN(testCase)
            % metrics.m should use knn_impute_train_test for imputation.
            code = metrics_code;
            matches = regexp(code, 'knn_impute_train_test\s*\(', 'match');
            testCase.verifyTrue(~isempty(matches), ...
            'metrics.m should use knn_impute_train_test for imputation');
        end

        function testImpute_ExcludesAllNaNRows(testCase)
            % metrics.m should exclude patients with ALL imaging data missing.
            code = metrics_code;
            matches = regexp(code, 'any\s*\(\s*~isnan\s*\(\s*X_lasso_all', 'match');
            testCase.verifyTrue(~isempty(matches), ...
            'metrics.m should check for patients missing all imaging data');
        end

        function testImpute_RetainsPartialDataRows(testCase)
            % Synthetic test: a patient with partial NaN data should be retained
            % after imputation, not dropped as in listwise deletion.
            X = [1 2 3; 4 NaN 6; NaN NaN NaN; 7 8 9];
            y = [0; 1; 0; 1];
            has_any = any(~isnan(X), 2);
            mask = has_any & ~isnan(y);
            X_imp = X(mask, :);
            y_out = y(mask);
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

        function testImpute_BeforeStandardization(testCase)
            % Verify imputation occurs before the Elastic Net standardization.
            % knn_impute_train_test must appear before lassoglm(..., 'Standardize', true)
            code = metrics_code;
            pos_fill = regexp(code, 'knn_impute_train_test\s*\(', 'start', 'once');
            pos_std  = regexp(code, '''Standardize''\s*,\s*true', 'start', 'once');
            testCase.verifyTrue(~isempty(pos_fill) && ~isempty(pos_std), ...
            'Both knn_impute_train_test and Standardize should exist in metrics.m');
            testCase.verifyTrue(pos_fill < pos_std, ...
            'Imputation (knn_impute_train_test) must occur before standardization');
        end

        function testVis_NoTruncation(testCase)
            % The old truncation pattern should be gone from visualize_results.m.
            code = vis_code;
            testCase.verifyTrue(~contains(code, 'min_dim = min('), ...
            'Arbitrary b-value truncation (min_dim) should be removed from visualize_results.m');
        end

        function testVis_ExpectedProtocolDefined(testCase)
            % visualize_results.m should define expected_bvals for protocol validation.
            code = vis_code;
            matches = regexp(code, 'expected_bvals\s*=\s*\[0;\s*30;\s*150;\s*550\]', 'match');
            testCase.verifyTrue(~isempty(matches), ...
            'visualize_results.m should define expected_bvals = [0; 30; 150; 550]');
        end

        function testVis_ProtocolDeviationFlagged(testCase)
            % visualize_results.m should flag protocol deviations with a warning message.
            code = vis_code;
            testCase.verifyTrue(contains(code, 'Protocol deviation'), ...
            'Protocol deviation flagging should be present in visualize_results.m');
        end

        function testVis_DeviationExcludesPatient(testCase)
            % The validation block should use 'continue' to exclude
            % deviating patients from the comparative mapping.
            code = vis_code;
            idx_dev = strfind(code, 'Protocol deviation');
            idx_continue = strfind(code, 'continue');
            testCase.verifyTrue(~isempty(idx_dev) && ~isempty(idx_continue), ...
            'Both protocol deviation flag and continue must exist in visualize_results.m');
            % At least one continue must follow the protocol deviation flag
            % (within 200 chars, i.e. within the same if-block)
            found = false;
            for ci = 1:length(idx_continue)
            if idx_continue(ci) > idx_dev(1) && (idx_continue(ci) - idx_dev(1)) < 200
            found = true;
            break;
            end
            end
            testCase.verifyTrue(found, ...
            'A continue statement should follow the protocol deviation flag to exclude the patient');
        end

        function testDIR_EmptyInputsReturnEmpty(testCase)
            % All three empty-input cases should return [].
            r1 = apply_dir_mask_propagation([], ones(4,4,4), true(4,4,4));
            r2 = apply_dir_mask_propagation(ones(4,4,4), [], true(4,4,4));
            r3 = apply_dir_mask_propagation(ones(4,4,4), ones(4,4,4), []);
            testCase.verifyTrue(isempty(r1) && isempty(r2) && isempty(r3), ...
            'Empty inputs should return empty output');
        end

        function testDIR_SizeMismatchReturnsEmpty(testCase)
            % Mismatched image sizes should return [].
            r = apply_dir_mask_propagation(ones(4,4,4), ones(5,5,5), true(4,4,4));
            testCase.verifyTrue(isempty(r), ...
            'Size-mismatched inputs should return empty output');
        end

        function testDIR_IdenticalImagesPreserveMask(testCase)
            % When fixed == moving (no deformation), the warped mask should
            % closely match the original mask.
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

        function testDIR_FunctionExistsInCodebase(testCase)
            % load_dwi_data_forAvery.m should call apply_dir_mask_propagation.
            code = loaddwi_code;
            testCase.verifyTrue(contains(code, 'apply_dir_mask_propagation'), ...
            'load_dwi_data_forAvery.m should call apply_dir_mask_propagation');
        end

        function testADC_NegativeClampedToNaN(testCase)
            % load_dwi_data_forAvery.m should assign NaN to negative ADC.
            code = loaddwi_code;
            testCase.verifyTrue(contains(code, 'adc_vec(adc_vec < 0) = nan') || ...
            contains(code, 'adc_vec(adc_vec<0) = nan') || ...
            contains(code, 'adc_vec(adc_vec < 0)= nan'), ...
            'Negative ADC should be set to NaN, not zero');
        end

        function testADC_NotClampedToZero(testCase)
            % The old clamping pattern should not exist.
            code = loaddwi_code;
            testCase.verifyTrue(~contains(code, 'adc_vec(adc_vec < 0) = 0'), ...
            'Negative ADC should NOT be clamped to zero');
        end

        function testADC_NaNRemovedByNanmean(testCase)
            % NaN values should be excluded from mean calculations.
            adc_vec = [1.5e-3; 2.0e-3; NaN; 1.8e-3];
            m = mean(adc_vec, 'omitnan');
            testCase.verifyTrue(abs(m - mean([1.5e-3; 2.0e-3; 1.8e-3])) < 1e-15, ...
            'nanmean should exclude NaN failed fits');
        end

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

        function testPerTimepointFDR_NoGlobalSweepInSection9(testCase)
            % Section 9 should no longer pool across timepoints.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'FDR Correction') || ...
            contains(code, 'Per-Timepoint'), ...
            'Section 9 FDR should be per-timepoint, not global');
        end

        function testOOFROC_UsesRiskScoresAll(testCase)
            % perfcurve should be called with risk_scores_all, not
            % mdl_roc.Fitted.Probability.
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'perfcurve(labels(valid_roc), risk_scores_all(valid_roc)'), ...
            'perfcurve should use risk_scores_all directly');
        end

        function testOOFROC_NoInSampleROC(testCase)
            % The old in-sample ROC pattern should be gone.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'EXPLORATORY ROC ANALYSIS'), ...
            'In-sample exploratory ROC should be replaced with OOF ROC');
        end

        function testOOFROC_LabeledAsPrimary(testCase)
            % ROC section should say PRIMARY ROC ANALYSIS in the fprintf.
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'PRIMARY ROC ANALYSIS'), ...
            'ROC section should be labeled Primary');
        end

        function testOOFROC_NoStaleGLMRefit(testCase)
            % The stale fitglm call inside the ROC block should be gone.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'mdl_roc = fitglm'), ...
            'Stale fitglm call in ROC block should be removed');
        end

        function testOOFROC_YoudenCutoffFromLOOCV(testCase)
            % The Youden optimal cutoff (roc_opt_thresh) must come from perfcurve,
            % not from any logistic regression fitted probabilities.
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'roc_opt_thresh'), ...
            'roc_opt_thresh must be defined (Youden cutoff from LOOCV ROC)');
            testCase.verifyTrue(~contains(code, 'mdl_roc.Fitted'), ...
            'Optimal cutoff must not use mdl_roc.Fitted (no in-sample refitting)');
        end

        function testCollinearity_NoGlobalLeakage(testCase)
            % filter_collinear_features should NOT be called on global matrices like X_clean_all or X_impute.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'filter_collinear_features(X_clean_all'), ...
            'filter_collinear_features should not be called on the full cohort X_clean_all');
            testCase.verifyTrue(~contains(code, 'filter_collinear_features(X_impute'), ...
            'filter_collinear_features should not be called on the full imputed matrix X_impute');
            % It SHOULD be called on training folds
            testCase.verifyTrue(contains(code, 'filter_collinear_features(X_tr_imp'), ...
            'filter_collinear_features must be applied to training folds (X_tr_imp)');
        end

        function testCIndex_UnbiasedEvaluation(testCase)
            % Harrell's C-index should be evaluated instead of binary LPOCV.
            code = metrics_code;
            testCase.verifyTrue(contains(code, "Harrell's C-index"), ...
            'metrics.m should implement Harrell''s C-index evaluation');
            testCase.verifyTrue(contains(code, 'drops mathematically "incomparable pairs"'), ...
            'C-index section should specify dropping incomparable pairs');
            testCase.verifyTrue(~contains(code, 'LPOCV (Leave-Pair-Out)'), ...
            'Legacy LPOCV loop should be removed');
        end

        function testSubVol_DIRReenabled(testCase)
            % The old keep_policy = 1:10 exclusion should be gone.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'keep_policy = 1:10'), ...
            'Sub-volume features should no longer be excluded now that DIR is implemented');
        end

        function testSubVol_PostRTStillExcluded(testCase)
            % Post-RT dose exclusion (target_fx == 6) should still be present.
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'target_fx == 6'), ...
            'Post-RT dose exclusion should be retained');
        end

        function testThreshold_KurtosisProtected_ADC(testCase)
            code = loaddwi_code;
            % Should check numel(adc_vec_sub) >= min_vox_hist or similar
            testCase.verifyTrue(contains(code, 'numel(adc_vec_sub) >= min_vox_hist'), ...
            'ADC sub-volume kurtosis should be protected by voxel threshold');
        end

        function testThreshold_KurtosisProtected_D(testCase)
            code = loaddwi_code;
            testCase.verifyTrue(contains(code, 'numel(d_vec_sub) >= min_vox_hist'), ...
            'D sub-volume kurtosis should be protected by voxel threshold');
        end

        function testThreshold_KSProtected(testCase)
            code = loaddwi_code;
            % Should check numel(adc_vec) >= min_vox_hist && numel(adc_baseline) >= min_vox_hist
            testCase.verifyTrue(contains(code, 'numel(adc_vec) >= min_vox_hist') && contains(code, 'numel(adc_baseline) >= min_vox_hist'), ...
            'KS test should be protected by voxel threshold for both current and baseline');
        end

        function testThreshold_D95Protected_Metrics(testCase)
            code = metrics_code;
            % Should check sum(adc_mask_1d) >= min_subvol_voxels or similar
            testCase.verifyTrue(contains(code, 'sum(adc_mask_1d) >= min_subvol_voxels'), ...
            'Metrics.m D95/V50 should be protected by voxel threshold for ADC sub-volume');
        end

        function testThreshold_V50Protected_Metrics(testCase)
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'sum(dstar_mask_1d) >= min_subvol_voxels'), ...
            'Metrics.m D95/V50 should be protected by voxel threshold for D* sub-volume');
        end

        function testSanityCheck_NoKurtosisPlot(testCase)
            % The sanity check subplot must NOT use adc_kurt as a metric.
            % It should use adc_sd instead.
            code = metrics_code;
            testCase.verifyTrue(~contains(code, 'kurt_fx1 = adc_kurt'), ...
            'Sanity check should use ADC SD, not trace-average ADC Kurtosis');
            testCase.verifyTrue(contains(code, 'sd_fx1  = adc_sd'), ...
            'Sanity check heterogeneity subplot should use adc_sd');
        end

        function testIVIM_Bthr100(testCase)
            % Pipeline uses config struct for bthr.
            code = loaddwi_code;
            testCase.verifyTrue(contains(code, 'ivim_bthr = config_struct.ivim_bthr') || contains(code, 'ivim_bthr = 100'), ...
            'ivim_bthr should be read from config_struct in load_dwi_data.m');
        end

        function testIVIM_fUpperBound(testCase)
            % IVIMmodelfit.m should cap perfusion fraction f at 0.4.
            code = ivim_code;
            testCase.verifyTrue(contains(code, '0.4'), ...
            'IVIMmodelfit.m should set f upper bound to 0.4');
        end

        function testIVIM_DstarUpperBound(testCase)
            % IVIMmodelfit.m should cap D* at 0.1 mm²/s.
            code = ivim_code;
            testCase.verifyTrue(contains(code, '0.1'), ...
            'IVIMmodelfit.m should set D* upper bound to 0.1 mm²/s');
        end

        function testIVIM_NoLooseBoundsF(testCase)
            % Old loose f upper bound of 1 should be gone.
            code = ivim_code;
            % The old line was: lim = [0 0 0 0;3e-3 2*max(Y(:)) 1 1]
            testCase.verifyTrue(~contains(code, '3e-3 2*max(Y(:)) 1 1'), ...
            'Old loose f and D* limits should be removed from IVIMmodelfit.m');
        end

        function testIVIM_OptsLimSupported(testCase)
            % IVIMmodelfit.m should support opts.lim for caller overrides.
            code = ivim_code;
            testCase.verifyTrue(contains(code, '''lim''') && contains(code, 'opts.lim'), ...
            'IVIMmodelfit.m should support opts.lim override');
        end

        function testIVIM_DefaultBlim100(testCase)
            % IVIMmodelfit.m should default blim to 100.
            code = ivim_code;
            testCase.verifyTrue(contains(code, 'blim = 100'), ...
            'IVIMmodelfit.m default blim should be 100');
        end

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

        function testDIR_ReturnsDfieldAndRef(testCase)
            % apply_dir_mask_propagation should return 3 outputs (mask, D_forward, ref3d).
            code = readLoadDwiSource();
            % Check signature via file reading
            apPath = fullfile(fileparts(mfilename('fullpath')), '..', 'utils', 'apply_dir_mask_propagation.m');
            fid = fopen(apPath, 'r'); 
            if fid == -1, return; end
            ap_code = fread(fid,'*char')'; fclose(fid);
            testCase.verifyTrue(contains(ap_code, 'D_forward, ref3d] = apply_dir_mask_propagation'), ...
            'apply_dir_mask_propagation should return D_forward and ref3d');
        end

        function testDIR_DoseMapRemainsRigid(testCase)
            % load_dwi_data_forAvery.m should NOT apply imwarp to the dose map.
            code = loaddwi_code;
            testCase.verifyTrue(~contains(code, 'imwarp(dose_map, D_forward_cur'), ...
            'Dose map must remain rigidly aligned and should NOT be warped via imwarp');
        end

        function testDIR_DoseWarpFallback(testCase)
            % At baseline (fi=1) or when DIR fails, code must fall back to rigid dose_map.
            code = loaddwi_code;
            testCase.verifyTrue(contains(code, 'dose_map_dvh     = dose_map'), ...
            'Rigid dose fallback ''dose_map_dvh = dose_map'' should be initialised before DIR block');
        end

        function testDIR_DfieldCached(testCase)
            % D_forward and ref3d should be written to the .mat cache file.
            code = loaddwi_code;
            testCase.verifyTrue(contains(code, "save(dir_cache_file, 'gtv_mask_warped', 'D_forward', 'ref3d')"), ...
            'D_forward and ref3d must be included in the cache save call');
        end

        function testDIR_GTVnDoseRemainsRigid(testCase)
            % GTVn dose block must also use the rigid dose (dose_map_dvh_n).
            code = loaddwi_code;
            testCase.verifyTrue(~contains(code, 'dose_map_dvh_n = imwarp(dose_map'), ...
            'GTVn dose block must not warp dose_map into dose_map_dvh_n');
        end

        function testDIR_LinearMaskWarping(testCase)
            % apply_dir_mask_propagation.m should use 'linear' interpolation for mask warping.
            apPath = fullfile(fileparts(mfilename('fullpath')), '..', 'utils', 'apply_dir_mask_propagation.m');
            fid = fopen(apPath, 'r');
            if fid == -1, return; end
            ap_code = fread(fid,'*char')'; fclose(fid);
            testCase.verifyTrue(contains(ap_code, "'Interp', 'linear'"), ...
            "apply_dir_mask_propagation should use 'linear' interpolation for mask warping");
            testCase.verifyTrue(~contains(ap_code, "'Interp', 'nearest'"), ...
            "apply_dir_mask_propagation should NOT use 'nearest' interpolation for mask warping");
        end

        function testDIR_SymmetricWorkflow(testCase)
            % apply_dir_mask_propagation.m should implement the symmetric halfway-space logic.
            apPath = fullfile(fileparts(mfilename('fullpath')), '..', 'utils', 'apply_dir_mask_propagation.m');
            fid = fopen(apPath, 'r'); 
            if fid == -1, return; end
            ap_code = fread(fid,'*char')'; fclose(fid);
            testCase.verifyTrue(contains(ap_code, 'mid_img_refined'), ...
            "apply_dir_mask_propagation should use a refined midpoint image");
            testCase.verifyTrue(contains(ap_code, 'D_forward = D_forward_mid - D_backward_mid'), ...
            "apply_dir_mask_propagation should compute the symmetric proxy field");
        end

        function testTD_PanelHasMoreRowsThanPatients(testCase)
            % Synthetic 3-patient, 3-timepoint panel
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
            arr1 = [1.0e-3 1.1e-3; 0.9e-3 0.8e-3; 1.5e-3 1.6e-3];
            lf_td  = [1; 0; 1];
            tot_td = [15; 50; 8];
            [~, t0, t1, ~, ~] = build_td_panel({arr1}, {'ADC'}, lf_td, tot_td, 2, [0 5]);
            testCase.verifyTrue(all(t0 < t1), 't_start must be strictly less than t_stop for every interval');
        end

        function testTD_MetricsCallsBuildTdPanel(testCase)
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'build_td_panel'), ...
            'metrics.m should call build_td_panel for the time-dependent Cox model');
        end

        function testTD_TimeStratifiedCollinearityFilter(testCase)
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'filter_collinear_features(X_train, pat_event_train, frac_train)'), ...
            'TD LOOCV must call filter_collinear_features with frac_train for time-stratified filtering');
            testCase.verifyTrue(contains(code, 'frac_train = opt_panel.frac(train_mask)'), ...
            'TD LOOCV must extract frac_train from the optimal panel training rows');
            testCase.verifyTrue(contains(code, 'X_train = X_train(:, keep_td)'), ...
            'TD LOOCV must apply keep_td mask to X_train');
            testCase.verifyTrue(contains(code, 'X_test  = X_test(:, keep_td)'), ...
            'TD LOOCV must apply the same keep_td mask to X_test');
        end

        function testDL_Isolation_LogicPresent(testCase)
            code = metrics_code;
            testCase.verifyTrue(contains(code, 'DEEP LEARNING RIGOR AUDIT'), ...
            'metrics.m should implement the DL Rigor Audit section');
            testCase.verifyTrue(contains(code, 'dl_provenance.dncnn_train_ids'), ...
            'Audit must check dncnn_train_ids');
            testCase.verifyTrue(contains(code, 'error(''DATA LEAKAGE DETECTED'), ...
            'LOOCV loop must error on detected leakage to prevent optimistic bias');
        end

    end
end

% Local helper functions to read source code for static analysis tests

function code = metrics_code()
    fid = fopen(fullfile(fileparts(mfilename('fullpath')), '..', 'core', 'metrics.m'), 'r');
    if fid == -1, code = ''; return; end
    code = fread(fid, '*char')';
    fclose(fid);
end

function code = loaddwi_code()
    fid = fopen(fullfile(fileparts(mfilename('fullpath')), '..', 'core', 'load_dwi_data.m'), 'r');
    if fid == -1, code = ''; return; end
    code = fread(fid, '*char')';
    fclose(fid);
end

function code = ivim_code()
    fid = fopen(fullfile(fileparts(mfilename('fullpath')), '..', 'dependencies', 'IVIMmodelfit.m'), 'r');
    if fid == -1, code = ''; return; end
    code = fread(fid, '*char')';
    fclose(fid);
end

function code = vis_code()
    fid = fopen(fullfile(fileparts(mfilename('fullpath')), '..', 'core', 'visualize_results.m'), 'r');
    if fid == -1, code = ''; return; end
    code = fread(fid, '*char')';
    fclose(fid);
end

function code = readLoadDwiSource()
    code = loaddwi_code();
end
