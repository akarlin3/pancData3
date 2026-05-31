classdef test_optimize_adc_threshold_nested_cv < matlab.unittest.TestCase
    % TEST_OPTIMIZE_ADC_THRESHOLD_NESTED_CV
    %
    % Correctness of the cross-validated Tactic-3 correction (CP1/CP3):
    %   (2) REMOVES OPTIMISM: on NULL data the out-of-fold association is
    %       LESS optimistic than the in-sample min-p selection (mean p shifts
    %       toward 0.5; fewer false "significant" hits).
    %   (3) RECOVERS REAL SIGNAL: on data with a TRUE threshold effect the
    %       out-of-fold association is significant and the recommended
    %       threshold lands in a plausible window around the truth.
    %   (4) LEAKAGE: perturbing a held-out patient must not change its own
    %       fold's selected threshold (same discipline as
    %       test_knn_temporal_leakage).
    %   (5) Inestimable folds NaN-propagate and are reported, not dropped.
    %
    % Margins are deliberately generous: this container's local execution is
    % not trusted, so the assertions test ROBUST directional properties that
    % CI's real MATLAB will confirm, not finely-tuned magnitudes.

    properties
        Thresholds
    end

    methods(TestMethodSetup)
        function setupPaths(testCase)
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'utils'));
            addpath(fullfile(repoRoot, 'dependencies'));
            if exist('OCTAVE_VERSION', 'builtin')
                addpath(fullfile(repoRoot, '.octave_compat'));
            end
            testCase.Thresholds = 0.8e-3 : 0.1e-3 : 2.0e-3;  % 13 values
        end
    end

    methods(Test)

        function testStructuralContract(testCase)
            [vf, lf] = makeNullCohort(40, 13, testCase.Thresholds, 11);
            cv = optimize_adc_threshold_nested_cv(testCase.Thresholds, vf, lf);
            fields = {'oof_vol_frac', 'oof_idx', 'oof_label', 'n_folds', ...
                'n_repeats', 'cv_scheme', 'n_inestimable_folds', 'oof_pvalue', ...
                'oof_auc', 'oof_n_lc', 'oof_n_lf', 'selected_threshold_counts', ...
                'recommended_thresh', 'recommended_idx', 'modal_fraction', ...
                'n_unique_selected', 'threshold_stability_std', ...
                'threshold_stability_range'};
            for i = 1:numel(fields)
                testCase.verifyTrue(isfield(cv, fields{i}), ...
                    sprintf('cv_results must contain field %s.', fields{i}));
            end
            testCase.verifyEqual(cv.n_folds, 40);
            testCase.verifyEqual(numel(cv.oof_vol_frac), 40);
            testCase.verifyEqual(numel(cv.selected_threshold_counts), 13);
            testCase.verifyEqual(cv.cv_scheme, 'leave-one-patient-out');
        end

        function testRepeatedKfoldRuns(testCase)
            % Repeated grouped K-fold path returns a populated, valid struct.
            [vf, lf] = makeSignalCohort(48, 24, testCase.Thresholds, 1.3e-3, 21);
            cv = optimize_adc_threshold_nested_cv(testCase.Thresholds, vf, lf, ...
                3, 5, 10, 99);
            testCase.verifyEqual(cv.n_repeats, 10);
            testCase.verifyTrue(any(~isnan(cv.oof_vol_frac)), ...
                'Repeated K-fold should produce out-of-fold values.');
            testCase.verifyGreaterThanOrEqual(sum(cv.selected_threshold_counts), 1);
        end

        function testRemovesOptimismOnNull(testCase)
            % CORE ASSERTION (CP3 test 2): out-of-fold is LESS optimistic than
            % the in-sample min-p selection on NULL data.
            n_cohorts = 25;
            insample_sig = 0; oof_sig = 0;
            insample_p = nan(n_cohorts, 1);
            oof_p = nan(n_cohorts, 1);
            for c = 1:n_cohorts
                [vf, lf] = makeNullCohort(42, 14, testCase.Thresholds, 2000 + c);
                [~, ~, ip] = select_significance_threshold(testCase.Thresholds, vf, lf);
                insample_p(c) = ip;
                if ~isnan(ip) && ip < 0.05, insample_sig = insample_sig + 1; end
                cv = optimize_adc_threshold_nested_cv(testCase.Thresholds, vf, lf);
                oof_p(c) = cv.oof_pvalue;
                if ~isnan(cv.oof_pvalue) && cv.oof_pvalue < 0.05
                    oof_sig = oof_sig + 1;
                end
            end
            mean_insample = mean(insample_p(~isnan(insample_p)));
            mean_oof = mean(oof_p(~isnan(oof_p)));

            % Directional, robust: out-of-fold p is larger (less optimistic).
            testCase.verifyGreaterThan(mean_oof, mean_insample, ...
                sprintf('Out-of-fold mean p (%.3f) should exceed in-sample min-p mean (%.3f).', ...
                mean_oof, mean_insample));
            % And it does not over-fire on the null.
            testCase.verifyGreaterThan(mean_oof, 0.25, ...
                sprintf('Out-of-fold mean p should move toward 0.5 (was %.3f).', mean_oof));
            testCase.verifyGreaterThanOrEqual(insample_sig, oof_sig, ...
                'In-sample selection should be "significant" at least as often as out-of-fold.');
        end

        function testRecoversTrueSignal(testCase)
            % CORE ASSERTION (CP3 test 3): repeated K-fold (the stronger
            % estimator) recovers a TRUE effect and the recommended threshold
            % lands in a plausible window around the truth.
            true_thresh = 1.3e-3;
            n_trials = 12;
            n_sig = 0; near_true = 0;
            for c = 1:n_trials
                [vf, lf] = makeSignalCohort(48, 24, testCase.Thresholds, ...
                    true_thresh, 3000 + c);
                cv = optimize_adc_threshold_nested_cv(testCase.Thresholds, vf, lf, ...
                    5, 10, 7, 4000 + c);
                if ~isnan(cv.oof_pvalue) && cv.oof_pvalue < 0.05
                    n_sig = n_sig + 1;
                end
                if ~isnan(cv.recommended_thresh) && ...
                        abs(cv.recommended_thresh - true_thresh) <= 0.4e-3
                    near_true = near_true + 1;
                end
            end
            % Real signal is recovered out-of-fold in the large majority.
            testCase.verifyGreaterThanOrEqual(n_sig, ceil(0.6 * n_trials), ...
                sprintf('Out-of-fold should detect true signal in most trials (%d/%d).', ...
                n_sig, n_trials));
            % Recommended threshold clusters in a plausible window around truth
            % (a broad window: threshold identification is honestly noisy at N).
            testCase.verifyGreaterThanOrEqual(near_true, ceil(0.5 * n_trials), ...
                sprintf('Recommended threshold should land near truth in >=half (%d/%d).', ...
                near_true, n_trials));
        end

        function testNoLeakageOfHeldOutPatient(testCase)
            % CORE ASSERTION (CP3 test 4): in LOO, perturbing a held-out
            % patient's data must NOT change the threshold selected for that
            % patient's own fold -- proving the inner selection never sees it.
            [vf, lf] = makeSignalCohort(48, 24, testCase.Thresholds, 1.3e-3, 77);
            cv1 = optimize_adc_threshold_nested_cv(testCase.Thresholds, vf, lf);

            i = find(~isnan(cv1.oof_idx), 1, 'first');
            testCase.assumeNotEmpty(i, 'Need at least one estimable fold.');
            vf2 = vf;
            vf2(i, :) = rand(1, numel(testCase.Thresholds));  % garbage row
            cv2 = optimize_adc_threshold_nested_cv(testCase.Thresholds, vf2, lf);

            testCase.verifyEqual(cv2.oof_idx(i), cv1.oof_idx(i), ...
                ['Held-out patient''s fold selection must be invariant to that ' ...
                 'patient''s own data (no leakage).']);
        end

        function testInestimableFoldsReported(testCase)
            % CP3 test 5: with only 3 LF patients, leaving any LF out drops a
            % group below the >=3 floor -> those folds are inestimable and
            % must be NaN-propagated and counted, not silently dropped.
            [vf, lf] = makeNullCohort(20, 3, testCase.Thresholds, 5);
            cv = optimize_adc_threshold_nested_cv(testCase.Thresholds, vf, lf);
            testCase.verifyGreaterThan(cv.n_inestimable_folds, 0, ...
                'With 3 LF patients some LOO folds must be inestimable.');
            % Folds left inestimable leave NaN in the per-patient index vector.
            testCase.verifyTrue(any(isnan(cv.oof_idx)), ...
                'Inestimable folds must leave NaN, not be silently dropped.');
        end

    end
end


function [vf, lf] = makeNullCohort(n_pts, n_lf, thresholds, seed)
% Synthetic NULL cohort: per-patient 2-component ADC voxel histogram with
% patient-specific random params; LF labels INDEPENDENT of histograms.
    rng(seed);
    n_vox = 150;
    n_t = numel(thresholds);
    vf = nan(n_pts, n_t);
    for j = 1:n_pts
        w  = 0.2 + 0.6 * rand();
        m1 = 1.0e-3 + 0.18e-3 * randn();
        m2 = 1.9e-3 + 0.22e-3 * randn();
        comp = rand(n_vox, 1) < w;
        adc = nan(n_vox, 1);
        adc(comp)  = m1 + 0.15e-3 * randn(sum(comp), 1);
        adc(~comp) = m2 + 0.30e-3 * randn(sum(~comp), 1);
        for ti = 1:n_t
            vf(j, ti) = mean(adc < thresholds(ti));
        end
    end
    lf = zeros(n_pts, 1);
    perm = randperm(n_pts);
    lf(perm(1:n_lf)) = 1;
end


function [vf, lf] = makeSignalCohort(n_pts, n_lf, thresholds, true_thresh, seed)
% Synthetic cohort with a TRUE effect: LF patients carry a larger cellular
% (low-ADC) core fraction, so their sub-volume fraction below `true_thresh`
% is genuinely higher; the LC/LF gap concentrates near `true_thresh`.
    rng(seed);
    n_vox = 150;
    n_t = numel(thresholds);
    lf = zeros(n_pts, 1);
    perm = randperm(n_pts);
    lf(perm(1:n_lf)) = 1;
    vf = nan(n_pts, n_t);
    for j = 1:n_pts
        if lf(j) == 1
            w = 0.55 + 0.10 * rand();
        else
            w = 0.25 + 0.10 * rand();
        end
        m1 = (true_thresh - 0.2e-3) + 0.05e-3 * randn();
        m2 = 2.0e-3 + 0.15e-3 * randn();
        comp = rand(n_vox, 1) < w;
        adc = nan(n_vox, 1);
        adc(comp)  = m1 + 0.10e-3 * randn(sum(comp), 1);
        adc(~comp) = m2 + 0.20e-3 * randn(sum(~comp), 1);
        for ti = 1:n_t
            vf(j, ti) = mean(adc < thresholds(ti));
        end
    end
end
