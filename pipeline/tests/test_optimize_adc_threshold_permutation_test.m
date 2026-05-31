classdef test_optimize_adc_threshold_permutation_test < matlab.unittest.TestCase
    % TEST_OPTIMIZE_ADC_THRESHOLD_PERMUTATION_TEST
    %
    % Proves the bias-measurement instrument is correct:
    %   (1) Structural contract of perm_results.
    %   (2) BIAS DETECTION: on NULL synthetic data (sub-volume fraction
    %       independent of LC/LF label) the naive min-p is biased LOW (the
    %       selection-on-outcome inflation) while the permutation-adjusted p
    %       is restored to ~uniform/near 0.5 -- i.e. the instrument actually
    %       detects and removes the bias.
    %   (3) The instrument NaN-propagates honestly when nothing is estimable.
    %
    % Mirrors the fixture/teardown style of test_optimize_adc_threshold.m
    % and the statistical-assertion style of test_bootstrap_ci.m.

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
            [vf, lf] = makeNullCohort(42, 14, testCase.Thresholds, 7);
            pr = optimize_adc_threshold_permutation_test( ...
                testCase.Thresholds, vf, lf, 500, 13);

            fields = {'observed_min_p', 'observed_thresh', 'observed_idx', ...
                'observed_pvalues', 'perm_adjusted_min_p', ...
                'per_threshold_adjusted_p', 'n_notable_naive', ...
                'n_notable_adjusted', 'null_min_p_dist', 'n_perm_requested', ...
                'n_perm_valid', 'n_labeled_patients', 'n_lc_total', ...
                'n_lf_total', 'rng_seed', 'method'};
            for i = 1:numel(fields)
                testCase.verifyTrue(isfield(pr, fields{i}), ...
                    sprintf('perm_results must contain field %s.', fields{i}));
            end
            testCase.verifyEqual(numel(pr.observed_pvalues), 13);
            testCase.verifyEqual(numel(pr.per_threshold_adjusted_p), 13);
            testCase.verifyEqual(numel(pr.null_min_p_dist), pr.n_perm_valid);
            % Adjusted p is a valid probability strictly in (0, 1].
            testCase.verifyGreaterThan(pr.perm_adjusted_min_p, 0);
            testCase.verifyLessThanOrEqual(pr.perm_adjusted_min_p, 1);
            testCase.verifyEqual(pr.n_labeled_patients, 42);
            testCase.verifyEqual(pr.n_lf_total, 14);
            testCase.verifyEqual(pr.n_lc_total, 28);
        end

        function testReproducibleSeed(testCase)
            [vf, lf] = makeNullCohort(40, 13, testCase.Thresholds, 3);
            pr1 = optimize_adc_threshold_permutation_test(testCase.Thresholds, vf, lf, 400, 99);
            pr2 = optimize_adc_threshold_permutation_test(testCase.Thresholds, vf, lf, 400, 99);
            testCase.verifyEqual(pr1.perm_adjusted_min_p, pr2.perm_adjusted_min_p, ...
                'Same seed must give identical adjusted p.');
        end

        function testDetectsBiasOnNullData(testCase)
            % CORE ASSERTION (CP0b / CP3 test 1).
            % Across many NULL cohorts (vol_frac independent of label):
            %   - the naive min-p is biased LOW (mean well below 0.5), and
            %   - the permutation-adjusted p is restored toward uniform
            %     (mean near 0.5), and is, on average, much larger than the
            %     naive min-p.
            % Means (not tail rates) are asserted so the test is not flaky.
            n_cohorts = 25;
            naive = nan(n_cohorts, 1);
            adj   = nan(n_cohorts, 1);
            for c = 1:n_cohorts
                [vf, lf] = makeNullCohort(42, 14, testCase.Thresholds, 1000 + c);
                pr = optimize_adc_threshold_permutation_test( ...
                    testCase.Thresholds, vf, lf, 400, 7 * c);
                naive(c) = pr.observed_min_p;
                adj(c)   = pr.perm_adjusted_min_p;
            end
            naive = naive(~isnan(naive));
            adj   = adj(~isnan(adj));

            mean_naive = mean(naive);
            mean_adj   = mean(adj);

            % The naive selection statistic is inflated: its mean sits well
            % below the 0.5 a calibrated p-value would average to.
            testCase.verifyLessThan(mean_naive, 0.40, ...
                sprintf(['Naive min-p should be biased low under the null ' ...
                'but mean was %.3f.'], mean_naive));
            % The adjustment restores calibration: mean adjusted p ~ 0.5.
            testCase.verifyGreaterThan(mean_adj, 0.38, ...
                sprintf('Adjusted p should be ~uniform (mean ~0.5) but was %.3f.', mean_adj));
            testCase.verifyLessThan(mean_adj, 0.62, ...
                sprintf('Adjusted p should be ~uniform (mean ~0.5) but was %.3f.', mean_adj));
            % And the adjustment moves the number in the de-biasing direction.
            testCase.verifyGreaterThan(mean_adj, mean_naive + 0.05, ...
                'Permutation adjustment must increase the p-value vs naive min-p.');
        end

        function testInestimableReturnsNaN(testCase)
            % Only 2 LF patients -> the >=3-per-group floor is never met, so
            % nothing is estimable: NaN-propagate, do not fabricate.
            [vf, lf] = makeNullCohort(20, 2, testCase.Thresholds, 5);
            pr = optimize_adc_threshold_permutation_test(testCase.Thresholds, vf, lf, 200, 1);
            testCase.verifyTrue(isnan(pr.observed_min_p));
            testCase.verifyTrue(isnan(pr.perm_adjusted_min_p));
            testCase.verifyEqual(pr.n_perm_valid, 0);
        end

    end
end


function [vf, lf] = makeNullCohort(n_pts, n_lf, thresholds, seed)
% Synthetic NULL cohort: each patient has an ADC voxel histogram (2-component
% mixture with patient-specific random parameters), so the rank order of
% patients by sub-volume fraction REORDERS across thresholds -- the per-
% threshold tests genuinely decorrelate.  LF labels are assigned at random
% and are INDEPENDENT of the histograms, so there is NO true LC/LF effect.
    rng(seed);
    n_vox = 150;
    n_t = numel(thresholds);
    vf = nan(n_pts, n_t);
    for j = 1:n_pts
        w  = 0.2 + 0.6 * rand();                 % cellular-core fraction
        m1 = 1.0e-3 + 0.18e-3 * randn();         % restricted component mean
        m2 = 1.9e-3 + 0.22e-3 * randn();         % free-water component mean
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
    lf(perm(1:n_lf)) = 1;   % labels independent of the histograms
end
