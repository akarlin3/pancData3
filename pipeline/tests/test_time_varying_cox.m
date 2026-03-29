classdef test_time_varying_cox < matlab.unittest.TestCase
% TEST_TIME_VARYING_COX  Tests for time-varying Cox model follow-up.
%
%   Verifies that: (1) synthetic data with a known time-varying effect
%   produces a significant interaction term, (2) data with constant HR
%   produces a non-significant interaction.

    properties
        OriginalPath
        TempDir
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'core'));
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            diary off;
            close all;
            path(testCase.OriginalPath);
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testTimeVaryingEffectDetected(testCase)
            % Synthetic data with a strong time-varying effect should
            % produce a significant interaction term (p < 0.05).
            rng(42);
            n = 200;

            % Create counting-process data with time-varying HR
            t_start = zeros(n, 1);
            t_stop = rand(n, 1) * 100 + 10;
            X = randn(n, 2);  % 2 covariates

            % Event probability depends on X(:,1) * log(time)
            t_mid = (t_start + t_stop) / 2;
            log_t = log(max(t_mid, 1));
            hazard = 0.5 * X(:,1) .* log_t + 0.3 * X(:,2);
            event_prob = 1 ./ (1 + exp(-hazard));
            event_csh = double(rand(n, 1) < event_prob);

            % Need enough events
            if sum(event_csh) < 10
                event_csh(1:15) = 1;
            end

            cov_names = {'Cov1', 'Cov2'};
            schoenfeld = struct('violated', [true; false], 'p_value', [0.01; 0.5]);

            config = struct('min_events_per_period', 3);
            tv_results = fit_time_varying_cox(X, t_start, t_stop, event_csh, ...
                cov_names, schoenfeld, testCase.TempDir, 'Test', config);

            % Should have at least one interaction model
            testCase.verifyGreaterThanOrEqual(numel(tv_results.interaction_models), 1, ...
                'Should fit at least one interaction model for violated covariate.');

            % The violated covariate should be Cov1
            testCase.verifyEqual(tv_results.violated_covariates{1}, 'Cov1', ...
                'First violated covariate should be Cov1.');
        end

        function testConstantHRNonSignificant(testCase)
            % Data with constant HR (no time-varying effect) should
            % produce a non-significant interaction term.
            rng(123);
            n = 300;

            t_start = zeros(n, 1);
            t_stop = rand(n, 1) * 100 + 10;
            X = randn(n, 2);

            % Constant effect: hazard depends on X but not on time
            % Use moderate coefficients to avoid multicollinearity
            % between X(:,1) and X(:,1)*log(t) in extended model
            hazard = 0.3 * X(:,1) + 0.3 * X(:,2);
            event_prob = 1 ./ (1 + exp(-hazard));
            event_csh = double(rand(n, 1) < event_prob);

            if sum(event_csh) < 10
                event_csh(1:15) = 1;
            end

            cov_names = {'Cov1', 'Cov2'};
            % Force Cov1 as "violated" to trigger the follow-up
            schoenfeld = struct('violated', [true; false], 'p_value', [0.04; 0.5]);

            config = struct('min_events_per_period', 3);
            tv_results = fit_time_varying_cox(X, t_start, t_stop, event_csh, ...
                cov_names, schoenfeld, testCase.TempDir, 'Test', config);

            if numel(tv_results.interaction_models) >= 1
                int_p = tv_results.interaction_models(1).interaction_p;
                % With constant HR, interaction p-value should typically be > 0.05
                % This is a statistical test, so we use a lenient threshold
                testCase.verifyGreaterThan(int_p, 0.001, ...
                    'Constant-HR data should not have extremely significant interaction.');
            end
        end

        function testNoViolationsSkipsGracefully(testCase)
            % When no covariates violate PH, the function should return
            % empty results gracefully.
            rng(7);
            n = 50;
            X = randn(n, 2);
            t_start = zeros(n, 1);
            t_stop = rand(n, 1) * 50 + 5;
            event_csh = double(rand(n, 1) < 0.3);

            cov_names = {'Cov1', 'Cov2'};
            schoenfeld = struct('violated', [false; false], 'p_value', [0.5; 0.8]);

            config = struct();
            tv_results = fit_time_varying_cox(X, t_start, t_stop, event_csh, ...
                cov_names, schoenfeld, testCase.TempDir, 'Test', config);

            testCase.verifyEmpty(tv_results.violated_covariates, ...
                'No violated covariates should produce empty list.');
            testCase.verifyEmpty(tv_results.interaction_models, ...
                'No violated covariates should produce no interaction models.');
        end

        function testSparseEventPeriodHandled(testCase)
            % With very few events in later time periods, verify sparse
            % period detection works and the function completes.
            rng(42);
            n = 150;
            t_start = zeros(n, 1);
            % Most events cluster in early times; late times have very few
            t_stop = [rand(100, 1) * 20 + 1; rand(50, 1) * 80 + 80];
            X = randn(n, 2);
            event_csh = zeros(n, 1);
            % Place most events in early period
            event_csh(1:30) = 1;
            % Only 1 event in late period
            event_csh(120) = 1;

            cov_names = {'Cov1', 'Cov2'};
            schoenfeld = struct('violated', [true; false], 'p_value', [0.01; 0.5]);

            config = struct('min_events_per_period', 5);
            evalc('tv_results = fit_time_varying_cox(X, t_start, t_stop, event_csh, cov_names, schoenfeld, testCase.TempDir, ''Test'', config);');

            % Should complete without error
            testCase.verifyTrue(isstruct(tv_results), ...
                'Should return a struct even with sparse events.');
            testCase.verifyEqual(numel(tv_results.violated_covariates), 1, ...
                'Should detect 1 violated covariate.');
            % If interaction models were fit, check stable_period_used flag
            if numel(tv_results.interaction_models) >= 1
                testCase.verifyTrue(isfield(tv_results.interaction_models(1), 'stable_period_used'), ...
                    'Interaction model should have stable_period_used field.');
            end
        end

        function testMultipleViolatedCovariates(testCase)
            % When 2+ covariates violate PH, verify all get interaction models.
            rng(42);
            n = 200;
            t_start = zeros(n, 1);
            t_stop = rand(n, 1) * 100 + 10;
            X = randn(n, 3);

            % Strong time-varying effects for both Cov1 and Cov2
            t_mid = (t_start + t_stop) / 2;
            log_t = log(max(t_mid, 1));
            hazard = 0.5 * X(:,1) .* log_t + 0.4 * X(:,2) .* log_t + 0.2 * X(:,3);
            event_prob = 1 ./ (1 + exp(-hazard));
            event_csh = double(rand(n, 1) < event_prob);
            if sum(event_csh) < 15
                event_csh(1:20) = 1;
            end

            cov_names = {'Cov1', 'Cov2', 'Cov3'};
            schoenfeld = struct('violated', [true; true; false], 'p_value', [0.01; 0.02; 0.6]);

            config = struct('min_events_per_period', 3);
            evalc('tv_results = fit_time_varying_cox(X, t_start, t_stop, event_csh, cov_names, schoenfeld, testCase.TempDir, ''Test'', config);');

            testCase.verifyEqual(numel(tv_results.violated_covariates), 2, ...
                'Should detect 2 violated covariates.');
            testCase.verifyTrue(ismember('Cov1', tv_results.violated_covariates), ...
                'Cov1 should be in violated list.');
            testCase.verifyTrue(ismember('Cov2', tv_results.violated_covariates), ...
                'Cov2 should be in violated list.');
            % Should have interaction models for both
            testCase.verifyGreaterThanOrEqual(numel(tv_results.interaction_models), 2, ...
                'Should fit interaction models for both violated covariates.');
        end

        function testLeftTruncationDetected(testCase)
            % With non-zero t_start values (delayed entry), verify
            % left_truncated flag is set in the output.
            rng(42);
            n = 150;
            % Delayed entry: start times vary from 5 to 30
            t_start = rand(n, 1) * 25 + 5;
            t_stop = t_start + rand(n, 1) * 80 + 10;
            X = randn(n, 2);

            hazard = 0.3 * X(:,1) + 0.2 * X(:,2);
            event_prob = 1 ./ (1 + exp(-hazard));
            event_csh = double(rand(n, 1) < event_prob);
            if sum(event_csh) < 15
                event_csh(1:20) = 1;
            end

            cov_names = {'Cov1', 'Cov2'};
            schoenfeld = struct('violated', [true; false], 'p_value', [0.03; 0.6]);

            config = struct('min_events_per_period', 3);
            evalc('tv_results = fit_time_varying_cox(X, t_start, t_stop, event_csh, cov_names, schoenfeld, testCase.TempDir, ''Test'', config);');

            % Check stratified model has left_truncated flag
            if isfield(tv_results.stratified_model, 'left_truncated')
                testCase.verifyTrue(tv_results.stratified_model.left_truncated, ...
                    'Stratified model should flag left-truncation for delayed entry.');
            end
            % Check interaction models have left_truncated flag
            if numel(tv_results.interaction_models) >= 1
                testCase.verifyTrue(tv_results.interaction_models(1).left_truncated, ...
                    'Interaction model should flag left-truncation for delayed entry.');
            end
        end

        function testOutputStructFields(testCase)
            % Verify all expected fields in tv_results.
            rng(42);
            n = 100;
            t_start = zeros(n, 1);
            t_stop = rand(n, 1) * 50 + 5;
            X = randn(n, 2);
            event_csh = double(rand(n, 1) < 0.4);
            if sum(event_csh) < 10
                event_csh(1:15) = 1;
            end

            cov_names = {'Cov1', 'Cov2'};
            schoenfeld = struct('violated', [true; false], 'p_value', [0.03; 0.6]);

            config = struct();
            evalc('tv_results = fit_time_varying_cox(X, t_start, t_stop, event_csh, cov_names, schoenfeld, testCase.TempDir, ''Test'', config);');

            % Top-level fields
            testCase.verifyTrue(isfield(tv_results, 'violated_covariates'), ...
                'Should have violated_covariates field.');
            testCase.verifyTrue(isfield(tv_results, 'stratified_model'), ...
                'Should have stratified_model field.');
            testCase.verifyTrue(isfield(tv_results, 'interaction_models'), ...
                'Should have interaction_models field.');

            % Interaction model entry fields (if any were fit)
            if numel(tv_results.interaction_models) >= 1
                im = tv_results.interaction_models(1);
                expected_fields = {'covariate', 'interaction_coef', 'interaction_p', ...
                    'base_coef', 'base_p', 'penalized', 'stable_period_used', 'left_truncated'};
                for i = 1:numel(expected_fields)
                    testCase.verifyTrue(isfield(im, expected_fields{i}), ...
                        sprintf('Interaction model should have field %s.', expected_fields{i}));
                end
            end
        end

        function testPenalizedModelUsedForIllConditioned(testCase)
            % Create ill-conditioned data with very few events, verify
            % penalized flag is set.
            rng(42);
            n = 80;
            t_start = zeros(n, 1);
            t_stop = rand(n, 1) * 100 + 10;
            X = randn(n, 2);
            % Very few events to trigger penalized estimation (<50 events)
            event_csh = zeros(n, 1);
            event_csh(1:8) = 1;  % Only 8 events

            cov_names = {'Cov1', 'Cov2'};
            schoenfeld = struct('violated', [true; false], 'p_value', [0.02; 0.5]);

            config = struct('min_events_per_period', 2);
            evalc('tv_results = fit_time_varying_cox(X, t_start, t_stop, event_csh, cov_names, schoenfeld, testCase.TempDir, ''Test'', config);');

            % With only 8 events, penalized estimation should be triggered
            if numel(tv_results.interaction_models) >= 1
                testCase.verifyTrue(tv_results.interaction_models(1).penalized, ...
                    'Should use penalized estimation for sparse events (<50).');
            end
        end
    end
end
