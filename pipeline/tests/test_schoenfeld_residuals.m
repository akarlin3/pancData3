classdef test_schoenfeld_residuals < matlab.unittest.TestCase
    % TEST_SCHOENFELD_RESIDUALS  Tests for PH assumption via Schoenfeld residuals.
    %
    % Covers:
    %   (a) Synthetic dataset where PH holds
    %   (b) Synthetic dataset with a known time-varying effect
    %   (c) Edge case with <5 events

    properties
        OriginalPath
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            close all;
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testPHHolds(testCase)
            % Synthetic data where PH holds: constant effect, no time interaction.
            rng(42);
            n = 50;
            p = 2;
            X = randn(n, p);
            beta = [0.5; -0.3];

            % Generate event times from exponential model (PH holds)
            t_start = zeros(n, 1);
            t_stop = exprnd(10, n, 1) .* exp(-X * beta);
            event = ones(n, 1);  % all events
            % Censor some
            cens_time = exprnd(15, n, 1);
            censored = t_stop > cens_time;
            t_stop(censored) = cens_time(censored);
            event(censored) = 0;

            feat_names = {'X1', 'X2'};
            results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);

            % PH should hold: p-values should be > 0.05
            testCase.verifyFalse(any(results.violated), ...
                'PH should not be violated for constant-effect data.');
            testCase.verifyEqual(length(results.p_value), p);
            testCase.verifyGreaterThan(min(results.p_value), 0);
        end

        function testTimeVaryingEffect(testCase)
            % Known time-varying effect: covariate effect increases with time.
            rng(123);
            n = 100;
            X = randn(n, 1);

            % Time-varying hazard: h(t) = h0(t) * exp(beta(t) * X)
            % where beta(t) = 0.1 * t (increasing with time)
            t_start = zeros(n, 1);
            t_stop = abs(randn(n, 1) * 5) + 1;
            % Make events correlate with X*t (time-varying effect)
            hazard = exp(0.3 * X .* t_stop);
            event = double(rand(n, 1) < 1 - exp(-hazard * 0.1));

            % Need enough events
            if sum(event) < 10
                event(1:10) = 1;
            end

            beta = 0.3;
            feat_names = {'X1'};
            results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);

            % Should return valid results
            testCase.verifyTrue(isfinite(results.rho(1)), ...
                'Should compute finite correlation.');
            testCase.verifyTrue(~isempty(results.residuals), ...
                'Should have residuals.');
        end

        function testFewEvents(testCase)
            % Edge case with <5 events: should skip gracefully.
            rng(7);
            n = 20;
            X = randn(n, 2);
            t_start = zeros(n, 1);
            t_stop = rand(n, 1) * 10 + 1;
            event = zeros(n, 1);
            event(1:2) = 1;  % only 2 events

            beta = [0.1; -0.1];
            feat_names = {'X1', 'X2'};

            output = evalc('results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);');

            % Should not crash and should return NaN p-values
            testCase.verifyTrue(all(isnan(results.p_value)), ...
                'Should return NaN p-values with <3 events.');
        end

        function testDiagnosticFigureSaved(testCase)
            % Verify that the diagnostic figure is saved when output_folder is provided.
            rng(42);
            n = 30;
            X = randn(n, 2);
            beta = [0.3; -0.2];
            t_start = zeros(n, 1);
            t_stop = exprnd(8, n, 1);
            event = ones(n, 1);
            event(end-5:end) = 0;

            tmp_dir = tempname;
            mkdir(tmp_dir);
            cleanup = onCleanup(@() rmdir(tmp_dir, 's'));

            feat_names = {'X1', 'X2'};
            evalc('results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names, tmp_dir, ''Test'');');

            fig_path = fullfile(tmp_dir, 'schoenfeld_residuals_Test.png');
            testCase.verifyTrue(exist(fig_path, 'file') == 2, ...
                'Diagnostic figure should be saved.');
        end

        function testResidualDimensions(testCase)
            % Verify residual matrix dimensions match n_events x p.
            rng(99);
            n = 40;
            p = 3;
            X = randn(n, p);
            beta = randn(p, 1) * 0.2;
            t_start = zeros(n, 1);
            t_stop = exprnd(5, n, 1);
            event = ones(n, 1);
            event(end-10:end) = 0;

            feat_names = arrayfun(@(i) sprintf('X%d', i), 1:p, 'UniformOutput', false);
            results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);

            n_events = sum(event == 1);
            testCase.verifyEqual(size(results.residuals, 1), n_events);
            testCase.verifyEqual(size(results.residuals, 2), p);
        end

    end
end
