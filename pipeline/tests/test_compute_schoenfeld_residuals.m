classdef test_compute_schoenfeld_residuals < matlab.unittest.TestCase
% TEST_COMPUTE_SCHOENFELD_RESIDUALS — Unit tests for compute_schoenfeld_residuals.m
%
% Validates Schoenfeld residual computation and PH assumption testing:
%   - Output struct field completeness
%   - Insufficient events / non-finite beta skip behavior
%   - Constant covariate handling
%   - Known PH violation detection
%   - Correct residual dimensions and event times
%   - Figure file generation
%   - Feature name usage

    properties
        origPath
        tmpDir
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            testCase.origPath = path;
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'core'));
            testCase.tmpDir = tempname;
            mkdir(testCase.tmpDir);
        end
    end

    methods (TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.origPath);
            diary off;
            if isfolder(testCase.tmpDir)
                rmdir(testCase.tmpDir, 's');
            end
        end
    end

    methods (Test)
        function testOutputStructFields(testCase)
            % Verify all expected fields exist in output struct
            rng(42);
            n = 30;
            X = randn(n, 2);
            t_start = zeros(n, 1);
            t_stop = (1:n)';
            event = zeros(n, 1);
            event(5) = 1; event(10) = 1; event(15) = 1; event(20) = 1;
            beta = [0.5; -0.3];
            feat_names = {'feat1', 'feat2'};
            evalc('results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);');
            testCase.verifyTrue(isfield(results, 'rho'), 'Missing field: rho');
            testCase.verifyTrue(isfield(results, 'chi2'), 'Missing field: chi2');
            testCase.verifyTrue(isfield(results, 'p_value'), 'Missing field: p_value');
            testCase.verifyTrue(isfield(results, 'violated'), 'Missing field: violated');
            testCase.verifyTrue(isfield(results, 'residuals'), 'Missing field: residuals');
            testCase.verifyTrue(isfield(results, 'event_times'), 'Missing field: event_times');
        end

        function testInsufficientEventsSkips(testCase)
            % With <3 events, should skip and return NaN results
            rng(42);
            n = 20;
            X = randn(n, 2);
            t_start = zeros(n, 1);
            t_stop = (1:n)';
            event = zeros(n, 1);
            event(5) = 1; event(10) = 1;  % Only 2 events
            beta = [0.5; -0.3];
            feat_names = {'feat1', 'feat2'};
            evalc('results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);');
            testCase.verifyTrue(all(isnan(results.rho)), ...
                'rho should be NaN when <3 events.');
            testCase.verifyTrue(all(isnan(results.p_value)), ...
                'p_value should be NaN when <3 events.');
            testCase.verifyEmpty(results.residuals, ...
                'residuals should be empty when <3 events.');
        end

        function testNonFiniteBetaSkips(testCase)
            % With Inf/NaN beta, should skip and return NaN results
            rng(42);
            n = 20;
            X = randn(n, 2);
            t_start = zeros(n, 1);
            t_stop = (1:n)';
            event = zeros(n, 1);
            event(3) = 1; event(7) = 1; event(12) = 1; event(16) = 1;
            beta = [Inf; 0.5];
            feat_names = {'feat1', 'feat2'};
            evalc('results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);');
            testCase.verifyTrue(all(isnan(results.rho)), ...
                'rho should be NaN when beta contains Inf.');
            testCase.verifyEmpty(results.residuals, ...
                'residuals should be empty when beta is non-finite.');

            % Also test with NaN beta
            beta2 = [NaN; -0.3];
            evalc('results2 = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta2, feat_names);');
            testCase.verifyTrue(all(isnan(results2.rho)), ...
                'rho should be NaN when beta contains NaN.');
        end

        function testConstantCovariateNoViolation(testCase)
            % A constant covariate should not show PH violation
            rng(42);
            n = 40;
            X = [ones(n, 1), randn(n, 1)];  % Col 1 is constant
            t_start = zeros(n, 1);
            t_stop = (1:n)';
            event = zeros(n, 1);
            event(5) = 1; event(10) = 1; event(15) = 1; event(20) = 1;
            event(25) = 1; event(30) = 1; event(35) = 1;
            beta = [0.1; 0.2];
            feat_names = {'constant', 'varying'};
            evalc('results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);');
            testCase.verifyEqual(results.rho(1), 0, 'AbsTol', 1e-12, ...
                'Constant covariate should have rho = 0.');
            testCase.verifyEqual(results.p_value(1), 1, 'AbsTol', 1e-12, ...
                'Constant covariate should have p = 1.');
            testCase.verifyFalse(results.violated(1), ...
                'Constant covariate should not be flagged as PH violation.');
        end

        function testKnownViolation(testCase)
            % Create data where covariate effect changes over time
            % (multiply by log(t)), which should cause PH violation
            rng(42);
            n = 200;
            t_stop = sort(rand(n, 1) * 100 + 1);
            t_start = [0; t_stop(1:end-1)];
            % Create a covariate whose effect explicitly depends on time
            X = randn(n, 1);
            % Generate events where hazard depends on X * log(t) — time-varying effect
            log_t = log(t_stop);
            hazard = exp(X .* log_t * 0.5);
            event = double(rand(n, 1) < (hazard ./ (1 + hazard)));
            % Ensure at least 10 events
            if sum(event) < 10
                idx = find(event == 0);
                event(idx(1:10)) = 1;
            end
            beta = [0.3];
            feat_names = {'time_varying_effect'};
            evalc('results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);');
            % With a strong time-varying effect, we expect violation
            % Use a lenient check: rho should be noticeably non-zero
            testCase.verifyGreaterThan(abs(results.rho(1)), 0.05, ...
                'Time-varying covariate effect should produce non-trivial correlation.');
        end

        function testNoViolationWithConstantHR(testCase)
            % Data with constant HR should NOT violate PH assumption
            rng(42);
            n = 100;
            t_stop = sort(rand(n, 1) * 50 + 1);
            t_start = [0; t_stop(1:end-1)];
            X = randn(n, 1);
            % Generate events with constant hazard ratio (no time dependence)
            beta_true = 0.5;
            hazard = exp(X * beta_true) * 0.3;
            event = double(rand(n, 1) < min(hazard, 0.8));
            % Ensure enough events
            if sum(event) < 10
                idx = find(event == 0);
                event(idx(1:min(10, length(idx)))) = 1;
            end
            beta = [beta_true];
            feat_names = {'constant_hr'};
            evalc('results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);');
            testCase.verifyGreaterThan(results.p_value(1), 0.01, ...
                'Constant HR data should generally not violate PH (p > 0.01).');
        end

        function testResidualsDimensions(testCase)
            % Verify residuals matrix has dimensions [n_events x p]
            rng(42);
            n = 30;
            p = 3;
            X = randn(n, p);
            t_start = zeros(n, 1);
            t_stop = (1:n)';
            event = zeros(n, 1);
            event(5) = 1; event(10) = 1; event(15) = 1;
            event(20) = 1; event(25) = 1;
            beta = randn(p, 1) * 0.3;
            feat_names = {'f1', 'f2', 'f3'};
            evalc('results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);');
            n_events = sum(event);
            testCase.verifyEqual(size(results.residuals), [n_events, p], ...
                'Residuals should be [n_events x p].');
        end

        function testEventTimesCorrect(testCase)
            % Verify event_times matches the stop times of event rows
            rng(42);
            n = 20;
            X = randn(n, 2);
            t_start = zeros(n, 1);
            t_stop = (10:10:200)';
            event = zeros(n, 1);
            event_rows = [3, 7, 12, 16];
            event(event_rows) = 1;
            beta = [0.2; -0.1];
            feat_names = {'a', 'b'};
            evalc('results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);');
            expected_times = t_stop(event_rows);
            testCase.verifyEqual(results.event_times, expected_times, 'AbsTol', 1e-12, ...
                'Event times should match t_stop values at event rows.');
        end

        function testFigureSaved(testCase)
            % When output_folder provided, verify PNG file is created
            rng(42);
            n = 30;
            X = randn(n, 2);
            t_start = zeros(n, 1);
            t_stop = (1:n)';
            event = zeros(n, 1);
            event(5) = 1; event(10) = 1; event(15) = 1;
            event(20) = 1; event(25) = 1;
            beta = [0.3; -0.2];
            feat_names = {'x1', 'x2'};
            dtype_label = 'Standard';
            evalc('results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names, testCase.tmpDir, dtype_label);');
            expected_file = fullfile(testCase.tmpDir, 'schoenfeld_residuals_Standard.png');
            testCase.verifyTrue(isfile(expected_file), ...
                'Schoenfeld diagnostic PNG should be saved to output folder.');
        end

        function testFeatNamesUsed(testCase)
            % Verify feature names appear correctly in output dimensions
            rng(42);
            n = 30;
            feat_names = {'adc_mean', 'f_median', 'dstar_kurt'};
            p = length(feat_names);
            X = randn(n, p);
            t_start = zeros(n, 1);
            t_stop = (1:n)';
            event = zeros(n, 1);
            event(5) = 1; event(10) = 1; event(15) = 1;
            event(20) = 1; event(25) = 1;
            beta = randn(p, 1) * 0.2;
            evalc('results = compute_schoenfeld_residuals(X, t_start, t_stop, event, beta, feat_names);');
            % Verify output vectors have correct length matching feature count
            testCase.verifyEqual(length(results.rho), p, ...
                'rho length should match number of features.');
            testCase.verifyEqual(length(results.p_value), p, ...
                'p_value length should match number of features.');
            testCase.verifyEqual(length(results.violated), p, ...
                'violated length should match number of features.');
            testCase.verifyEqual(size(results.residuals, 2), p, ...
                'Residuals column count should match number of features.');
        end
    end
end
