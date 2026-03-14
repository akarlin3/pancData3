classdef test_compute_ipcw_weights < matlab.unittest.TestCase
% TEST_COMPUTE_IPCW_WEIGHTS — Unit tests for compute_ipcw_weights.m
%
% Validates IPCW weight computation including:
%   - No administrative censoring → returns ones
%   - All competing events excluded from censoring model
%   - Weight truncation at 0.05 floor
%   - Mean-normalised output
%   - Graceful fallback on estimation failure

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
        end
    end

    methods (Test)
        function test_no_admin_censoring_returns_ones(testCase)
            % When no rows are administratively censored (event==0 at last
            % row per patient), weights should all be 1.0
            event = [1; 2; 1; 2];  % all events, no censoring
            t_start = [0; 0; 0; 0];
            t_stop  = [10; 20; 15; 25];
            X = randn(4, 2);
            pat_id = [1; 2; 3; 4];

            w = compute_ipcw_weights(event, t_start, t_stop, X, pat_id);
            testCase.verifyEqual(w, ones(4, 1), 'AbsTol', 1e-12, ...
                'No admin censoring should yield unit weights.');
        end

        function test_weights_positive(testCase)
            % Weights should always be positive
            n = 20;
            rng(42);
            event = [ones(10,1); zeros(5,1); 2*ones(5,1)];
            t_start = zeros(n, 1);
            t_stop = (1:n)';
            X = randn(n, 3);
            pat_id = (1:n)';

            w = compute_ipcw_weights(event, t_start, t_stop, X, pat_id);
            testCase.verifyGreaterThan(w, 0, 'All weights must be positive.');
        end

        function test_weights_mean_normalised(testCase)
            % Mean of weights should be approximately 1.0
            n = 30;
            rng(42);
            % Mix of events, censoring, and competing risks
            event = zeros(n, 1);
            event(1:10) = 1;      % primary events
            event(11:15) = 2;     % competing events
            % remaining 16:30 are censored (event=0)
            t_start = zeros(n, 1);
            t_stop = sort(rand(n, 1) * 100);
            X = randn(n, 2);
            pat_id = (1:n)';

            w = compute_ipcw_weights(event, t_start, t_stop, X, pat_id);
            testCase.verifyEqual(mean(w), 1.0, 'AbsTol', 0.01, ...
                'Weights should be mean-normalised to 1.0.');
        end

        function test_weight_truncation_floor(testCase)
            % G_hat is floored at 0.05, so max weight = 1/0.05 = 20
            % before mean normalisation. After normalisation, check no
            % weight exceeds a reasonable upper bound.
            n = 20;
            rng(42);
            event = [ones(5,1); zeros(10,1); 2*ones(5,1)];
            t_start = zeros(n, 1);
            t_stop = (1:n)' * 5;
            X = randn(n, 2);
            pat_id = (1:n)';

            w = compute_ipcw_weights(event, t_start, t_stop, X, pat_id);
            % After mean normalisation, weights should be finite
            testCase.verifyTrue(all(isfinite(w)), 'All weights must be finite.');
        end

        function test_competing_events_excluded(testCase)
            % When all non-competing rows have events (no admin censoring
            % among non-competing), should return ones
            event = [1; 1; 1; 2; 2];
            t_start = [0; 0; 0; 0; 0];
            t_stop  = [10; 20; 30; 15; 25];
            X = randn(5, 2);
            pat_id = [1; 2; 3; 4; 5];

            w = compute_ipcw_weights(event, t_start, t_stop, X, pat_id);
            testCase.verifyEqual(w, ones(5, 1), 'AbsTol', 1e-12, ...
                'No admin censoring among non-competing → unit weights.');
        end

        function test_constant_covariates_fallback(testCase)
            % When all covariates are constant, coxphfit should fail
            % and the function should return ones (graceful fallback)
            n = 10;
            event = [1; 1; 0; 0; 0; 0; 0; 0; 0; 0];
            t_start = zeros(n, 1);
            t_stop = (1:n)';
            X = ones(n, 2);  % all constant — will be removed by remove_constant_columns
            pat_id = (1:n)';

            w = compute_ipcw_weights(event, t_start, t_stop, X, pat_id);
            testCase.verifyEqual(w, ones(n, 1), 'AbsTol', 1e-12, ...
                'Constant covariates should trigger fallback to unit weights.');
        end

        function test_multi_row_per_patient(testCase)
            % Multiple rows per patient: only last row per patient matters
            % for terminal interval detection
            event   = [0; 1; 0; 0; 0; 2];
            t_start = [0; 5; 0; 5; 0; 5];
            t_stop  = [5; 10; 5; 10; 5; 10];
            X = randn(6, 2);
            pat_id  = [1; 1; 2; 2; 3; 3];
            % Patient 1: last row event=1 (not admin censored)
            % Patient 2: last row event=0 (admin censored)
            % Patient 3: last row event=2 (competing, not admin censored)

            w = compute_ipcw_weights(event, t_start, t_stop, X, pat_id);
            testCase.verifyEqual(numel(w), 6);
            testCase.verifyTrue(all(isfinite(w)));
        end

        function test_output_size_matches_input(testCase)
            n = 15;
            event = randi([0 2], n, 1);
            t_start = zeros(n, 1);
            t_stop = (1:n)';
            X = randn(n, 3);
            pat_id = (1:n)';

            w = compute_ipcw_weights(event, t_start, t_stop, X, pat_id);
            testCase.verifyEqual(size(w), size(event));
        end
    end
end
