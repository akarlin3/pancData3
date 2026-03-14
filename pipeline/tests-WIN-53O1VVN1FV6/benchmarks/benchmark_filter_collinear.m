classdef benchmark_filter_collinear < matlab.unittest.TestCase
% BENCHMARK_FILTER_COLLINEAR  Performance benchmark for filter_collinear_features.
%   Creates a large synthetic feature matrix (200 samples x 2000 features) with
%   50 deliberately collinear feature pairs, then times how long the collinearity
%   pruning takes. Verifies that the function runs to completion and returns a
%   valid subset of feature indices.

    methods (Test)
        function testFilterCollinearRuns(testCase)
        %TESTFILTERCOLLINEARRUNS Benchmark filter_collinear_features on a
        %   large matrix with known collinear pairs.

            % Define benchmark dimensions
            n_samples = 200;
            n_features = 2000;
            n_collinear_pairs = 50;

            % Generate random feature matrix with fixed seed for reproducibility
            rng(123);
            X = randn(n_samples, n_features);
            y = randi([0, 1], n_samples, 1);

            % Inject collinear pairs: make the last 50 features near-copies of
            % the first 50, with small Gaussian noise to simulate real-world
            % collinearity (r ~ 0.99)
            for i = 1:n_collinear_pairs
                X(:, n_features - i + 1) = X(:, i) + 0.1 * randn(n_samples, 1);
            end

            % Warmup call to ensure JIT compilation does not skew timing
            filter_collinear_features(X(1:20, 1:20), y(1:20));

            % Timed benchmark run
            tic;
            keep_idx = filter_collinear_features(X, y);
            elapsed = toc;

            fprintf('Time taken: %.6f seconds\n', elapsed);
            fprintf('Features kept: %d / %d\n', length(keep_idx), n_features);

            % Verify the function returned a non-empty, valid subset
            testCase.verifyGreaterThan(length(keep_idx), 0);
            testCase.verifyLessThanOrEqual(length(keep_idx), n_features);
        end
    end
end
