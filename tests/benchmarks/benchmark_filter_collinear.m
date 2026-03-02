classdef benchmark_filter_collinear < matlab.unittest.TestCase
    % Benchmark for filter_collinear_features performance

    methods (Test)
        function testFilterCollinearRuns(testCase)
            n_samples = 200;
            n_features = 2000;
            n_collinear_pairs = 50;

            rng(123);
            X = randn(n_samples, n_features);
            y = randi([0, 1], n_samples, 1);

            for i = 1:n_collinear_pairs
                X(:, n_features - i + 1) = X(:, i) + 0.1 * randn(n_samples, 1);
            end

            % Warmup
            filter_collinear_features(X(1:20, 1:20), y(1:20));

            tic;
            keep_idx = filter_collinear_features(X, y);
            elapsed = toc;

            fprintf('Time taken: %.6f seconds\n', elapsed);
            fprintf('Features kept: %d / %d\n', length(keep_idx), n_features);

            testCase.verifyGreaterThan(length(keep_idx), 0);
            testCase.verifyLessThanOrEqual(length(keep_idx), n_features);
        end
    end
end
