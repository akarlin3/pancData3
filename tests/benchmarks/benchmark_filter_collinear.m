function benchmark_filter_collinear()
    % Benchmark for filter_collinear_features

    fprintf('Benchmarking filter_collinear_features...\n');

    % Parameters
    n_samples = 200;
    n_features = 2000;
    n_collinear_pairs = 50; % Only a few collinear pairs

    rng(123);
    X = randn(n_samples, n_features);
    y = randi([0, 1], n_samples, 1);

    % Create some collinear features
    % For the first n_collinear_pairs features, make the next n_collinear_pairs correlated
    for i = 1:n_collinear_pairs
        X(:, n_features - i + 1) = X(:, i) + 0.1 * randn(n_samples, 1);
    end

    % Warmup
    filter_collinear_features(X(1:20, 1:20), y(1:20));

    % Measure time
    tic;
    keep_idx = filter_collinear_features(X, y);
    elapsed = toc;

    fprintf('Time taken: %.6f seconds\n', elapsed);
    fprintf('Features kept: %d / %d\n', length(keep_idx), n_features);

end
