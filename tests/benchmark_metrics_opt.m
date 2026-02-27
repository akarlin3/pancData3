% Benchmark script for optimizing dynamic array growth in metrics.m
% This script simulates the data collection loop in Section 8 of metrics.m
% and compares the performance of dynamic array growth vs. pre-allocation.

function benchmark_metrics_opt()

    fprintf('--- Benchmark: Dynamic Array Growth vs. Pre-allocation ---\n');

    % 1. Setup Mock Data
    n_patients = 100;
    n_timepoints = 6;
    n_sets = 4;
    n_metrics_per_set = 5;

    % Simulate metric sets
    metric_sets = cell(1, n_sets);
    set_names = cell(1, n_sets);
    time_labels = {'Fx1', 'Fx2', 'Fx3', 'Fx4', 'Fx5', 'Post'};

    for s = 1:n_sets
        current_metrics = cell(1, n_metrics_per_set);
        current_names = cell(1, n_metrics_per_set);
        for m = 1:n_metrics_per_set
            % Random data [n_patients x n_timepoints]
            current_metrics{m} = rand(n_patients, n_timepoints);
            current_names{m} = sprintf('Metric_S%d_M%d', s, m);
        end
        metric_sets{s} = current_metrics;
        set_names{s} = current_names;
    end

    % Simulate patient outcomes (50% LC, 50% LF)
    lf_group = [zeros(n_patients/2, 1); ones(n_patients/2, 1)];
    valid_pts = true(n_patients, 1);

    % Number of iterations to make the benchmark measurable
    n_iterations = 100;

    fprintf('Iterations: %d\n', n_iterations);
    fprintf('Mock Data: %d Sets, %d Metrics/Set, %d Timepoints\n', n_sets, n_metrics_per_set, n_timepoints);

    % --- Approach 1: Dynamic Array Growth (Original) ---
    tic;
    for iter = 1:n_iterations
        sig_metric = {};
        sig_fraction = {};
        sig_pval = [];
        sig_mean_LC = [];
        sig_mean_LF = [];

        for s = 1:length(metric_sets)
            current_metrics = metric_sets{s};
            current_names = set_names{s};
            num_metrics = length(current_metrics);

            for m = 1:num_metrics
                metric_data = current_metrics{m};
                cols_to_plot = min(size(metric_data, 2), length(time_labels));

                for tp = 1:cols_to_plot
                    % Simulate statistical test (assume always significant for stress test)
                    p = 0.01;
                    mean_LC = 0.5;
                    mean_LF = 0.6;

                    if p < 0.05
                         % Append to storage arrays (The bottleneck)
                        sig_metric{end+1, 1} = current_names{m};
                        sig_fraction{end+1, 1} = time_labels{tp};
                        sig_pval(end+1, 1) = p;
                        sig_mean_LC(end+1, 1) = mean_LC;
                        sig_mean_LF(end+1, 1) = mean_LF;
                    end
                end
            end
        end
    end
    time_dynamic = toc;
    fprintf('Dynamic Growth Time: %.4f seconds\n', time_dynamic);

    % --- Approach 2: Pre-allocation (Optimized) ---
    tic;
    for iter = 1:n_iterations
        % Pre-calculation of maximum possible size
        % In the real code, we iterate to count, or upper bound it.
        % Here we know the exact count: n_sets * n_metrics_per_set * n_timepoints
        % But let's simulate the counting loop as part of the cost.

        total_checks = 0;
        for s = 1:length(metric_sets)
             total_checks = total_checks + length(metric_sets{s}) * length(time_labels);
        end

        % Pre-allocate
        sig_metric = cell(total_checks, 1);
        sig_fraction = cell(total_checks, 1);
        sig_pval = nan(total_checks, 1);
        sig_mean_LC = nan(total_checks, 1);
        sig_mean_LF = nan(total_checks, 1);

        sig_count = 0;

        for s = 1:length(metric_sets)
            current_metrics = metric_sets{s};
            current_names = set_names{s};
            num_metrics = length(current_metrics);

            for m = 1:num_metrics
                metric_data = current_metrics{m};
                cols_to_plot = min(size(metric_data, 2), length(time_labels));

                for tp = 1:cols_to_plot
                    % Simulate statistical test
                    p = 0.01;
                    mean_LC = 0.5;
                    mean_LF = 0.6;

                    if p < 0.05
                        sig_count = sig_count + 1;
                        sig_metric{sig_count, 1} = current_names{m};
                        sig_fraction{sig_count, 1} = time_labels{tp};
                        sig_pval(sig_count, 1) = p;
                        sig_mean_LC(sig_count, 1) = mean_LC;
                        sig_mean_LF(sig_count, 1) = mean_LF;
                    end
                end
            end
        end

        % Truncate
        if sig_count < total_checks
            sig_metric = sig_metric(1:sig_count, :);
            sig_fraction = sig_fraction(1:sig_count, :);
            sig_pval = sig_pval(1:sig_count, :);
            sig_mean_LC = sig_mean_LC(1:sig_count, :);
            sig_mean_LF = sig_mean_LF(1:sig_count, :);
        end

    end
    time_prealloc = toc;
    fprintf('Pre-allocation Time: %.4f seconds\n', time_prealloc);

    fprintf('Speedup: %.2fx\n', time_dynamic / time_prealloc);

end
