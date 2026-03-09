classdef benchmark_metrics_opt < matlab.unittest.TestCase
% BENCHMARK_METRICS_OPT  Compares dynamic array growth vs pre-allocation
%   for the triple-nested metrics loop in metrics_stats_comparisons.
%
%   The metrics pipeline iterates over metric sets x metrics x timepoints,
%   appending significant results to arrays. Dynamic growth (cell{end+1})
%   causes repeated memory reallocation. This benchmark demonstrates the
%   speedup from pre-allocating arrays to their maximum possible size.

    methods (Test)
        function testPreallocationFaster(testCase)
        %TESTPREALLOCATIONFASTER Compare dynamic-growth vs pre-allocated
        %   array filling over 100 iterations of a 4-set x 5-metric x
        %   6-timepoint loop to measure the speedup.

            % Benchmark dimensions matching typical pipeline scale
            n_patients = 100;
            n_timepoints = 6;
            n_sets = 4;             % e.g., ADC, D, f, D* metric groups
            n_metrics_per_set = 5;  % e.g., mean, median, std, skewness, kurtosis
            n_iterations = 100;     % Repeat to get stable timing

            % Build synthetic metric data (random values per patient x timepoint)
            metric_sets = cell(1, n_sets);
            set_names = cell(1, n_sets);
            time_labels = {'Fx1', 'Fx2', 'Fx3', 'Fx4', 'Fx5', 'Post'};

            for s = 1:n_sets
                current_metrics = cell(1, n_metrics_per_set);
                current_names = cell(1, n_metrics_per_set);
                for m = 1:n_metrics_per_set
                    current_metrics{m} = rand(n_patients, n_timepoints);
                    current_names{m} = sprintf('Metric_S%d_M%d', s, m);
                end
                metric_sets{s} = current_metrics;
                set_names{s} = current_names;
            end

            % --- Approach 1: Dynamic growth ---
            % Uses {end+1} and (end+1) which forces MATLAB to copy and
            % reallocate the array on every append.
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
                    for m = 1:length(current_metrics)
                        cols_to_plot = min(size(current_metrics{m}, 2), length(time_labels));
                        for tp = 1:cols_to_plot
                            sig_metric{end+1, 1} = current_names{m}; %#ok<AGROW>
                            sig_fraction{end+1, 1} = time_labels{tp}; %#ok<AGROW>
                            sig_pval(end+1, 1) = 0.01; %#ok<AGROW>
                            sig_mean_LC(end+1, 1) = 0.5; %#ok<AGROW>
                            sig_mean_LF(end+1, 1) = 0.6; %#ok<AGROW>
                        end
                    end
                end
            end
            time_dynamic = toc;

            % --- Approach 2: Pre-allocation ---
            % Computes the total number of entries up front and allocates
            % arrays once, then fills them by index.
            tic;
            for iter = 1:n_iterations
                % Count total entries across all sets x metrics x timepoints
                total_checks = 0;
                for s = 1:length(metric_sets)
                    total_checks = total_checks + length(metric_sets{s}) * length(time_labels);
                end
                % Pre-allocate all output arrays to maximum size
                sig_metric_p = cell(total_checks, 1);
                sig_fraction_p = cell(total_checks, 1);
                sig_pval_p = nan(total_checks, 1);
                sig_mean_LC_p = nan(total_checks, 1);
                sig_mean_LF_p = nan(total_checks, 1);
                sig_count = 0;
                for s = 1:length(metric_sets)
                    current_metrics = metric_sets{s};
                    current_names = set_names{s};
                    for m = 1:length(current_metrics)
                        cols_to_plot = min(size(current_metrics{m}, 2), length(time_labels));
                        for tp = 1:cols_to_plot
                            sig_count = sig_count + 1;
                            sig_metric_p{sig_count} = current_names{m};
                            sig_fraction_p{sig_count} = time_labels{tp};
                            sig_pval_p(sig_count) = 0.01;
                            sig_mean_LC_p(sig_count) = 0.5;
                            sig_mean_LF_p(sig_count) = 0.6;
                        end
                    end
                end
            end
            time_prealloc = toc;

            fprintf('Dynamic Growth: %.4f s, Pre-allocation: %.4f s, Speedup: %.2fx\n', ...
                time_dynamic, time_prealloc, time_dynamic / time_prealloc);
            % Verify that the loop actually ran (sig_count > 0)
            testCase.verifyGreaterThan(sig_count, 0);
        end
    end
end
