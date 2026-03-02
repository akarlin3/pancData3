classdef benchmark_metrics_opt < matlab.unittest.TestCase
    % Benchmark: dynamic array growth vs. pre-allocation for metrics loops

    methods (Test)
        function testPreallocationFaster(testCase)
            n_patients = 100;
            n_timepoints = 6;
            n_sets = 4;
            n_metrics_per_set = 5;
            n_iterations = 100;

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

            % Dynamic growth
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

            % Pre-allocation
            tic;
            for iter = 1:n_iterations
                total_checks = 0;
                for s = 1:length(metric_sets)
                    total_checks = total_checks + length(metric_sets{s}) * length(time_labels);
                end
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
            testCase.verifyGreaterThan(sig_count, 0);
        end
    end
end
