classdef benchmark_scale_td < matlab.unittest.TestCase
    % Benchmark for temporal panel scaling approaches

    methods (Test)
        function testScalingApproaches(testCase)
            num_patients = 100;
            num_weeks = 50;
            num_features = 20;
            num_rows = num_patients * num_weeks * 2;

            X_td_raw = rand(num_rows, num_features);
            pat_id_td = repmat((1:num_patients)', num_rows/num_patients, 1);
            t_start_td = randi([0, 50*7], num_rows, 1);
            train_pat_ids = 1:50;

            temporal_week_td = zeros(size(t_start_td));
            pos_mask = t_start_td > 0;
            temporal_week_td(pos_mask) = ceil(t_start_td(pos_mask) / 7);
            unique_weeks = unique(temporal_week_td(~isnan(temporal_week_td)));
            is_train_row = ismember(pat_id_td, train_pat_ids);

            num_runs = 5;

            % Current approach
            tic;
            for iter = 1:num_runs
                benchmark_scale_td.runCurrent(X_td_raw, pat_id_td, temporal_week_td, unique_weeks, is_train_row, num_features);
            end
            t_curr = toc;

            % Optimized approach
            tic;
            for iter = 1:num_runs
                benchmark_scale_td.runOptimized(X_td_raw, pat_id_td, temporal_week_td, unique_weeks, is_train_row, num_features);
            end
            t_opt = toc;

            fprintf('Current: %f s, Optimized: %f s, Speedup: %.2fx\n', t_curr, t_opt, t_curr / t_opt);
            testCase.verifyGreaterThan(t_curr, 0);
            testCase.verifyGreaterThan(t_opt, 0);
        end
    end

    methods (Static)
        function runCurrent(X_td_raw, pat_id_td, temporal_week_td, unique_weeks, is_train_row, num_features)
            for fi = 1:num_features
                for fn = 1:length(unique_weeks)
                    week_val = unique_weeks(fn);
                    week_mask = (temporal_week_td == week_val);
                    train_week_mask = week_mask & is_train_row;
                    week_pat_ids = pat_id_td(train_week_mask);
                    [~, unique_idx] = unique(week_pat_ids, 'first');
                    train_week_indices = find(train_week_mask);
                    first_occurrence_indices = train_week_indices(unique_idx);
                    vals = X_td_raw(first_occurrence_indices, fi);
                    unique_vals = vals(~isnan(vals));
                    valid_cnt = length(unique_vals);
                    if valid_cnt > 1
                        mu_col = mean(unique_vals); %#ok<NASGU>
                        sd_col = std(unique_vals); %#ok<NASGU>
                    end
                end
            end
        end

        function runOptimized(X_td_raw, pat_id_td, temporal_week_td, unique_weeks, is_train_row, num_features)
            n_weeks = length(unique_weeks);
            first_occ_indices_cell = cell(n_weeks, 1);
            for fn = 1:n_weeks
                week_val = unique_weeks(fn);
                week_mask = (temporal_week_td == week_val);
                train_week_mask = week_mask & is_train_row;
                [~, unique_idx] = unique(pat_id_td(train_week_mask), 'first');
                train_week_indices = find(train_week_mask);
                first_occ_indices_cell{fn} = train_week_indices(unique_idx);
            end
            for fi = 1:num_features
                for fn = 1:n_weeks
                    first_occurrence_indices = first_occ_indices_cell{fn};
                    vals = X_td_raw(first_occurrence_indices, fi);
                    unique_vals = vals(~isnan(vals));
                    valid_cnt = length(unique_vals);
                    if valid_cnt > 1
                        mu_col = mean(unique_vals); %#ok<NASGU>
                        sd_col = std(unique_vals); %#ok<NASGU>
                    end
                end
            end
        end
    end
end
