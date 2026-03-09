classdef benchmark_scale_td < matlab.unittest.TestCase
% BENCHMARK_SCALE_TD  Compares two approaches for per-week feature scaling
%   in the time-dependent (TD) panel used by scale_td_panel.
%
%   The "current" approach recomputes first-occurrence indices for every
%   feature x week combination. The "optimized" approach pre-computes the
%   first-occurrence indices once per week and reuses them across features,
%   avoiding redundant unique()/find() calls.

    methods (Test)
        function testScalingApproaches(testCase)
        %TESTSCALINGAPPROACHES Generate a synthetic TD panel and compare
        %   wall-clock time for current vs optimized scaling loops.

            % Synthetic panel dimensions
            num_patients = 100;
            num_weeks = 50;
            num_features = 20;
            num_rows = num_patients * num_weeks * 2;  % ~2 rows per patient-week

            % Generate random feature matrix and patient/time metadata
            X_td_raw = rand(num_rows, num_features);
            pat_id_td = repmat((1:num_patients)', num_rows/num_patients, 1);
            t_start_td = randi([0, 50*7], num_rows, 1);
            train_pat_ids = 1:50;  % First 50 patients are training set

            % Convert days to week numbers for temporal grouping
            temporal_week_td = zeros(size(t_start_td));
            pos_mask = t_start_td > 0;
            temporal_week_td(pos_mask) = ceil(t_start_td(pos_mask) / 7);
            unique_weeks = unique(temporal_week_td(~isnan(temporal_week_td)));
            is_train_row = ismember(pat_id_td, train_pat_ids);

            num_runs = 5;  % Repeat for stable timing

            % Time the current (naive) approach: recomputes indices per feature x week
            tic;
            for iter = 1:num_runs
                benchmark_scale_td.runCurrent(X_td_raw, pat_id_td, temporal_week_td, unique_weeks, is_train_row, num_features);
            end
            t_curr = toc;

            % Time the optimized approach: pre-computes indices per week, reuses across features
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
        %RUNCURRENT Naive scaling loop: recomputes first-occurrence row
        %   indices for every (feature, week) pair.  The unique() and find()
        %   calls inside the inner loop are the bottleneck.
            for fi = 1:num_features
                for fn = 1:length(unique_weeks)
                    week_val = unique_weeks(fn);
                    week_mask = (temporal_week_td == week_val);
                    train_week_mask = week_mask & is_train_row;
                    % Find first row per patient within this week (training set only)
                    week_pat_ids = pat_id_td(train_week_mask);
                    [~, unique_idx] = unique(week_pat_ids, 'first');
                    train_week_indices = find(train_week_mask);
                    first_occurrence_indices = train_week_indices(unique_idx);
                    % Compute mean/std from one value per patient (no duplicates)
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
        %RUNOPTIMIZED Optimized scaling loop: pre-computes the first-
        %   occurrence row indices for each week once, then reuses them
        %   across all features.  Eliminates redundant unique()/find() calls.

            % Phase 1: Pre-compute per-week first-occurrence indices
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

            % Phase 2: Iterate features using cached indices
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
