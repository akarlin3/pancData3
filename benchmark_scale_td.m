function benchmark_scale_td()
    % Set up synthetic data
    num_patients = 100;
    num_weeks = 50;
    num_features = 20;
    num_rows = num_patients * num_weeks * 2; % Each patient has 2 entries per week

    X_td_raw = rand(num_rows, num_features);

    % Mock the behavior locally without contains
    pat_id_td = repmat((1:num_patients)', num_rows/num_patients, 1);
    t_start_td = randi([0, 50*7], num_rows, 1);
    train_pat_ids = 1:50;

    temporal_week_td = zeros(size(t_start_td));
    pos_mask = t_start_td > 0;
    temporal_week_td(pos_mask) = ceil(t_start_td(pos_mask) / 7);
    unique_weeks = unique(temporal_week_td(~isnan(temporal_week_td)));
    is_train_row = ismember(pat_id_td, train_pat_ids);

    disp('Warming up Original...');
    tic;
    for iter = 1:5
        run_original(X_td_raw, pat_id_td, temporal_week_td, unique_weeks, is_train_row, num_features);
    end
    toc;

    disp('Warming up Optimized...');
    tic;
    for iter = 1:5
        run_optimized(X_td_raw, pat_id_td, temporal_week_td, unique_weeks, is_train_row, num_features);
    end
    toc;

    num_runs = 20;

    disp('Running Original...');
    tic;
    for iter = 1:num_runs
        run_original(X_td_raw, pat_id_td, temporal_week_td, unique_weeks, is_train_row, num_features);
    end
    t_orig = toc;

    disp('Running Optimized...');
    tic;
    for iter = 1:num_runs
        run_optimized(X_td_raw, pat_id_td, temporal_week_td, unique_weeks, is_train_row, num_features);
    end
    t_opt = toc;

    fprintf('Original: %f seconds\n', t_orig);
    fprintf('Optimized: %f seconds\n', t_opt);
    fprintf('Speedup: %fx\n', t_orig / t_opt);

end

function run_original(X_td_raw, pat_id_td, temporal_week_td, unique_weeks, is_train_row, num_features)
    X_td_scaled = X_td_raw;
    for fi = 1:num_features
        for fn = 1:length(unique_weeks)
            week_val = unique_weeks(fn);
            week_mask = (temporal_week_td == week_val);
            train_week_mask = week_mask & is_train_row;

            train_pats_in_week = unique(pat_id_td(train_week_mask));
            unique_vals = zeros(length(train_pats_in_week), 1);

            valid_cnt = 0;
            for p_idx = 1:length(train_pats_in_week)
                pid = train_pats_in_week(p_idx);
                % First row for this patient at this week
                idx = find(train_week_mask & (pat_id_td == pid), 1, 'first');
                val = X_td_raw(idx, fi);
                if ~isnan(val)
                    valid_cnt = valid_cnt + 1;
                    unique_vals(valid_cnt) = val;
                end
            end

            unique_vals = unique_vals(1:valid_cnt);
            if valid_cnt > 1
                mu_col = mean(unique_vals);
                sd_col = std(unique_vals);
            else
                mu_col = 0; sd_col = 1;
            end
        end
    end
end

function run_optimized(X_td_raw, pat_id_td, temporal_week_td, unique_weeks, is_train_row, num_features)
    X_td_scaled = X_td_raw;
    for fi = 1:num_features
        for fn = 1:length(unique_weeks)
            week_val = unique_weeks(fn);
            week_mask = (temporal_week_td == week_val);
            train_week_mask = week_mask & is_train_row;

            % Get patient IDs for training rows in this week
            week_pat_ids = pat_id_td(train_week_mask);

            % Find the first index of each unique patient ID in this subset
            [~, unique_idx] = unique(week_pat_ids, 'first');

            % Get the actual indices in the main arrays corresponding to these subset indices
            train_week_indices = find(train_week_mask);
            first_occurrence_indices = train_week_indices(unique_idx);

            % Extract the values using these indices
            vals = X_td_raw(first_occurrence_indices, fi);
            unique_vals = vals(~isnan(vals));
            valid_cnt = length(unique_vals);

            if valid_cnt > 1
                mu_col = mean(unique_vals);
                sd_col = std(unique_vals);
            else
                mu_col = 0; sd_col = 1;
            end
        end
    end
end
