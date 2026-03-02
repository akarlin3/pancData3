function [X_td_scaled] = scale_td_panel(X_td_raw, feat_names, pat_id_td, t_start_td, train_pat_ids)
% scale_td_panel  Applies timepoint-specific standard scaling to TD features.
%
%   [X_td_scaled] = scale_td_panel(X_td_raw, feat_names, pat_id_td, t_start_td, train_pat_ids)
%
%   Computes independent mu and sigma for each feature and each temporal week 
%   strictly from the rows belonging to train_pat_ids, and scales ALL rows 
%   in X_td_raw using those specific parameters.
%
%   Inputs:
%       X_td_raw      - Raw feature matrix from start-stop Cox counting process
%       feat_names    - Cell array of strings representing feature names
%       pat_id_td     - Array of patient IDs corresponding to each row
%       t_start_td    - Array of start times corresponding to each row
%       train_pat_ids - Array of patient IDs strictly belonging to the training set
%
%   Outputs:
%       X_td_scaled   - Scaled feature matrix where training set sets the mean/stdev
%
    
    n_feat = length(feat_names);
    X_td_scaled = X_td_raw;
    
    % Compute temporal week: Day 0 -> Week 0; Days 1-7 -> Week 1; etc.
    temporal_week_td = zeros(size(t_start_td));
    pos_mask = t_start_td > 0;
    temporal_week_td(pos_mask) = ceil(t_start_td(pos_mask) / 7);
    
    unique_weeks = unique(temporal_week_td(~isnan(temporal_week_td)));
    
    % Identify which rows belong to the training set
    is_train_row = ismember(pat_id_td, train_pat_ids);

    for fi = 1:n_feat
        name_fi = feat_names{fi};
        is_derivative = contains(name_fi, 'Delta', 'IgnoreCase', true) || ...
                        contains(name_fi, 'Change', 'IgnoreCase', true) || ...
                        contains(name_fi, 'pct', 'IgnoreCase', true) || ...
                        contains(name_fi, 'diff', 'IgnoreCase', true);

        for fn = 1:length(unique_weeks)
            week_val = unique_weeks(fn);
            
            % Rows for this specific temporal week
            week_mask = (temporal_week_td == week_val);
            
            % Training rows for this week
            train_week_mask = week_mask & is_train_row;
            
            mu_col = 0;
            sd_col = 1;

            if is_derivative && week_val == 0
                % Derivative features at Baseline (Week 0) are identically 0
                mu_col = 0;
                sd_col = 1;
            else
                % Extract training values for this specific week.
                % To prevent row-weighted bias, extract the first occurrence per patient.
                % Find the first index of each unique patient ID in this week's training subset
                [~, unique_idx] = unique(pat_id_td(train_week_mask), 'first');
                
                % Get the actual row indices in X_td_raw
                train_week_indices = find(train_week_mask);
                first_occurrence_indices = train_week_indices(unique_idx);
                
                % Extract the values and remove NaNs
                vals = X_td_raw(first_occurrence_indices, fi);
                unique_vals = vals(~isnan(vals));
                valid_cnt = length(unique_vals);
                
                if valid_cnt > 1
                    mu_col = mean(unique_vals);
                    sd_col = std(unique_vals);
                    if sd_col == 0
                        sd_col = 1; % Prevent division by zero
                    end
                elseif valid_cnt == 1
                    mu_col = unique_vals(1);
                    sd_col = 1;
                end
            end
            
            % Apply scaling to ALL rows (train + test) for this week
            cols_to_scale = X_td_raw(week_mask, fi);
            X_td_scaled(week_mask, fi) = (cols_to_scale - mu_col) / sd_col;
        end
    end
end
