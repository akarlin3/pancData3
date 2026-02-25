function [X_td_scaled] = scale_td_panel(X_td_raw, feat_names, pat_id_td, frac_td, train_pat_ids)
% scale_td_panel  Applies timepoint-specific standard scaling to TD features.
%
%   [X_td_scaled] = scale_td_panel(X_td_raw, feat_names, pat_id_td, frac_td, train_pat_ids)
%
%   Computes independent mu and sigma for each feature and each fraction 
%   strictly from the rows belonging to train_pat_ids, and scales ALL rows 
%   in X_td_raw using those specific parameters.
    
    n_feat = length(feat_names);
    n_rows = size(X_td_raw, 1);
    X_td_scaled = X_td_raw;
    unique_fracs = unique(frac_td(~isnan(frac_td)));
    
    % Identify which rows belong to the training set
    is_train_row = ismember(pat_id_td, train_pat_ids);

    for fi = 1:n_feat
        name_fi = feat_names{fi};
        is_derivative = contains(name_fi, 'Delta', 'IgnoreCase', true) || ...
                        contains(name_fi, 'Change', 'IgnoreCase', true) || ...
                        contains(name_fi, 'pct', 'IgnoreCase', true) || ...
                        contains(name_fi, 'diff', 'IgnoreCase', true);

        for fn = 1:length(unique_fracs)
            tp_val = unique_fracs(fn);
            
            % Rows for this specific fraction
            frac_mask = (frac_td == tp_val);
            
            % Training rows for this fraction
            train_frac_mask = frac_mask & is_train_row;
            
            mu_col = 0;
            sd_col = 1;

            if is_derivative && tp_val == 1
                % Derivative features at Fx1 are identically 0
                mu_col = 0;
                sd_col = 1;
            else
                % Extract training values for this specific fraction.
                % To prevent row-weighted bias (e.g. from 90-day decay splitting), 
                % we extract the first occurrence per patient.
                train_pats_in_frac = unique(pat_id_td(train_frac_mask));
                unique_vals = zeros(length(train_pats_in_frac), 1);
                
                valid_cnt = 0;
                for p_idx = 1:length(train_pats_in_frac)
                    pid = train_pats_in_frac(p_idx);
                    % First row for this patient at this fraction
                    idx = find(train_frac_mask & (pat_id_td == pid), 1, 'first');
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
                    if sd_col == 0
                        sd_col = 1; % Prevent division by zero
                    end
                elseif valid_cnt == 1
                    mu_col = unique_vals(1);
                    sd_col = 1;
                end
            end
            
            % Apply scaling to ALL rows (train + test) for this fraction
            cols_to_scale = X_td_raw(frac_mask, fi);
            X_td_scaled(frac_mask, fi) = (cols_to_scale - mu_col) / sd_col;
        end
    end
end
