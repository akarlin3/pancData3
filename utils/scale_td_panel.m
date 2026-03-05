function [X_td_scaled] = scale_td_panel(X_td_raw, feat_names, pat_id_td, t_start_td, train_pat_ids, scaling_mode)
% scale_td_panel  Applies standard scaling to TD features.
%
%   [X_td_scaled] = scale_td_panel(X_td_raw, feat_names, pat_id_td, t_start_td, train_pat_ids)
%   [X_td_scaled] = scale_td_panel(X_td_raw, feat_names, pat_id_td, t_start_td, train_pat_ids, scaling_mode)
%
%   Computes mu and sigma for each feature strictly from rows belonging to
%   train_pat_ids, and scales ALL rows in X_td_raw using those parameters.
%
%   Inputs:
%       X_td_raw      - Raw feature matrix from start-stop Cox counting process
%       feat_names    - Cell array of strings representing feature names
%       pat_id_td     - Array of patient IDs corresponding to each row
%       t_start_td    - Array of start times corresponding to each row
%       train_pat_ids - Array of patient IDs strictly belonging to the training set
%       scaling_mode  - (Optional) 'per_week' (default) computes mu/sigma per
%                       temporal week; 'baseline' uses only t_start==0 rows.
%
%   Outputs:
%       X_td_scaled   - Scaled feature matrix where training set sets the mean/stdev
%

    if nargin < 6 || isempty(scaling_mode)
        scaling_mode = 'per_week';
    end

    n_feat = length(feat_names);
    X_td_scaled = X_td_raw;

    % Identify which rows belong to the training set
    is_train_row = ismember(pat_id_td, train_pat_ids);

    if strcmp(scaling_mode, 'baseline')
        % --- Baseline-only scaling: use t_start==0 training rows ---
        is_train_base_mask = is_train_row & (t_start_td == 0);

        for fi = 1:n_feat
            is_train_base = is_train_base_mask & ~isnan(X_td_raw(:, fi));
            base_vals = X_td_raw(is_train_base, fi);

            if isempty(base_vals)
                warning('scale_td_panel:noBaseline', ...
                    'No baseline (t_start==0) training data for feature %s. Using all training rows.', feat_names{fi});
                mu  = mean(X_td_raw(is_train_row, fi), 'omitnan');
                sig = std(X_td_raw(is_train_row, fi), 0, 'omitnan');
            else
                mu  = mean(base_vals);
                sig = std(base_vals);
            end

            if sig == 0 || isnan(sig)
                % Clamp sigma to 1 so that values are centered (mean-
                % subtracted) without rescaling.  This preserves distance
                % information for test-set observations that may differ
                % from the constant training value.
                sig = 1;
            end
            X_td_scaled(:, fi) = (X_td_raw(:, fi) - mu) ./ sig;
        end
    else
        % --- Per-week scaling (default): independent mu/sigma per temporal week ---
        % Compute temporal week: Day 0 -> Week 0; Days 1-7 -> Week 1; etc.
        temporal_week_td = zeros(size(t_start_td));
        pos_mask = t_start_td > 0;
        temporal_week_td(pos_mask) = ceil(t_start_td(pos_mask) / 7);

        % Discover weeks from training rows only to avoid information
        % leakage (knowing which weeks exist in test data).  Test rows
        % with weeks not seen in training are handled after the loop.
        unique_weeks = unique(temporal_week_td(is_train_row & ~isnan(temporal_week_td)));

        % Pre-compute week masks and first occurrence indices outside the feature loop
        n_weeks = length(unique_weeks);
        week_masks_cell = cell(n_weeks, 1);
        first_occ_indices_cell = cell(n_weeks, 1);

        for fn = 1:n_weeks
            week_val = unique_weeks(fn);
            week_mask = (temporal_week_td == week_val);
            train_week_mask = week_mask & is_train_row;
            week_masks_cell{fn} = week_mask;

            [~, unique_idx] = unique(pat_id_td(train_week_mask), 'first');
            train_week_indices = find(train_week_mask);
            first_occ_indices_cell{fn} = train_week_indices(unique_idx);
        end

        for fi = 1:n_feat
            name_fi = feat_names{fi};
            % Use word-boundary matching (\b) for 'diff' to avoid false
            % positives on feature names like 'diffusion_coefficient'.
            is_derivative = contains(name_fi, 'Delta', 'IgnoreCase', true) || ...
                            contains(name_fi, 'Change', 'IgnoreCase', true) || ...
                            contains(name_fi, 'pct', 'IgnoreCase', true) || ...
                            ~isempty(regexp(name_fi, '\bdiff\b', 'once', 'ignorecase'));

            for fn = 1:n_weeks
                week_val = unique_weeks(fn);
                week_mask = week_masks_cell{fn};

                mu_col = 0;
                sd_col = 0;  % default: zero out (no valid training data for this week)

                if is_derivative && week_val == 0
                    % Derivatives at baseline are structurally zero.  Use
                    % first-occurrence-per-patient training-set statistics
                    % to avoid over-weighting patients with duplicate rows.
                    first_occ = first_occ_indices_cell{fn};
                    bl_vals = X_td_raw(first_occ, fi);
                    bl_vals = bl_vals(~isnan(bl_vals));
                    if length(bl_vals) > 1
                        mu_col = mean(bl_vals);
                        sd_col = std(bl_vals);
                        % sd_col == 0 → zero-variance; clamped to 1 by guard below
                    end
                    % else: keep defaults mu_col=0, sd_col=0 (clamped to 1 below)
                else
                    first_occurrence_indices = first_occ_indices_cell{fn};
                    vals = X_td_raw(first_occurrence_indices, fi);
                    unique_vals = vals(~isnan(vals));
                    valid_cnt = length(unique_vals);

                    if valid_cnt > 1
                        mu_col = mean(unique_vals);
                        sd_col = std(unique_vals);
                        % sd_col == 0 → zero-variance; clamped to 1
                        % by guard below, consistent with baseline mode
                    elseif valid_cnt == 1
                        mu_col = unique_vals(1);
                        sd_col = 1;  % single value → clamp sigma to 1 (center only)
                    end
                end

                cols_to_scale = X_td_raw(week_mask, fi);
                if sd_col == 0
                    % Clamp sigma to 1: center without rescaling to
                    % preserve distance information for test-set
                    % observations.
                    sd_col = 1;
                end
                X_td_scaled(week_mask, fi) = (cols_to_scale - mu_col) ./ sd_col;
            end

            % Handle test-only weeks not seen in training: use the
            % nearest training week's statistics (already computed above
            % in the per-week loop, stored in mu/sd arrays below).
        end

        % Scale any rows whose temporal week was not in the training set
        % (these remain at their raw values after the loop above).
        all_weeks = unique(temporal_week_td(~isnan(temporal_week_td)));
        unseen_weeks = setdiff(all_weeks, unique_weeks);
        if ~isempty(unseen_weeks) && ~isempty(unique_weeks)
            for fi = 1:n_feat
                for uw = 1:length(unseen_weeks)
                    wk = unseen_weeks(uw);
                    % Find nearest training week
                    [~, nearest_idx] = min(abs(unique_weeks - wk));
                    nearest_wk = unique_weeks(nearest_idx);
                    nearest_mask = (temporal_week_td == nearest_wk) & is_train_row;

                    % Compute mu/sd from nearest training week
                    [~, u_idx] = unique(pat_id_td(nearest_mask), 'first');
                    near_indices = find(nearest_mask);
                    near_first = near_indices(u_idx);
                    vals = X_td_raw(near_first, fi);
                    vals = vals(~isnan(vals));
                    if length(vals) > 1
                        mu_nn = mean(vals);
                        sd_nn = std(vals);
                        if sd_nn == 0, sd_nn = 1; end
                    elseif length(vals) == 1
                        mu_nn = vals(1);
                        sd_nn = 1;
                    else
                        mu_nn = 0;
                        sd_nn = 1;
                    end

                    unseen_mask = (temporal_week_td == wk);
                    X_td_scaled(unseen_mask, fi) = (X_td_raw(unseen_mask, fi) - mu_nn) ./ sd_nn;
                end
            end
        end
    end
end
