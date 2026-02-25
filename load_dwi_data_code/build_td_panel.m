function [X_td, t_start, t_stop, event_td, pat_id_td] = build_td_panel( ...
    feat_arrays, ~, lf_vec, total_time_vec, nTp, scan_days)
% build_td_panel  Converts longitudinal DWI arrays into a counting-process
%   (start–stop) panel for fitting time-dependent Cox PH models.
%
%   [X_td, t_start, t_stop, event_td, pat_id_td] = build_td_panel( ...
%       feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days)
%
%   Inputs
%   ------
%   feat_arrays  - cell array of [n_pts × nTp] matrices, one per covariate.
%                  NaN values are permitted; rows with all-NaN covariates are
%                  silently excluded.
%   feat_names   - cell array of strings naming each covariate (for display).
%   lf_vec       - [n_pts × 1] binary event indicator (1 = local failure).
%   total_time_vec - [n_pts × 1] time-to-event / censoring (days).
%   nTp          - number of timepoints (columns in feat_arrays).
%   scan_days    - [1 × nTp] approximate day of each scan relative to Fx1.
%                  Default: [0 5 10 15 20 90] (Fx1 through Post-RT).
%
%   Outputs
%   -------
%   X_td      - [n_intervals × n_feats] covariate matrix (Z-scored within
%               panel for comparable HRs across features).
%   t_start   - [n_intervals × 1] interval start times (days).
%   t_stop    - [n_intervals × 1] interval stop times (days).
%   event_td  - [n_intervals × 1] logical: 1 on the terminal row of an LF
%               patient when the scan falls before the event time.
%   pat_id_td - [n_intervals × 1] patient index for cluster-robust SE.
%
%   Algorithm (counting-process representation)
%   -------------------------------------------
%   For patient j, for each scan tp where scan_days(tp) < total_time_vec(j):
%     t_start = scan_days(tp)
%     t_stop  = min(scan_days(tp+1), total_time_vec(j))  [or total_time for last scan]
%     event   = lf_vec(j) AND (this is the last valid interval)
%
%   Requires: Statistics and Machine Learning Toolbox (for z-scoring only;
%   the core construction uses only base MATLAB).

    n_pts  = length(lf_vec);
    n_feat = length(feat_arrays);

    % --- Default scan-day schedule (fraction days from RT start) ---
    if nargin < 6 || isempty(scan_days)
        scan_days = [0, 5, 10, 15, 20, 90];
    end
    scan_days = scan_days(1:nTp);  % trim to available timepoints

    % Pre-allocate with a generous upper bound (n_pts * nTp * 2 rows due to splitting)
    max_rows = n_pts * nTp * 2;
    X_buf       = nan(max_rows, n_feat);
    t_start_buf = nan(max_rows, 1);
    t_stop_buf  = nan(max_rows, 1);
    event_buf   = false(max_rows, 1);
    pat_buf     = zeros(max_rows, 1);

    row = 0;
    for j = 1:n_pts
        T_j   = total_time_vec(j);   % event/censoring time for this patient
        lf_j  = lf_vec(j);           % 1 = local failure

        if isnan(T_j) || isnan(lf_j)
            continue;
        end

        % Collect valid scan intervals that START before the event time
        last_valid_row = -1;
        for tp = 1:nTp
            day_tp = scan_days(tp);

            if day_tp >= T_j
                break;  % scan is at or after the event/censoring — stop
            end

            % Covariate values at this timepoint
            cov_row = nan(1, n_feat);
            for fi = 1:n_feat
                arr = feat_arrays{fi};
                if tp <= size(arr, 2)
                    cov_row(fi) = arr(j, tp);
                end
            end

            % Skip this interval if ALL covariates are missing
            if all(isnan(cov_row))
                continue;
            end

            % Interval stop: next scan or event/censoring time, whichever is earlier
            if tp < nTp
                day_next = min(scan_days(tp + 1), T_j);
            else
                day_next = T_j;
            end

            % Safety: start must be strictly less than stop
            if day_next <= day_tp
                continue;
            end

            % Implement 90-day biological validity window
            if day_next - day_tp > 90
                % First interval: current fraction values for 90 days
                row = row + 1;
                X_buf(row, :)       = cov_row;
                t_start_buf(row)    = day_tp;
                t_stop_buf(row)     = day_tp + 90;
                event_buf(row)      = false;
                pat_buf(row)        = j;
                last_valid_row      = row;
                
                % Second interval: remainder of time reverting to baseline (Fraction 1)
                baseline_cov_row = nan(1, n_feat);
                for fi = 1:n_feat
                    arr = feat_arrays{fi};
                    if size(arr, 2) >= 1
                        baseline_cov_row(fi) = arr(j, 1);
                    end
                end
                
                row = row + 1;
                X_buf(row, :)       = baseline_cov_row;
                t_start_buf(row)    = day_tp + 90;
                t_stop_buf(row)     = day_next;
                event_buf(row)      = false;
                pat_buf(row)        = j;
                last_valid_row      = row;
            else
                % Normal interval
                row = row + 1;
                X_buf(row, :)       = cov_row;
                t_start_buf(row)    = day_tp;
                t_stop_buf(row)     = day_next;
                event_buf(row)      = false;   % will set on last valid row below
                pat_buf(row)        = j;
                last_valid_row      = row;
            end
        end

        % Mark the last interval for this patient as the event row (if LF)
        if last_valid_row > 0 && lf_j == 1
            event_buf(last_valid_row) = true;
        end
    end

    % Trim to actual data
    X_td      = X_buf(1:row, :);
    t_start   = t_start_buf(1:row);
    t_stop    = t_stop_buf(1:row);
    event_td  = event_buf(1:row);
    pat_id_td = pat_buf(1:row);

    if row == 0
        warning('build_td_panel:noData', 'No valid intervals were generated. Check inputs.');
        return;
    end

    % Z-score each covariate column (using available, non-NaN data)
    % This makes HR coefficients comparable across features with different scales.
    for fi = 1:n_feat
        col     = X_td(:, fi);
        mu_col  = mean(col, 'omitnan');
        sd_col  = std(col,  'omitnan');
        if sd_col > 0
            X_td(:, fi) = (col - mu_col) / sd_col;
        end
    end

    fprintf('  [TD Panel] %d patients → %d intervals (%d events)\n', ...
        n_pts, row, sum(event_td));
end
