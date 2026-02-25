function [X_td, t_start, t_stop, event_td, pat_id_td, frac_td] = build_td_panel( ...
    feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, decay_half_life_months)
% build_td_panel  Converts longitudinal DWI arrays into a counting-process
%   (start–stop) panel for fitting time-dependent Cox PH models.
%
%   [X_td, t_start, t_stop, event_td, pat_id_td, frac_td] = build_td_panel( ...
%       feat_arrays, feat_names, lf_vec, total_time_vec, nTp, scan_days, decay_half_life_months)
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
%   decay_half_life_months - [scalar] biological half-life for washout (months).
%                            Default: 18.
%
%   Outputs
%   -------
%   X_td      - [n_intervals × n_feats] covariate matrix (Z-scored within
%               panel for comparable HRs across features).
%   t_start   - [n_intervals × 1] interval start times (days).
%   t_stop    - [n_intervals × 1] interval stop times (days).
%   event_td  - [n_intervals × 1] categorical event indicator:
%               0 = censored, 1 = disease progression (event of interest),
%               2 = competing risk (e.g., non-cancer death).
%               The event code is assigned to the terminal row of each patient.
%   pat_id_td - [n_intervals × 1] patient index for cluster-robust SE.
%   frac_td   - [n_intervals × 1] Measurement fraction index (1 to nTp) generating row.
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
    if nargin < 7 || isempty(decay_half_life_months)
        decay_half_life_months = 18;  % Default half-life of 18 months
    end
    scan_days = scan_days(1:nTp);  % trim to available timepoints

    % --- Collect Global Event Times for Exact Risk Set Splitting ---
    % Since Cox Partial Likelihood is only evaluated at failure times, splitting
    % at these times provides mathematically exact continuous evaluation.
    global_event_times = unique(total_time_vec(lf_vec == 1));
    global_event_times(isnan(global_event_times)) = [];
    global_event_times = sort(global_event_times);

    % Pre-allocate with a generous upper bound (n_pts * (nTp + total_events))
    max_rows = n_pts * (nTp + length(global_event_times));
    X_buf       = nan(max_rows, n_feat);
    t_start_buf = nan(max_rows, 1);
    t_stop_buf  = nan(max_rows, 1);
    event_buf   = zeros(max_rows, 1);
    pat_buf     = zeros(max_rows, 1);
    frac_buf    = nan(max_rows, 1);

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
                % 1. Acute phase: current fraction values for 90 days
                row = row + 1;
                X_buf(row, :)       = cov_row;
                t_start_buf(row)    = day_tp;
                t_stop_buf(row)     = day_tp + 90;
                event_buf(row)      = false;
                pat_buf(row)        = j;
                frac_buf(row)       = tp;
                last_valid_row      = row;
                
                % 2. Decay phase: asymptotic return to baseline
                % Define lambda for 50% decay (e.g. at 18 months (~547.5 days))
                % Average month is ~30.416 days (365/12).
                total_decay_days = decay_half_life_months * 30.416;
                decay_half_life = total_decay_days - 90; 
                if decay_half_life <= 0 
                    decay_half_life = 1; % guard against immediate washout
                end
                lambda = -log(0.5) / decay_half_life;
                
                % Implement exact risk-set splitting at global event times
                % Find all event times that fall within this decay window
                t_curr = day_tp + 90;
                interval_events = global_event_times(global_event_times > t_curr & global_event_times < day_next);
                
                % Define all unique stop-times (chunk boundaries) for this interval
                chunk_boundaries = unique([interval_events; day_next]);
                chunk_boundaries = sort(chunk_boundaries);
                
                for b_idx = 1:length(chunk_boundaries)
                    t_chunk_end = chunk_boundaries(b_idx);
                    
                    % MATLAB evaluates the covariate exactly at event time t.
                    % Therefore, we evaluate the continuous decay function at t_chunk_end.
                    exact_eval_time = t_chunk_end;
                    decay_factor = exp(-lambda * (exact_eval_time - day_tp - 90));
                    
                    decay_cov_row = nan(1, n_feat);
                    for fi = 1:n_feat
                        arr = feat_arrays{fi};
                        if size(arr, 2) >= 1
                            x_base = arr(j, 1);
                            x_frac = cov_row(fi);
                            % Asymptotic decay towards baseline. 
                            decay_cov_row(fi) = x_base + (x_frac - x_base) * decay_factor;
                        end
                    end
                    
                    row = row + 1;
                    X_buf(row, :)       = decay_cov_row;
                    t_start_buf(row)    = t_curr;
                    t_stop_buf(row)     = t_chunk_end;
                    event_buf(row)      = false;
                    pat_buf(row)        = j;
                    frac_buf(row)       = tp;
                    last_valid_row      = row;
                    
                    t_curr = t_chunk_end;
                end
            else
                % Normal interval
                row = row + 1;
                X_buf(row, :)       = cov_row;
                t_start_buf(row)    = day_tp;
                t_stop_buf(row)     = day_next;
                event_buf(row)      = false;   % will set on last valid row below
                pat_buf(row)        = j;
                frac_buf(row)       = tp;
                last_valid_row      = row;
            end
        end

        % Mark the last interval for this patient with the event code (if not censored)
        if last_valid_row > 0 && lf_j > 0
            event_buf(last_valid_row) = lf_j;
        end
    end

    % Trim to actual data
    X_td      = X_buf(1:row, :);
    t_start   = t_start_buf(1:row);
    t_stop    = t_stop_buf(1:row);
    event_td  = event_buf(1:row);
    pat_id_td = pat_buf(1:row);
    frac_td   = frac_buf(1:row);

    if row == 0
        warning('build_td_panel:noData', 'No valid intervals were generated. Check inputs.');
        return;
    end

    fprintf('  [TD Panel] %d patients → %d intervals (%d events)\n', ...
        n_pts, row, sum(event_td));
end
