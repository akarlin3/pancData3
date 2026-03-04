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
%   X_td      - [n_intervals × n_feats] raw covariate matrix (NOT scaled).
%               Callers MUST apply scale_td_panel() with explicit train/test
%               patient lists to avoid cross-patient leakage.  Do NOT Z-score
%               inside this function — all standardisation is delegated to
%               scale_td_panel().
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

    % Validate scan_days is strictly increasing
    assert(all(diff(scan_days) > 0), 'build_td_panel:scanDays', 'scan_days must be strictly increasing');

    % --- Vectorized Generation of General Intervals ---
    % Preallocate arrays safely overestimating bounds
    max_rows = n_pts * nTp;
    X_buf       = nan(max_rows, n_feat);
    t_start_buf = nan(max_rows, 1);
    t_stop_buf  = nan(max_rows, 1);
    event_buf   = zeros(max_rows, 1);
    pat_buf     = zeros(max_rows, 1);
    frac_buf    = nan(max_rows, 1);

    % Exponential decay imputation parameters
    if decay_half_life_months <= 0
        warning('build_td_panel:invalidHalfLife', ...
            'decay_half_life_months = %.2f is not positive. Decay imputation will be disabled.', ...
            decay_half_life_months);
    end
    use_decay = (decay_half_life_months > 0);
    if use_decay
        hl_days      = decay_half_life_months * 30.44;
        lambda_decay = log(2) / hl_days;
    end

    % Find valid patients (T_j and lf_j not NaN)
    valid_pts = find(~isnan(total_time_vec) & ~isnan(lf_vec));

    row = 0;
    for p_idx = 1:length(valid_pts)
        j = valid_pts(p_idx);
        T_j   = total_time_vec(j);
        lf_j  = lf_vec(j);

        last_valid_row = -1;
        last_X = nan(1, n_feat);
        last_t = 0;
        % Track original (non-imputed) observations and their times for
        % decay imputation.  Using the already-imputed last_X would cause
        % compounding decay across consecutive missing scans, making values
        % decay faster than the intended biological half-life.
        orig_X = nan(1, n_feat);   % last observed (non-NaN) value per feature
        orig_t = zeros(1, n_feat); % day of that observation

        for tp = 1:nTp
            day_tp = scan_days(tp);
            if day_tp >= T_j
                break;
            end

            % Extract covariates for this patient at this timepoint
            cov_row = nan(1, n_feat);
            for fi = 1:n_feat
                arr = feat_arrays{fi};
                if tp <= size(arr, 2)
                    cov_row(fi) = arr(j, tp);
                end
            end

            % Save the raw (pre-imputation) covariate row so we can track
            % which features were genuinely observed at this timepoint.
            cov_row_raw = cov_row;

            % Apply exponential decay imputation for missing covariates.
            % Decay from the *original* observed value and its timestamp to
            % avoid compounding decay across consecutive missing scans.
            if use_decay && any(isnan(cov_row))
                decay_mask = isnan(cov_row) & ~isnan(orig_X);
                dt_per_feat = day_tp - orig_t;
                if any(dt_per_feat(decay_mask) < 0)
                    warning('build_td_panel:nonMonotonicTime', ...
                        'Negative time delta detected (patient %d, tp %d). Clamping to zero.', j, tp);
                    dt_per_feat = max(dt_per_feat, 0);
                end
                cov_row(decay_mask) = orig_X(decay_mask) .* exp(-lambda_decay * dt_per_feat(decay_mask));
            end

            if all(isnan(cov_row))
                continue;
            end

            % Compute interval end: next observed scan or event/censoring time
            if tp < nTp
                % Look ahead for the next non-missing timepoint
                day_next = T_j;
                for t_next = (tp+1):nTp
                    next_row = nan(1, n_feat);
                    for fi = 1:n_feat
                        arr = feat_arrays{fi};
                        if t_next <= size(arr, 2)
                            next_row(fi) = arr(j, t_next);
                        end
                    end
                    if ~all(isnan(next_row))
                        day_next = min(scan_days(t_next), T_j);
                        break;
                    end
                end
            else
                day_next = T_j;
            end

            if day_next <= day_tp
                continue;
            end

            is_terminal = (day_next == T_j);

            row = row + 1;
            X_buf(row, :)       = cov_row;
            t_start_buf(row)    = day_tp;
            t_stop_buf(row)     = day_next;
            pat_buf(row)        = j;
            frac_buf(row)       = tp;

            if is_terminal
                event_buf(row) = lf_j;
            else
                event_buf(row) = 0;
            end

            last_valid_row      = row;
            last_X = cov_row;
            last_t = day_tp;
            % Update per-feature original observation tracking: only for
            % features that were actually observed (not decay-imputed).
            observed = ~isnan(cov_row_raw);
            orig_X(observed) = cov_row_raw(observed);
            orig_t(observed) = day_tp;

            if is_terminal, break; end
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

    fprintf('  [TD Panel] %d patients → %d intervals (%d events of interest, %d competing)\n', ...
        n_pts, row, sum(event_td == 1), sum(event_td == 2));
end
