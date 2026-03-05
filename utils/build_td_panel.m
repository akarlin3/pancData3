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
%   lf_vec       - [n_pts × 1] event indicator (0 = censored, 1 = local
%                  failure, 2 = competing risk).
%   total_time_vec - [n_pts × 1] time-to-event / censoring (days).
%   nTp          - number of timepoints (columns in feat_arrays).
%   scan_days    - [1 × nTp] approximate day of each scan relative to Fx1.
%                  Default: [0 5 10 15 20 90] (Fx1 through Post-RT).
%   decay_half_life_months - [scalar or 1×n_feat vector] biological half-life
%                            for washout (months).  When a vector is provided,
%                            each covariate uses its own half-life (e.g.,
%                            faster decay for perfusion parameters).
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
%   frac_td   - [n_intervals × 1] Original measurement fraction index
%               (column index in feat_arrays before NaN scan-day removal).
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
%
%   Analytical Background:
%   ----------------------
%   Standard Cox proportional hazards regression assumes covariates are
%   fixed at baseline.  In pancreatic RT, DWI parameters (ADC, D, f, D*)
%   change throughout the 5-week treatment course as radiation induces
%   tumor cell death, edema, and vascular disruption.  Ignoring these
%   intra-treatment changes would bias hazard ratio estimates toward
%   baseline tumor biology and miss the prognostic signal carried by
%   treatment-induced diffusion changes.
%
%   The counting-process (start-stop) representation splits each patient's
%   follow-up into intervals aligned with scan dates.  Within each
%   interval, the covariates are held constant at the most recent scan
%   values — a piecewise-constant approximation to the continuously
%   evolving tumor microenvironment.  This is the standard Anderson-Gill
%   extension for time-dependent covariates in survival analysis.
%
%   The event indicator uses a competing-risk encoding (0/1/2) rather
%   than simple binary, because pancreatic cancer patients face
%   substantial non-cancer mortality that must be modeled separately
%   to avoid upward bias in cancer-specific hazard estimates.

    n_pts  = length(lf_vec);
    n_feat = length(feat_arrays);

    % --- Default scan-day schedule (fraction days from RT start) ---
    % The default [0, 5, 10, 15, 20, 90] corresponds to weekly scans
    % during a standard 5-fraction SBRT course (Fx1 on day 0 through Fx5
    % on day 20) plus a ~3-month post-RT follow-up scan (day 90).  This
    % schedule captures the acute radiation response (first 3 weeks) and
    % the early subacute recovery phase.
    if nargin < 6 || isempty(scan_days)
        scan_days = [0, 5, 10, 15, 20, 90];
    end
    if nargin < 7 || isempty(decay_half_life_months)
        % 18-month half-life: empirical estimate for how quickly the
        % radiation-induced diffusion signal reverts toward pre-treatment
        % baseline.  Used for exponential decay imputation when scans are
        % missing — models the biological washout of treatment effect on
        % tissue microstructure (e.g., fibrosis replacing acute edema).
        decay_half_life_months = 18;
    end
    % Trim scan_days to the number of timepoints actually available in the
    % data.  Some cohorts may have fewer scans than the full schedule.
    scan_days = scan_days(1:nTp);

    % Remove NaN entries (fractions with no valid scan dates) while
    % preserving the mapping between scan_days indices and feature
    % array columns via tp_map.  NaN scan days arise when a fraction
    % was scheduled but the MRI scan was not acquired (e.g., patient
    % too sick, scanner unavailable).  The tp_map allows us to index
    % back into the original feat_arrays columns after NaN removal.
    valid_sd = ~isnan(scan_days);
    tp_map   = find(valid_sd);       % tp_map(k) = original column in feat_arrays
    scan_days = scan_days(valid_sd);
    nTp = length(scan_days);

    % Validate scan_days is strictly increasing — the counting-process
    % representation requires non-overlapping intervals with t_start <
    % t_stop.  Non-monotonic scan days would create degenerate intervals
    % that violate the Cox model likelihood assumptions.
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

    % Exponential decay imputation parameters.
    % Support per-feature half-lives: expand scalar to vector if needed.
    if isscalar(decay_half_life_months)
        decay_half_life_months = repmat(decay_half_life_months, 1, n_feat);
    end
    if length(decay_half_life_months) ~= n_feat
        warning('build_td_panel:halfLifeLength', ...
            'decay_half_life_months has %d elements but %d features. Using first element for all.', ...
            length(decay_half_life_months), n_feat);
        decay_half_life_months = repmat(decay_half_life_months(1), 1, n_feat);
    end
    use_decay = any(decay_half_life_months > 0);
    if use_decay
        hl_days      = decay_half_life_months * 30.44;
        lambda_decay = log(2) ./ hl_days;
        lambda_decay(decay_half_life_months <= 0) = 0;  % disable decay for non-positive half-lives
    end

    % Find valid patients (T_j and lf_j not NaN)
    valid_pts = find(~isnan(total_time_vec) & ~isnan(lf_vec));

    row = 0;
    for p_idx = 1:length(valid_pts)
        j = valid_pts(p_idx);
        T_j   = total_time_vec(j);
        lf_j  = lf_vec(j);

        % Track original (non-imputed) observations and their times for
        % decay imputation.  Using the already-imputed last_X would cause
        % compounding decay across consecutive missing scans, making values
        % decay faster than the intended biological half-life.
        orig_X = nan(1, n_feat);   % last observed (non-NaN) value per feature
        orig_t = zeros(1, n_feat); % day of that observation
        baseline_X = nan(1, n_feat);  % patient's baseline (tp=1) values for decay target

        for tp = 1:nTp
            day_tp = scan_days(tp);
            if day_tp >= T_j
                break;
            end

            % Extract covariates for this patient at this timepoint
            % tp_map translates the scan_days index back to the
            % original feature-array column (handles NaN-dropped fractions).
            orig_tp = tp_map(tp);
            cov_row = nan(1, n_feat);
            for fi = 1:n_feat
                arr = feat_arrays{fi};
                if orig_tp <= size(arr, 2)
                    cov_row(fi) = arr(j, orig_tp);
                end
            end

            % Save the raw (pre-imputation) covariate row so we can track
            % which features were genuinely observed at this timepoint.
            cov_row_raw = cov_row;

            % Record baseline values at the first timepoint for use as
            % the decay asymptote.
            if tp == 1
                baseline_X = cov_row;
            end

            % Apply exponential decay imputation for missing covariates.
            % Decay from the *original* observed value toward the patient's
            % baseline (not zero) to avoid artificially low imputed values
            % for parameters with physiological floor values (ADC, D).
            if use_decay && any(isnan(cov_row))
                decay_mask = isnan(cov_row) & ~isnan(orig_X);
                dt_per_feat = day_tp - orig_t;
                if any(dt_per_feat(decay_mask) < 0)
                    warning('build_td_panel:nonMonotonicTime', ...
                        'Negative time delta detected (patient %d, tp %d). Clamping affected features to zero.', j, tp);
                    dt_per_feat(decay_mask) = max(dt_per_feat(decay_mask), 0);
                end
                % Decay toward baseline: baseline + (last_obs - baseline) * exp(-lambda * dt)
                bl = baseline_X(decay_mask);
                % When baseline is NaN (patient missing first scan for this
                % feature), fall back to the last observed value itself
                % (LOCF — no decay).  The previous fallback to 0 produced
                % non-physiological imputed values for parameters like ADC
                % and D which have positive physiological floor values.
                bl(isnan(bl)) = orig_X(decay_mask & isnan(baseline_X));
                decay_factor = exp(-lambda_decay(decay_mask) .* dt_per_feat(decay_mask));
                cov_row(decay_mask) = bl + (orig_X(decay_mask) - bl) .* decay_factor;
            end

            if all(isnan(cov_row))
                continue;
            end

            % Compute interval end: next observed scan or event/censoring time
            if tp < nTp
                % Look ahead for the next non-missing timepoint
                day_next = T_j;
                for t_next = (tp+1):nTp
                    orig_t_next = tp_map(t_next);
                    next_row = nan(1, n_feat);
                    for fi = 1:n_feat
                        arr = feat_arrays{fi};
                        if orig_t_next <= size(arr, 2)
                            next_row(fi) = arr(j, orig_t_next);
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
            frac_buf(row)       = orig_tp;

            if is_terminal
                event_buf(row) = lf_j;
            else
                event_buf(row) = 0;
            end

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
        numel(unique(pat_id_td)), row, sum(event_td == 1), sum(event_td == 2));
end
