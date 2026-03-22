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

    % ================================================================
    % PANEL CONSTRUCTION: Anderson-Gill counting-process representation
    %
    % The key insight is that DWI parameters (ADC, D, f, D*) measured at
    % each MRI scan reflect the instantaneous state of the tumor
    % microenvironment at that timepoint.  Between scans, we assume these
    % values remain constant (piecewise-constant hazard model).  This is
    % a reasonable approximation because radiation-induced changes in
    % tissue cellularity and vascularity evolve over days-to-weeks, which
    % is comparable to our inter-scan interval.
    %
    % Each patient contributes multiple rows (one per scan interval),
    % enabling the Cox model to estimate how changes in diffusion
    % parameters during treatment predict subsequent clinical outcomes
    % (local failure vs. competing risks like distant metastasis).
    % ================================================================
    n_pts  = length(lf_vec);
    n_feat = length(feat_arrays);

    % Validate that all feature arrays have consistent row counts matching
    % the number of patients.  Mismatched dimensions would cause silent
    % data corruption or cryptic indexing errors deep in the loop.
    for fi = 1:n_feat
        if size(feat_arrays{fi}, 1) ~= n_pts
            error('build_td_panel:rowMismatch', ...
                'feat_arrays{%d} has %d rows but lf_vec has %d entries.', ...
                fi, size(feat_arrays{fi}, 1), n_pts);
        end
        if size(feat_arrays{fi}, 2) ~= nTp
            error('build_td_panel:colMismatch', ...
                'feat_arrays{%d} has %d columns but nTp=%d.', ...
                fi, size(feat_arrays{fi}, 2), nTp);
        end
    end

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
    if length(scan_days) > nTp
        fprintf('  💡 build_td_panel: trimming scan_days from %d to %d entries (nTp).\n', ...
            length(scan_days), nTp);
    end
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

    % ================================================================
    % EXPONENTIAL DECAY IMPUTATION
    %
    % When a scan is missing (patient too ill for MRI, scanner down), we
    % impute the covariate value using exponential decay from the last
    % observed value back toward the patient's own baseline.  This models
    % the biological reality that acute radiation effects (edema causing
    % elevated ADC, vascular disruption reducing f) gradually resolve
    % post-treatment as fibrosis replaces acute inflammation.
    %
    % Per-feature half-lives allow different decay rates for different
    % IVIM parameters: perfusion fraction (f) and pseudo-diffusion (D*)
    % may recover faster than true diffusion (D) because vascular
    % remodeling precedes structural tissue reorganization.
    %
    % The decay target is the patient's own baseline (not a population
    % mean) to respect inter-patient heterogeneity in tumor biology.
    % ================================================================
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
        % Convert months to days (30.44 days/month average) for consistent
        % time units with scan_days (which are in days from RT start).
        hl_days      = decay_half_life_months * 30.44;
        % Exponential decay rate: lambda = ln(2) / half-life, so the
        % imputed value at time dt after last observation is:
        %   baseline + (last_obs - baseline) * exp(-lambda * dt)
        lambda_decay = log(2) ./ hl_days;
        lambda_decay(decay_half_life_months <= 0) = 0;  % disable decay for non-positive half-lives
    end

    % Find valid patients: those with known event time AND known event
    % status.  Patients with NaN event time or NaN event indicator cannot
    % contribute to the Cox partial likelihood — including them would
    % introduce undefined terms in the risk set computation.
    valid_pts = find(~isnan(total_time_vec) & ~isnan(lf_vec));

    % --- Vectorized Generation of Intervals ---
    % Pre-compute all patient-timepoint combinations and filter valid intervals
    n_valid_pts = length(valid_pts);
    
    % Generate all patient-timepoint combinations
    [tp_grid, pat_grid] = meshgrid(1:nTp, 1:n_valid_pts);
    tp_vec = tp_grid(:);
    pat_idx_vec = pat_grid(:);
    
    % Map to actual patient indices (force column to avoid orientation
    % mismatches when valid_pts, scan_days, or total_time_vec are row vectors)
    actual_pat_vec = valid_pts(pat_idx_vec);
    actual_pat_vec = actual_pat_vec(:);

    % Get scan days for all combinations
    scan_day_vec = scan_days(tp_vec);
    scan_day_vec = scan_day_vec(:);

    % Get event times for all combinations
    event_time_vec = total_time_vec(actual_pat_vec);
    event_time_vec = event_time_vec(:);
    
    % Filter for valid intervals (scan_day < event_time)
    valid_intervals = scan_day_vec < event_time_vec;
    
    tp_vec = tp_vec(valid_intervals);
    actual_pat_vec = actual_pat_vec(valid_intervals);
    scan_day_vec = scan_day_vec(valid_intervals);
    event_time_vec = event_time_vec(valid_intervals);
    
    n_intervals = length(tp_vec);
    
    if n_intervals == 0
        warning('build_td_panel:noData', 'No valid intervals were generated. Check inputs.');
        X_td = []; t_start = []; t_stop = []; event_td = []; pat_id_td = []; frac_td = [];
        return;
    end
    
    % Pre-allocate output arrays (all column vectors)
    X_td = nan(n_intervals, n_feat);
    t_start = scan_day_vec(:);
    t_stop = nan(n_intervals, 1);
    event_td = zeros(n_intervals, 1);
    pat_id_td = actual_pat_vec(:);
    frac_td = tp_map(tp_vec);
    frac_td = frac_td(:);
    
    % --- Vectorized Covariate Extraction ---
    % Extract all covariates at once using linear indexing.
    % Force tp_map(tp_vec) to column to match actual_pat_vec orientation.
    tp_col = tp_map(tp_vec);
    tp_col = tp_col(:);
    for fi = 1:n_feat
        arr = feat_arrays{fi};
        linear_idx = sub2ind(size(arr), actual_pat_vec, tp_col);
        X_td(:, fi) = arr(linear_idx);
    end
    
    % Store original (non-imputed) values for decay computation
    X_td_raw = X_td;
    
    % --- Vectorized Decay Imputation ---
    if use_decay
        % Get unique patients and their indices
        [unique_pats, ~, pat_group] = unique(actual_pat_vec);
        
        for p_idx = 1:length(unique_pats)
            j = unique_pats(p_idx);
            pat_mask = (pat_group == p_idx);
            pat_intervals = find(pat_mask);
            pat_tps = tp_vec(pat_mask);
            pat_days = scan_day_vec(pat_mask);
            
            % Sort by timepoint to ensure proper temporal ordering
            [pat_tps_sorted, sort_idx] = sort(pat_tps);
            pat_intervals_sorted = pat_intervals(sort_idx);
            pat_days_sorted = pat_days(sort_idx);
            
            % Get baseline values (first timepoint) for this patient
            baseline_X = nan(1, n_feat);
            if ~isempty(pat_tps_sorted)
                first_tp = pat_tps_sorted(1);
                for fi = 1:n_feat
                    arr = feat_arrays{fi};
                    baseline_X(fi) = arr(j, tp_map(first_tp));
                end
            end
            
            % Track last observed values and times for decay
            orig_X = nan(1, n_feat);
            orig_t = zeros(1, n_feat);
            
            for k = 1:length(pat_intervals_sorted)
                interval_idx = pat_intervals_sorted(k);
                day_tp = pat_days_sorted(k);
                cov_row = X_td_raw(interval_idx, :);
                cov_row_imputed = cov_row;
                
                % Update baseline if this is first timepoint
                if k == 1
                    baseline_X = cov_row;
                end
                
                % Apply decay imputation for missing values
                decay_mask = isnan(cov_row) & ~isnan(orig_X);
                if any(decay_mask)
                    dt_per_feat = day_tp - orig_t;
                    dt_per_feat = max(dt_per_feat, 0); % Clamp negative times
                    
                    bl = baseline_X(decay_mask);
                    bl(isnan(bl)) = orig_X(decay_mask & isnan(baseline_X));
                    
                    decay_factor = exp(-lambda_decay(decay_mask) .* dt_per_feat(decay_mask));
                    cov_row_imputed(decay_mask) = bl + (orig_X(decay_mask) - bl) .* decay_factor;
                end
                
                % Update imputed values
                X_td(interval_idx, :) = cov_row_imputed;
                
                % Update tracking for next iteration (only for observed values)
                observed = ~isnan(cov_row);
                orig_X(observed) = cov_row(observed);
                orig_t(observed) = day_tp;
            end
        end
    end
    
    % --- Vectorized Interval End Times and Event Assignment ---
    % For each interval, find the next valid timepoint or event time
    for i = 1:n_intervals
        j = actual_pat_vec(i);
        tp = tp_vec(i);
        T_j = event_time_vec(i);
        
        % Find next valid timepoint
        if tp < nTp
            % Look for next non-missing timepoint
            next_day = T_j;
            for t_next = (tp+1):nTp
                orig_t_next = tp_map(t_next);
                % Check if any feature has valid data at next timepoint
                has_data = false;
                for fi = 1:n_feat
                    arr = feat_arrays{fi};
                    if orig_t_next <= size(arr, 2) && ~isnan(arr(j, orig_t_next))
                        has_data = true;
                        break;
                    end
                end
                if has_data
                    next_day = min(scan_days(t_next), T_j);
                    break;
                end
            end
            t_stop(i) = next_day;
        else
            t_stop(i) = T_j;
        end
        
        % Assign event for terminal intervals
        if t_stop(i) == T_j
            event_td(i) = lf_vec(j);
        end
    end
    
    % Remove intervals with invalid stop times
    valid_intervals = (t_stop > t_start) & ~isnan(t_stop);
    X_td = X_td(valid_intervals, :);
    t_start = t_start(valid_intervals);
    t_stop = t_stop(valid_intervals);
    event_td = event_td(valid_intervals);
    pat_id_td = pat_id_td(valid_intervals);
    frac_td = frac_td(valid_intervals);
    
    % Remove intervals with all-NaN covariates
    valid_cov = ~all(isnan(X_td), 2);
    X_td = X_td(valid_cov, :);
    t_start = t_start(valid_cov);
    t_stop = t_stop(valid_cov);
    event_td = event_td(valid_cov);
    pat_id_td = pat_id_td(valid_cov);
    frac_td = frac_td(valid_cov);
    
    n_final = length(t_start);
    
    if n_final == 0
        warning('build_td_panel:noData', 'No valid intervals were generated after filtering. Check inputs.');
        return;
    end

    fprintf('  [TD Panel] %d patients → %d intervals (%d events of interest, %d competing)\n', ...
        numel(unique(pat_id_td)), n_final, sum(event_td == 1), sum(event_td == 2));
end