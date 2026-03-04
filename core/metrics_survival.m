function metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label, m_gtv_vol, output_folder)
% METRICS_SURVIVAL — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 5/5 of the metrics step. Fits a Time-Dependent Cox Proportional Hazards
% model with dynamic covariate updating.
%
% Inputs:
%   valid_pts         - Logical mask of patients mapped to LF/LC groups
%   ADC_abs, D_abs, f_abs, Dstar_abs - Baseline covariate values matrices
%   m_lf              - Local failure clinical grouping statuses
%   m_total_time      - Time-to-local-failure (events)
%   m_total_follow_up_time - Time-to-censoring for event-free subjects
%   nTp               - Counter of number of timepoints max
%   fx_label          - Fraction labels used in logging
%   dtype_label       - DWI type name used in output
%   m_gtv_vol         - (Optional) GTV volume matrix (patients x fractions)
%   output_folder     - (Optional) Directory for diary output
%
% Outputs:
%   None. Outputs printed to console (HR and p-value tables).
%

% Handle optional arguments for backward compatibility
if nargin < 12, m_gtv_vol = []; end
if nargin < 13, output_folder = ''; end

fprintf('\n--- TIME-DEPENDENT COX PH MODEL (Counting Process) ---\n');

% Diary: capture console output to output_folder
if ~isempty(output_folder)
    diary_file = fullfile(output_folder, ['metrics_survival_output_' dtype_label '.txt']);
    if exist(diary_file, 'file'), delete(diary_file); end
    diary(diary_file);
end

td_scan_days = [0, 5, 10, 15, 20, 90];   % update if exact scan dates are available

% Covariates: all four IVIM parameters (absolute, all fractions)
td_feat_arrays = { ADC_abs(valid_pts,:), D_abs(valid_pts,:), ...
                   f_abs(valid_pts,:),   Dstar_abs(valid_pts,:) };
td_feat_names  = {'ADC', 'D', 'f', 'D*'};

% Include baseline GTV volume as a time-constant confounder when available.
% Tumor volume is a known prognostic factor; omitting it risks confounding
% imaging biomarker effect estimates.
has_vol = ~isempty(m_gtv_vol) && any(isfinite(m_gtv_vol(valid_pts, 1)));
if has_vol
    vol_baseline = m_gtv_vol(valid_pts, 1);
    % Replicate baseline volume across timepoints (time-constant covariate)
    vol_rep = repmat(vol_baseline, 1, nTp);
    td_feat_arrays{end+1} = vol_rep;
    td_feat_names{end+1}  = 'GTVvol';
    fprintf('  Including baseline GTV volume as a covariate.\n');
end

td_n_feat      = numel(td_feat_arrays);

td_lf       = m_lf(valid_pts);
td_tot_time = m_total_time(valid_pts);

% Censored patients use follow-up time; events use time-to-event
cens_mask_td = (td_lf == 0) & ~isnan(m_total_follow_up_time(valid_pts));
td_tot_time(cens_mask_td) = m_total_follow_up_time(valid_pts & (m_lf(:)==0) & ~isnan(m_total_follow_up_time(:)));

[X_td_def, t_start_td_def, t_stop_td_def, event_td_def, pat_id_td_def, frac_td_def] = ...
    build_td_panel(td_feat_arrays, td_feat_names, td_lf, td_tot_time, nTp, td_scan_days, 18);

% NOTE: Scaling is deferred until after landmark subsetting (below) to
% prevent pre-landmark covariate distributions from contaminating the
% statistics used to standardise the post-landmark analysis set.
td_ok = (sum(event_td_def) >= 3) && (size(X_td_def, 1) > td_n_feat + 1);

half_life_grid = [3, 6, 12, 18, 24];
n_halflife = length(half_life_grid);
td_panels = cell(n_halflife, 1);
for hl_idx = 1:n_halflife
    text_progress_bar(hl_idx, n_halflife, 'Half-life sensitivity');
    [X_i, t_start_i, t_stop_i, ev_i, pid_i, frac_i] = ...
        build_td_panel(td_feat_arrays, td_feat_names, td_lf, td_tot_time, nTp, td_scan_days, half_life_grid(hl_idx));
    td_panels{hl_idx} = struct('X', X_i, 't_start', t_start_i, 't_stop', t_stop_i, 'event', ev_i, 'pat_id', pid_i, 'frac', frac_i);
end

if ~td_ok
    fprintf('  Insufficient events (%d) or intervals for time-dependent Cox model.\n', sum(event_td_def));
    return;
end

X_td = X_td_def; t_start_td = t_start_td_def; t_stop_td = t_stop_td_def; event_td = event_td_def; pat_id_td = pat_id_td_def; frac_td = frac_td_def;

% ---- Landmark analysis: discard intervals before end-of-RT -----------
% Using covariates from scans that occurred AFTER the event they predict
% violates the time-dependent Cox assumption.  Landmark subsetting keeps
% only patients still at risk after the last treatment fraction, using
% covariates measured at or before the landmark.
landmark_day = td_scan_days(min(5, length(td_scan_days)));  % end of RT (Fx5, day 20)
lm_keep = (t_start_td >= landmark_day);
if any(lm_keep)
    n_before = length(t_start_td);
    X_td        = X_td(lm_keep, :);
    t_start_td  = t_start_td(lm_keep);
    t_stop_td   = t_stop_td(lm_keep);
    event_td    = event_td(lm_keep);
    pat_id_td   = pat_id_td(lm_keep);
    frac_td     = frac_td(lm_keep);
    fprintf('  Landmark at day %d: %d → %d intervals (%d patients, %d events)\n', ...
        landmark_day, n_before, sum(lm_keep), numel(unique(pat_id_td)), sum(event_td > 0));
end

% Scale AFTER landmark subsetting so that standardisation statistics are
% computed exclusively on post-landmark intervals, preventing pre-landmark
% covariate distributions from leaking into the primary analysis.
X_td_global = scale_td_panel(X_td, td_feat_names, pat_id_td, t_start_td, unique(pat_id_td));

% ---- Compute IPCW weights for informative censoring correction --------
% Cause-Specific Hazards: competing events (event=2) are treated as
% censored when modelling the cause of interest (event=1).
event_td_csh = event_td;
event_td_csh(event_td_csh == 2) = 0;  % CSH: competing risks → censored

% IPCW corrects for informative *administrative* censoring only.
% Competing events are NOT uninformative censoring — including them in
% the censoring model violates the IPCW independence assumption and can
% bias hazard ratios in either direction (Robins & Finkelstein 2000).
% We therefore exclude competing-event rows from the IPCW censoring
% model and only model the probability of true administrative censoring.
ipcw_weights = ones(size(event_td));   % default: unweighted
has_admin_cens = any(event_td == 0);
if has_admin_cens
    try
        % Restrict IPCW model to rows that are NOT competing events.
        % In this subset: event_td==0 is admin censored, event_td==1 is the event.
        not_competing = (event_td ~= 2);
        is_admin_cens_subset = double(event_td(not_competing) == 0);
        T_cens = [t_start_td(not_competing), t_stop_td(not_competing)];
        X_cens_subset = X_td_global(not_competing, :);
        w_ipcw = warning('off', 'all');
        [b_cens, ~, ~, stats_cens] = coxphfit(X_cens_subset, T_cens, ...
            'Censoring', (is_admin_cens_subset == 0), 'Ties', 'breslow');
        warning(w_ipcw);

        % Compute survival function of censoring (Kaplan-Meier-like via Cox)
        % Use the censoring model coefficients (fit on non-competing subset)
        % to compute censoring probabilities for ALL rows.
        lp_cens = X_td_global * b_cens;           % linear predictor (all rows)
        % Baseline hazard estimated from the non-competing subset only
        lp_cens_sub = X_cens_subset * b_cens;
        t_start_sub = t_start_td(not_competing);
        t_stop_sub  = t_stop_td(not_competing);
        uniq_times = unique(t_stop_sub);
        G_hat = ones(size(t_stop_td));             % P(C > t | X)
        for ui = 1:length(uniq_times)
            t_u = uniq_times(ui);
            at_risk_sub  = (t_start_sub < t_u) & (t_stop_sub >= t_u);
            events_u_sub = at_risk_sub & (t_stop_sub == t_u) & (is_admin_cens_subset == 1);
            if any(events_u_sub) && any(at_risk_sub)
                h0 = sum(events_u_sub) / sum(exp(lp_cens_sub(at_risk_sub)));
                rows_after = (t_stop_td >= t_u);
                G_hat(rows_after) = G_hat(rows_after) .* exp(-h0 * exp(lp_cens(rows_after)));
            end
        end

        % Stabilised weights: G_marginal / G_conditional (truncated)
        G_hat = max(G_hat, 0.05);                  % truncate to avoid extreme weights
        ipcw_weights = 1 ./ G_hat;
        ipcw_weights = ipcw_weights / mean(ipcw_weights);  % stabilise
        fprintf('  IPCW weights applied (admin censoring model, %d competing events excluded). Range: [%.2f, %.2f]\n', ...
            sum(event_td == 2), min(ipcw_weights), max(ipcw_weights));
    catch ME_ipcw
        fprintf('  ⚠️  IPCW weight estimation failed (%s). Proceeding unweighted.\n', ME_ipcw.message);
        ipcw_weights = ones(size(event_td));
    end
end

% ---- Fit the time-dependent Cox model --------------------------------
% coxphfit accepts a two-column [t_start t_stop] matrix for start-stop data.
warning('error', 'stats:coxphfit:FitWarning');
warning('error', 'stats:coxphfit:IterationLimit');
try
    T_td = [t_start_td, t_stop_td];
    is_censored = (event_td_csh == 0);

    % Suppress trivial warnings like Constant Term, but KEEP error triggers for convergence failures
    w_temp = warning('off', 'all');
    warning('error', 'stats:coxphfit:FitWarning');
    warning('error', 'stats:coxphfit:IterationLimit');
    [b_td, logl_td, ~, stats_td_raw] = coxphfit(X_td_global, T_td, ...
        'Censoring', is_censored, 'Ties', 'breslow', ...
        'Frequency', round(ipcw_weights * 100));   % integer pseudo-counts
    warning(w_temp);
    stats_td.se = stats_td_raw.se;
    stats_td.p  = stats_td_raw.p;
catch ME_td
    if exist('w_temp', 'var'), warning(w_temp); end
    if ~isempty(strfind(ME_td.identifier, 'FitWarning')) || ...
       ~isempty(strfind(ME_td.identifier, 'IterationLimit')) || ...
       ~isempty(strfind(lower(ME_td.message), 'perfect'))
        % Cox model failed to converge (e.g., perfect separation).
        % Set outputs to NaN rather than falling back to a different model
        % to avoid mixing Hazard Ratios with Odds Ratios.
        fprintf('  ⚠️  Cox model did not converge: %s\n', ME_td.message);
        b_td        = nan(td_n_feat, 1);
        stats_td.se = nan(td_n_feat, 1);
        stats_td.p  = nan(td_n_feat, 1);
        logl_td     = NaN;
    else
        warning('on', 'stats:coxphfit:FitWarning');
        warning('on', 'stats:coxphfit:IterationLimit');
        rethrow(ME_td);
    end
end
warning('on', 'stats:coxphfit:FitWarning');
warning('on', 'stats:coxphfit:IterationLimit');

% ---- Likelihood-ratio test vs. null model (no covariates) -----------
try
    is_censored_null = (event_td == 0);
    w_temp_null = warning('off', 'all');
    cleanupObj = onCleanup(@() warning(w_temp_null));
    [~, logl_null_td] = coxphfit(zeros(size(X_td_global,1),1), [t_start_td, t_stop_td], ...
        'Censoring', is_censored_null, 'Ties', 'breslow', ...
        'Options', statset('MaxIter', 100));
    LRT_stat = 2 * (logl_td - logl_null_td);
    LRT_p    = 1 - chi2cdf(LRT_stat, td_n_feat);
catch
    LRT_stat = NaN; LRT_p = NaN;
end

% ---- Print results table -------------------------------------------
fprintf('  %-10s  %6s  %6s  %6s  %6s\n', 'Covariate', 'HR ', 'CI_lo', 'CI_hi', 'p');
fprintf('  %s\n', repmat('-', 1, 52));
for fi = 1:td_n_feat
    hr_i  = exp(b_td(fi));
    ci_lo = exp(b_td(fi) - 1.96*stats_td.se(fi));
    ci_hi = exp(b_td(fi) + 1.96*stats_td.se(fi));
    fprintf('  %-10s  %6.3f  %6.3f  %6.3f  %6.4f\n', ...
        td_feat_names{fi}, hr_i, ci_lo, ci_hi, stats_td.p(fi));
end
if ~isnan(LRT_p)
    fprintf('  Global LRT: chi2(%d) = %.2f, p = %.4f\n', td_n_feat, LRT_stat, LRT_p);
end

% ---- Imputation half-life sensitivity analysis -----------------------
% Report how HRs change across the half-life grid so users can assess
% robustness of the exponential decay imputation assumption.
fprintf('\n  --- Imputation Sensitivity (half-life months → HR) ---\n');
fprintf('  %-6s', 'HL');
for fi = 1:td_n_feat
    fprintf('  %10s', td_feat_names{fi});
end
fprintf('\n');
for hl_idx = 1:length(half_life_grid)
    pnl = td_panels{hl_idx};
    try
        ev_csh = pnl.event; ev_csh(ev_csh == 2) = 0;
        X_hl = scale_td_panel(pnl.X, td_feat_names, pnl.pat_id, pnl.t_start, unique(pnl.pat_id));
        w_hl = warning('off', 'all');
        [b_hl] = coxphfit(X_hl, [pnl.t_start, pnl.t_stop], 'Censoring', ev_csh==0, 'Ties', 'breslow');
        warning(w_hl);
        fprintf('  %-6d', half_life_grid(hl_idx));
        for fi = 1:td_n_feat
            fprintf('  %10.3f', exp(b_hl(fi)));
        end
        fprintf('\n');
    catch
        fprintf('  %-6d  (model did not converge)\n', half_life_grid(hl_idx));
    end
end

fprintf('\nMetrics module sequence completed successfully.\n');

if ~isempty(output_folder)
    diary off;
end
end

%% ========================================================================
%  Helper Functions for Time-Dependent Cox and Elastic Net (Copied from metrics.m)
% =========================================================================

function [X_panel, t_start_out, t_stop_out, event_out, pat_id_out, frac_out] = build_td_panel(feat_arrays, feat_names, lf, tot_time, nTp, scan_days, hl_months)
    % Unpacks subject-level data into a start-stop panel for time-dependent Cox.
    n_pts   = size(feat_arrays{1}, 1);
    n_feat  = numel(feat_arrays);
    n_rows_est = n_pts * nTp;
    
    X_panel_pre = nan(n_rows_est, n_feat);
    t_start_pre = nan(n_rows_est, 1);
    t_stop_pre  = nan(n_rows_est, 1);
    event_pre   = zeros(n_rows_est, 1);
    pat_id_pre  = zeros(n_rows_est, 1);
    frac_pre    = zeros(n_rows_est, 1);
    
    hl_days     = hl_months * 30.44;
    lambda_decay= log(2) / hl_days;
    
    row_count = 0;
    
    for pi = 1:n_pts
        surv_T  = tot_time(pi);
        ev_stat = lf(pi);
        
        if isnan(surv_T), continue; end
        
        last_t  = 0;
        last_X  = nan(1, n_feat);
        
        for ti = 1:nTp
            curr_X = nan(1, n_feat);
            for fi = 1:n_feat
                curr_X(fi) = feat_arrays{fi}(pi, ti);
            end
            
            % Exponential decay imputation if missing
            if any(isnan(curr_X))
                dt = scan_days(ti) - last_t;
                curr_X(isnan(curr_X)) = last_X(isnan(curr_X)) * exp(-lambda_decay * dt);
            end
            
            if any(isnan(curr_X))
                continue; % Still missing (e.g., Fx1 was NA)
            end
            
            t_int_start = scan_days(ti);
            
            t_int_stop = surv_T;
            for t_next = (ti+1):nTp
                next_X = nan(1, n_feat);
                for fi = 1:n_feat
                    next_X(fi) = feat_arrays{fi}(pi, t_next);
                end
                if ~any(isnan(next_X))
                    t_int_stop = scan_days(t_next);
                    break;
                end
            end
            
            if t_int_start >= surv_T
                break;
            end
            
            t_int_stop = min(t_int_stop, surv_T);
            is_terminal = (t_int_stop == surv_T);
            
            row_count = row_count + 1;
            X_panel_pre(row_count, :) = curr_X;
            t_start_pre(row_count)    = t_int_start;
            t_stop_pre(row_count)     = t_int_stop;
            pat_id_pre(row_count)     = pi;
            frac_pre(row_count)       = ti;
            
            if is_terminal
                event_pre(row_count) = ev_stat;
            else
                event_pre(row_count) = 0;
            end
            
            last_X = curr_X;
            last_t = t_int_start;
            
            if is_terminal, break; end 
        end
    end
    
    X_panel     = X_panel_pre(1:row_count, :);
    t_start_out = t_start_pre(1:row_count);
    t_stop_out  = t_stop_pre(1:row_count);
    event_out   = event_pre(1:row_count);
    pat_id_out  = pat_id_pre(1:row_count);
    frac_out    = frac_pre(1:row_count);
end

function X_scaled = scale_td_panel(X_td, feat_names, pat_id, t_start, train_pids)
    % Z-scores the panel data using the baseline (t_start==0) mean and std
    % strictly from the training patients.
    
    X_scaled = X_td;
    n_feat   = size(X_td, 2);
    
    for fi = 1:n_feat
        is_train_base = ismember(pat_id, train_pids) & (t_start == 0) & ~isnan(X_td(:, fi));
        
        base_vals = X_td(is_train_base, fi);
        
        if isempty(base_vals)
            mu  = mean(X_td(ismember(pat_id, train_pids), fi), 'omitnan');
            sig = std(X_td(ismember(pat_id, train_pids), fi), 0, 'omitnan');
        else
            mu  = mean(base_vals);
            sig = std(base_vals);
        end
        
        if sig == 0 || isnan(sig), sig = 1; end
        
        X_scaled(:, fi) = (X_td(:, fi) - mu) / sig;
    end
end
