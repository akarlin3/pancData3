function metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label)
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
%
% Outputs:
%   None. Outputs printed to console (HR and p-value tables).
%

fprintf('\n--- TIME-DEPENDENT COX PH MODEL (Counting Process) ---\n');

td_scan_days = [0, 5, 10, 15, 20, 90];   % update if exact scan dates are available

% Covariates: all four IVIM parameters (absolute, all fractions)
td_feat_arrays = { ADC_abs(valid_pts,:), D_abs(valid_pts,:), ...
                   f_abs(valid_pts,:),   Dstar_abs(valid_pts,:) };
td_feat_names  = {'ADC', 'D', 'f', 'D*'};
td_n_feat      = numel(td_feat_arrays);

td_lf       = m_lf(valid_pts);
td_tot_time = m_total_time(valid_pts);

% Censored patients use follow-up time; events use time-to-event
cens_mask_td = (td_lf == 0) & ~isnan(m_total_follow_up_time(valid_pts));
td_tot_time(cens_mask_td) = m_total_follow_up_time(valid_pts & (m_lf(:)==0) & ~isnan(m_total_follow_up_time(:)));

[X_td_def, t_start_td_def, t_stop_td_def, event_td_def, pat_id_td_def, frac_td_def] = ...
    build_td_panel(td_feat_arrays, td_feat_names, td_lf, td_tot_time, nTp, td_scan_days, 18);

% Global exploratory scaling using all patients
X_td_global_def = scale_td_panel(X_td_def, td_feat_names, pat_id_td_def, t_start_td_def, unique(pat_id_td_def));

td_ok = (sum(event_td_def) >= 3) && (size(X_td_global_def, 1) > td_n_feat + 1);

half_life_grid = [3, 6, 12, 18, 24];
td_panels = cell(length(half_life_grid), 1);
for hl_idx = 1:length(half_life_grid)
    [X_i, t_start_i, t_stop_i, ev_i, pid_i, frac_i] = ...
        build_td_panel(td_feat_arrays, td_feat_names, td_lf, td_tot_time, nTp, td_scan_days, half_life_grid(hl_idx));
    td_panels{hl_idx} = struct('X', X_i, 't_start', t_start_i, 't_stop', t_stop_i, 'event', ev_i, 'pat_id', pid_i, 'frac', frac_i);
end

if ~td_ok
    fprintf('  Insufficient events (%d) or intervals for time-dependent Cox model.\n', sum(event_td_def));
    return;
end

X_td = X_td_def; t_start_td = t_start_td_def; t_stop_td = t_stop_td_def; event_td = event_td_def; pat_id_td = pat_id_td_def; frac_td = frac_td_def;
X_td_global = X_td_global_def;

% ---- Fit the time-dependent Cox model --------------------------------
% coxphfit accepts a two-column [t_start t_stop] matrix for start-stop data.
warning('error', 'stats:coxphfit:FitWarning');
warning('error', 'stats:coxphfit:IterationLimit');
is_firth_td = false;
try
    % Treat competing risks (2) as right-censored (0) for Cause-Specific Hazard.
    event_td_csh = event_td;
    event_td_csh(event_td_csh == 2) = 0;
    T_td = [t_start_td, t_stop_td];
    is_censored = (event_td_csh == 0);
    
    % Suppress trivial warnings like Constant Term, but KEEP error triggers for Firth fallback
    w_temp = warning('off', 'all');
    warning('error', 'stats:coxphfit:FitWarning');
    warning('error', 'stats:coxphfit:IterationLimit');
    mdl_td = fitcox(X_td_global, T_td, 'Censoring', is_censored, 'TieBreakMethod', 'breslow');
    warning(w_temp);
    b_td        = mdl_td.Coefficients.Beta;
    logl_td     = mdl_td.LogLikelihood;
    stats_td.se = mdl_td.Coefficients.SE;
    stats_td.p  = mdl_td.Coefficients.pValue;
    
    % Extract scaled Schoenfeld residuals for PH check
    stats_td.sschres = mdl_td.Residuals.ScaledSchoenfeld;
catch ME_td
    if ~isempty(strfind(ME_td.identifier, 'FitWarning')) || ...
       ~isempty(strfind(ME_td.identifier, 'IterationLimit')) || ...
       ~isempty(strfind(lower(ME_td.message), 'perfect'))
        % Firth penalised logistic fallback: use last observed covariate per patient
        n_vp = sum(valid_pts);
        td_last_X = nan(n_vp, td_n_feat);
        for ii = 1:n_vp
            rows_ii = (pat_id_td == ii);
            if any(rows_ii)
                td_last_X(ii, :) = X_td_global(find(rows_ii, 1, 'last'), :);
            end
        end
        valid_firth_td = ~any(isnan(td_last_X), 2) & ~isnan(td_lf);
        mdl_firth_td = fitglm(td_last_X(valid_firth_td,:), td_lf(valid_firth_td), ...
            'Distribution', 'binomial', 'LikelihoodPenalty', 'jeffreys-prior', ...
            'Options', statset('MaxIter', 10000));
        b_td           = mdl_firth_td.Coefficients.Estimate(2:end);
        stats_td.se    = mdl_firth_td.Coefficients.SE(2:end);
        stats_td.p     = mdl_firth_td.Coefficients.pValue(2:end);
        stats_td.sschres = zeros(size(X_td_global,1), td_n_feat);
        logl_td        = mdl_firth_td.LogLikelihood;
        is_firth_td    = true;
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
    w_temp_null = warning('off', 'stats:coxphfit:ConstantTerm');
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
if is_firth_td
    fprintf('  [Fallback to Firth logistic due to separation]\n');
    hr_label = 'OR ';
else
    hr_label = 'HR ';
end
fprintf('  %-10s  %6s  %6s  %6s  %6s\n', 'Covariate', hr_label, 'CI_lo', 'CI_hi', 'p');
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

fprintf('\nMetrics module sequence completed successfully.\n');
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
            mu  = nanmean(X_td(ismember(pat_id, train_pids), fi));
            sig = nanstd(X_td(ismember(pat_id, train_pids), fi));
        else
            mu  = mean(base_vals);
            sig = std(base_vals);
        end
        
        if sig == 0 || isnan(sig), sig = 1; end
        
        X_scaled(:, fi) = (X_td(:, fi) - mu) / sig;
    end
end
