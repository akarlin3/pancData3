function metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label, m_gtv_vol, output_folder, actual_scan_days)
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
%   actual_scan_days  - (Optional) Vector of actual scan days from DICOM
%                       headers or clinical records.  When provided,
%                       replaces the default [0,5,10,15,20,90] to avoid
%                       immortal time bias.
%
% Outputs:
%   None. Outputs printed to console (HR and p-value tables).
%

% Handle optional arguments for backward compatibility
if nargin < 12, m_gtv_vol = []; end
if nargin < 13, output_folder = ''; end
if nargin < 14, actual_scan_days = []; end

fprintf('\n--- TIME-DEPENDENT COX PH MODEL (Counting Process) ---\n');

% Diary: capture console output to output_folder
if ~isempty(output_folder)
    diary_file = fullfile(output_folder, ['metrics_survival_output_' dtype_label '.txt']);
    if exist(diary_file, 'file'), delete(diary_file); end
    diary(diary_file);
end

% Scan day timing for the time-dependent Cox model.
if ~isempty(actual_scan_days)
    td_scan_days = actual_scan_days;
    fprintf('  Using provided scan days [%s].\n', num2str(td_scan_days));
else
    % Default scan days assume 5 on-treatment fractions + 1 post-treatment scan.
    td_scan_days = [0, 5, 10, 15, 20, 90];
    warning('metrics_survival:defaultScanDays', ...
        'Using default scan days [%s]. Replace with actual timing to avoid immortal time bias.', ...
        num2str(td_scan_days));
    fprintf('  ⚠️  CAUTION: Using default scan days [%s].\n', num2str(td_scan_days));
    fprintf('      Pass actual DICOM-derived scan days via config.json td_scan_days field\n');
    fprintf('      or as the 14th argument to avoid immortal time bias.\n');
end

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
follow_up_valid = m_total_follow_up_time(valid_pts);
cens_mask_td = (td_lf == 0) & ~isnan(follow_up_valid);
td_tot_time(cens_mask_td) = follow_up_valid(cens_mask_td);

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
% The last element of td_scan_days is the post-treatment scan.
% Select the last on-treatment fraction as the landmark.  This adapts
% automatically to protocols with fewer or more than 5 fractions.
n_on_tx = max(length(td_scan_days) - 1, 1);  % exclude post-treatment element
landmark_idx = n_on_tx;
landmark_day = td_scan_days(landmark_idx);  % end of RT
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

    % Re-validate event count after landmark subsetting — the pre-landmark
    % check (line 76) may have passed, but removing early intervals can
    % reduce event count below the minimum for reliable Cox estimation.
    n_events_post_lm = sum(event_td > 0);
    if n_events_post_lm < 3 || size(X_td, 1) <= td_n_feat + 1
        fprintf('  Insufficient events (%d) or intervals after landmark subsetting for Cox model.\n', n_events_post_lm);
        if ~isempty(output_folder), diary off; end
        return;
    end
end

% Scale AFTER landmark subsetting so that standardisation statistics are
% computed exclusively on post-landmark intervals, preventing pre-landmark
% covariate distributions from leaking into the primary analysis.
X_td_global = scale_td_panel(X_td, td_feat_names, pat_id_td, t_start_td, unique(pat_id_td), 'baseline');

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
        [X_cens_clean, keep_ipcw] = remove_constant_columns(X_cens_subset);
        if size(X_cens_clean, 2) == 0
            error('IPCW:NoVariableColumns', 'All covariate columns are constant.');
        end
        w_ipcw = warning('off', 'all');
        [b_cens_short, ~, ~, stats_cens] = coxphfit(X_cens_clean, T_cens, ...
            'Censoring', (is_admin_cens_subset == 0), 'Ties', 'breslow');
        warning(w_ipcw);
        b_cens = zeros(td_n_feat, 1);
        b_cens(keep_ipcw) = b_cens_short;

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
                at_risk_full = (t_start_td < t_u) & (t_stop_td >= t_u);
                G_hat(at_risk_full) = G_hat(at_risk_full) .* exp(-h0 * exp(lp_cens(at_risk_full)));
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

    % Remove constant columns to prevent DisallowedConstantTerm warning
    [X_td_clean, keep_main] = remove_constant_columns(X_td_global);
    if sum(~keep_main) > 0
        fprintf('  💡 Removed %d constant column(s) before Cox fit: %s\n', ...
            sum(~keep_main), strjoin(td_feat_names(~keep_main), ', '));
    end
    if size(X_td_clean, 2) == 0
        error('Cox:NoVariableColumns', 'All covariate columns are constant after scaling.');
    end

    % Suppress trivial warnings, but KEEP error triggers for convergence failures
    w_temp = warning('off', 'all');
    warning('error', 'stats:coxphfit:FitWarning');
    warning('error', 'stats:coxphfit:IterationLimit');
    % IPCW weighting via coxphfit's 'Frequency' parameter.  'Frequency'
    % expects integer counts, so we scale weights by 100 before rounding to
    % preserve relative differences.  A scale factor of 1 (i.e., simple
    % rounding of mean-stabilised ~1 weights) collapses most weights to 1,
    % losing the IPCW correction.  Scaling by 100 retains two decimal
    % places of precision.  The ideal solution would be a custom weighted
    % Cox partial likelihood, but MATLAB's coxphfit does not support
    % continuous probability weights natively.
    ipcw_scale = 100;
    ipcw_freq = max(1, round(ipcw_weights * ipcw_scale));
    fprintf('  IPCW→Frequency: scale=%d, effective N=%d (actual N=%d)\n', ...
        ipcw_scale, sum(ipcw_freq), length(ipcw_weights));
    [b_td_short, logl_td, ~, stats_td_short] = coxphfit(X_td_clean, T_td, ...
        'Censoring', is_censored, 'Ties', 'breslow', ...
        'Frequency', ipcw_freq);
    warning(w_temp);

    % Correct variance inflation from the Frequency workaround.
    % coxphfit's Frequency parameter treats each "replicate" as an
    % independent observation, inflating effective N by ipcw_scale and
    % deflating SEs by sqrt(ipcw_scale).  Rescale SEs and recompute
    % p-values to recover correct inference.
    stats_td_short.se = stats_td_short.se * sqrt(ipcw_scale);
    stats_td_short.p  = 2 * (1 - normcdf(abs(b_td_short ./ stats_td_short.se)));

    % Map back to full feature space (removed columns get coef=0, SE/p=NaN)
    b_td = zeros(td_n_feat, 1);
    b_td(keep_main) = b_td_short;
    stats_td.se = nan(td_n_feat, 1);
    stats_td.p  = nan(td_n_feat, 1);
    stats_td.se(keep_main) = stats_td_short.se;
    stats_td.p(keep_main)  = stats_td_short.p;
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
    is_censored_null = (event_td_csh == 0);
    w_temp_null = warning('off', 'all');
    cleanupObj = onCleanup(@() warning(w_temp_null));
    [~, logl_null_td] = coxphfit(zeros(size(X_td_global,1),1), [t_start_td, t_stop_td], ...
        'Censoring', is_censored_null, 'Ties', 'breslow', ...
        'Frequency', ipcw_freq, ...
        'Options', statset('MaxIter', 100));
    LRT_stat = 2 * (logl_td - logl_null_td);
    LRT_df   = sum(keep_main);  % degrees of freedom = number of non-constant features actually fit
    LRT_p    = 1 - chi2cdf(LRT_stat, LRT_df);
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
    fprintf('  Global LRT: chi2(%d) = %.2f, p = %.4f\n', LRT_df, LRT_stat, LRT_p);
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
        % Apply landmark subsetting to match the primary analysis
        lm_keep_hl = (pnl.t_start >= landmark_day);
        if ~any(lm_keep_hl), error('sensitivity:noData', 'No post-landmark intervals.'); end
        pnl_X = pnl.X(lm_keep_hl, :);
        pnl_tstart = pnl.t_start(lm_keep_hl);
        pnl_tstop  = pnl.t_stop(lm_keep_hl);
        pnl_event  = pnl.event(lm_keep_hl);
        pnl_pid    = pnl.pat_id(lm_keep_hl);
        ev_csh = pnl_event; ev_csh(ev_csh == 2) = 0;
        X_hl = scale_td_panel(pnl_X, td_feat_names, pnl_pid, pnl_tstart, unique(pnl_pid), 'baseline');
        [X_hl_clean, keep_hl] = remove_constant_columns(X_hl);
        % Recompute IPCW weights for this sensitivity panel independently.
        % The primary ipcw_weights are indexed by lm_keep (from the default
        % half-life panel), NOT lm_keep_hl.  Slicing ipcw_weights(lm_keep)
        % then taking the first sum(lm_keep_hl) elements silently assigns
        % wrong weights when the two panels select different rows.
        ipcw_hl_sens = ones(sum(lm_keep_hl), 1);
        if has_admin_cens
            try
                not_comp_hl = (pnl_event ~= 2);
                is_cens_hl = double(pnl_event(not_comp_hl) == 0);
                T_cens_hl = [pnl_tstart(not_comp_hl), pnl_tstop(not_comp_hl)];
                X_cens_hl = X_hl(not_comp_hl, :);
                [X_cens_hl_clean, keep_ipcw_hl] = remove_constant_columns(X_cens_hl);
                if size(X_cens_hl_clean, 2) > 0
                    w_ipcw_hl = warning('off', 'all');
                    [b_cens_hl_short] = coxphfit(X_cens_hl_clean, T_cens_hl, ...
                        'Censoring', (is_cens_hl == 0), 'Ties', 'breslow');
                    warning(w_ipcw_hl);
                    b_cens_hl_full = zeros(td_n_feat, 1);
                    b_cens_hl_full(keep_ipcw_hl) = b_cens_hl_short;
                    lp_cens_hl = X_hl * b_cens_hl_full;
                    lp_cens_hl_sub = X_cens_hl * b_cens_hl_full;
                    t_start_hl_sub = pnl_tstart(not_comp_hl);
                    t_stop_hl_sub = pnl_tstop(not_comp_hl);
                    uniq_t_hl = unique(t_stop_hl_sub);
                    G_hl = ones(size(pnl_tstop));
                    for ui_hl = 1:length(uniq_t_hl)
                        t_u_hl = uniq_t_hl(ui_hl);
                        ar_sub = (t_start_hl_sub < t_u_hl) & (t_stop_hl_sub >= t_u_hl);
                        ev_sub = ar_sub & (t_stop_hl_sub == t_u_hl) & (is_cens_hl == 1);
                        if any(ev_sub) && any(ar_sub)
                            h0_hl = sum(ev_sub) / sum(exp(lp_cens_hl_sub(ar_sub)));
                            ar_full = (pnl_tstart < t_u_hl) & (pnl_tstop >= t_u_hl);
                            G_hl(ar_full) = G_hl(ar_full) .* exp(-h0_hl * exp(lp_cens_hl(ar_full)));
                        end
                    end
                    G_hl = max(G_hl, 0.05);
                    ipcw_hl_sens = 1 ./ G_hl;
                    ipcw_hl_sens = ipcw_hl_sens / mean(ipcw_hl_sens);
                end
            catch
                ipcw_hl_sens = ones(sum(lm_keep_hl), 1);
            end
        end
        w_hl = warning('off', 'all');
        ipcw_freq_hl = max(1, round(ipcw_hl_sens * ipcw_scale));
        [b_hl_short] = coxphfit(X_hl_clean, [pnl_tstart, pnl_tstop], ...
            'Censoring', ev_csh==0, 'Ties', 'breslow', ...
            'Frequency', ipcw_freq_hl);
        warning(w_hl);
        b_hl = zeros(td_n_feat, 1);
        b_hl(keep_hl) = b_hl_short;
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

% ---- Local helper functions --------------------------------------------

function [X_clean, keep] = remove_constant_columns(X)
%REMOVE_CONSTANT_COLUMNS  Drop columns with zero variance.
%   [X_clean, keep] = remove_constant_columns(X) returns the submatrix of
%   columns whose variance exceeds eps, and a logical mask of kept columns.
    col_var = var(X, 0, 1);
    keep = col_var > eps;
    X_clean = X(:, keep);
end
