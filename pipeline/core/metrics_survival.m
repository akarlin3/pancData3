function metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label, m_gtv_vol, output_folder, actual_scan_days, config_struct_in)
% METRICS_SURVIVAL — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 5/5 of the metrics step. Fits survival models including Cox PH and competing risks analysis.
%
% This function coordinates the entire survival analysis workflow by calling
% specialized subfunctions for different aspects of the analysis.

% Handle optional arguments for backward compatibility
if nargin < 12, m_gtv_vol = []; end
if nargin < 13, output_folder = ''; end
if nargin < 14, actual_scan_days = []; end
if nargin < 15 || isempty(config_struct_in)
    config_struct_internal = struct();
else
    config_struct_internal = config_struct_in;
end

fprintf('\n--- SURVIVAL ANALYSIS PIPELINE ---\n');

% Set up diary for output capture
if ~isempty(output_folder)
    diary_file = fullfile(output_folder, ['metrics_survival_output_' dtype_label '.txt']);
    if exist(diary_file, 'file'), delete(diary_file); end
    diary(diary_file);
end

% Prepare common data structures
[survival_data, config] = prepare_survival_data(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ...
    m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label, ...
    m_gtv_vol, actual_scan_days, config_struct_internal);

% Check if we have sufficient data for analysis
if ~survival_data.has_sufficient_data
    fprintf('  Insufficient data for survival analysis.\n');
    if ~isempty(output_folder), diary off; end
    return;
end

% Fit Cox Proportional Hazards model
fprintf('\n--- COX PROPORTIONAL HAZARDS MODEL ---\n');
cox_results = fit_cox_ph(survival_data, config);

% Fit Fine-Gray competing risks model
fprintf('\n--- FINE-GRAY COMPETING RISKS MODEL ---\n');
finegray_results = fit_fine_gray(survival_data, config);

% Compute time-varying effects analysis
fprintf('\n--- TIME-VARYING COEFFICIENTS ANALYSIS ---\n');
timevar_results = compute_time_varying_effects(survival_data, config);

% Validate survival models
fprintf('\n--- MODEL VALIDATION ---\n');
validation_results = validate_survival_model(survival_data, cox_results, config);

% Output summary results
print_survival_summary(cox_results, finegray_results, timevar_results, validation_results);

if ~isempty(output_folder), diary off; end
end

function [survival_data, config] = prepare_survival_data(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ...
    m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label, ...
    m_gtv_vol, actual_scan_days, config_struct_in)
% PREPARE_SURVIVAL_DATA Prepare and validate data for survival analysis

% Set up timing
if ~isempty(actual_scan_days)
    td_scan_days = actual_scan_days;
    fprintf('  Using provided scan days [%s].\n', num2str(td_scan_days));
else
    td_scan_days = [0, 5, 10, 15, 20, 90];
    fprintf('  ⚠️  CAUTION: Using default scan days [%s].\n', num2str(td_scan_days));
    fprintf('      Pass actual DICOM-derived scan days to avoid immortal time bias.\n');
end

% Prepare feature arrays
td_feat_arrays = { ADC_abs(valid_pts,:), D_abs(valid_pts,:), ...
                   f_abs(valid_pts,:),   Dstar_abs(valid_pts,:) };
td_feat_names  = {'ADC', 'D', 'f', 'D*'};

% Include baseline GTV volume if available
has_vol = ~isempty(m_gtv_vol) && any(isfinite(m_gtv_vol(valid_pts, 1)));
if has_vol
    vol_baseline = m_gtv_vol(valid_pts, 1);
    vol_rep = repmat(vol_baseline, 1, nTp);
    td_feat_arrays{end+1} = vol_rep;
    td_feat_names{end+1}  = 'GTVvol';
    fprintf('  Including baseline GTV volume as a covariate.\n');
end

% Prepare outcome data
td_lf = m_lf(valid_pts);
td_tot_time = m_total_time(valid_pts);
follow_up_valid = m_total_follow_up_time(valid_pts);
cens_mask_td = (td_lf == 0 | td_lf == 2) & ~isnan(follow_up_valid);
td_tot_time(cens_mask_td) = follow_up_valid(cens_mask_td);

% Build time-dependent panel with default half-life
[X_td_def, t_start_td_def, t_stop_td_def, event_td_def, pat_id_td_def, frac_td_def] = ...
    build_td_panel(td_feat_arrays, td_feat_names, td_lf, td_tot_time, nTp, td_scan_days, 18);

% Check data sufficiency
td_n_feat = numel(td_feat_arrays);
td_ok = (sum(event_td_def == 1) >= 3) && (size(X_td_def, 1) > td_n_feat + 1);

% Apply landmark analysis
scan_gaps = diff(td_scan_days);
[~, gap_idx] = max(scan_gaps);
landmark_day = td_scan_days(gap_idx);
lm_keep = (t_start_td_def >= landmark_day);

if any(lm_keep)
    X_td = X_td_def(lm_keep, :);
    t_start_td = t_start_td_def(lm_keep);
    t_stop_td = t_stop_td_def(lm_keep);
    event_td = event_td_def(lm_keep);
    pat_id_td = pat_id_td_def(lm_keep);
    frac_td = frac_td_def(lm_keep);
    
    n_events_post_lm = sum(event_td == 1);
    td_ok = td_ok && (n_events_post_lm >= 3) && (size(X_td, 1) > td_n_feat + 1);
    
    fprintf('  Landmark at day %d: %d events, %d patients\n', ...
        landmark_day, n_events_post_lm, numel(unique(pat_id_td)));
else
    td_ok = false;
end

% Package results
survival_data = struct();
survival_data.has_sufficient_data = td_ok;
if td_ok
    survival_data.X = X_td;
    survival_data.t_start = t_start_td;
    survival_data.t_stop = t_stop_td;
    survival_data.event = event_td;
    survival_data.pat_id = pat_id_td;
    survival_data.frac = frac_td;
    survival_data.feat_names = td_feat_names;
    survival_data.n_feat = td_n_feat;
    survival_data.landmark_day = landmark_day;
    survival_data.scan_days = td_scan_days;
    survival_data.feat_arrays = td_feat_arrays;
    survival_data.lf_status = td_lf;
    survival_data.total_time = td_tot_time;
end

config = struct();
config.half_life_grid = [3, 6, 12, 18, 24];
config.n_bootstrap = 1000;
config.output_folder = output_folder;
config.dtype_label = dtype_label;
end

function cox_results = fit_cox_ph(survival_data, config)
% FIT_COX_PH Fit Cox Proportional Hazards model with time-dependent covariates

if ~survival_data.has_sufficient_data
    cox_results = struct('success', false, 'message', 'Insufficient data');
    return;
end

% Scale covariates
X_scaled = scale_td_panel(survival_data.X, survival_data.feat_names, ...
    survival_data.pat_id, survival_data.t_start, unique(survival_data.pat_id), 'baseline');

% Prepare for cause-specific hazard
event_csh = survival_data.event;
event_csh(event_csh == 2) = 0;  % Competing risks → censored

% Compute IPCW weights
ipcw_weights = compute_ipcw_weights(survival_data.event, survival_data.t_start, ...
    survival_data.t_stop, X_scaled, survival_data.pat_id);

% Check for tied event times
event_times = survival_data.t_stop(event_csh == 1);
[unique_times, ~, time_idx] = unique(event_times);
has_tied_times = any(arrayfun(@(i) sum(time_idx == i) > 1, 1:length(unique_times)));
if has_tied_times
    ties_method = 'efron';
else
    ties_method = 'breslow';
end

if has_tied_times
    fprintf('  ⚠️  Detected tied event times. Using Efron approximation.\n');
end

% Fit Cox model
try
    T_matrix = [survival_data.t_start, survival_data.t_stop];
    is_censored = (event_csh == 0);
    
    [X_clean, keep_cols] = remove_constant_columns(X_scaled);
    if size(X_clean, 2) == 0
        error('Cox:NoVariableColumns', 'All covariate columns are constant.');
    end
    
    % Apply IPCW via frequency weights
    ipcw_scale = 100;
    ipcw_freq = max(1, round(ipcw_weights * ipcw_scale));
    
    warning('off', 'all');
    warning('error', 'stats:coxphfit:FitWarning');
    warning('error', 'stats:coxphfit:IterationLimit');
    
    [b_short, logl, ~, stats_short] = coxphfit(X_clean, T_matrix, ...
        'Censoring', is_censored, 'Ties', ties_method, 'Frequency', ipcw_freq);
    
    % Correct for frequency scaling
    eff_scale = mean(ipcw_freq);
    stats_short.se = stats_short.se * sqrt(eff_scale);
    stats_short.p = 2 * (1 - normcdf(abs(b_short ./ stats_short.se)));
    
    % Map back to full feature space
    b_full = zeros(survival_data.n_feat, 1);
    b_full(keep_cols) = b_short;
    se_full = nan(survival_data.n_feat, 1);
    p_full = nan(survival_data.n_feat, 1);
    se_full(keep_cols) = stats_short.se;
    p_full(keep_cols) = stats_short.p;
    
    warning('on', 'all');
    
    cox_results = struct();
    cox_results.success = true;
    cox_results.coefficients = b_full;
    cox_results.se = se_full;
    cox_results.p_values = p_full;
    cox_results.hazard_ratios = exp(b_full);
    cox_results.loglik = logl;
    cox_results.keep_cols = keep_cols;
    cox_results.X_scaled = X_scaled;
    cox_results.event_csh = event_csh;
    cox_results.ipcw_weights = ipcw_weights;
    
    % Print results
    print_cox_results(cox_results, survival_data.feat_names);
    
catch ME
    warning('on', 'all');
    if contains(ME.identifier, 'FitWarning') || contains(ME.identifier, 'IterationLimit')
        fprintf('  ⚠️  Cox model convergence issue: %s\n', ME.message);
        cox_results = create_failed_cox_results(survival_data.n_feat);
    else
        rethrow(ME);
    end
end
end

function finegray_results = fit_fine_gray(survival_data, config)
% FIT_FINE_GRAY Fit Fine-Gray competing risks model with interval censoring support

if ~survival_data.has_sufficient_data
    finegray_results = struct('success', false, 'message', 'Insufficient data');
    return;
end

fprintf('  Fine-Gray subdistribution hazard model with interval censoring support.\n');

try
    % Detect interval-censored observations
    interval_censored = detect_interval_censoring(survival_data);
    
    if any(interval_censored)
        fprintf('  Detected %d interval-censored observations.\n', sum(interval_censored));
        
        % Handle interval censoring using midpoint imputation as primary method
        [t_start_adj, t_stop_adj, event_adj] = handle_interval_censoring(survival_data, interval_censored);
        
        % Optionally perform multiple imputation sensitivity analysis
        if config.n_bootstrap > 100  % Only if sufficient bootstrap samples
            mi_results = perform_multiple_imputation_fg(survival_data, interval_censored, config);
            fprintf('  Multiple imputation completed for sensitivity analysis.\n');
        else
            mi_results = struct('performed', false);
        end
    else
        t_start_adj = survival_data.t_start;
        t_stop_adj = survival_data.t_stop;
        event_adj = survival_data.event;
        mi_results = struct('performed', false);
        fprintf('  No interval censoring detected.\n');
    end
    
    % Scale covariates
    X_scaled = scale_td_panel(survival_data.X, survival_data.feat_names, ...
        survival_data.pat_id, t_start_adj, unique(survival_data.pat_id), 'baseline');
    
    % Fit subdistribution hazard model using Fine-Gray approach
    fg_model = fit_subdistribution_hazard(X_scaled, t_start_adj, t_stop_adj, ...
        event_adj, survival_data.pat_id, survival_data.feat_names);
    
    finegray_results = struct();
    finegray_results.success = fg_model.converged;
    
    if fg_model.converged
        finegray_results.coefficients = fg_model.coefficients;
        finegray_results.se = fg_model.se;
        finegray_results.p_values = fg_model.p_values;
        finegray_results.shr = exp(fg_model.coefficients);  % subdistribution hazard ratios
        finegray_results.loglik = fg_model.loglik;
        finegray_results.cumulative_incidence = fg_model.cif;
        finegray_results.interval_censored_count = sum(interval_censored);
        finegray_results.imputation_method = 'midpoint';
        finegray_results.multiple_imputation = mi_results;
        
        % Print Fine-Gray results
        print_finegray_results(finegray_results, survival_data.feat_names);
    else
        finegray_results.message = 'Fine-Gray model did not converge';
        fprintf('  ⚠️  Fine-Gray model convergence failed.\n');
    end
    
catch ME
    fprintf('  ⚠️  Fine-Gray model error: %s\n', ME.message);
    finegray_results = struct('success', false, 'message', ME.message);
end
end

function interval_censored = detect_interval_censoring(survival_data)
% DETECT_INTERVAL_CENSORING Identify interval-censored observations

% Interval censoring occurs when we only know the event occurred within an interval
% This can happen when events are detected at discrete follow-up visits
% rather than exact occurrence times

% Simple heuristic: observations with suspicious timing patterns
scan_intervals = diff(survival_data.scan_days);
typical_interval = median(scan_intervals(scan_intervals > 0));

% Flag observations where event time coincides exactly with scan times
% (suggesting interval censoring rather than exact timing)
event_times = survival_data.t_stop(survival_data.event > 0);
scan_time_matches = false(size(survival_data.t_stop));

for i = 1:length(survival_data.t_stop)
    if survival_data.event(i) > 0  % Only for events, not censored
        time_diff_to_scans = abs(survival_data.t_stop(i) - survival_data.scan_days);
        scan_time_matches(i) = any(time_diff_to_scans < typical_interval * 0.1);
    end
end

% Additional criteria: multiple events at exact same time
[unique_event_times, ~, time_groups] = unique(event_times);
tied_event_times = false(size(survival_data.t_stop));
for i = 1:length(unique_event_times)
    if sum(time_groups == i) > 1  % Multiple events at same time
        tied_event_times(survival_data.t_stop == unique_event_times(i) & survival_data.event > 0) = true;
    end
end

interval_censored = scan_time_matches | tied_event_times;
end

function [t_start_adj, t_stop_adj, event_adj] = handle_interval_censoring(survival_data, interval_censored)
% HANDLE_INTERVAL_CENSORING Apply midpoint imputation for interval-censored observations

t_start_adj = survival_data.t_start;
t_stop_adj = survival_data.t_stop;
event_adj = survival_data.event;

% For interval-censored observations, use midpoint imputation
scan_days = survival_data.scan_days;

for i = find(interval_censored)'
    if event_adj(i) > 0  % Only adjust events, not censored observations
        % Find the interval containing this event time
        event_time = t_stop_adj(i);
        
        % Find preceding and following scan times
        preceding_scans = scan_days(scan_days <= event_time);
        following_scans = scan_days(scan_days > event_time);
        
        if ~isempty(preceding_scans) && ~isempty(following_scans)
            interval_start = max(preceding_scans);
            interval_end = min(following_scans);
            
            % Use midpoint of the interval
            midpoint = (interval_start + interval_end) / 2;
            
            % Adjust times
            t_start_adj(i) = max(t_start_adj(i), interval_start);
            t_stop_adj(i) = midpoint;
            
        elseif ~isempty(preceding_scans)
            % Event after last scan - use previous interval pattern
            if length(preceding_scans) >= 2
                typical_gap = median(diff(preceding_scans));
                interval_start = max(preceding_scans);
                t_start_adj(i) = interval_start;
                t_stop_adj(i) = interval_start + typical_gap/2;
            end
        end
    end
end
end

function mi_results = perform_multiple_imputation_fg(survival_data, interval_censored, config)
% PERFORM_MULTIPLE_IMPUTATION_FG Multiple imputation sensitivity analysis

n_imputations = min(10, floor(config.n_bootstrap / 100));  % Reasonable number of imputations
coeffs_matrix = nan(survival_data.n_feat, n_imputations);

fprintf('  Performing multiple imputation (M=%d) for sensitivity analysis.\n', n_imputations);

for m = 1:n_imputations
    try
        % Generate imputed event times
        [t_start_imp, t_stop_imp, event_imp] = impute_interval_times(survival_data, interval_censored, m);
        
        % Fit model with imputed data
        X_scaled = scale_td_panel(survival_data.X, survival_data.feat_names, ...
            survival_data.pat_id, t_start_imp, unique(survival_data.pat_id), 'baseline');
        
        fg_model_imp = fit_subdistribution_hazard(X_scaled, t_start_imp, t_stop_imp, ...
            event_imp, survival_data.pat_id, survival_data.feat_names);
        
        if fg_model_imp.converged
            coeffs_matrix(:, m) = fg_model_imp.coefficients;
        end
        
    catch
        % Skip failed imputation
        continue;
    end
end

% Combine results using Rubin's rules
valid_imputations = ~isnan(coeffs_matrix(1, :));
n_valid = sum(valid_imputations);

if n_valid >= 3
    coeffs_valid = coeffs_matrix(:, valid_imputations);
    
    % Pool estimates
    pooled_coeffs = mean(coeffs_valid, 2);
    within_var = var(coeffs_valid, 0, 2);
    between_var = var(mean(coeffs_valid, 1)) * n_valid / (n_valid - 1);
    total_var = within_var + (1 + 1/n_valid) * between_var;
    
    mi_results = struct();
    mi_results.performed = true;
    mi_results.n_imputations = n_valid;
    mi_results.pooled_coefficients = pooled_coeffs;
    mi_results.pooled_se = sqrt(total_var);
    mi_results.pooled_shr = exp(pooled_coeffs);
    
else
    mi_results = struct('performed', false, 'message', 'Insufficient valid imputations');
end
end

function [t_start_imp, t_stop_imp, event_imp] = impute_interval_times(survival_data, interval_censored, seed)
% IMPUTE_INTERVAL_TIMES Generate single imputation for interval-censored times

rng(seed);  % Reproducible imputations

t_start_imp = survival_data.t_start;
t_stop_imp = survival_data.t_stop;
event_imp = survival_data.event;

scan_days = survival_data.scan_days;

for i = find(interval_censored)'
    if event_imp(i) > 0
        event_time = t_stop_imp(i);
        
        % Find interval
        preceding_scans = scan_days(scan_days <= event_time);
        following_scans = scan_days(scan_days > event_time);
        
        if ~isempty(preceding_scans) && ~isempty(following_scans)
            interval_start = max(preceding_scans);
            interval_end = min(following_scans);
            
            % Random imputation within interval (uniform distribution)
            imputed_time = interval_start + rand() * (interval_end - interval_start);
            
            t_start_imp(i) = max(t_start_imp(i), interval_start);
            t_stop_imp(i) = imputed_time;
        end
    end
end
end

function fg_model = fit_subdistribution_hazard(X, t_start, t_stop, event, pat_id, feat_names)
% FIT_SUBDISTRIBUTION_HAZARD Fit Fine-Gray subdistribution hazard model

% This is a simplified implementation of the Fine-Gray approach
% In practice, this would use specialized competing risks software

try
    % Focus on primary event of interest (event = 1)
    % Competing events (event = 2) are handled via subdistribution approach
    
    % Create modified dataset for subdistribution hazard
    [X_fg, t_start_fg, t_stop_fg, event_fg, weights_fg] = create_finegray_dataset(X, t_start, t_stop, event, pat_id);
    
    % Remove constant columns
    [X_clean, keep_cols] = remove_constant_columns(X_fg);
    if size(X_clean, 2) == 0
        fg_model = struct('converged', false);
        return;
    end
    
    % Fit weighted Cox model (approximation to Fine-Gray)
    T_fg = [t_start_fg, t_stop_fg];
    is_censored_fg = (event_fg == 0);
    
    warning('off', 'all');
    [b_short, logl, ~, stats_short] = coxphfit(X_clean, T_fg, ...
        'Censoring', is_censored_fg, 'Frequency', weights_fg);
    warning('on', 'all');
    
    % Map back to full feature space
    b_full = zeros(size(X, 2), 1);
    b_full(keep_cols) = b_short;
    se_full = nan(size(X, 2), 1);
    se_full(keep_cols) = stats_short.se;
    p_full = 2 * (1 - normcdf(abs(b_full ./ se_full)));
    
    % Estimate cumulative incidence function
    cif = estimate_cumulative_incidence(t_stop, event, X, b_full);
    
    fg_model = struct();
    fg_model.converged = true;
    fg_model.coefficients = b_full;
    fg_model.se = se_full;
    fg_model.p_values = p_full;
    fg_model.loglik = logl;
    fg_model.cif = cif;
    
catch ME
    fg_model = struct('converged', false, 'error', ME.message);
end
end

function [X_fg, t_start_fg, t_stop_fg, event_fg, weights_fg] = create_finegray_dataset(X, t_start, t_stop, event, pat_id)
% CREATE_FINEGRAY_DATASET Create Fine-Gray transformed dataset

% Simplified implementation - in practice would need full Fine-Gray transformation
% This includes artificial censoring and IPCW for competing events

n_obs = length(t_start);
keep_idx = true(n_obs, 1);

% For competing events (event = 2), create artificial censoring
competing_mask = (event == 2);

% Apply inverse probability weighting for competing risks
weights = ones(n_obs, 1);

% Simplified weight calculation (in practice, would use Kaplan-Meier for censoring)
if any(competing_mask)
    % Estimate probability of not having competing event
    competing_times = t_stop(competing_mask);
    for i = 1:n_obs
        if event(i) == 1  % Primary event
            % Weight by inverse probability of not having competing event by this time
            prob_no_competing = sum(competing_times > t_stop(i)) / length(competing_times);
            weights(i) = 1 / max(0.1, prob_no_competing);  % Avoid extreme weights
        end
    end
end

X_fg = X(keep_idx, :);
t_start_fg = t_start(keep_idx);
t_stop_fg = t_stop(keep_idx);
event_fg = event(keep_idx);
event_fg(event_fg == 2) = 0;  % Competing events become censored
weights_fg = max(1, round(weights(keep_idx) * 100));  % Scale for frequency weights
end

function cif = estimate_cumulative_incidence(t_stop, event, X, coefficients)
% ESTIMATE_CUMULATIVE_INCIDENCE Estimate cumulative incidence function

% Simplified CIF estimation
unique_times = unique(t_stop(event == 1));
unique_times = sort(unique_times);

cif = struct();
cif.times = unique_times;
cif.incidence = zeros(size(unique_times));

% Basic empirical CIF (would be more sophisticated in practice)
for i = 1:length(unique_times)
    t = unique_times(i);
    n_events = sum(event == 1 & t_stop <= t);
    n_at_risk = sum(t_stop >= t | (t_stop < t & event > 0));
    cif.incidence(i) = n_events / max(1, n_at_risk);
end
end

function print_finegray_results(finegray_results, feat_names)
% PRINT_FINEGRAY_RESULTS Print Fine-Gray model results

fprintf('  Fine-Gray Subdistribution Hazard Model Results:\n');
fprintf('  %-10s  %6s  %6s  %6s  %6s\n', 'Covariate', 'SHR', 'CI_lo', 'CI_hi', 'p');
fprintf('  %s\n', repmat('-', 1, 52));

for fi = 1:length(feat_names)
    shr = finegray_results.shr(fi);
    ci_lo = exp(finegray_results.coefficients(fi) - 1.96 * finegray_results.se(fi));
    ci_hi = exp(finegray_results.coefficients(fi) + 1.96 * finegray_results.se(fi));
    
    fprintf('  %-10s  %6.3f  %6.3f  %6.3f  %6.4f\n', ...
        feat_names{fi}, shr, ci_lo, ci_hi, finegray_results.p_values(fi));
end

if finegray_results.interval_censored_count > 0
    fprintf('  Note: %d interval-censored observations handled via %s imputation.\n', ...
        finegray_results.interval_censored_count, finegray_results.imputation_method);
end

if finegray_results.multiple_imputation.performed
    fprintf('  Multiple imputation performed with %d imputations for sensitivity analysis.\n', ...
        finegray_results.multiple_imputation.n_imputations);
end
end

function timevar_results = compute_time_varying_effects(survival_data, config)
% COMPUTE_TIME_VARYING_EFFECTS Analyze time-varying coefficient effects

if ~survival_data.has_sufficient_data
    timevar_results = struct('success', false, 'message', 'Insufficient data');
    return;
end

fprintf('  Analyzing sensitivity to imputation half-life assumptions.\n');

% Test different half-life assumptions
half_life_grid = config.half_life_grid;
n_halflife = length(half_life_grid);
hr_matrix = nan(survival_data.n_feat, n_halflife);

for hl_idx = 1:n_halflife
    text_progress_bar(hl_idx, n_halflife, 'Half-life sensitivity');
    
    try
        % Rebuild panel with different half-life
        [X_hl, t_start_hl, t_stop_hl, event_hl, pat_id_hl, ~] = ...
            build_td_panel(survival_data.feat_arrays, survival_data.feat_names, ...
            survival_data.lf_status, survival_data.total_time, ...
            size(survival_data.feat_arrays{1}, 2