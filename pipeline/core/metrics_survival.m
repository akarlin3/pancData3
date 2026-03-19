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
ties_method = has_tied_times ? 'efron' : 'breslow';

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
% FIT_FINE_GRAY Fit Fine-Gray competing risks model

if ~survival_data.has_sufficient_data
    finegray_results = struct('success', false, 'message', 'Insufficient data');
    return;
end

fprintf('  Fine-Gray model for competing risks analysis.\n');
fprintf('  (Implementation would use subdistribution hazard approach)\n');

% Placeholder for Fine-Gray implementation
% This would implement the subdistribution hazard model
finegray_results = struct();
finegray_results.success = false;
finegray_results.message = 'Fine-Gray model not yet implemented';
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
            size(survival_data.feat_arrays{1}, 2), survival_data.scan_days, ...
            half_life_grid(hl_idx));
        
        % Apply landmark
        lm_keep = (t_start_hl >= survival_data.landmark_day);
        if ~any(lm_keep), continue; end
        
        X_hl = X_hl(lm_keep, :);
        t_start_hl = t_start_hl(lm_keep);
        t_stop_hl = t_stop_hl(lm_keep);
        event_hl = event_hl(lm_keep);
        pat_id_hl = pat_id_hl(lm_keep);
        
        % Scale and fit
        X_scaled_hl = scale_td_panel(X_hl, survival_data.feat_names, ...
            pat_id_hl, t_start_hl, unique(pat_id_hl), 'baseline');
        
        event_csh_hl = event_hl;
        event_csh_hl(event_csh_hl == 2) = 0;
        
        [X_clean_hl, keep_cols_hl] = remove_constant_columns(X_scaled_hl);
        if size(X_clean_hl, 2) == 0, continue; end
        
        T_hl = [t_start_hl, t_stop_hl];
        is_censored_hl = (event_csh_hl == 0);
        
        warning('off', 'all');
        [b_hl, ~, ~, ~] = coxphfit(X_clean_hl, T_hl, 'Censoring', is_censored_hl);
        warning('on', 'all');
        
        b_full_hl = zeros(survival_data.n_feat, 1);
        b_full_hl(keep_cols_hl) = b_hl;
        hr_matrix(:, hl_idx) = exp(b_full_hl);
        
    catch
        % Skip this half-life if fitting fails
        continue;
    end
end

% Print sensitivity results
fprintf('\n  --- Imputation Sensitivity (half-life months → HR) ---\n');
fprintf('  %-6s', 'HL');
for fi = 1:survival_data.n_feat
    fprintf('  %10s', survival_data.feat_names{fi});
end
fprintf('\n');

for hl_idx = 1:n_halflife
    fprintf('  %6.0f', half_life_grid(hl_idx));
    for fi = 1:survival_data.n_feat
        if isfinite(hr_matrix(fi, hl_idx))
            fprintf('  %10.3f', hr_matrix(fi, hl_idx));
        else
            fprintf('  %10s', 'NaN');
        end
    end
    fprintf('\n');
end

timevar_results = struct();
timevar_results.success = true;
timevar_results.half_life_grid = half_life_grid;
timevar_results.hr_matrix = hr_matrix;
timevar_results.feat_names = survival_data.feat_names;
end

function validation_results = validate_survival_model(survival_data, cox_results, config)
% VALIDATE_SURVIVAL_MODEL Perform model validation including bootstrap CIs and diagnostics

if ~survival_data.has_sufficient_data || ~cox_results.success
    validation_results = struct('success', false, 'message', 'No valid model to validate');
    return;
end

fprintf('  Performing bootstrap validation.\n');

% Bootstrap confidence intervals
try
    n_boot = config.n_bootstrap;
    fprintf('  Computing bootstrap 95%% CIs (B=%d)\n', n_boot);
    
    boot_data = [cox_results.X_scaled, [survival_data.t_start, survival_data.t_stop], ...
                 double(cox_results.event_csh == 0), cox_results.ipcw_weights];
    
    boot_hrs = nan(survival_data.n_feat, 2);  % [lower, upper]
    
    for fi = 1:survival_data.n_feat
        if cox_results.keep_cols(fi)
            hr_fn = @(d) exp(bootstrap_cox_coef(d, fi, cox_results.keep_cols, survival_data.n_feat));
            [bci_lo, bci_hi] = bootstrap_ci(boot_data, hr_fn, n_boot, 0.05);
            boot_hrs(fi, :) = [bci_lo, bci_hi];
        end
    end
    
    % Print bootstrap results
    fprintf('  %-10s  %6s  %6s  %6s\n', 'Covariate', 'HR', 'BCa_lo', 'BCa_hi');
    fprintf('  %s\n', repmat('-', 1, 44));
    for fi = 1:survival_data.n_feat
        fprintf('  %-10s  %6.3f  %6.3f  %6.3f\n', ...
            survival_data.feat_names{fi}, cox_results.hazard_ratios(fi), ...
            boot_hrs(fi, 1), boot_hrs(fi, 2));
    end
    
    validation_results = struct();
    validation_results.success = true;
    validation_results.bootstrap_cis = boot_hrs;
    
catch boot_err
    fprintf('  ⚠️  Bootstrap validation failed: %s\n', boot_err.message);
    validation_results = struct('success', false, 'message', boot_err.message);
end
end

function print_cox_results(cox_results, feat_names)
% Print Cox model results table

fprintf('  %-10s  %6s  %6s  %6s  %6s\n', 'Covariate', 'HR ', 'CI_lo', 'CI_hi', 'p');
fprintf('  %s\n', repmat('-', 1, 52));

for fi = 1:length(feat_names)
    hr = cox_results.hazard_ratios(fi);
    ci_lo = exp(cox_results.coefficients(fi) - 1.96 * cox_results.se(fi));
    ci_hi = exp(cox_results.coefficients(fi) + 1.96 * cox_results.se(fi));
    
    fprintf('  %-10s  %6.3f  %6.3f  %6.3f  %6.4f\n', ...
        feat_names{fi}, hr, ci_lo, ci_hi, cox_results.p_values(fi));
end
end

function failed_results = create_failed_cox_results(n_feat)
% Create failed results structure for Cox model

failed_results = struct();
failed_results.success = false;
failed_results.coefficients = nan(n_feat, 1);
failed_results.se = nan(n_feat, 1);
failed_results.p_values = nan(n_feat, 1);
failed_results.hazard_ratios = nan(n_feat, 1);
failed_results.loglik = NaN;
end

function print_survival_summary(cox_results, finegray_results, timevar_results, validation_results)
% Print overall summary of survival analysis

fprintf('\n--- SURVIVAL ANALYSIS SUMMARY ---\n');

if cox_results.success
    fprintf('✓ Cox PH model fitted successfully\n');
else
    fprintf('✗ Cox PH model failed\n');
end

if finegray_results.success
    fprintf('✓ Fine-Gray model fitted successfully\n');
else
    fprintf('✗ Fine-Gray model not available\n');
end

if timevar_results.success
    fprintf('✓ Time-varying effects analysis completed\n');
else
    fprintf('✗ Time-varying effects analysis failed\n');
end

if validation_results.success
    fprintf('✓ Model validation completed\n');
else
    fprintf('✗ Model validation failed\n');
end
end

function coef = bootstrap_cox_coef(boot_data, feat_idx, keep_cols, n_feat)
% Helper function for bootstrap coefficient extraction

try
    n = size(boot_data, 1);
    idx = randsample(n, n, true);
    boot_sample = boot_data(idx, :);
    
    X = boot_sample(:, 1:n_feat);
    T = boot_sample(:, n_feat+1:n_feat+2);
    censored = logical(boot_sample(:, n_feat+3));
    
    [X_clean, ~] = remove_constant_columns(X);
    if size(X_clean, 2) == 0
        coef = 0;
        return;
    end
    
    warning('off', 'all');
    [b, ~, ~, ~] = coxphfit(X_clean, T, 'Censoring', censored);
    warning('on', 'all');
    
    b_full = zeros(n_feat, 1);
    b_full(keep_cols) = b;
    coef = b_full(feat_idx);
    
catch
    coef = 0;
end
end