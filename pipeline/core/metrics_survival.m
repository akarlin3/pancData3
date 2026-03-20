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

% Determine landmark day
% Priority 1: Use explicitly configured landmark_day (clinically motivated)
% Priority 2: Fall back to gap-based heuristic with a warning
if isfield(config_struct_in, 'landmark_day') && ~isempty(config_struct_in.landmark_day)
    landmark_day = config_struct_in.landmark_day;
    fprintf('  Using configured landmark_day = %d (from config).\n', landmark_day);
else
    % Gap-based heuristic fallback: choose the scan day preceding the largest
    % gap between consecutive scans. This is a fragile proxy for "end of
    % treatment" and should be replaced with a clinically motivated value.
    scan_gaps = diff(td_scan_days);
    [~, gap_idx] = max(scan_gaps);
    landmark_day = td_scan_days(gap_idx);
    fprintf('  ⚠️  WARNING: config.landmark_day not set. Falling back to gap-based heuristic (landmark_day = %d).\n', landmark_day);
    fprintf('      This heuristic is fragile and may not reflect a clinically meaningful timepoint.\n');
    fprintf('      Set config.landmark_day to the last intra-treatment scan day for robust results.\n');
end

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
% Uses robust sandwich (Huber-White) variance estimator for IPCW-weighted models.

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

% Fit Cox model with robust sandwich variance estimation
try
    T_matrix = [survival_data.t_start, survival_data.t_stop];
    is_censored = (event_csh == 0);
    
    [X_clean, keep_cols] = remove_constant_columns(X_scaled);
    if size(X_clean, 2) == 0
        error('Cox:NoVariableColumns', 'All covariate columns are constant.');
    end
    
    % --- Step 1: Fit unweighted Cox model to get initial estimates ---
    % We first fit the unweighted model, then compute the sandwich variance
    % using the IPCW weights and score residuals.
    
    warning('off', 'all');
    warning('error', 'stats:coxphfit:FitWarning');
    warning('error', 'stats:coxphfit:IterationLimit');
    
    % Apply IPCW via frequency weights for point estimation
    ipcw_scale = 100;
    ipcw_freq = max(1, round(ipcw_weights * ipcw_scale));
    
    [b_short, logl, H_out, stats_short] = coxphfit(X_clean, T_matrix, ...
        'Censoring', is_censored, 'Ties', ties_method, 'Frequency', ipcw_freq);
    
    % --- Step 2: Compute robust sandwich (Huber-White) standard errors ---
    % The sandwich variance is: V_robust = A^{-1} B A^{-1}
    % where A = negative Hessian (observed information), B = meat matrix
    % B = sum over clusters (patients) of (score_i * w_i)' * (score_i * w_i)
    
    fprintf('  Computing robust sandwich variance estimator for IPCW weights.\n');
    
    sandwich_se = compute_sandwich_se_cox(b_short, X_clean, T_matrix, ...
        is_censored, ipcw_weights, survival_data.pat_id, ties_method);
    
    if ~isempty(sandwich_se) && all(isfinite(sandwich_se)) && all(sandwich_se > 0)
        % Use sandwich SE
        stats_short.se = sandwich_se;
        fprintf('  Using robust sandwich SE estimator.\n');
    else
        % Fallback: use clustered bootstrap SE stratified by patient ID
        fprintf('  Sandwich SE computation failed; falling back to clustered bootstrap.\n');
        bootstrap_se = compute_bootstrap_se_cox(X_clean, T_matrix, is_censored, ...
            ipcw_weights, survival_data.pat_id, ties_method, config.n_bootstrap);
        
        if ~isempty(bootstrap_se) && all(isfinite(bootstrap_se)) && all(bootstrap_se > 0)
            stats_short.se = bootstrap_se;
            fprintf('  Using clustered bootstrap SE (%d replicates).\n', config.n_bootstrap);
        else
            % Last resort: keep model-based SE with a warning
            fprintf('  ⚠️  Both sandwich and bootstrap SE failed; using model-based SE (interpret with caution).\n');
        end
    end
    
    % Recompute p-values with corrected SE
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

function sandwich_se = compute_sandwich_se_cox(beta, X, T_matrix, is_censored, ipcw_weights, pat_id, ties_method)
% COMPUTE_SANDWICH_SE_COX Compute robust sandwich (Huber-White) SE for IPCW-weighted Cox model.
%
% The sandwich estimator accounts for the IPCW weighting by computing:
%   V_sandwich = A^{-1} * B * A^{-1}
% where A is the observed information (negative Hessian) and B is the "meat"
% matrix formed from weighted score residuals clustered by patient.

try
    n = size(X, 1);
    p = size(X, 2);
    t_stop = T_matrix(:, 2);
    event = double(~is_censored);  % 1 = event, 0 = censored
    
    % Compute linear predictor and risk scores
    eta = X * beta;
    risk = exp(eta);
    w = ipcw_weights(:);  % per-observation IPCW weights
    
    % Get unique ordered event times
    event_times = sort(unique(t_stop(event == 1)));
    n_events = length(event_times);
    
    if n_events == 0
        sandwich_se = [];
        return;
    end
    
    % --- Compute score residuals for each observation ---
    % Score residual for observation i:
    %   U_i = delta_i * w_i * (x_i - xbar(t_i)) 
    %         - w_i * sum_{j: t_j <= t_i, delta_j=1} [ w_j * risk_i * (x_i - xbar(t_j)) / S0(t_j) ] * I(i in risk set at t_j)
    % where S0(t) = sum_{k in R(t)} w_k * risk_k
    %       xbar(t) = sum_{k in R(t)} w_k * risk_k * x_k / S0(t)
    
    score_resid = zeros(n, p);
    
    % Precompute risk set quantities at each event time
    S0 = zeros(n_events, 1);
    S1 = zeros(n_events, p);
    
    for j = 1:n_events
        tj = event_times(j);
        in_risk = (T_matrix(:, 1) < tj) & (t_stop >= tj);
        if ~any(in_risk), continue; end
        
        w_risk = w(in_risk) .* risk(in_risk);
        S0(j) = sum(w_risk);
        S1(j, :) = sum(bsxfun(@times, X(in_risk, :), w_risk), 1);
    end
    
    % Avoid division by zero
    S0(S0 == 0) = eps;
    xbar = bsxfun(@rdivide, S1, S0);  % n_events x p
    
    % For each observation, accumulate score residual contributions
    for i = 1:n
        ti = t_stop(i);
        
        % Event contribution (if this observation had an event)
        if event(i) == 1
            % Find which event time this corresponds to
            eidx = find(event_times == ti, 1);
            if ~isempty(eidx)
                score_resid(i, :) = score_resid(i, :) + w(i) * (X(i, :) - xbar(eidx, :));
            end
        end
        
        % At-risk contribution: for each event time t_j where i is in the risk set
        for j = 1:n_events
            tj = event_times(j);
            if T_matrix(i, 1) < tj && t_stop(i) >= tj
                % Number of events at this time (for ties)
                n_events_at_tj = sum(event == 1 & t_stop == tj);
                
                % Contribution from being in risk set
                contrib = w(i) * risk(i) * n_events_at_tj * (X(i, :) - xbar(j, :)) / S0(j);
                score_resid(i, :) = score_resid(i, :) - w(i) * contrib;
            end
        end
    end
    
    % --- Compute the Hessian (observed information matrix A) ---
    % A = sum over event times of [ S2(t)/S0(t) - (S1(t)/S0(t))' * (S1(t)/S0(t)) ]
    % weighted by number of events at each time
    A = zeros(p, p);
    for j = 1:n_events
        tj = event_times(j);
        in_risk = (T_matrix(:, 1) < tj) & (t_stop >= tj);
        if ~any(in_risk), continue; end
        
        w_risk = w(in_risk) .* risk(in_risk);
        X_risk = X(in_risk, :);
        
        % S2 matrix
        S2 = (bsxfun(@times, X_risk, w_risk))' * X_risk;
        
        n_events_at_tj = sum(event == 1 & t_stop == tj);
        
        A = A + n_events_at_tj * (S2 / S0(j) - xbar(j, :)' * xbar(j, :));
    end
    
    % --- Compute meat matrix B, clustered by patient ---
    unique_pats = unique(pat_id);
    B = zeros(p, p);
    for k = 1:length(unique_pats)
        idx_k = (pat_id == unique_pats(k));
        U_k = sum(score_resid(idx_k, :), 1);  % 1 x p: sum score residuals within cluster
        B = B + U_k' * U_k;
    end
    
    % --- Sandwich covariance: A^{-1} * B * A^{-1} ---
    A_inv = pinv(A);  % Use pseudoinverse for numerical stability
    V_sandwich = A_inv * B * A_inv;
    
    % Extract SE from diagonal
    sandwich_var = diag(V_sandwich);
    
    if any(sandwich_var < 0)
        fprintf('  ⚠️  Negative sandwich variance detected; attempting nearest PSD correction.\n');
        % Force positive semi-definiteness
        [V_eig, D_eig] = eig(V_sandwich);
        D_eig = diag(max(diag(D_eig), 0));
        V_sandwich = V_eig * D_eig * V_eig';
        sandwich_var = diag(V_sandwich);
    end
    
    sandwich_se = sqrt(sandwich_var);
    
    if any(~isfinite(sandwich_se))
        sandwich_se = [];
    end
    
catch ME
    fprintf('  Sandwich SE computation error: %s\n', ME.message);
    sandwich_se = [];
end
end

function bootstrap_se = compute_bootstrap_se_cox(X, T_matrix, is_censored, ipcw_weights, pat_id, ties_method, n_boot)
% COMPUTE_BOOTSTRAP_SE_COX Clustered bootstrap SE estimation stratified by patient ID.
%
% Resamples patients (with replacement), refits the IPCW-weighted Cox model,
% and computes the empirical SE of the bootstrap coefficient distribution.

try
    unique_pats = unique(pat_id);
    n_pats = length(unique_pats);
    p = size(X, 2);
    
    % Limit bootstrap replicates for computational feasibility
    n_boot_actual = min(n_boot, 500);
    boot_coefs = nan(p, n_boot_actual);
    
    ipcw_scale = 100;
    
    for b = 1:n_boot_actual
        % Resample patients with replacement
        rng_state = rng;  %#ok<NASGU>
        boot_pat_idx = randsample(n_pats, n_pats, true);
        
        % Build bootstrap dataset
        boot_rows = [];
        for k = 1:n_pats
            rows_k = find(pat_id == unique_pats(boot_pat_idx(k)));
            boot_rows = [boot_rows; rows_k]; %#ok<AGROW>
        end
        
        X_b = X(boot_rows, :);
        T_b = T_matrix(boot_rows, :);
        cens_b = is_censored(boot_rows);
        w_b = ipcw_weights(boot_rows);
        freq_b = max(1, round(w_b * ipcw_scale));
        
        % Remove constant columns in bootstrap sample
        col_var = var(X_b, 0, 1);
        var_cols = col_var > eps;
        if sum(var_cols) < p
            % Skip this replicate if we lose columns
            continue;
        end
        
        try
            warning('off', 'all');
            [b_boot, ~, ~, ~] = coxphfit(X_b, T_b, ...
                'Censoring', cens_b, 'Ties', ties_method, 'Frequency', freq_b);
            boot_coefs(:, b) = b_boot;
            warning('on', 'all');
        catch
            warning('on', 'all');
            continue;
        end
    end
    
    % Compute SE from valid bootstrap replicates
    valid_boots = all(isfinite(boot_coefs), 1);
    n_valid = sum(valid_boots);
    
    if n_valid >= 50
        bootstrap_se = std(boot_coefs(:, valid_boots), 0, 2);
        fprintf('  Bootstrap: %d/%d valid replicates.\n', n_valid, n_boot_actual);
    else
        fprintf('  Bootstrap: only %d valid replicates (need >= 50).\n', n_valid);
        bootstrap_se = [];
    end
    
catch ME
    fprintf('  Bootstrap SE error: %s\n', ME.message);
    bootstrap_se = [];
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

%