function metrics_survival(varargin)
% METRICS_SURVIVAL — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 5/5 of the metrics step. Fits survival models including Cox PH and competing risks analysis.
%
% This function coordinates the entire survival analysis workflow by calling
% specialized subfunctions for different aspects of the analysis.
%
% NEW INTERFACE (recommended):
%   metrics_survival(data_struct, opts)
%     data_struct — struct with fields:
%       .valid_pts, .ADC_abs, .D_abs, .f_abs, .Dstar_abs,
%       .m_lf, .m_total_time, .m_total_follow_up_time, .nTp,
%       .fx_label, .dtype_label
%     Optional fields:
%       .m_gtv_vol, .actual_scan_days
%     opts — struct with configuration options (passed through as config).
%            Any fields not present receive defaults.
%            Recognized config fields:
%              .td_halflife_days — Exponential decay half-life (in days) used
%                  by build_td_panel to weight imaging features from prior
%                  fractions. Controls how quickly older measurements lose
%                  influence when constructing the time-dependent covariate
%                  panel. Default: 18 days.
%                  Clinical guidance: For a standard 5-fraction SBRT course
%                  spanning ~20 days, 18 days means Fx1 features retain ~46%
%                  weight at Fx5. Shorter half-lives (e.g., 6–12 days) are
%                  appropriate when rapid biological changes are expected;
%                  longer half-lives (e.g., 24–36 days) suit slowly evolving
%                  responses. This parameter can be tuned via cross-validation
%                  over the half_life_grid in the time-varying effects analysis.
%
% LEGACY INTERFACE (deprecated, preserved for backward compatibility):
%   metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ...
%       m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, ...
%       dtype_label, m_gtv_vol, output_folder, actual_scan_days, config_struct_in)

% --- Argument parsing: detect new vs legacy interface ---
if nargin >= 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'valid_pts')
    % ---- New struct-based interface ----
    data_struct = varargin{1};
    if nargin >= 2 && isstruct(varargin{2})
        opts = varargin{2};
    else
        opts = struct();
    end

    % Unpack required fields
    valid_pts              = data_struct.valid_pts;
    ADC_abs                = data_struct.ADC_abs;
    D_abs                  = data_struct.D_abs;
    f_abs                  = data_struct.f_abs;
    Dstar_abs              = data_struct.Dstar_abs;
    m_lf                   = data_struct.m_lf;
    m_total_time           = data_struct.m_total_time;
    m_total_follow_up_time = data_struct.m_total_follow_up_time;
    nTp                    = data_struct.nTp;
    fx_label               = data_struct.fx_label;
    dtype_label            = data_struct.dtype_label;

    % Optional data fields
    if isfield(data_struct, 'm_gtv_vol')
        m_gtv_vol = data_struct.m_gtv_vol;
    else
        m_gtv_vol = [];
    end
    if isfield(data_struct, 'actual_scan_days')
        actual_scan_days = data_struct.actual_scan_days;
    else
        actual_scan_days = [];
    end

    % Config: opts acts as the config struct; extract output_folder from it
    if isfield(opts, 'output_folder')
        output_folder = opts.output_folder;
    else
        output_folder = '';
    end
    config_struct_internal = opts;

else
    % ---- Legacy positional interface (deprecated) ----
    warning('metrics_survival:DeprecatedSignature', ...
        ['The 15-argument positional interface is deprecated and will be removed in a future release. ', ...
         'Use metrics_survival(data_struct, opts) instead.']);

    if nargin < 11
        error('metrics_survival:InsufficientArguments', ...
            'At least 11 arguments are required for the legacy interface.');
    end

    valid_pts              = varargin{1};
    ADC_abs                = varargin{2};
    D_abs                  = varargin{3};
    f_abs                  = varargin{4};
    Dstar_abs              = varargin{5};
    m_lf                   = varargin{6};
    m_total_time           = varargin{7};
    m_total_follow_up_time = varargin{8};
    nTp                    = varargin{9};
    fx_label               = varargin{10};
    dtype_label            = varargin{11};

    if nargin >= 12, m_gtv_vol = varargin{12}; else, m_gtv_vol = []; end
    if nargin >= 13, output_folder = varargin{13}; else, output_folder = ''; end
    if nargin >= 14, actual_scan_days = varargin{14}; else, actual_scan_days = []; end
    if nargin >= 15 && ~isempty(varargin{15})
        config_struct_internal = varargin{15};
    else
        config_struct_internal = struct();
    end
end

fprintf('\n--- SURVIVAL ANALYSIS PIPELINE ---\n');

% Set up diary for output capture.
% We use a two-phase approach: first perform all filesystem operations that
% could fail, then start the diary and create the cleanup guard. This
% ensures that if file deletion or diary activation throws, the diary is
% not left in an inconsistent state.
cleanupDiary = []; %#ok<NASGU> — will hold onCleanup handle if diary is started
if ~isempty(output_folder)
    diary_file = fullfile(output_folder, ['metrics_survival_output_' dtype_label '.txt']);
    % Phase 1: filesystem checks (may throw on read-only filesystem)
    try
        if exist(diary_file, 'file')
            delete(diary_file);
        end
    catch ME_del
        warning('metrics_survival:DiaryDeleteFailed', ...
            'Could not delete existing diary file: %s. Proceeding without diary.', ME_del.message);
        diary_file = '';
    end
    % Phase 2: start diary only if filesystem prep succeeded
    if ~isempty(diary_file)
        try
            diary(diary_file);
            cleanupDiary = onCleanup(@() diary('off')); %#ok<NASGU>
        catch ME_diary
            % diary() itself failed — ensure it is off and warn
            try
                diary('off');
            catch
                % diary('off') may also fail; nothing more we can do
            end
            warning('metrics_survival:DiaryStartFailed', ...
                'Could not start diary: %s. Proceeding without diary.', ME_diary.message);
        end
    end
end

% Prepare common data structures using the extracted utility functions.
% Previously this was a single local function (prepare_survival_data);
% now the logic is split across independently testable functions in
% pipeline/utils/:
%   prepare_outcome_data  — censoring logic
%   select_landmark_day   — landmark day selection
%   build_survival_features — feature assembly and panel construction
[survival_data, config] = prepare_survival_data_dispatch(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ...
    m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label, ...
    m_gtv_vol, actual_scan_days, config_struct_internal, output_folder);

% Check if we have sufficient data for analysis
if ~survival_data.has_sufficient_data
    fprintf('  Insufficient data for survival analysis.\n');
    return;
end

% Initialize success/failure flags for each analysis
cox_success = false;
finegray_success = false;
timevar_success = false;
validation_success = false;

% Initialize all result variables before the try-catch blocks so that they
% are always defined in the outer scope regardless of which analyses
% succeed or fail. This prevents 'Undefined function or variable' errors
% in print_survival_summary.
cox_results        = struct('success', false);
finegray_results   = struct('success', false);
timevar_results    = struct('success', false);
validation_results = struct('success', false);

% Fit Cox Proportional Hazards model
fprintf('\n--- COX PROPORTIONAL HAZARDS MODEL ---\n');
try
    cox_results = fit_cox_ph(survival_data, config);
    cox_success = true;
catch ME_cox
    fprintf('  ⚠️  Cox PH analysis failed: %s\n', ME_cox.message);
    fprintf('      Identifier: %s\n', ME_cox.identifier);
    cox_results = struct('success', false, 'message', ...
        sprintf('Exception: %s', ME_cox.message));
end

% Validate survival models (Schoenfeld residuals for PH assumption)
% Run validation before time-varying analysis because the Schoenfeld
% results identify which covariates violate PH — needed by the
% time-varying Cox model.
fprintf('\n--- MODEL VALIDATION ---\n');
try
    validation_results = validate_survival_model(survival_data, cox_results, config);
    validation_success = validation_results.success;
catch ME_val
    fprintf('  ⚠️  Model validation failed: %s\n', ME_val.message);
    fprintf('      Identifier: %s\n', ME_val.identifier);
    validation_results = struct('success', false, 'message', ...
        sprintf('Exception: %s', ME_val.message));
end

% Compute time-varying effects analysis (requires Schoenfeld results
% from validation to know which covariates violate PH)
fprintf('\n--- TIME-VARYING COEFFICIENTS ANALYSIS ---\n');
try
    timevar_results = compute_time_varying_effects(survival_data, cox_results, validation_results, config);
    timevar_success = timevar_results.success;
catch ME_tv
    fprintf('  ⚠️  Time-varying effects analysis failed: %s\n', ME_tv.message);
    fprintf('      Identifier: %s\n', ME_tv.identifier);
    timevar_results = struct('success', false, 'message', ...
        sprintf('Exception: %s', ME_tv.message));
end

% Fit Fine-Gray competing risks model
fprintf('\n--- FINE-GRAY COMPETING RISKS MODEL ---\n');
try
    finegray_results = fit_fine_gray(survival_data, config);
    finegray_success = finegray_results.success;
catch ME_fg
    fprintf('  ⚠️  Fine-Gray analysis failed: %s\n', ME_fg.message);
    fprintf('      Identifier: %s\n', ME_fg.identifier);
    finegray_results = struct('success', false, 'message', ...
        sprintf('Exception: %s', ME_fg.message));
end

% Output summary results, including success/failure status for each analysis
analysis_status = struct();
analysis_status.cox_success = cox_success;
analysis_status.finegray_success = finegray_success;
analysis_status.timevar_success = timevar_success;
analysis_status.validation_success = validation_success;

print_survival_summary(cox_results, finegray_results, timevar_results, ...
    validation_results, analysis_status);

end

function [survival_data, config] = prepare_survival_data_dispatch(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ...
    m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label, ...
    m_gtv_vol, actual_scan_days, config_struct_in, output_folder)
% PREPARE_SURVIVAL_DATA_DISPATCH Orchestrate data preparation by delegating
% to independently testable utility functions:
%   build_survival_features  — feature array assembly
%   prepare_outcome_data     — censoring/competing risk logic
%   select_landmark_day      — landmark day selection

% --- 1. Scan days ---
if ~isempty(actual_scan_days)
    td_scan_days = actual_scan_days;
    fprintf('  Using provided scan days [%s].\n', num2str(td_scan_days));
else
    td_scan_days = [0, 5, 10, 15, 20, 90];
    fprintf('  ⚠️  CAUTION: Using default scan days [%s].\n', num2str(td_scan_days));
    fprintf('      Pass actual DICOM-derived scan days to avoid immortal time bias.\n');
end

% --- 2. Build feature arrays (delegated to build_survival_features) ---
[td_feat_arrays, td_feat_names] = build_survival_features(valid_pts, ...
    ADC_abs, D_abs, f_abs, Dstar_abs, m_gtv_vol, nTp);

% --- 3. Prepare outcome data with censoring logic (delegated to prepare_outcome_data) ---
[td_lf, td_tot_time] = prepare_outcome_data(valid_pts, m_lf, m_total_time, m_total_follow_up_time);

% --- 4. Parse half-life config ---
if isfield(config_struct_in, 'td_halflife_days') && ~isempty(config_struct_in.td_halflife_days)
    td_halflife_days = config_struct_in.td_halflife_days;
    fprintf('  Using configured td_halflife_days = %.1f (from config).\n', td_halflife_days);
else
    td_halflife_days = 18;
    fprintf('  Using default td_halflife_days = %.1f.\n', td_halflife_days);
    fprintf('      Set config.td_halflife_days to customize the decay half-life for time-dependent feature weighting.\n');
end

% --- 5. Build time-dependent panel ---
[X_td_def, t_start_td_def, t_stop_td_def, event_td_def, pat_id_td_def, frac_td_def] = ...
    build_td_panel(td_feat_arrays, td_feat_names, td_lf, td_tot_time, nTp, td_scan_days, td_halflife_days);

% --- 6. Check initial data sufficiency ---
td_n_feat = numel(td_feat_arrays);
td_ok = (sum(event_td_def == 1) >= 3) && (size(X_td_def, 1) > td_n_feat + 1);

% --- 7. Select landmark day (delegated to select_landmark_day) ---
landmark_day = select_landmark_day(td_scan_days, config_struct_in);

% --- 8. Apply landmark filtering ---
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

% --- 9. Package results ---
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
config.td_halflife_days = td_halflife_days;
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
    
    % Suppress only the specific coxphfit warnings, with guaranteed restoration
    ws = warning('off', 'stats:coxphfit:FitWarning');
    cleanupWarn = onCleanup(@() warning(ws));
    
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
    boot_coeffs = nan(p, n_boot_actual);

    for b = 1:n_boot_actual
        % Resample patients with replacement (clustered bootstrap)
        boot_pat_idx = randsample(n_pats, n_pats, true);
        boot_rows = [];
        for bp = 1:n_pats
            boot_rows = [boot_rows; find(pat_id == unique_pats(boot_pat_idx(bp)))]; %#ok<AGROW>
        end

        X_boot = X(boot_rows, :);
        T_boot = T_matrix(boot_rows, :);
        cens_boot = is_censored(boot_rows);
        w_boot = ipcw_weights(boot_rows);

        % Skip degenerate bootstrap samples
        if length(unique(cens_boot)) < 2
            continue;
        end

        try
            ipcw_freq_boot = max(1, round(w_boot * 100));
            [beta_boot, ~, ~, ~] = coxphfit(X_boot, T_boot, ...
                'Censoring', cens_boot, 'Ties', ties_method, 'Frequency', ipcw_freq_boot);
            boot_coeffs(:, b) = beta_boot;
        catch
            % Skip failed bootstrap replicate
        end
    end

    % Compute SE from non-NaN bootstrap replicates
    valid_boots = ~any(isnan(boot_coeffs), 1);
    if sum(valid_boots) < 10
        bootstrap_se = [];
    else
        bootstrap_se = std(boot_coeffs(:, valid_boots), 0, 2);
    end

catch ME
    fprintf('  Bootstrap SE computation error: %s\n', ME.message);
    bootstrap_se = [];
end
end

function [td_feat_arrays, td_feat_names] = build_survival_features(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_gtv_vol, nTp)
% BUILD_SURVIVAL_FEATURES Assemble feature arrays for time-dependent Cox model.
td_feat_arrays = { ADC_abs(valid_pts,:), D_abs(valid_pts,:), ...
                   f_abs(valid_pts,:),   Dstar_abs(valid_pts,:) };
td_feat_names  = {'ADC', 'D', 'f', 'D*'};

% Include baseline GTV volume as a time-constant confounder when available.
has_vol = ~isempty(m_gtv_vol) && any(isfinite(m_gtv_vol(valid_pts, 1)));
if has_vol
    vol_baseline = m_gtv_vol(valid_pts, 1);
    vol_rep = repmat(vol_baseline, 1, nTp);
    td_feat_arrays{end+1} = vol_rep;
    td_feat_names{end+1}  = 'GTVvol';
    fprintf('  Including baseline GTV volume as a covariate.\n');
end
end

function [td_lf, td_tot_time] = prepare_outcome_data(valid_pts, m_lf, m_total_time, m_total_follow_up_time)
% PREPARE_OUTCOME_DATA Prepare outcome data with censoring logic.
td_lf       = m_lf(valid_pts);
td_tot_time = m_total_time(valid_pts);

% Censored and competing-risk patients use follow-up time
follow_up_valid = m_total_follow_up_time(valid_pts);
cens_mask_td = (td_lf == 0 | td_lf == 2) & ~isnan(follow_up_valid);
td_tot_time(cens_mask_td) = follow_up_valid(cens_mask_td);
end

function landmark_day = select_landmark_day(td_scan_days, config_struct_in)
% SELECT_LANDMARK_DAY Select the landmark day for left-truncation.
if isfield(config_struct_in, 'landmark_day') && ~isempty(config_struct_in.landmark_day)
    landmark_day = config_struct_in.landmark_day;
else
    % Default: use the second scan day (first post-baseline)
    if numel(td_scan_days) >= 2
        landmark_day = td_scan_days(2);
    else
        landmark_day = 0;
    end
end
end

function results = fit_fine_gray(survival_data, config)
% FIT_FINE_GRAY  Fine-Gray competing risks sub-distribution hazard model.
%   Not yet implemented — compute_fine_gray does not exist in the codebase.
%   Returns a clean skip result instead of silently failing.
results = struct('success', false, 'message', 'Fine-Gray model not yet implemented');
fprintf('  💡 Fine-Gray competing risks model is not yet implemented — skipping.\n');
end

function results = compute_time_varying_effects(survival_data, cox_results, validation_results, config)
% COMPUTE_TIME_VARYING_EFFECTS  Assess time-varying coefficients.
%   Requires successful Cox PH fit and Schoenfeld residual results from
%   model validation to identify which covariates violate PH.
results = struct('success', false, 'message', 'Time-varying effects not computed');
if ~survival_data.has_sufficient_data
    return;
end
if ~cox_results.success
    results.message = 'Skipped: Cox PH model did not succeed';
    fprintf('  💡 Cox PH model did not succeed — skipping time-varying analysis.\n');
    return;
end
if ~validation_results.success || ~isfield(validation_results, 'schoenfeld_results')
    results.message = 'Skipped: Schoenfeld residuals not available';
    fprintf('  💡 Schoenfeld residuals not available — skipping time-varying analysis.\n');
    return;
end
try
    results = fit_time_varying_cox(cox_results.X_scaled, survival_data.t_start, ...
        survival_data.t_stop, cox_results.event_csh, survival_data.feat_names, ...
        validation_results.schoenfeld_results, config.output_folder, ...
        config.dtype_label, config);
    results.success = true;
catch ME
    results = struct('success', false, 'message', ME.message);
end
end

function results = validate_survival_model(survival_data, cox_results, config)
% VALIDATE_SURVIVAL_MODEL  Model validation (calibration, discrimination).
%   Computes Schoenfeld residuals to test the PH assumption and stores
%   the full results struct for downstream time-varying Cox analysis.
results = struct('success', false, 'message', 'Model validation not computed');
if ~survival_data.has_sufficient_data
    results.message = 'Skipped: insufficient data';
    return;
end
if ~cox_results.success
    results.message = 'Skipped: Cox PH model did not succeed';
    fprintf('  💡 Cox PH model did not succeed — skipping model validation.\n');
    return;
end
try
    % Schoenfeld residuals for PH assumption
    if isfield(cox_results, 'coefficients') && isfield(cox_results, 'X_scaled')
        schoenfeld = compute_schoenfeld_residuals( ...
            cox_results.X_scaled, survival_data.t_start, survival_data.t_stop, ...
            cox_results.event_csh, cox_results.coefficients, ...
            survival_data.feat_names, config.output_folder, config.dtype_label);
        results.schoenfeld_results = schoenfeld;
        results.ph_pvals = schoenfeld.p_value;
    end
    results.success = true;
catch ME
    results = struct('success', false, 'message', ME.message);
end
end

function print_survival_summary(cox_results, finegray_results, timevar_results, ...
    validation_results, analysis_status)
% PRINT_SURVIVAL_SUMMARY  Print a compact summary table of all survival analyses.
fprintf('\n===== SURVIVAL ANALYSIS SUMMARY =====\n');
if analysis_status.cox_success && cox_results.success
    fprintf('  Cox PH model:           SUCCESS\n');
else
    fprintf('  Cox PH model:           FAILED\n');
end
if analysis_status.finegray_success && isfield(finegray_results, 'success') && finegray_results.success
    fprintf('  Fine-Gray model:        SUCCESS\n');
else
    fprintf('  Fine-Gray model:        FAILED/SKIPPED\n');
end
if analysis_status.timevar_success && isfield(timevar_results, 'success') && timevar_results.success
    fprintf('  Time-varying effects:   SUCCESS\n');
else
    fprintf('  Time-varying effects:   FAILED/SKIPPED\n');
end
if analysis_status.validation_success && isfield(validation_results, 'success') && validation_results.success
    fprintf('  Model validation:       SUCCESS\n');
else
    fprintf('  Model validation:       FAILED/SKIPPED\n');
end
fprintf('=====================================\n');
end

function print_cox_results(cox_results, feat_names)
% PRINT_COX_RESULTS  Print Cox PH regression results table.
if ~cox_results.success
    fprintf('  Cox model did not converge.\n');
    return;
end
fprintf('\n  %-10s %8s %8s %12s %8s\n', 'Feature', 'Coeff', 'HR', '95% CI', 'p-value');
fprintf('  %s\n', repmat('-', 1, 52));
for fi = 1:numel(feat_names)
    b = cox_results.coefficients(fi);
    hr = cox_results.hazard_ratios(fi);
    se = cox_results.se(fi);
    p = cox_results.p_values(fi);
    ci_lo = exp(b - 1.96 * se);
    ci_hi = exp(b + 1.96 * se);
    fprintf('  %-10s %8.3f %8.3f [%5.2f-%5.2f] %8.4f\n', ...
        feat_names{fi}, b, hr, ci_lo, ci_hi, p);
end
end

function cox_results = create_failed_cox_results(n_feat)
% CREATE_FAILED_COX_RESULTS  Return a failure struct with correct dimensions.
cox_results = struct();
cox_results.success = false;
cox_results.message = 'Cox model failed to converge';
cox_results.coefficients = nan(n_feat, 1);
cox_results.se = nan(n_feat, 1);
cox_results.p_values = nan(n_feat, 1);
cox_results.hazard_ratios = nan(n_feat, 1);
end