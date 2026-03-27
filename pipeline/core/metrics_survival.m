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
if isfield(config, 'fit_time_varying_cox') && ~config.fit_time_varying_cox
    fprintf('  Time-varying Cox analysis disabled by config (fit_time_varying_cox=false).\n');
    timevar_results = struct('success', false, 'message', 'Disabled by config');
else
    try
        timevar_results = compute_time_varying_effects(survival_data, cox_results, validation_results, config);
        timevar_success = timevar_results.success;
    catch ME_tv
        fprintf('  ⚠️  Time-varying effects analysis failed: %s\n', ME_tv.message);
        fprintf('      Identifier: %s\n', ME_tv.identifier);
        timevar_results = struct('success', false, 'message', ...
            sprintf('Exception: %s', ME_tv.message));
    end
end

% Fit Fine-Gray competing risks model
fprintf('\n--- FINE-GRAY COMPETING RISKS MODEL ---\n');
if isfield(config, 'compute_fine_gray') && ~config.compute_fine_gray
    fprintf('  Fine-Gray model disabled by config (compute_fine_gray=false).\n');
    finegray_results = struct('success', false, 'message', 'Disabled by config');
else
    try
        finegray_results = fit_fine_gray(survival_data, cox_results, config);
        finegray_success = finegray_results.success;
    catch ME_fg
        fprintf('  ⚠️  Fine-Gray analysis failed: %s\n', ME_fg.message);
        fprintf('      Identifier: %s\n', ME_fg.identifier);
        finegray_results = struct('success', false, 'message', ...
            sprintf('Exception: %s', ME_fg.message));
    end
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

% Propagate user config flags for optional analyses
if isfield(config_struct_in, 'compute_fine_gray')
    config.compute_fine_gray = config_struct_in.compute_fine_gray;
else
    config.compute_fine_gray = true;
end
if isfield(config_struct_in, 'fit_time_varying_cox')
    config.fit_time_varying_cox = config_struct_in.fit_time_varying_cox;
else
    config.fit_time_varying_cox = true;
end
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

function results = fit_fine_gray(survival_data, cox_results, config)
% FIT_FINE_GRAY  Fine-Gray subdistribution hazard model for competing risks.
%   Estimates the subdistribution hazard ratio (sHR) which accounts for
%   competing events by keeping competing-event subjects in the risk set
%   with decreasing weights w(t) = G(t)/G(t_comp), where G is the KM
%   estimate of the censoring distribution.
%
%   References
%   ----------
%   Fine JP, Gray RJ. A proportional hazards model for the subdistribution
%   of a competing risk. JASA. 1999;94(446):496-509.

results = struct('success', false, 'message', '');

if ~survival_data.has_sufficient_data
    results.message = 'Insufficient data';
    fprintf('  Insufficient data for Fine-Gray model.\n');
    return;
end

% Count event types
n_competing = sum(survival_data.event == 2);
n_primary = sum(survival_data.event == 1);

if n_primary < 3
    results.message = sprintf('Insufficient primary events (%d)', n_primary);
    fprintf('  Insufficient primary events (%d) for Fine-Gray model.\n', n_primary);
    return;
end

if n_competing == 0
    results.message = 'No competing events: Fine-Gray equivalent to CSH — skipped';
    fprintf('  No competing events detected (0 event=2). Fine-Gray reduces to CSH — skipping.\n');
    return;
end

fprintf('  Fitting Fine-Gray model: %d primary events, %d competing events.\n', ...
    n_primary, n_competing);

% --- Stage 1: Compute per-patient terminal outcomes ---
[unique_pats, ~, pat_group] = unique(survival_data.pat_id);
n_pats_fg = length(unique_pats);
pat_terminal_time = zeros(n_pats_fg, 1);
pat_terminal_event = zeros(n_pats_fg, 1);
for pi = 1:n_pats_fg
    pat_mask = (pat_group == pi);
    pat_rows = find(pat_mask);
    [~, max_idx] = max(survival_data.t_stop(pat_rows));
    terminal_row = pat_rows(max_idx);
    pat_terminal_time(pi) = survival_data.t_stop(terminal_row);
    pat_terminal_event(pi) = survival_data.event(terminal_row);
end

% --- Stage 2: Compute censoring KM G(t) ---
% For censoring KM: "event" = administrative censoring (event==0),
% all actual events (1 or 2) are "censored" in this reversed view.
cens_indicator = double(pat_terminal_event == 0);
[sorted_times, sort_idx] = sort(pat_terminal_time);
sorted_cens = cens_indicator(sort_idx);

uniq_times_fg = unique(sorted_times);
G_values = ones(length(uniq_times_fg), 1);
G_current = 1.0;
for k = 1:length(uniq_times_fg)
    tk = uniq_times_fg(k);
    n_at_risk = sum(sorted_times >= tk);
    n_censored = sum(sorted_times == tk & sorted_cens == 1);
    if n_at_risk > 0
        G_current = G_current * (1 - n_censored / n_at_risk);
    end
    G_values(k) = max(G_current, 0.01);  % floor to avoid div-by-zero
end

% --- Stage 3: Compute subdistribution weights ---
fg_weights = ones(size(survival_data.event));

for pi = 1:n_pats_fg
    if pat_terminal_event(pi) ~= 2
        continue;  % Only competing-event patients need modified weights
    end

    t_comp = pat_terminal_time(pi);
    G_t_comp = interp_G(t_comp, uniq_times_fg, G_values);

    pat_rows = find(pat_group == pi);
    for ri = 1:length(pat_rows)
        row_idx = pat_rows(ri);
        t_row = survival_data.t_stop(row_idx);
        if t_row > t_comp
            G_t = interp_G(t_row, uniq_times_fg, G_values);
            fg_weights(row_idx) = G_t / G_t_comp;
        end
    end
end

% --- Stage 4: Fit weighted Cox model ---
% Reuse CSH scaling if available
if cox_results.success && isfield(cox_results, 'X_scaled')
    X_fg = cox_results.X_scaled;
else
    X_fg = scale_td_panel(survival_data.X, survival_data.feat_names, ...
        survival_data.pat_id, survival_data.t_start, unique(survival_data.pat_id), 'baseline');
end

% Primary event (1) stays as event; competing (2) and censored (0)
% are both "censored" for the subdistribution hazard. Competing-event
% subjects remain in the risk set via the subdistribution weights.
event_fg = survival_data.event;
event_fg(event_fg == 2) = 0;

[X_fg_clean, keep_cols_fg] = remove_constant_columns(X_fg);
if size(X_fg_clean, 2) == 0
    results.message = 'All covariate columns are constant';
    fprintf('  All covariate columns constant — cannot fit Fine-Gray model.\n');
    return;
end

% Apply weights via Frequency (discretized, same pattern as IPCW)
fg_scale = 100;
fg_freq = max(1, round(fg_weights * fg_scale));

ws = warning('off', 'stats:coxphfit:FitWarning');
cleanupWarn = onCleanup(@() warning(ws));

try
    T_fg = [survival_data.t_start, survival_data.t_stop];
    is_censored_fg = (event_fg == 0);

    [b_fg, logl_fg, ~, stats_fg] = coxphfit(X_fg_clean, T_fg, ...
        'Censoring', is_censored_fg, 'Ties', 'efron', 'Frequency', fg_freq);

    % Sandwich SE for weighted model
    sandwich_se_fg = compute_sandwich_se_cox(b_fg, X_fg_clean, T_fg, ...
        is_censored_fg, fg_weights, survival_data.pat_id, 'efron');

    if ~isempty(sandwich_se_fg) && all(isfinite(sandwich_se_fg)) && all(sandwich_se_fg > 0)
        stats_fg.se = sandwich_se_fg;
        fprintf('  Using robust sandwich SE for Fine-Gray model.\n');
    end

    stats_fg.p = 2 * (1 - normcdf(abs(b_fg ./ stats_fg.se)));

    % Map to full feature space
    b_full_fg = zeros(survival_data.n_feat, 1);
    b_full_fg(keep_cols_fg) = b_fg;
    se_full_fg = nan(survival_data.n_feat, 1);
    p_full_fg = nan(survival_data.n_feat, 1);
    se_full_fg(keep_cols_fg) = stats_fg.se;
    p_full_fg(keep_cols_fg) = stats_fg.p;

    results.success = true;
    results.coefficients = b_full_fg;
    results.se = se_full_fg;
    results.p_values = p_full_fg;
    results.hazard_ratios = exp(b_full_fg);  % subdistribution HR (sHR)
    results.loglik = logl_fg;
    results.n_competing = n_competing;
    results.n_primary = n_primary;
    results.fg_weights = fg_weights;
catch ME_fit
    results.message = sprintf('Fine-Gray fit failed: %s', ME_fit.message);
    fprintf('  Fine-Gray model did not converge: %s\n', ME_fit.message);
    return;
end

% --- Stage 5: Print comparison table and CIF plot ---
print_fine_gray_comparison(cox_results, results, survival_data.feat_names);

if ~isempty(config.output_folder)
    try
        plot_cumulative_incidence(survival_data, config);
    catch ME_cif
        fprintf('  ⚠️  CIF plot failed: %s\n', ME_cif.message);
    end
end

fprintf('  Fine-Gray model completed: %d competing events, %d primary events.\n', ...
    n_competing, n_primary);
end

function results = compute_time_varying_effects(survival_data, cox_results, validation_results, config)
% COMPUTE_TIME_VARYING_EFFECTS  Assess time-varying coefficients.
%   Requires successful Cox PH fit and Schoenfeld residual results from
%   model validation to identify which covariates violate PH.
results = struct('success', false, 'message', 'Time-varying effects not computed');
if ~survival_data.has_sufficient_data
    results.message = 'Skipped: insufficient data';
    fprintf('  Time-varying analysis skipped: insufficient data.\n');
    return;
end
if ~cox_results.success
    results.message = sprintf('Skipped: Cox PH model did not succeed (%s)', ...
        get_cox_failure_reason(cox_results));
    fprintf('  💡 Cox PH model did not succeed — skipping time-varying analysis.\n');
    return;
end
if ~validation_results.success || ~isfield(validation_results, 'schoenfeld_results')
    val_reason = 'unknown';
    if isfield(validation_results, 'message')
        val_reason = validation_results.message;
    end
    results.message = sprintf('Skipped: Schoenfeld residuals not available (%s)', val_reason);
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
    fprintf('  ⚠️  Time-varying Cox failed: %s\n', ME.message);
end
end

function results = validate_survival_model(survival_data, cox_results, config)
% VALIDATE_SURVIVAL_MODEL  Model validation (calibration, discrimination).
%   Computes Schoenfeld residuals to test the PH assumption and stores
%   the full results struct for downstream time-varying Cox analysis.
results = struct('success', false, 'message', 'Model validation not computed');
if ~survival_data.has_sufficient_data
    results.message = 'Skipped: insufficient data for survival analysis';
    fprintf('  Model validation skipped: insufficient data.\n');
    return;
end
if ~cox_results.success
    results.message = sprintf('Skipped: Cox PH model did not succeed (%s)', ...
        get_cox_failure_reason(cox_results));
    fprintf('  💡 Cox PH model did not succeed — skipping model validation.\n');
    return;
end
try
    % Schoenfeld residuals for PH assumption
    if isfield(cox_results, 'coefficients') && isfield(cox_results, 'X_scaled')
        % Pre-check: verify kept coefficients are finite
        if isfield(cox_results, 'keep_cols') && ...
                any(~isfinite(cox_results.coefficients(cox_results.keep_cols)))
            results.message = 'Skipped: non-finite Cox coefficients for active covariates';
            fprintf('  Model validation skipped: non-finite Cox coefficients.\n');
            return;
        end

        schoenfeld = compute_schoenfeld_residuals( ...
            cox_results.X_scaled, survival_data.t_start, survival_data.t_stop, ...
            cox_results.event_csh, cox_results.coefficients, ...
            survival_data.feat_names, config.output_folder, config.dtype_label);
        results.schoenfeld_results = schoenfeld;
        results.ph_pvals = schoenfeld.p_value;
    else
        results.message = 'Skipped: Cox results missing coefficients or X_scaled';
        fprintf('  Model validation skipped: Cox results incomplete.\n');
        return;
    end
    results.success = true;
catch ME
    results = struct('success', false, 'message', ...
        sprintf('Schoenfeld residuals failed: %s', ME.message));
    fprintf('  ⚠️  Model validation failed: %s\n', ME.message);
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

function print_fine_gray_comparison(cox_results, fg_results, feat_names)
% PRINT_FINE_GRAY_COMPARISON  CSH HR vs Fine-Gray sHR comparison table.
fprintf('\n  --- CSH vs Fine-Gray Subdistribution Hazard Comparison ---\n');
fprintf('  %-10s %8s %8s %8s %8s  |  %8s %8s %8s %8s\n', ...
    'Feature', 'CSH_HR', 'CI_lo', 'CI_hi', 'p', 'sHR', 'CI_lo', 'CI_hi', 'p');
fprintf('  %s\n', repmat('-', 1, 86));

for fi = 1:numel(feat_names)
    % CSH results
    if cox_results.success
        b_csh = cox_results.coefficients(fi);
        se_csh = cox_results.se(fi);
        hr_csh = cox_results.hazard_ratios(fi);
        p_csh = cox_results.p_values(fi);
        ci_lo_csh = exp(b_csh - 1.96 * se_csh);
        ci_hi_csh = exp(b_csh + 1.96 * se_csh);
    else
        hr_csh = NaN; ci_lo_csh = NaN; ci_hi_csh = NaN; p_csh = NaN;
    end

    % Fine-Gray results
    b_fg = fg_results.coefficients(fi);
    se_fg = fg_results.se(fi);
    hr_fg = fg_results.hazard_ratios(fi);
    p_fg = fg_results.p_values(fi);
    ci_lo_fg = exp(b_fg - 1.96 * se_fg);
    ci_hi_fg = exp(b_fg + 1.96 * se_fg);

    fprintf('  %-10s %8.3f %8.3f %8.3f %8.4f  |  %8.3f %8.3f %8.3f %8.4f\n', ...
        feat_names{fi}, hr_csh, ci_lo_csh, ci_hi_csh, p_csh, ...
        hr_fg, ci_lo_fg, ci_hi_fg, p_fg);
end
end

function plot_cumulative_incidence(survival_data, config)
% PLOT_CUMULATIVE_INCIDENCE  Non-parametric CIF plot via Aalen-Johansen.

% Get per-patient terminal outcomes
[unique_pats, ~, pat_group] = unique(survival_data.pat_id);
n_pats = length(unique_pats);
pat_time = zeros(n_pats, 1);
pat_event = zeros(n_pats, 1);
for pi = 1:n_pats
    rows = find(pat_group == pi);
    [~, max_idx] = max(survival_data.t_stop(rows));
    pat_time(pi) = survival_data.t_stop(rows(max_idx));
    pat_event(pi) = survival_data.event(rows(max_idx));
end

[sorted_t, si] = sort(pat_time);
sorted_ev = pat_event(si);
uniq_t = unique(sorted_t);

% Aalen-Johansen estimator
S_prev = 1.0;
CIF_primary = zeros(length(uniq_t), 1);
CIF_competing = zeros(length(uniq_t), 1);
cum_primary = 0;
cum_competing = 0;

for k = 1:length(uniq_t)
    tk = uniq_t(k);
    n_risk = sum(sorted_t >= tk);
    d_primary = sum(sorted_t == tk & sorted_ev == 1);
    d_competing = sum(sorted_t == tk & sorted_ev == 2);

    if n_risk > 0
        h_primary = d_primary / n_risk;
        h_competing = d_competing / n_risk;
        cum_primary = cum_primary + S_prev * h_primary;
        cum_competing = cum_competing + S_prev * h_competing;
        S_prev = S_prev * (1 - (d_primary + d_competing) / n_risk);
    end

    CIF_primary(k) = cum_primary;
    CIF_competing(k) = cum_competing;
end

fig = figure('Visible', 'off', 'Position', [100 100 700 500]);
stairs(uniq_t, CIF_primary, 'b-', 'LineWidth', 2);
hold on;
stairs(uniq_t, CIF_competing, 'r--', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Cumulative Incidence');
title(sprintf('Cumulative Incidence Functions (%s)', config.dtype_label));
legend('Local Failure', 'Competing Events', 'Location', 'northwest');
grid on;
ylim([0, 1]);

fig_path = fullfile(config.output_folder, sprintf('cif_plot_%s.png', config.dtype_label));
saveas(fig, fig_path);
close(fig);
fprintf('  📁 CIF plot saved: %s\n', fig_path);
end

function G_val = interp_G(t, uniq_times, G_values)
% INTERP_G  Step-function lookup for censoring survival G(t).
    idx = find(uniq_times <= t, 1, 'last');
    if isempty(idx)
        G_val = 1.0;
    else
        G_val = G_values(idx);
    end
end

function reason = get_cox_failure_reason(cox_results)
% GET_COX_FAILURE_REASON  Extract failure reason from Cox results.
    if isfield(cox_results, 'message') && ~isempty(cox_results.message)
        reason = cox_results.message;
    else
        reason = 'unknown';
    end
end