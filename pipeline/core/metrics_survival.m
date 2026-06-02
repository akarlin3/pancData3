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

% ============================================================
% SECTION INDEX
% ------------------------------------------------------------
%  Main function:
%    Argument parsing (new struct vs legacy positional interface)
%    Diary setup (two-phase: filesystem ops, then activation)
%    Prepare survival data (dispatch to subfunctions)
%    Fit Cox PH model (with IPCW and robust sandwich SE)
%    Validate Cox model (Schoenfeld residual diagnostics)
%    Time-varying effects analysis
%    Fit Fine-Gray competing risks model
%    Output summary results
%  Local functions (in this file):
%    build_survival_features
%    prepare_outcome_data
%    select_landmark_day
%    compute_time_varying_effects
%    validate_survival_model
%    print_survival_summary
%    print_cox_results
%    create_failed_cox_results
%    get_cox_failure_reason
%  Subfunctions extracted to their own files in pipeline/core/:
%    prepare_survival_data_dispatch, fit_cox_ph, compute_sandwich_se_cox,
%    compute_bootstrap_se_cox, fit_fine_gray, print_fine_gray_comparison,
%    plot_cumulative_incidence, interp_G
% ============================================================
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

function reason = get_cox_failure_reason(cox_results)
% GET_COX_FAILURE_REASON  Extract failure reason from Cox results.
    if isfield(cox_results, 'message') && ~isempty(cox_results.message)
        reason = cox_results.message;
    else
        reason = 'unknown';
    end
end