function [risk_scores_all, is_high_risk, times_km, events_km] = metrics_stats_predictive(varargin)
% METRICS_STATS_PREDICTIVE — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 4b/5 of the metrics step: Predictive Modeling (Elastic Net & Cox prep).
% Performs multivariate feature selection via Elastic Net logistic regression
% and generates out-of-fold risk scores via nested Leave-One-Out Cross-Validation.
%
% INTERFACE:
%   Two calling conventions are supported:
%
%   (1) STRUCT-BASED (recommended):
%       results = metrics_stats_predictive(data_struct)
%       where data_struct is a struct with named fields (see below).
%
%   (2) LEGACY POSITIONAL (deprecated, backward-compatible):
%       [risk_scores_all, is_high_risk, times_km, events_km] = ...
%           metrics_stats_predictive(valid_pts, lf_group, dtype_label, ...
%               output_folder, dataloc, nTp, m_gtv_vol, adc_sd, ...
%               ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, ...
%               f_delta, Dstar_pct, m_d95_gtvp, m_v50gy_gtvp, ...
%               d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, ...
%               d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub, ...
%               id_list, dtype, dl_provenance, x_labels, m_lf, ...
%               m_total_time, m_total_follow_up_time, config_struct)
%       A deprecation warning is issued when this form is used.
%
% STRUCT FIELDS (for struct-based interface):
%   valid_pts, lf_group, dtype_label, output_folder, dataloc, nTp,
%   m_gtv_vol, adc_sd, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct,
%   D_pct, f_delta, Dstar_pct, m_d95_gtvp, m_v50gy_gtvp,
%   d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub,
%   v50_f_sub, d95_dstar_sub, v50_dstar_sub, id_list, dtype,
%   dl_provenance, x_labels, m_lf, m_total_time,
%   m_total_follow_up_time, config_struct (optional)
%
% ANALYTICAL OVERVIEW:
%   This module builds a multivariate predictive model to identify patients
%   at high risk of local failure using combined imaging and dosimetric
%   features.  The analysis pipeline is:
%
%   1. FEATURE ASSEMBLY — Combines 22 candidate features per timepoint:
%      - 4 baseline covariates (Fx1 absolute ADC, D, f, D*)
%      - 4 absolute values at the target fraction
%      - 4 change metrics (percent delta or absolute delta)
%      - 2 whole-GTV dose metrics (D95, V50)
%      - 8 sub-volume dose metrics (D95 and V50 for each of 4 diffusion-
%        defined resistant sub-volumes)
%      Baseline covariates are always included to adjust for pre-existing
%      differences in tumour characteristics.
%
%   2. ELASTIC NET FEATURE SELECTION (alpha=0.5) — Elastic net combines
%      L1 (lasso, promotes sparsity) and L2 (ridge, handles collinearity)
%      penalties.  Alpha=0.5 is a balanced mixture appropriate for correlated
%      imaging features.  The optimal regularisation lambda is selected via
%      5-fold grouped CV (patient-stratified to prevent leakage).
%
%   3. NESTED LOOCV FOR UNBIASED RISK SCORES — Each patient is held out
%      in turn; a complete inner 5-fold CV selects lambda, fits elastic net,
%      and predicts the held-out patient's risk score.  This nested design
%      prevents optimistic bias from using the same data for both feature
%      selection and performance estimation.
%
%   4. ROC ANALYSIS — Out-of-fold risk scores are evaluated via ROC curve
%      and AUC to quantify discriminative ability.  The Youden index
%      identifies the optimal operating point for clinical decision-making.
%
%   5. SANITY CHECKS — Volume change, ADC heterogeneity (SD), and
%      signal-vs-noise floor scatter plots verify that selected biomarkers
%      are not confounded by tumour size changes or measurement noise.
%
%   Data leakage prevention is enforced at multiple levels:
%     - Patient-stratified CV folds (no intra-patient leakage)
%     - KNN imputation with same-patient distance blocking
%     - DL provenance checks (dnCNN/IVIMnet training set exclusion)
%     - Collinearity filtering computed per-fold on training data only
%
% Outputs:
%   risk_scores_all   - LOOCV out-of-fold elastic-net computed risk scores.
%                        If no timepoint yields significant predictive features,
%                        this is returned as a NaN-filled column vector of size
%                        [sum(valid_pts), 1] and a warning is emitted.
%   is_high_risk      - Binarised array demarcating patients > median training risk.
%                        If no timepoint is predictive, returned as NaN-filled
%                        column vector of size [sum(valid_pts), 1].
%   times_km          - Times used for subsequent Kaplan-Meier models.
%                        If no timepoint is predictive, returned as NaN-filled
%                        column vector of size [sum(valid_pts), 1].
%   events_km         - Events matched to times_km.
%                        If no timepoint is predictive, returned as NaN-filled
%                        column vector of size [sum(valid_pts), 1].
%
% NOTE ON NO-SIGNIFICANT-FEATURES CASE:
%   When no timepoint yields significant predictive features (plausible with
%   small cohorts or non-predictive biomarkers), the function returns NaN-filled
%   arrays for all outputs and emits a warning. Downstream code (e.g.,
%   metrics_survival.m) should check for all-NaN risk_scores_all before
%   attempting Kaplan-Meier stratification or Cox modeling.

%% ====================================================================
%  ARGUMENT PARSING: struct-based vs legacy positional
%  ====================================================================
if nargin == 1 && isstruct(varargin{1})
    % ---- Struct-based interface (recommended) ----
    S = varargin{1};
    required_fields = {'valid_pts', 'lf_group', 'dtype_label', 'output_folder', ...
        'dataloc', 'nTp', 'm_gtv_vol', 'adc_sd', 'ADC_abs', 'D_abs', ...
        'f_abs', 'Dstar_abs', 'ADC_pct', 'D_pct', 'f_delta', 'Dstar_pct', ...
        'm_d95_gtvp', 'm_v50gy_gtvp', 'd95_adc_sub', 'v50_adc_sub', ...
        'd95_d_sub', 'v50_d_sub', 'd95_f_sub', 'v50_f_sub', ...
        'd95_dstar_sub', 'v50_dstar_sub', 'id_list', 'dtype', ...
        'dl_provenance', 'x_labels', 'm_lf', 'm_total_time', ...
        'm_total_follow_up_time'};
    missing_fields = setdiff(required_fields, fieldnames(S));
    if ~isempty(missing_fields)
        error('metrics_stats_predictive:missingFields', ...
            'The following required fields are missing from data_struct: %s', ...
            strjoin(missing_fields, ', '));
    end
    valid_pts              = S.valid_pts;
    lf_group               = S.lf_group;
    dtype_label            = S.dtype_label;
    output_folder          = S.output_folder;
    dataloc                = S.dataloc;
    nTp                    = S.nTp;
    m_gtv_vol              = S.m_gtv_vol;
    adc_sd                 = S.adc_sd;
    ADC_abs                = S.ADC_abs;
    D_abs                  = S.D_abs;
    f_abs                  = S.f_abs;
    Dstar_abs              = S.Dstar_abs;
    ADC_pct                = S.ADC_pct;
    D_pct                  = S.D_pct;
    f_delta                = S.f_delta;
    Dstar_pct              = S.Dstar_pct;
    m_d95_gtvp             = S.m_d95_gtvp;
    m_v50gy_gtvp           = S.m_v50gy_gtvp;
    d95_adc_sub            = S.d95_adc_sub;
    v50_adc_sub            = S.v50_adc_sub;
    d95_d_sub              = S.d95_d_sub;
    v50_d_sub              = S.v50_d_sub;
    d95_f_sub              = S.d95_f_sub;
    v50_f_sub              = S.v50_f_sub;
    d95_dstar_sub          = S.d95_dstar_sub;
    v50_dstar_sub          = S.v50_dstar_sub;
    id_list                = S.id_list;
    dtype                  = S.dtype;
    dl_provenance          = S.dl_provenance;
    x_labels               = S.x_labels;
    m_lf                   = S.m_lf;
    m_total_time           = S.m_total_time;
    m_total_follow_up_time = S.m_total_follow_up_time;
    if isfield(S, 'config_struct') && ~isempty(S.config_struct)
        config_struct = S.config_struct;
    else
        config_struct = struct();
    end

    % ---- Dimension consistency checks for struct-based interface ----
    % Verify that all patient-indexed arrays have consistent first dimension
    % (equal to numel(valid_pts)) and all timepoint-indexed arrays have
    % second dimension equal to nTp. This catches mismatched dimensions at
    % the interface boundary rather than deep inside assemble_predictive_features.
    n_pts_expected = numel(valid_pts);
    nTp_expected   = nTp;

    % Define arrays that must be [n_pts_expected x nTp_expected]
    patient_tp_arrays = {ADC_abs, D_abs, f_abs, Dstar_abs, ...
                         ADC_pct, D_pct, f_delta, Dstar_pct, ...
                         m_d95_gtvp, m_v50gy_gtvp, m_gtv_vol, adc_sd};
    patient_tp_names  = {'ADC_abs', 'D_abs', 'f_abs', 'Dstar_abs', ...
                         'ADC_pct', 'D_pct', 'f_delta', 'Dstar_pct', ...
                         'm_d95_gtvp', 'm_v50gy_gtvp', 'm_gtv_vol', 'adc_sd'};

    % Define sub-volume arrays that must also be [n_pts_expected x nTp_expected]
    subvol_arrays = {d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, ...
                     d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub};
    subvol_names  = {'d95_adc_sub', 'v50_adc_sub', 'd95_d_sub', 'v50_d_sub', ...
                     'd95_f_sub', 'v50_f_sub', 'd95_dstar_sub', 'v50_dstar_sub'};

    all_2d_arrays = [patient_tp_arrays, subvol_arrays];
    all_2d_names  = [patient_tp_names, subvol_names];

    % Collect mismatches for a single comprehensive error message
    dim_errors = {};
    for di = 1:numel(all_2d_arrays)
        arr_i = all_2d_arrays{di};
        if isempty(arr_i)
            continue;  % allow legitimately empty arrays
        end
        [nr_i, nc_i] = size(arr_i);
        if nr_i ~= n_pts_expected
            dim_errors{end+1} = sprintf('  %s: %d rows (expected %d = numel(valid_pts))', ...
                all_2d_names{di}, nr_i, n_pts_expected); %#ok<AGROW>
        end
        if nc_i ~= nTp_expected
            dim_errors{end+1} = sprintf('  %s: %d columns (expected %d = nTp)', ...
                all_2d_names{di}, nc_i, nTp_expected); %#ok<AGROW>
        end
    end

    % Check 1-D patient-indexed arrays
    oneD_arrays = {m_lf, m_total_time, m_total_follow_up_time};
    oneD_names  = {'m_lf', 'm_total_time', 'm_total_follow_up_time'};
    for di = 1:numel(oneD_arrays)
        if numel(oneD_arrays{di}) ~= n_pts_expected
            dim_errors{end+1} = sprintf('  %s: %d elements (expected %d = numel(valid_pts))', ...
                oneD_names{di}, numel(oneD_arrays{di}), n_pts_expected); %#ok<AGROW>
        end
    end

    % Check lf_group: should be sum(valid_pts) or numel(valid_pts)
    n_valid_expected = sum(valid_pts);
    if numel(lf_group) ~= n_valid_expected && numel(lf_group) ~= n_pts_expected
        dim_errors{end+1} = sprintf('  lf_group: %d elements (expected %d = sum(valid_pts) or %d = numel(valid_pts))', ...
            numel(lf_group), n_valid_expected, n_pts_expected);
    end

    % Check x_labels
    if numel(x_labels) < nTp_expected
        dim_errors{end+1} = sprintf('  x_labels: %d elements (expected at least %d = nTp)', ...
            numel(x_labels), nTp_expected);
    end

    if ~isempty(dim_errors)
        error('metrics_stats_predictive:dimensionMismatch', ...
            ['Dimension consistency check failed for struct-based interface. ', ...
             'The following arrays have mismatched sizes:\n%s\n', ...
             'All patient-indexed arrays must have numel(valid_pts)=%d rows ', ...
             'and all timepoint-indexed arrays must have nTp=%d columns.'], ...
            strjoin(dim_errors, '\n'), n_pts_expected, nTp_expected);
    end

elseif nargin >= 33
    % ---- Legacy 34-positional-argument interface (deprecated) ----
    warning('metrics_stats_predictive:deprecatedPositionalArgs', ...
        ['The 34-positional-argument calling convention is deprecated and will be ', ...
         'removed in a future release. Please switch to the struct-based interface: ', ...
         'metrics_stats_predictive(data_struct). See help for field names.']);
    valid_pts              = varargin{1};
    lf_group               = varargin{2};
    dtype_label            = varargin{3};
    output_folder          = varargin{4};
    dataloc                = varargin{5};
    nTp                    = varargin{6};
    m_gtv_vol              = varargin{7};
    adc_sd                 = varargin{8};
    ADC_abs                = varargin{9};
    D_abs                  = varargin{10};
    f_abs                  = varargin{11};
    Dstar_abs              = varargin{12};
    ADC_pct                = varargin{13};
    D_pct                  = varargin{14};
    f_delta                = varargin{15};
    Dstar_pct              = varargin{16};
    m_d95_gtvp             = varargin{17};
    m_v50gy_gtvp           = varargin{18};
    d95_adc_sub            = varargin{19};
    v50_adc_sub            = varargin{20};
    d95_d_sub              = varargin{21};
    v50_d_sub              = varargin{22};
    d95_f_sub              = varargin{23};
    v50_f_sub              = varargin{24};
    d95_dstar_sub          = varargin{25};
    v50_dstar_sub          = varargin{26};
    id_list                = varargin{27};
    dtype                  = varargin{28};
    dl_provenance          = varargin{29};
    x_labels               = varargin{30};
    m_lf                   = varargin{31};
    m_total_time           = varargin{32};
    m_total_follow_up_time = varargin{33};
    if nargin >= 34 && ~isempty(varargin{34})
        config_struct = varargin{34};
    else
        config_struct = struct();
    end
else
    error('metrics_stats_predictive:invalidArgs', ...
        ['Invalid calling convention. Use either:\n', ...
         '  (1) metrics_stats_predictive(data_struct)  [struct with named fields]\n', ...
         '  (2) metrics_stats_predictive(valid_pts, lf_group, ..., config_struct)  [34 positional args, deprecated]\n', ...
         'Received %d arguments.'], nargin);
end

%% ====================================================================
%  CONFIG DEFAULTS
%  ====================================================================
if ~isfield(config_struct, 'use_firth_refit')
    config_struct.use_firth_refit = true;
end
use_firth = config_struct.use_firth_refit;

%% ====================================================================
%  DIMENSION VALIDATION
%  ====================================================================
% Determine the number of patients (n_patients) from valid_pts and nTp
% from the scalar. Validate that all [n_patients × nTp] arrays are
% consistently sized. This catches silent errors from swapped arguments
% (e.g., d95_adc_sub vs v50_adc_sub) when they differ in row/column count.
n_patients_total = numel(valid_pts);

% Helper: validate that an array has size [n_patients_total, nTp]
    function validate_array_dims(arr, arr_name)
        if isempty(arr)
            return;  % allow empty arrays (some may be legitimately empty)
        end
        [nr, nc] = size(arr);
        if nr ~= n_patients_total
            error('metrics_stats_predictive:dimensionMismatch', ...
                'Array ''%s'' has %d rows but expected %d (numel(valid_pts)).', ...
                arr_name, nr, n_patients_total);
        end
        if nc ~= nTp
            error('metrics_stats_predictive:dimensionMismatch', ...
                'Array ''%s'' has %d columns but expected %d (nTp).', ...
                arr_name, nc, nTp);
        end
    end

% Validate all [n_patients × nTp] arrays
arrays_to_check = {ADC_abs, D_abs, f_abs, Dstar_abs, ...
                   ADC_pct, D_pct, f_delta, Dstar_pct, ...
                   m_d95_gtvp, m_v50gy_gtvp, ...
                   d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, ...
                   d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub, ...
                   m_gtv_vol, adc_sd};
array_names = {'ADC_abs', 'D_abs', 'f_abs', 'Dstar_abs', ...
               'ADC_pct', 'D_pct', 'f_delta', 'Dstar_pct', ...
               'm_d95_gtvp', 'm_v50gy_gtvp', ...
               'd95_adc_sub', 'v50_adc_sub', 'd95_d_sub', 'v50_d_sub', ...
               'd95_f_sub', 'v50_f_sub', 'd95_dstar_sub', 'v50_dstar_sub', ...
               'm_gtv_vol', 'adc_sd'};
for vi = 1:numel(arrays_to_check)
    validate_array_dims(arrays_to_check{vi}, array_names{vi});
end

% Validate 1-D arrays that should have n_patients_total elements
assert(numel(lf_group) == sum(valid_pts) || numel(lf_group) == n_patients_total, ...
    'metrics_stats_predictive:dimensionMismatch', ...
    'lf_group has %d elements; expected %d (sum(valid_pts)) or %d (numel(valid_pts)).', ...
    numel(lf_group), sum(valid_pts), n_patients_total);
assert(numel(m_lf) == n_patients_total, ...
    'metrics_stats_predictive:dimensionMismatch', ...
    'm_lf has %d elements; expected %d.', numel(m_lf), n_patients_total);
assert(numel(m_total_time) == n_patients_total, ...
    'metrics_stats_predictive:dimensionMismatch', ...
    'm_total_time has %d elements; expected %d.', numel(m_total_time), n_patients_total);
assert(numel(m_total_follow_up_time) == n_patients_total, ...
    'metrics_stats_predictive:dimensionMismatch', ...
    'm_total_follow_up_time has %d elements; expected %d.', numel(m_total_follow_up_time), n_patients_total);
assert(numel(x_labels) >= nTp, ...
    'metrics_stats_predictive:dimensionMismatch', ...
    'x_labels has %d elements; expected at least %d (nTp).', numel(x_labels), nTp);

%% ====================================================================
%  MAIN ANALYSIS
%  ====================================================================
fprintf('  --- SECTION 10: Per-Timepoint Analysis Loop ---\n');

% Diary: capture console output to output_folder
diary_file = fullfile(output_folder, ['metrics_stats_predictive_output_' dtype_label '.txt']);
if exist(diary_file, 'file'), delete(diary_file); end
diary(diary_file);
cleanupDiary = onCleanup(@() diary('off'));

% Initialize risk outputs to be returned and used by metrics_survival.m
% for Kaplan-Meier stratification and Cox proportional hazards modeling.
n_valid = sum(valid_pts);
risk_scores_all = [];  % continuous LOOCV out-of-fold predicted risk scores
is_high_risk = [];     % binary high/low risk classification (threshold = median training risk)
times_km = [];         % time-to-event (days from RT end) for KM plotting
events_km = [];        % event indicator (0=censored, 1=LF, 2=competing risk)
% Track the earliest timepoint with significant features.  Earlier timepoints
% are clinically more actionable — if treatment resistance can be detected
% at Fx2 (after 1 week) rather than Fx5 (after 5 weeks), there is more
% time to adapt the treatment plan (e.g., dose escalation, change in
% chemotherapy regimen, or surgical intervention).
best_risk_fx = Inf;

% --- Pre-index all parameter arrays with valid_pts ONCE before the loop ---
% This avoids redundant re-indexing of ~30 large arrays at each timepoint
% iteration (4 iterations × ~30 arrays = ~120 redundant indexing operations).
ADC_abs_valid      = ADC_abs(valid_pts, :);
D_abs_valid        = D_abs(valid_pts, :);
f_abs_valid        = f_abs(valid_pts, :);
Dstar_abs_valid    = Dstar_abs(valid_pts, :);
ADC_pct_valid      = ADC_pct(valid_pts, :);
D_pct_valid        = D_pct(valid_pts, :);
f_delta_valid      = f_delta(valid_pts, :);
Dstar_pct_valid    = Dstar_pct(valid_pts, :);
m_d95_gtvp_valid   = m_d95_gtvp(valid_pts, :);
m_v50gy_gtvp_valid = m_v50gy_gtvp(valid_pts, :);
d95_adc_sub_valid  = d95_adc_sub(valid_pts, :);
v50_adc_sub_valid  = v50_adc_sub(valid_pts, :);
d95_d_sub_valid    = d95_d_sub(valid_pts, :);
v50_d_sub_valid    = v50_d_sub(valid_pts, :);
d95_f_sub_valid    = d95_f_sub(valid_pts, :);
v50_f_sub_valid    = v50_f_sub(valid_pts, :);
d95_dstar_sub_valid = d95_dstar_sub(valid_pts, :);
v50_dstar_sub_valid = v50_dstar_sub(valid_pts, :);
m_gtv_vol_valid    = m_gtv_vol(valid_pts, :);
adc_sd_valid       = adc_sd(valid_pts, :);

% --- Validate that all *_valid arrays have at least nTp columns ---
% After pre-indexing with valid_pts, the row count is guaranteed to be
% n_valid = sum(valid_pts).  However, if any source array had fewer columns
% than nTp (e.g., sub-volume dosimetry computed for a subset of timepoints,
% or column-indexed differently), the downstream indexing inside
% assemble_predictive_features would fail with a cryptic dimension mismatch.
% Catch this early with descriptive error messages.
valid_arrays_to_check = {ADC_abs_valid, D_abs_valid, f_abs_valid, Dstar_abs_valid, ...
                         ADC_pct_valid, D_pct_valid, f_delta_valid, Dstar_pct_valid, ...
                         m_d95_gtvp_valid, m_v50gy_gtvp_valid, ...
                         d95_adc_sub_valid, v50_adc_sub_valid, ...
                         d95_d_sub_valid, v50_d_sub_valid, ...
                         d95_f_sub_valid, v50_f_sub_valid, ...
                         d95_dstar_sub_valid, v50_dstar_sub_valid, ...
                         m_gtv_vol_valid, adc_sd_valid};
valid_array_names = {'ADC_abs_valid', 'D_abs_valid', 'f_abs_valid', 'Dstar_abs_valid', ...
                     'ADC_pct_valid', 'D_pct_valid', 'f_delta_valid', 'Dstar_pct_valid', ...
                     'm_d95_gtvp_valid', 'm_v50gy_gtvp_valid', ...
                     'd95_adc_sub_valid', 'v50_adc_sub_valid', ...
                     'd95_d_sub_valid', 'v50_d_sub_valid', ...
                     'd95_f_sub_valid', 'v50_f_sub_valid', ...
                     'd95_dstar_sub_valid', 'v50_dstar_sub_valid', ...
                     'm_gtv_vol_valid', 'adc_sd_valid'};
for vi = 1:numel(valid_arrays_to_check)
    arr_vi = valid_arrays_to_check{vi};
    if ~isempty(arr_vi)
        assert(size(arr_vi, 2) >= nTp, ...
            'metrics_stats_predictive:validArrayColumnMismatch', ...
            ['Array ''%s'' has %d columns after valid_pts indexing, but nTp = %d. ', ...
             'All pre-indexed arrays must have at least nTp columns to support ', ...
             'timepoint indexing in assemble_predictive_features.'], ...
            valid_array_names{vi}, size(arr_vi, 2), nTp);
    end
end

% Iterate from Fx2 onwards (Fx1 is baseline — no change to analyse).
% Each timepoint is analysed independently to identify the earliest
% fraction at which treatment response prediction becomes feasible.
for target_fx = 2:nTp
    fx_label = x_labels{target_fx};
    fprintf('\n=== Analyzing %s ===\n', fx_label);

    %% --- Feature Assembly ---
    % Pass pre-indexed (valid_pts subset) arrays to avoid redundant indexing.
    % assemble_predictive_features receives data already subset to valid
    % patients, so it should use a trivial index (true mask or 1:size) internally.
    trivial_mask = true(size(ADC_abs_valid, 1), 1);
    [X_lasso_all, feat_names_lasso, original_feature_indices, feat_names_lasso_full] = assemble_predictive_features( ...
        trivial_mask, target_fx, nTp, fx_label, output_folder, ...
        ADC_abs_valid, D_abs_valid, f_abs_valid, Dstar_abs_valid, ...
        ADC_pct_valid, D_pct_valid, f_delta_valid, Dstar_pct_valid, ...
        m_d95_gtvp_valid, m_v50gy_gtvp_valid, ...
        d95_adc_sub_valid, v50_adc_sub_valid, d95_d_sub_valid, v50_d_sub_valid, ...
        d95_f_sub_valid, v50_f_sub_valid, d95_dstar_sub_valid, v50_dstar_sub_valid);

    y_lasso_all = lf_group;
    % Exclude competing risk patients (lf==2) from the binomial model.
    % Previously, competing events were recoded as LC (lf==0), which
    % misclassifies patients who died of non-cancer causes before LF could
    % be observed.  This biased the elastic net by diluting the LC group.
    competing_mask = (y_lasso_all == 2);
    y_lasso_all(competing_mask) = NaN;  % mark for exclusion below

    % Filter to patients with at least SOME imaging data in the first 8
    % columns (baseline + absolute parameters).
    base_cols = min(8, size(X_lasso_all, 2));
    has_any_imaging = any(~isnan(X_lasso_all(:, 1:base_cols)), 2);
    impute_mask = has_any_imaging & ~isnan(y_lasso_all);
    X_impute = X_lasso_all(impute_mask, :);
    y_clean  = y_lasso_all(impute_mask);

    assert(all(y_clean == 0 | y_clean == 1), ...
        'metrics_stats_predictive:unexpectedOutcome', ...
        'Unexpected outcome values in lf_group after competing risk exclusion.');

    id_list_valid = id_list(valid_pts);
    id_list_impute = id_list_valid(impute_mask);

    if isempty(X_impute)
        fprintf('  No patients with any imaging data at %s. Skipping predictive modeling.\n', fx_label);
        continue;
    end

    fprintf('  ⚠️  Single imputation: CIs may be anti-conservative (imputation uncertainty not propagated).\n');

    % --- DL provenance leakage check ---
    if isfield(dl_provenance, 'manifest_loaded') && ~dl_provenance.manifest_loaded && (dtype == 2 || dtype == 3)
        fprintf('  ⚠️  DL provenance manifest not loaded for %s — leakage guard inactive. Skipping predictive modeling.\n', dtype_label);
        continue;
    end
    if dtype == 2
        for chk_i = 1:numel(id_list_impute)
            if any(strcmp(dl_provenance.dncnn_train_ids, id_list_impute{chk_i}))
                error('metrics_stats_predictive:dataLeakage', 'DATA LEAKAGE DETECTED: Patient %s was used to train the DnCNN model.', id_list_impute{chk_i});
            end
        end
    elseif dtype == 3
        for chk_i = 1:numel(id_list_impute)
            if any(strcmp(dl_provenance.ivimnet_train_ids, id_list_impute{chk_i}))
                error('metrics_stats_predictive:dataLeakage', 'DATA LEAKAGE DETECTED: Patient %s was used to train the IVIMnet model.', id_list_impute{chk_i});
            end
        end
    end

    %% --- Elastic Net Feature Selection ---
    [selected_indices, opt_lambda, common_Lambda, cv_failed, keep_fold_counts, coefs_en, final_feature_indices] = run_elastic_net_cv( ...
        X_impute, y_clean, id_list_impute, 5, use_firth, ...
        original_feature_indices, feat_names_lasso_full, fx_label);

    %% --- LOOCV For Risk Scores ---
    [risk_scores_oof, is_high_risk_oof] = run_loocv_risk_scores( ...
        X_impute, y_clean, id_list_impute, dl_provenance, dtype, dtype_label, use_firth);

    risk_scores_all_target = nan(sum(valid_pts), 1);
    risk_scores_all_target(impute_mask) = risk_scores_oof;

    is_high_risk_target = nan(sum(valid_pts), 1);
    is_high_risk_target(impute_mask) = is_high_risk_oof;

    %% --- Build feature metadata for diagnostics ---
    all_feat_data  = {ADC_abs_valid,       D_abs_valid,       f_abs_valid,       Dstar_abs_valid, ...
                      ADC_abs_valid,       D_abs_valid,       f_abs_valid,       Dstar_abs_valid, ...
                      ADC_pct_valid,       D_pct_valid,       f_delta_valid,       Dstar_pct_valid, ...
                      m_d95_gtvp_valid,    m_v50gy_gtvp_valid, ...
                      d95_adc_sub_valid,   v50_adc_sub_valid, ...
                      d95_d_sub_valid,     v50_d_sub_valid, ...
                      d95_f_sub_valid,     v50_f_sub_valid, ...
                      d95_dstar_sub_valid, v50_dstar_sub_valid};

    all_feat_names = {'ADC BL',      'D BL',      'f BL',      'D* BL', ...
                      'ADC',         'D',         'f',         'D*', ...
                      'ADC',         'D',         'f',         'D*', ...
                      'D95 GTVp',    'V50 GTVp', ...
                      'D95 Sub(ADC)','V50 Sub(ADC)', ...
                      'D95 Sub(D)',  'V50 Sub(D)', ...
                      'D95 Sub(f)',  'V50 Sub(f)', ...
                      'D95 Sub(D*)', 'V50 Sub(D*)'};

    all_feat_is_abs = [true          true         true         true  ...
                       true          true         true         true  ...
                       false         false        false        false ...
                       true          false        true         false ...
                       true          false        true         false ...
                       true          false];

    all_feat_disp  = {'BL ADC',      'BL D',      'BL f',      'BL D*', ...
                      'Abs ADC',     'Abs D',     'Abs f',     'Abs D*', ...
                      '\Delta ADC',  '\Delta D',  '\Delta f',  '\Delta D*', ...
                      'D95 GTVp',    'V50 GTVp', ...
                      'D95 Sub(ADC)','V50 Sub(ADC)', ...
                      'D95 Sub(D)',  'V50 Sub(D)', ...
                      'D95 Sub(f)',  'V50 Sub(f)', ...
                      'D95 Sub(D*)', 'V50 Sub(D*)'};

    all_feat_units = {'mm^2/s',      'mm^2/s',    'Fraction',  'mm^2/s', ...
                      'mm^2/s',      'mm^2/s',    'Fraction',  'mm^2/s', ...
                      '%',           '%',         'Fraction',  '%', ...
                      'Gy',          '%', ...
                      'Gy',          '%', ...
                      'Gy',          '%', ...
                      'Gy',          '%'};

    n_sig = length(selected_indices);

    sig_data_selected = cell(1, n_sig);
    sig_abs_data      = cell(1, n_sig);
    sig_pct_data      = cell(1, n_sig);
    sig_names         = cell(1, n_sig);
    sig_is_abs        = false(1, n_sig);
    sig_is_pct_imaging = false(1, n_sig);
    sig_disp_names    = cell(1, n_sig);
    sig_units         = cell(1, n_sig);
    sig_col_idx       = zeros(1, n_sig);

    for si = 1:n_sig
        fi = selected_indices(si);
        sig_data_selected{si} = all_feat_data{fi};
        sig_names{si}         = all_feat_names{fi};
        sig_is_abs(si)        = all_feat_is_abs(fi);
        sig_is_pct_imaging(si) = (fi >= 9 && fi <= 12);
        sig_disp_names{si}    = all_feat_disp{fi};
        sig_units{si}         = all_feat_units{fi};
        if fi <= 4
            sig_col_idx(si) = 1;
        else
            sig_col_idx(si) = target_fx;
        end

        if fi <= 4
            sig_abs_data{si} = all_feat_data{fi};
            sig_pct_data{si} = all_feat_data{min(fi + 8, numel(all_feat_data))};
        elseif fi >= 5 && fi <= 8
            sig_abs_data{si} = all_feat_data{fi};
            sig_pct_data{si} = all_feat_data{fi + 4};
        elseif fi >= 9 && fi <= 12
            sig_abs_data{si} = all_feat_data{fi - 4};
            sig_pct_data{si} = all_feat_data{fi};
        else
            sig_abs_data{si} = all_feat_data{fi};
            sig_pct_data{si} = all_feat_data{fi};
        end
    end

    fprintf('Significant variables at %s: ', fx_label);
    if n_sig == 0
        fprintf('NONE. Skipping downstream analyses for %s.\n', fx_label);
        continue;
    else
        fprintf('%s\n', strjoin(sig_disp_names, ', '));
    end

    % Construct time-to-event arrays for Kaplan-Meier survival analysis.
    times_km = m_total_time;
    cens_or_cr = (m_lf == 0 | m_lf == 2) & ~isnan(m_total_follow_up_time);
    times_km(cens_or_cr) = m_total_follow_up_time(cens_or_cr);
    events_km = m_lf;

    times_km = times_km(valid_pts);
    events_km = events_km(valid_pts);

    %% --- Calibration Assessment ---
    cal_metrics = compute_calibration_metrics( ...
        risk_scores_oof, y_clean, 5, output_folder, dtype_label, fx_label);

    %% --- Diagnostic Plots (ROC, Sanity Checks, 2D Scatter) ---
    plot_predictive_diagnostics( ...
        selected_indices, n_sig, sig_data_selected, sig_names, sig_is_abs, ...
        sig_is_pct_imaging, sig_disp_names, sig_units, sig_col_idx, ...
        sig_abs_data, sig_pct_data, ...
        risk_scores_all_target, lf_group, valid_pts, ...
        m_gtv_vol_valid, adc_sd_valid, ADC_abs_valid, ...
        target_fx, fx_label, dtype_label, dtype, output_folder, use_firth);

    %% --- Decision Curve Analysis ---
    try
        dca_results = decision_curve_analysis( ...
            y_clean, risk_scores_oof, [], output_folder, dtype_label, fx_label);
    catch ME_dca
        fprintf('  ⚠️  Decision curve analysis failed: %s\n', ME_dca.message);
    end

    %% --- Net Reclassification Improvement (NRI) ---
    try
        baseline_cols = [];
        for bc = 1:numel(feat_names_lasso)
            if any(strcmp(feat_names_lasso{bc}, {'D95_GTVp', 'V50_GTVp'}))
                baseline_cols(end+1) = bc; %#ok<AGROW>
            end
        end
        if numel(baseline_cols) >= 1 && numel(y_clean) >= 10
            X_baseline = X_impute(:, baseline_cols);
            prob_baseline = nan(numel(y_clean), 1);
            for li = 1:numel(y_clean)
                train_mask = true(numel(y_clean), 1);
                train_mask(li) = false;
                X_tr = X_baseline(train_mask, :);
                y_tr = y_clean(train_mask);
                X_te = X_baseline(li, :);
                try
                    b_bl = glmfit(X_tr, y_tr, 'binomial', 'Link', 'logit');
                    prob_baseline(li) = glmval(b_bl, X_te, 'logit');
                catch
                    prob_baseline(li) = mean(y_tr);
                end
            end

            prob_full = 1 ./ (1 + exp(-risk_scores_oof));
            prob_full(isnan(prob_full)) = 0.5;

            fprintf('\n  --- Net Reclassification Improvement (Full vs Baseline) ---\n');
            nri_results = compute_nri(y_clean, prob_baseline, prob_full);
        end
    catch ME_nri
        fprintf('  ⚠️  NRI computation failed: %s\n', ME_nri.message);
    end

    %% --- Imputation Sensitivity Analysis (optional, expensive) ---
    if isfield(config_struct, 'run_imputation_sensitivity') && config_struct.run_imputation_sensitivity
        try
            fprintf('\n');
            imp_sens = imputation_sensitivity(X_impute, id_list_impute, ...
                feat_names_lasso, y_clean, id_list_impute, dl_provenance, ...
                dtype, dtype_label, use_firth, output_folder);
        catch ME_imp
            fprintf('  ⚠️  Imputation sensitivity analysis failed: %s\n', ME_imp.message);
        end
    end

    % Keep the earliest timepoint with significant features.
    if target_fx < best_risk_fx
        risk_scores_all = risk_scores_all_target;
        is_high_risk = is_high_risk_target;
        best_risk_fx = target_fx;
        fprintf('  Retaining risk scores from %s (earliest significant timepoint so far).\n', fx_label);
    end
end

%% --- External Validation Model Export (optional) ---
if isfield(config_struct, 'export_validation_model') && config_struct.export_validation_model
    if ~isempty(risk_scores_all) && exist('coefs_en', 'var') && exist('feat_names_lasso', 'var')
        try
            trained_model = struct();
            trained_model.coefficients = coefs_en;
            trained_model.selected_features = selected_indices;
            trained_model.feature_names = feat_names_lasso;
            trained_model.scaling_mu = mean(X_impute, 1, 'omitnan');
            trained_model.scaling_sigma = std(X_impute, 0, 1, 'omitnan');
            trained_model.imputation_ref = X_impute;
            valid_rs = ~isnan(risk_scores_oof) & ~isnan(y_clean);
            if sum(valid_rs) >= 5
                [~, ~, thresholds_roc, ~] = perfcurve(y_clean(valid_rs), risk_scores_oof(valid_rs), 1);
                trained_model.risk_threshold = median(thresholds_roc);
            else
                trained_model.risk_threshold = 0.5;
            end
            if sum(valid_rs) >= 5
                [~, ~, ~, auc_export] = perfcurve(y_clean(valid_rs), risk_scores_oof(valid_rs), 1);
            else
                auc_export = NaN;
            end
            trained_model.auc = auc_export;
            trained_model.n_patients = numel(y_clean);
            trained_model.event_rate = mean(y_clean == 1);

            val_path = fullfile(output_folder, sprintf('validation_model_%s.mat', dtype_label));
            prepare_external_validation(trained_model, config_struct, val_path);
        catch ME_exp
            fprintf('  ⚠️  Validation model export failed: %s\n', ME_exp.message);
        end
    end
end

if use_firth
    [~, last_warn_id] = lastwarn;
    if strcmp(last_warn_id, 'stats:lassoGlm:PerfectSeparation')
        lastwarn('');
    end
end

diary off;
end