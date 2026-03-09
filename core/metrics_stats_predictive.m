function [risk_scores_all, is_high_risk, times_km, events_km] = metrics_stats_predictive(valid_pts, lf_group, dtype_label, output_folder, dataloc, nTp, m_gtv_vol, adc_sd, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub, id_list, dtype, dl_provenance, x_labels, m_lf, m_total_time, m_total_follow_up_time, config_struct)
% METRICS_STATS_PREDICTIVE — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 4b/5 of the metrics step: Predictive Modeling (Elastic Net & Cox prep).
% Performs multivariate feature selection via Elastic Net logistic regression
% and generates out-of-fold risk scores via nested Leave-One-Out Cross-Validation.
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
% Inputs:
%   valid_pts         - Logical mask of patients mapped to LF/LC groups
%   lf_group          - Outcome variable arrays (0=LC, 1=LF)
%   dtype_label       - Used in naming figures
%   output_folder     - Where output figures should be saved
%   dataloc           - Directory path containing input data
%   nTp               - Total count of timepoints
%   m_*               - Multiple parameter arrays (GTV volumes, Dose, absolute parameters, etc)
%   d95_*, v50_*      - Sub-volume dose coverage arrays
%   dl_provenance     - Matrix or struct ensuring training leak prevention in cross-val
%   x_labels          - Labels for the time component
%   config_struct     - (Optional, 34th arg) Pipeline config struct.  When
%                       config_struct.use_firth_refit is true, the final model
%                       and each LOOCV fold are refitted using Firth penalized
%                       logistic regression (Jeffreys prior) on the features
%                       selected by elastic net.  This produces finite, bias-
%                       corrected coefficient estimates even under perfect or
%                       quasi-perfect separation.  Default: use_firth_refit=true.
%
% Outputs:
%   risk_scores_all   - LOOCV out-of-fold elastic-net computed risk scores
%   is_high_risk      - Binarised array demarcating patients > median training risk
%   times_km          - Times used for subsequent Kaplan-Meier models
%   events_km         - Events matched to times_km

% Backward-compatible: config_struct is optional (34th arg).
if nargin < 34 || isempty(config_struct)
    config_struct = struct();
end
if ~isfield(config_struct, 'use_firth_refit')
    config_struct.use_firth_refit = true;
end
use_firth = config_struct.use_firth_refit;

fprintf('  --- SECTION 10: Per-Timepoint Analysis Loop ---\n');

% Diary: capture console output to output_folder
diary_file = fullfile(output_folder, ['metrics_stats_predictive_output_' dtype_label '.txt']);
if exist(diary_file, 'file'), delete(diary_file); end
diary(diary_file);

% Initialize risk outputs to be returned and used by survival
risk_scores_all = [];
is_high_risk = [];
times_km = [];
events_km = [];
% Track the earliest timepoint with significant features.  Earlier timepoints
% are clinically more actionable — if treatment resistance can be detected
% at Fx2 (after 1 week) rather than Fx5 (after 5 weeks), there is more
% time to adapt the treatment plan (e.g., dose escalation, change in
% chemotherapy regimen, or surgical intervention).
best_risk_fx = Inf;

% Iterate from Fx2 onwards (Fx1 is baseline — no change to analyse).
% Each timepoint is analysed independently to identify the earliest
% fraction at which treatment response prediction becomes feasible.
for target_fx = 2:nTp
    fx_label = x_labels{target_fx};
    fprintf('\n=== Analyzing %s ===\n', fx_label);

    %% --- Feature Assembly ---
    [X_lasso_all, feat_names_lasso, original_feature_indices, feat_names_lasso_full] = assemble_predictive_features( ...
        valid_pts, target_fx, nTp, fx_label, output_folder, ...
        ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
        m_d95_gtvp, m_v50gy_gtvp, ...
        d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, ...
        d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub);

    y_lasso_all = lf_group;
    % Exclude competing risk patients (lf==2) from the binomial model.
    % Previously, competing events were recoded as LC (lf==0), which
    % misclassifies patients who died of non-cancer causes before LF could
    % be observed.  This biased the elastic net by diluting the LC group.
    % Consistent with the GLME approach in metrics_stats_comparisons.m
    % which also excludes competing risk patients.
    competing_mask = (y_lasso_all == 2);
    y_lasso_all(competing_mask) = NaN;  % mark for exclusion below

    base_cols = min(8, size(X_lasso_all, 2));
    has_any_imaging = any(~isnan(X_lasso_all(:, 1:base_cols)), 2);
    impute_mask = has_any_imaging & ~isnan(y_lasso_all);
    X_impute = X_lasso_all(impute_mask, :);
    y_clean  = y_lasso_all(impute_mask);

    id_list_valid = id_list(valid_pts);
    id_list_impute = id_list_valid(impute_mask);

    if isempty(X_impute)
        fprintf('  No patients with any imaging data at %s. Skipping predictive modeling.\n', fx_label);
        continue;
    end

    % LIMITATION: KNN imputation produces a single completed dataset.
    % Confidence intervals and p-values may be anti-conservative because
    % they do not account for imputation uncertainty.  For definitive
    % inference, consider multiple imputation (e.g., MICE) and pool
    % estimates via Rubin's rules.  Single imputation is retained here
    % because the small cohort size makes stable MI infeasible.
    fprintf('  ⚠️  Single imputation: CIs may be anti-conservative (imputation uncertainty not propagated).\n');

    % --- DL provenance leakage check (covers outer CV, not just LOOCV) ---
    if isfield(dl_provenance, 'manifest_loaded') && ~dl_provenance.manifest_loaded && (dtype == 2 || dtype == 3)
        fprintf('  ⚠️  DL provenance manifest not loaded for %s — leakage guard inactive. Skipping predictive modeling.\n', dtype_label);
        continue;
    end
    if dtype == 2
        for chk_i = 1:numel(id_list_impute)
            if any(strcmp(dl_provenance.dncnn_train_ids, id_list_impute{chk_i}))
                error('DATA LEAKAGE DETECTED: Patient %s was used to train the DnCNN model.', id_list_impute{chk_i});
            end
        end
    elseif dtype == 3
        for chk_i = 1:numel(id_list_impute)
            if any(strcmp(dl_provenance.ivimnet_train_ids, id_list_impute{chk_i}))
                error('DATA LEAKAGE DETECTED: Patient %s was used to train the IVIMnet model.', id_list_impute{chk_i});
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
    all_feat_data  = {ADC_abs,       D_abs,       f_abs,       Dstar_abs, ...   % 1-4: baseline covariates
                      ADC_abs,       D_abs,       f_abs,       Dstar_abs, ...   % 5-8: absolute at target_fx
                      ADC_pct,       D_pct,       f_delta,       Dstar_pct, ... % 9-12: percent change
                      m_d95_gtvp,    m_v50gy_gtvp, ...                          % 13-14: dose
                      d95_adc_sub,   v50_adc_sub, ...                           % 15-16
                      d95_d_sub,     v50_d_sub, ...                             % 17-18
                      d95_f_sub,     v50_f_sub, ...                             % 19-20
                      d95_dstar_sub, v50_dstar_sub};                            % 21-22

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
    sig_col_idx       = zeros(1, n_sig);  % which column to plot for each feature

    for si = 1:n_sig
        fi = selected_indices(si);
        sig_data_selected{si} = all_feat_data{fi};
        sig_names{si}         = all_feat_names{fi};
        sig_is_abs(si)        = all_feat_is_abs(fi);
        sig_is_pct_imaging(si) = (fi >= 9 && fi <= 12);
        sig_disp_names{si}    = all_feat_disp{fi};
        sig_units{si}         = all_feat_units{fi};
        % Baseline features (indices 1-4) used column 1 in the model;
        % all other features used column target_fx.
        if fi <= 4
            sig_col_idx(si) = 1;
        else
            sig_col_idx(si) = target_fx;
        end

        if fi <= 4
            % Baseline covariates: use baseline as abs, pair with pct at fi+8
            sig_abs_data{si} = all_feat_data{fi};
            sig_pct_data{si} = all_feat_data{min(fi + 8, numel(all_feat_data))};
        elseif fi >= 5 && fi <= 8
            % Absolute at target_fx: pair with pct at fi+4
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

    times_km = m_total_time;
    % Censored (lf==0) and competing risk (lf==2) patients use follow-up
    % time, consistent with the CSH approach in metrics_survival.m.
    cens_or_cr = (m_lf == 0 | m_lf == 2) & ~isnan(m_total_follow_up_time);
    times_km(cens_or_cr) = m_total_follow_up_time(cens_or_cr);
    events_km = m_lf;

    times_km = times_km(valid_pts);
    events_km = events_km(valid_pts);

    %% --- Diagnostic Plots (ROC, Sanity Checks, 2D Scatter) ---
    plot_predictive_diagnostics( ...
        selected_indices, n_sig, sig_data_selected, sig_names, sig_is_abs, ...
        sig_is_pct_imaging, sig_disp_names, sig_units, sig_col_idx, ...
        sig_abs_data, sig_pct_data, ...
        risk_scores_all_target, lf_group, valid_pts, ...
        m_gtv_vol, adc_sd, ADC_abs, ...
        target_fx, fx_label, dtype_label, dtype, output_folder, use_firth);

    % Keep the earliest timepoint with significant features so that
    % survival analysis uses the most clinically actionable (early) risk
    % scores rather than silently overwriting with the last timepoint.
    if target_fx < best_risk_fx
        risk_scores_all = risk_scores_all_target;
        is_high_risk = is_high_risk_target;
        best_risk_fx = target_fx;
        fprintf('  Retaining risk scores from %s (earliest significant timepoint so far).\n', fx_label);
    end
end

% Clear any PerfectSeparation warning from lastwarn buffer when Firth
% handled the separation.  The orchestrator (run_dwi_pipeline.m) checks
% lastwarn after each module and logs it to error.log — clearing it here
% prevents logging a warning that Firth has already resolved.
if use_firth
    [~, last_warn_id] = lastwarn;
    if strcmp(last_warn_id, 'stats:lassoGlm:PerfectSeparation')
        lastwarn('');
    end
end

diary off;
end
