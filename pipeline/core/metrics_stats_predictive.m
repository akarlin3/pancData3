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

% Initialize risk outputs to be returned and used by metrics_survival.m
% for Kaplan-Meier stratification and Cox proportional hazards modeling.
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

    % Filter to patients with at least SOME imaging data in the first 8
    % columns (baseline + absolute parameters).  Patients with entirely
    % missing imaging at this timepoint cannot contribute to the model.
    % The impute_mask also requires a valid (non-NaN) outcome label.
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
    % Map each of the 22 model features back to its source data array,
    % display name, units, and whether it is an absolute or change metric.
    % This metadata drives the diagnostic scatter plots and ROC annotations.
    % Feature indices 1-4 are baseline covariates (Fx1), 5-8 are absolute
    % values at target_fx, 9-12 are percent/absolute changes, 13-14 are
    % whole-GTV dose, and 15-22 are sub-volume dose coverage metrics.
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

    % Construct time-to-event arrays for Kaplan-Meier survival analysis.
    % LF patients (m_lf==1) use time-to-failure; censored (m_lf==0) and
    % competing risk (m_lf==2) patients use time-to-last-follow-up.
    % This is the Cause-Specific Hazard (CSH) convention: competing events
    % are treated as censored for the LF endpoint, while still recording
    % the event type for later competing risk analysis.
    times_km = m_total_time;
    cens_or_cr = (m_lf == 0 | m_lf == 2) & ~isnan(m_total_follow_up_time);
    times_km(cens_or_cr) = m_total_follow_up_time(cens_or_cr);
    events_km = m_lf;

    times_km = times_km(valid_pts);
    events_km = events_km(valid_pts);

    %% --- Calibration Assessment ---
    % Compute calibration metrics for the LOOCV risk scores at this timepoint
    cal_metrics = compute_calibration_metrics( ...
        risk_scores_oof, y_clean, 5, output_folder, dtype_label, fx_label);

    %% --- Diagnostic Plots (ROC, Sanity Checks, 2D Scatter) ---
    plot_predictive_diagnostics( ...
        selected_indices, n_sig, sig_data_selected, sig_names, sig_is_abs, ...
        sig_is_pct_imaging, sig_disp_names, sig_units, sig_col_idx, ...
        sig_abs_data, sig_pct_data, ...
        risk_scores_all_target, lf_group, valid_pts, ...
        m_gtv_vol, adc_sd, ADC_abs, ...
        target_fx, fx_label, dtype_label, dtype, output_folder, use_firth);

    %% --- Decision Curve Analysis ---
    try
        dca_results = decision_curve_analysis( ...
            y_clean, risk_scores_oof, [], output_folder, dtype_label, fx_label);
    catch ME_dca
        fprintf('  ⚠️  Decision curve analysis failed: %s\n', ME_dca.message);
    end

    %% --- Net Reclassification Improvement (NRI) ---
    % Compare full model vs baseline model (GTV volume + D95 only)
    try
        % Build baseline model features: GTV volume (col 13) + D95 (col 14) only
        % Map to available columns in imputed data
        baseline_cols = [];
        for bc = 1:numel(feat_names_lasso)
            if any(strcmp(feat_names_lasso{bc}, {'D95_GTVp', 'V50_GTVp'}))
                baseline_cols(end+1) = bc; %#ok<AGROW>
            end
        end
        if numel(baseline_cols) >= 1 && numel(y_clean) >= 10
            % Fit baseline logistic model via LOOCV
            X_baseline = X_impute(:, baseline_cols);
            prob_baseline = nan(numel(y_clean), 1);
            for li = 1:numel(y_clean)
                train_mask = true(numel(y_clean), 1);
                train_mask(li) = false;
                X_tr = X_baseline(train_mask, :);
                y_tr = y_clean(train_mask);
                X_te = X_baseline(li, :);
                % Simple logistic regression for baseline
                try
                    b_bl = glmfit(X_tr, y_tr, 'binomial', 'Link', 'logit');
                    prob_baseline(li) = glmval(b_bl, X_te, 'logit');
                catch
                    prob_baseline(li) = mean(y_tr);
                end
            end

            % Full model probabilities
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

%% --- External Validation Model Export (optional) ---
if isfield(config_struct, 'export_validation_model') && config_struct.export_validation_model
    if ~isempty(risk_scores_all) && exist('coefs_en', 'var') && exist('feat_names_lasso', 'var')
        try
            trained_model = struct();
            trained_model.coefficients = coefs_en;
            trained_model.selected_features = selected_indices;
            trained_model.feature_names = feat_names_lasso;
            % Compute scaling from training data
            trained_model.scaling_mu = mean(X_impute, 1, 'omitnan');
            trained_model.scaling_sigma = std(X_impute, 0, 1, 'omitnan');
            trained_model.imputation_ref = X_impute;
            % Youden threshold from ROC
            valid_rs = ~isnan(risk_scores_oof) & ~isnan(y_clean);
            if sum(valid_rs) >= 5
                [~, ~, thresholds_roc, ~] = perfcurve(y_clean(valid_rs), risk_scores_oof(valid_rs), 1);
                trained_model.risk_threshold = median(thresholds_roc);
            else
                trained_model.risk_threshold = 0.5;
            end
            % Compute AUC from LOOCV risk scores
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

% Clear any PerfectSeparation warning from lastwarn buffer when Firth
% handled the separation.  The orchestrator (run_dwi_pipeline.m) checks
% lastwarn after each module and logs it to error.log — clearing it here
% prevents logging a warning that Firth has already resolved.
% Perfect separation occurs when a linear combination of features perfectly
% predicts the outcome (common with small N and many features), causing
% standard MLE to diverge to infinity.  Firth's penalized likelihood
% (Jeffreys prior) guarantees finite estimates in this scenario.
if use_firth
    [~, last_warn_id] = lastwarn;
    if strcmp(last_warn_id, 'stats:lassoGlm:PerfectSeparation')
        lastwarn('');
    end
end

diary off;
end
