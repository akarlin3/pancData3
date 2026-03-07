function [risk_scores_all, is_high_risk, times_km, events_km] = metrics_stats_predictive(valid_pts, lf_group, dtype_label, output_folder, dataloc, nTp, m_gtv_vol, adc_sd, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub, id_list, dtype, dl_provenance, x_labels, m_lf, m_total_time, m_total_follow_up_time)
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
%
% Outputs:
%   risk_scores_all   - LOOCV out-of-fold elastic-net computed risk scores
%   is_high_risk      - Binarised array demarcating patients > median training risk
%   times_km          - Times used for subsequent Kaplan-Meier models
%   events_km         - Events matched to times_km
%

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

    %% --- Elastic Net Feature Selection ---
    % Assemble the feature matrix with 22 candidate predictors.
    % Include baseline (Fx1) absolute values as covariates to adjust for
    % pre-existing differences in tumor burden and diffusion properties.
    % Without baseline adjustment, observed differences in change metrics
    % could be confounded by baseline tumour characteristics (e.g., larger
    % tumours may show smaller percent changes due to dilution effects).
    X_lasso_all = [ADC_abs(valid_pts, 1), D_abs(valid_pts, 1), ...
                   f_abs(valid_pts, 1),   Dstar_abs(valid_pts, 1), ...
                   ADC_abs(valid_pts, target_fx), D_abs(valid_pts, target_fx), ...
                   f_abs(valid_pts, target_fx),   Dstar_abs(valid_pts, target_fx), ...
                   ADC_pct(valid_pts, target_fx), D_pct(valid_pts, target_fx), ...
                   f_delta(valid_pts, target_fx),   Dstar_pct(valid_pts, target_fx), ...
                   m_d95_gtvp(valid_pts, target_fx), m_v50gy_gtvp(valid_pts, target_fx), ...
                   d95_adc_sub(valid_pts, target_fx), v50_adc_sub(valid_pts, target_fx), ...
                   d95_d_sub(valid_pts, target_fx),   v50_d_sub(valid_pts, target_fx), ...
                   d95_f_sub(valid_pts, target_fx),   v50_f_sub(valid_pts, target_fx), ...
                   d95_dstar_sub(valid_pts, target_fx), v50_dstar_sub(valid_pts, target_fx)];


    feat_names_lasso = {'ADC_BL', 'D_BL', 'f_BL', 'Dstar_BL', ...
                        'ADC_Abs', 'D_Abs', 'f_Abs', 'Dstar_Abs', ...
                        'ADC_Pct', 'D_Pct', 'f_Delta', 'Dstar_Pct', ...
                        'D95_GTVp', 'V50_GTVp', ...
                        'D95_Sub_ADC', 'V50_Sub_ADC', ...
                        'D95_Sub_D', 'V50_Sub_D', ...
                        'D95_Sub_f', 'V50_Sub_f', ...
                        'D95_Sub_Dstar', 'V50_Sub_Dstar'};

    original_feature_indices = 1:22;

    if target_fx == nTp || target_fx == 6
        % Post-treatment timepoint: exclude dose features (columns 13-22).
        % Dose metrics are only meaningful during active treatment when the
        % dose is being delivered.  At the post-treatment scan (typically
        % 3 months after RT), the full dose has been delivered and there is
        % no additional dose to correlate with — the therapeutic window for
        % dose-response analysis has closed.
        X_lasso_all = X_lasso_all(:, 1:12);
        feat_names_lasso = feat_names_lasso(1:12);
        original_feature_indices = original_feature_indices(1:12);
    end
    
    if iscell(X_lasso_all)
        vars = {'ADC_abs_BL', 'D_abs_BL', 'f_abs_BL', 'Dstar_abs_BL', 'ADC_abs', 'D_abs', 'f_abs', 'Dstar_abs', 'ADC_pct', 'D_pct', 'f_delta', 'Dstar_pct', 'm_d95_gtvp', 'm_v50gy_gtvp', 'd95_adc_sub', 'v50_adc_sub', 'd95_d_sub', 'v50_d_sub', 'd95_f_sub', 'v50_f_sub', 'd95_dstar_sub', 'v50_dstar_sub'};
        vars_vals = {ADC_abs, D_abs, f_abs, Dstar_abs, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub};
        fp = fopen(fullfile(output_folder, 'debug_concat_error.txt'), 'a');
        fprintf(fp, '\n--- X_lasso_all is a cell array at fx_label=%s ---\n', fx_label);
        for i_v = 1:length(vars)
            tmp_v = vars_vals{i_v};
            fprintf(fp, '%s -> Size: %s, Class: %s\n', vars{i_v}, mat2str(size(tmp_v)), class(tmp_v));
        end
        fclose(fp);
        error('Invalid data type: cell array detected in X_lasso_all concatenation. See debug_concat_error.txt');
    end

    valid_cols = ~all(isnan(X_lasso_all), 1);
    feat_names_lasso_full = feat_names_lasso;  % keep unfiltered copy for display
    X_lasso_all = X_lasso_all(:, valid_cols);
    feat_names_lasso = feat_names_lasso(valid_cols);
    original_feature_indices = original_feature_indices(valid_cols);

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

    % Fix random seed for reproducibility of fold assignments.
    % make_grouped_folds ensures patient-stratified folds: all scans from
    % one patient are in the same fold to prevent intra-patient data leakage
    % (where the model could memorise patient-specific patterns from training
    % data and exploit them when the same patient appears in the test fold).
    rng(42);
    fold_id_en = make_grouped_folds(id_list_impute, y_clean, 5);
    n_folds_en = max(fold_id_en);
    n_lambdas = 10;  % lambda grid size for elastic net regularisation path

    common_Lambda = [];
    cv_failed = false;
    n_features_impute = size(X_impute, 2);
    keep_fold_counts = zeros(1, n_features_impute);  % track feature retention across folds

    w_state_cv = warning('off', 'all');
    for cv_i = 1:n_folds_en
        text_progress_bar(cv_i, n_folds_en, 'Elastic Net CV');
        tr_idx = (fold_id_en ~= cv_i);
        te_idx = (fold_id_en == cv_i);
        X_tr = X_impute(tr_idx, :); y_tr = y_clean(tr_idx);
        X_te = X_impute(te_idx, :); y_te = y_clean(te_idx);

        id_tr = id_list_impute(tr_idx);
        id_te = id_list_impute(te_idx);
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, 5, id_tr, id_te);

        % Compute collinearity mask per fold from training data only
        keep_fold = filter_collinear_features(X_tr_imp, y_tr);
        keep_mask = false(1, n_features_impute);
        keep_mask(keep_fold) = true;
        keep_fold_counts = keep_fold_counts + double(keep_mask);
        X_tr_kept = X_tr_imp(:, keep_fold);
        X_te_kept = X_te_imp(:, keep_fold);
        
        try
            if cv_i == 1
                [B_fold, FitInfo_fold] = lassoglm(X_tr_kept, y_tr, 'binomial', ...
                    'Alpha', 0.5, 'NumLambda', n_lambdas, 'Standardize', true, 'MaxIter', 1e7);
                common_Lambda = FitInfo_fold.Lambda;
                all_deviance = zeros(length(common_Lambda), n_folds_en);
            else
                [B_fold, FitInfo_fold] = lassoglm(X_tr_kept, y_tr, 'binomial', ...
                    'Alpha', 0.5, 'Lambda', common_Lambda, 'Standardize', true, 'MaxIter', 1e7);
            end
            
            % Compute test-fold predicted probabilities via logistic sigmoid
            eta = X_te_kept * B_fold + FitInfo_fold.Intercept;
            p = 1 ./ (1 + exp(-eta));       % logistic transform: log-odds → probability
            p = max(min(p, 1 - 1e-10), 1e-10);  % clamp to avoid log(0) in deviance
            y_te_mat = repmat(y_te, 1, length(common_Lambda));

            % Binomial deviance = -2 * log-likelihood.  This is the CV loss
            % function for selecting the optimal lambda: the lambda that
            % minimises mean deviance across folds provides the best
            % bias-variance tradeoff for prediction.
            dev_val = -2 * sum(y_te_mat .* log(p) + (1 - y_te_mat) .* log(1 - p), 1);
            all_deviance(:, cv_i) = dev_val';
        catch
            cv_failed = true;
            break;
        end
    end
    warning(w_state_cv);

    if ~cv_failed && ~isempty(common_Lambda)
        mean_deviance = mean(all_deviance, 2);
        [~, idx_min] = min(mean_deviance);
        opt_lambda = common_Lambda(idx_min);
        
        % NOTE: Final model uses in-sample KNN imputation (entire dataset as
        % its own reference).  This is intentional — the final model is NOT
        % used for performance estimation (LOOCV handles that).  In-sample
        % imputation here maximises the data available for coefficient
        % estimation in the deployed model.  See LOOCV section below for
        % the unbiased performance estimate which uses proper train/test
        % imputation splits.
        X_clean_all = knn_impute_train_test(X_impute, [], 5, id_list_impute);

        % Use consensus collinearity mask from CV folds (feature kept in
        % majority of folds) to avoid data leakage from filtering on the
        % full cohort.  See filter_collinear_features docstring.
        keep_final = (keep_fold_counts > n_folds_en / 2);
        X_clean_kept = X_clean_all(:, keep_final);
        final_feature_indices = original_feature_indices(keep_final);

        w_state_final = warning('off', 'stats:lassoGlm:IterationLimit');
        warning('off', 'stats:lassoGlm:PerfectSeparation');
        try
            [B_final, FitInfo_final] = lassoglm(X_clean_kept, y_clean, 'binomial', ...
                'Alpha', 0.5, 'Lambda', opt_lambda, 'Standardize', true, 'MaxIter', 1e7);

            coefs_en = B_final;

            % Map nonzero coefficient indices to original feature names.
            nz_in_kept = find(coefs_en ~= 0);
            selected_indices = final_feature_indices(nz_in_kept);

            fprintf('Elastic Net Selected Features for %s (Opt Lambda=%.4f): %s\n', ...
                fx_label, opt_lambda, strjoin(feat_names_lasso_full(selected_indices), ', '));
        catch
            cv_failed = true;
        end
        warning(w_state_final);
    end
    
    if cv_failed || isempty(common_Lambda)
        fprintf('Elastic Net failed to converge for %s. Fallback to empty selection.\n', fx_label);
        selected_indices = [];
    end

    %% --- LOOCV For Risk Scores ---
    % NESTED LEAVE-ONE-OUT CROSS-VALIDATION (LOOCV)
    % For each patient i:
    %   Outer loop: hold out patient i as test
    %   Inner loop: 5-fold CV on remaining N-1 patients to select lambda
    %   Fit elastic net with optimal lambda on N-1 patients
    %   Predict risk score for held-out patient i
    %
    % This produces truly out-of-fold risk scores: no patient's data
    % influences its own risk prediction.  LOOCV is used instead of k-fold
    % because the small cohort size (typical N=30-60 in pancreatic DWI
    % studies) means k-fold would have very small test sets, and LOOCV
    % maximises training set size for stable coefficient estimation.
    n_pts_impute = size(X_impute, 1);
    risk_scores_oof = nan(n_pts_impute, 1);
    is_high_risk_oof = false(n_pts_impute, 1);

    fprintf('  Generating unbiased out-of-fold risk scores via nested LOOCV...\n');
    for loo_i = 1:n_pts_impute
        text_progress_bar(loo_i, n_pts_impute, 'LOOCV risk scores');
        pat_id_i = id_list_impute{loo_i};
        is_leaky = false;
        if dtype == 2 
            if any(strcmp(dl_provenance.dncnn_train_ids, pat_id_i))
                is_leaky = true;
            end
        elseif dtype == 3 
            if any(strcmp(dl_provenance.ivimnet_train_ids, pat_id_i))
                is_leaky = true;
            end
        end
        
        if is_leaky
            error('DATA LEAKAGE DETECTED: Patient %s was used to train the %s model. Fundamental isolation broken.', ...
                pat_id_i, dtype_label);
        end

        train_mask = true(n_pts_impute, 1);
        train_mask(loo_i) = false;
        X_tr_fold = X_impute(train_mask, :);
        y_tr_fold = y_clean(train_mask);
        X_te_fold = X_impute(loo_i, :);
        
        id_tr_loo = id_list_impute(train_mask);
        id_te_loo = id_list_impute(loo_i);
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr_fold, X_te_fold, 5, id_tr_loo, id_te_loo);

        % Compute collinearity mask per LOO fold from training data only
        keep_fold = filter_collinear_features(X_tr_imp, y_tr_fold);
        X_tr_kept = X_tr_imp(:, keep_fold);
        X_te_kept = X_te_imp(:, keep_fold);
        
        w_state_loo = warning('off', 'all');
        try
            inn_fold_id = make_grouped_folds(id_tr_loo, y_tr_fold, 5);
            n_inn_folds = max(inn_fold_id);
            inn_Lambda = []; inn_deviance = [];
            for inn_i = 1:n_inn_folds
                inn_tr = (inn_fold_id ~= inn_i);
                inn_te = (inn_fold_id == inn_i);
                if inn_i == 1
                    [B_inn, FI_inn] = lassoglm(X_tr_kept(inn_tr,:), y_tr_fold(inn_tr), 'binomial', ...
                        'Alpha', 0.5, 'NumLambda', 10, 'Standardize', true, 'MaxIter', 1e7);
                    inn_Lambda = FI_inn.Lambda;
                    inn_deviance = zeros(numel(inn_Lambda), n_inn_folds);
                else
                    [B_inn, FI_inn] = lassoglm(X_tr_kept(inn_tr,:), y_tr_fold(inn_tr), 'binomial', ...
                        'Alpha', 0.5, 'Lambda', inn_Lambda, 'Standardize', true, 'MaxIter', 1e7);
                end
                p_inn = 1 ./ (1 + exp(-(X_tr_kept(inn_te,:) * B_inn + FI_inn.Intercept)));
                p_inn = max(min(p_inn, 1 - 1e-10), 1e-10);
                y_te_m = repmat(y_tr_fold(inn_te), 1, numel(inn_Lambda));
                inn_deviance(:, inn_i) = (-2 * sum(y_te_m .* log(p_inn) + (1 - y_te_m) .* log(1 - p_inn), 1))';
            end
            if isempty(inn_Lambda)
                error('grouped_cv:noLambda', 'No lambda path computed in inner CV.');
            end
            [~, best_idx] = min(mean(inn_deviance, 2));
            opt_lam_loo = inn_Lambda(best_idx);
            [B_loo, FitInfo_loo] = lassoglm(X_tr_kept, y_tr_fold, 'binomial', ...
                'Alpha', 0.5, 'Lambda', opt_lam_loo, 'Standardize', true, 'MaxIter', 1e7);
            coefs_loo = B_loo;
            intercept_loo = FitInfo_loo.Intercept;
        catch
            % Mark this patient's risk score as NaN rather than producing a
            % meaningless zero-coefficient prediction that could corrupt
            % downstream median-based stratification and ROC analysis.
            coefs_loo = nan(size(X_tr_kept, 2), 1);
            intercept_loo = NaN;
        end
        warning(w_state_loo);

        risk_scores_oof(loo_i) = X_te_kept * coefs_loo + intercept_loo;
    end
    
    % Stratify patients into high/low risk groups using the median of all
    % out-of-fold risk scores.  The median is computed AFTER all LOOCV
    % iterations complete — computing it per-fold would be circular
    % (each fold's threshold would depend on which patient was left out).
    % Median-based stratification creates balanced groups for Kaplan-Meier
    % analysis downstream, which maximises statistical power for detecting
    % survival differences.
    valid_oof = ~isnan(risk_scores_oof);
    oof_median = median(risk_scores_oof(valid_oof));
    is_high_risk_oof = zeros(size(risk_scores_oof));
    is_high_risk_oof(valid_oof) = risk_scores_oof(valid_oof) > oof_median;
    % Patients with NaN risk scores (failed LOOCV folds) are excluded from
    % stratification rather than silently assigned to the low-risk group.
    is_high_risk_oof(~valid_oof) = NaN;

    risk_scores_all_target = nan(sum(valid_pts), 1);
    risk_scores_all_target(impute_mask) = risk_scores_oof;

    is_high_risk_target = nan(sum(valid_pts), 1);
    is_high_risk_target(impute_mask) = is_high_risk_oof;

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
    
    labels = lf_group; % 0 = LC, 1 = LF
    % Exclude competing-risk patients (lf==2) from ROC analysis.
    % perfcurve treats non-positive-class as negatives, so lf==2 patients
    % would be lumped with genuine LC patients, inflating AUC.
    roc_eligible = (labels <= 1);
    valid_roc = roc_eligible & ~isnan(risk_scores_all_target) & ~isnan(labels);

    if sum(valid_roc) > 0
        % ROC curve from out-of-fold risk scores.  Because these scores
        % were generated via nested LOOCV (no patient influenced its own
        % prediction), the resulting AUC is an unbiased estimate of the
        % model's discriminative ability for new patients.
        [roc_X, roc_Y, roc_T, roc_AUC] = perfcurve(labels(valid_roc), risk_scores_all_target(valid_roc), 1);

        % Youden's J statistic = Sensitivity + Specificity - 1
        % The optimal threshold maximises the sum of sensitivity and
        % specificity, balancing the cost of missing true failures (false
        % negatives) against incorrectly flagging controlled patients
        % (false positives).
        [~, roc_opt_idx] = max(roc_Y - roc_X);
        roc_opt_thresh = roc_T(roc_opt_idx);
        
        figure('Name', ['ROC Analysis - ' fx_label ' — ' dtype_label], 'Position', [200, 200, 700, 600]);
        hold on;
        
        plot(roc_X, roc_Y, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2.5);
        leg_entries = {sprintf('LOOCV OOF Risk Score (AUC = %.3f)', roc_AUC)};
        
        plot(roc_X(roc_opt_idx), roc_Y(roc_opt_idx), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
        leg_entries{end+1} = sprintf('Youden Cutoff (score = %.3f)', roc_opt_thresh);
        
        plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);
        leg_entries{end+1} = 'Random Guess';
        
        xlabel('False Positive Rate (1 - Specificity)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('True Positive Rate (Sensitivity)', 'FontSize', 12, 'FontWeight', 'bold');
        title(['PRIMARY ROC Curve: LOOCV Out-of-Fold Risk Score (' fx_label ', ' dtype_label ')'], 'FontSize', 14);
        
        legend(leg_entries, 'Location', 'SouthEast', 'FontSize', 11);
        grid on; box on;
        hold off;
        set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);
        saveas(gcf, fullfile(output_folder, ['ROC_OOF_Risk_Score_' fx_label '_' dtype_label '.png']));
        close(gcf);
        
        fprintf('\n--- PRIMARY ROC ANALYSIS (LOOCV Out-of-Fold Risk Score) for %s ---\n', fx_label);
        fprintf('  AUC  = %.3f\n  Youden Optimal Score Cutoff = %.4f\n', roc_AUC, roc_opt_thresh);
        fprintf('  Sensitivity = %.1f%%  |  Specificity = %.1f%%\n\n', ...
            roc_Y(roc_opt_idx)*100, (1-roc_X(roc_opt_idx))*100);
    else
        fprintf('\n--- PRIMARY ROC ANALYSIS (LOOCV OOF) for %s ---\n', fx_label);
        fprintf('Insufficient data for out-of-fold ROC analysis.\n\n');
    end

    %% ---------- 5. Sanity Checks & Scatter Plots ----------
    % These plots verify that selected biomarkers are not confounded:
    %   Panel 1 (Volume): If GTV volume change differs between LC/LF,
    %     diffusion changes could be an artifact of partial-volume effects
    %     (smaller tumours have more edge voxels contaminated by normal tissue).
    %   Panel 2 (ADC SD): Changes in intra-tumour heterogeneity (texture)
    %     may provide independent prognostic information beyond mean values.
    %   Panel 3 (Signal vs Noise): Overlays the Coefficient of
    %     Reproducibility (CoR) band on the scatter plot.  Data points
    %     falling within the CoR band cannot be distinguished from
    %     measurement noise — only changes exceeding CoR represent real
    %     biological signal.
    for vi = 1:n_sig
        curr_sig_pct_full = sig_data_selected{vi};
        curr_sig_name = sig_names{vi};
        if sig_is_abs(vi), curr_sig_disp = ['Abs ' curr_sig_name]; curr_sig_file = ['Abs_' curr_sig_name];
        else, curr_sig_disp = ['\Delta ' curr_sig_name]; curr_sig_file = ['Delta_' curr_sig_name]; end
        
        figure('Name', ['Sanity Checks ' curr_sig_disp ' ' fx_label ' — ' dtype_label], 'Position', [100, 100, 1200, 500]);
        subplot(1, 3, 1);
        vol_fx1 = m_gtv_vol(valid_pts, 1);          
        vol_fx3 = m_gtv_vol(valid_pts, target_fx);  
        vol_pct = (vol_fx3 - vol_fx1) ./ vol_fx1 * 100;  
        
        % Exclude competing risk patients (lf==2) from sanity check plots
        non_competing = (lf_group <= 1);
        boxplot(vol_pct(non_competing), lf_group(non_competing), 'Labels', {'LC (0)', 'LF (1)'});
        ylabel(['% Change in GTV Volume (' fx_label ')']);
        title('Confounder Check: Volume', 'FontSize', 12, 'FontWeight', 'bold');
        p_vol = perform_statistical_test(vol_pct(non_competing), lf_group(non_competing), 'ranksum');
        
        y_lim = ylim;
        if numel(y_lim) >= 2 && all(isfinite(y_lim)) && y_lim(2) > y_lim(1)
            text(1.5, y_lim(1) + 0.9*(y_lim(2)-y_lim(1)), format_p_value(p_vol), ...
                'HorizontalAlignment', 'center', 'FontSize', 11);
        end
        grid on;
        if p_vol > 0.05, xlabel('Conclusion: No Volumetric Bias');
        else, xlabel('Warning: Volume is a Confounder'); end
        
        subplot(1, 3, 2);
        sd_fx1  = adc_sd(valid_pts, 1, dtype);
        sd_fxN  = adc_sd(valid_pts, target_fx, dtype);
        sd_delta = sd_fxN - sd_fx1;
        
        boxplot(sd_delta(non_competing), lf_group(non_competing), 'Labels', {'LC (0)', 'LF (1)'});
        ylabel(['\Delta ADC SD (' fx_label ') [mm^2/s]']);
        title('Heterogeneity: ADC SD Change', 'FontSize', 12, 'FontWeight', 'bold');
        p_sd = perform_statistical_test(sd_delta(non_competing), lf_group(non_competing), 'ranksum');
        
        y_lim = ylim;
        if numel(y_lim) >= 2 && all(isfinite(y_lim)) && y_lim(2) > y_lim(1)
            text(1.5, y_lim(1) + 0.9*(y_lim(2)-y_lim(1)), format_p_value(p_sd), ...
                'HorizontalAlignment', 'center', 'FontSize', 11);
        end
        grid on;
        
        subplot(1, 3, 3);
        hold on;
        % Derive wCV from baseline ADC SD and mean instead of hardcoding.
        % wCV = SD/mean; CoR for percent change = 1.96*sqrt(2)*wCV*100.
        baseline_sd  = adc_sd(valid_pts, 1, dtype);
        baseline_adc_vals = ADC_abs(valid_pts, 1);
        wcv_vals = baseline_sd ./ baseline_adc_vals;
        % Guard against Inf from zero/near-zero baseline ADC values
        wcv_vals(~isfinite(wcv_vals)) = NaN;
        wcv_est = median(wcv_vals, 'omitnan');          % median wCV as fraction
        cor_est = 1.96 * sqrt(2) * wcv_est * 100;       % CoR in percent
        
        % Exclude competing-risk patients (lf==2) from scatter, consistent
        % with the boxplot exclusion above (non_competing mask).
        scatter_mask = (lf_group <= 1);
        x_scatter = ones(sum(scatter_mask), 1);
        x_scatter(lf_group(scatter_mask)==1) = 2;
        x_scatter = x_scatter + (rand(size(x_scatter))-0.5)*0.2;
        scatter_vals = curr_sig_pct_full(valid_pts, sig_col_idx(vi));
        scatter(x_scatter, scatter_vals(scatter_mask), 50, 'filled', 'MarkerEdgeColor', 'k');
        
        base_idx = mod(selected_indices(vi)-1, 4) + 1;
        
        if sig_is_pct_imaging(vi) && base_idx == 1
            yfill = [-cor_est cor_est cor_est -cor_est];
            xfill = [0.5 0.5 2.5 2.5];
            fill(xfill, yfill, [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            yline(0, 'k-');
            yline(cor_est, 'k--', sprintf('CoR (+%.1f%%)', cor_est));
            yline(-cor_est, 'k--', sprintf('CoR (-%.1f%%)', cor_est));
        elseif ~sig_is_abs(vi) && ~isempty(strfind(sig_units{vi}, '%'))
            yline(0, 'k-', 'Alpha', 0.3);
        end
        
        xticks([1 2]); xticklabels({'LC', 'LF'});
        ylbl = sprintf('%s (%s)', sig_disp_names{vi}, sig_units{vi});
        ylabel(ylbl);
        title('Signal vs. Noise Floor', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        xlim([0.5 2.5]);
        
        sgtitle(['Validation (' curr_sig_disp '): Volume, Texture, and Noise (' fx_label ', ' dtype_label ')'], 'FontSize', 14, 'FontWeight', 'bold');
        allAx = findall(gcf, 'Type', 'Axes');
        for k = 1:numel(allAx)
            pos = get(allAx(k), 'Position');
            set(allAx(k), 'Position', [pos(1), pos(2) * 0.92, pos(3), pos(4) * 0.92]);
        end
        safe_name = strrep(curr_sig_file, '*', 'star');
        set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);
        saveas(gcf, fullfile(output_folder, ['Sanity_Checks_' safe_name '_' fx_label '_' dtype_label '.png']));
        close(gcf);
    end

    % 2D FEATURE SPACE SCATTER PLOTS
    % When 2+ features are selected, pairwise scatter plots show how
    % combinations of biomarkers separate LC from LF in feature space.
    % A logistic decision boundary is overlaid for visualisation (fitted
    % on the full dataset, NOT cross-validated — see legend disclaimer).
    % Clinically, if two features together provide better separation than
    % either alone, this suggests a multivariate signature that captures
    % complementary aspects of treatment resistance (e.g., cellularity
    % via D + vascular damage via f).
    if n_sig >= 2
        for fi = 1:(n_sig-1)
            for fj = (fi+1):n_sig
                figure('Name', sprintf('2D Feature Space %s vs %s %s — %s', sig_names{fi}, sig_names{fj}, fx_label, dtype_label), 'Position', [100, 100, 800, 600]);
                hold on;
                
                x_val = sig_data_selected{fi}(valid_pts, sig_col_idx(fi));
                y_val = sig_data_selected{fj}(valid_pts, sig_col_idx(fj));
                % Exclude competing risk patients (lf==2) from scatter plots
                group = lf_group;
                scatter_mask = (group <= 1);
                x_val = x_val(scatter_mask);
                y_val = y_val(scatter_mask);
                group = group(scatter_mask);

                scatter(x_val(group==0), y_val(group==0), 80, [0 0.4470 0.7410], 'filled', 'MarkerEdgeColor', 'k');
                scatter(x_val(group==1), y_val(group==1), 80, [0.8500 0.3250 0.0980], 'filled', 'MarkerEdgeColor', 'k');

                % NOTE: Decision boundary is fitted on the full displayed
                % dataset (not cross-validated) for visualization only.
                w_state = warning('off', 'all');
                mdl = fitglm([x_val, y_val], group, 'Distribution', 'binomial', 'Options', statset('MaxIter', 1e7));
                warning(w_state);
                coefs = mdl.Coefficients.Estimate;
                if numel(coefs) >= 3 && coefs(3) ~= 0
                    x_range = linspace(min(x_val), max(x_val), 100);
                    y_boundary = -(coefs(1) + coefs(2)*x_range) / coefs(3);
                    plot(x_range, y_boundary, 'k--', 'LineWidth', 2);
                end

                if sig_is_abs(fi), xl = sprintf('%s at %s (%s)', sig_names{fi}, fx_label, sig_units{fi}); else, xl = sprintf('\\Delta %s at %s (%s)', sig_names{fi}, fx_label, sig_units{fi}); end
                if sig_is_abs(fj), yl = sprintf('%s at %s (%s)', sig_names{fj}, fx_label, sig_units{fj}); else, yl = sprintf('\\Delta %s at %s (%s)', sig_names{fj}, fx_label, sig_units{fj}); end
                xlabel(xl, 'FontSize', 12, 'FontWeight', 'bold');
                ylabel(yl, 'FontSize', 12, 'FontWeight', 'bold');
                title(sprintf('Biomarker Interaction: Separation of LC vs LF (%s, %s)', fx_label, dtype_label), 'FontSize', 14);
                if numel(coefs) >= 3 && coefs(3) ~= 0
                    legend({'Local Control', 'Local Failure', 'Logistic Boundary (illustrative, not CV)'}, 'Location', 'NorthWest');
                else
                    legend({'Local Control', 'Local Failure'}, 'Location', 'NorthWest');
                end
                
                grid on; box on;
                xline(0, 'k-', 'Alpha', 0.2); yline(0, 'k-', 'Alpha', 0.2);
                
                safe_name1 = strrep(sig_names{fi}, '*', 'star');
                safe_name2 = strrep(sig_names{fj}, '*', 'star');
                set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);
                saveas(gcf, fullfile(output_folder, sprintf('2D_Space_%s_vs_%s_%s_%s.png', safe_name1, safe_name2, fx_label, dtype_label)));
                close(gcf);
            end
        end
    else
        fprintf('Skipping 2D scatter plots: requires at least 2 significant variables.\n');
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

diary off;
end
