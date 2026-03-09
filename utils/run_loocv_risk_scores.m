function [risk_scores_oof, is_high_risk_oof] = run_loocv_risk_scores( ...
    X_impute, y_clean, id_list_impute, dl_provenance, dtype, dtype_label, use_firth)
% RUN_LOOCV_RISK_SCORES  Nested LOOCV for unbiased out-of-fold risk scores.
%
%   For each patient i:
%     Outer loop: hold out patient i as test
%     Inner loop: 5-fold CV on remaining N-1 patients to select lambda
%     Fit elastic net with optimal lambda on N-1 patients
%     Predict risk score for held-out patient i
%
%   This produces truly out-of-fold risk scores: no patient's data
%   influences its own risk prediction.  LOOCV is used instead of k-fold
%   because the small cohort size (typical N=30-60 in pancreatic DWI
%   studies) means k-fold would have very small test sets, and LOOCV
%   maximises training set size for stable coefficient estimation.
%
%   DL provenance leakage checks are enforced per-patient to prevent
%   contamination from deep learning training sets.
%
% Inputs:
%   X_impute         - Feature matrix (patients x features), may contain NaN
%   y_clean          - Binary outcome vector (0/1)
%   id_list_impute   - Cell array of patient IDs for grouped folds
%   dl_provenance    - Struct with DL training provenance for leakage checks
%   dtype            - DWI type index (1=Standard, 2=dnCNN, 3=IVIMnet)
%   dtype_label      - DWI type label string for error messages
%   use_firth        - Whether to refit with Firth penalty
%
% Outputs:
%   risk_scores_oof  - Out-of-fold risk scores (n_pts x 1), NaN for failed folds
%   is_high_risk_oof - Binary stratification (0/1/NaN) based on median risk score

    n_pts_impute = size(X_impute, 1);
    risk_scores_oof = nan(n_pts_impute, 1);   % out-of-fold predicted risk scores
    is_high_risk_oof = false(n_pts_impute, 1); % binary high/low risk stratification

    fprintf('  Generating unbiased out-of-fold risk scores via nested LOOCV...\n');
    for loo_i = 1:n_pts_impute
        text_progress_bar(loo_i, n_pts_impute, 'LOOCV risk scores');
        pat_id_i = id_list_impute{loo_i};

        % --- Deep learning provenance leakage check ---
        % For DnCNN (dtype=2) and IVIMnet (dtype=3) pipelines, verify that
        % the held-out patient was NOT used to train the DL denoising model.
        % If a patient's raw images were used to train DnCNN/IVIMnet, their
        % denoised outputs would carry training-set bias, making any
        % downstream predictive model evaluation fundamentally invalid.
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

        % --- Outer loop: hold out patient i as test ---
        train_mask = true(n_pts_impute, 1);
        train_mask(loo_i) = false;
        X_tr_fold = X_impute(train_mask, :);
        y_tr_fold = y_clean(train_mask);
        X_te_fold = X_impute(loo_i, :);  % single held-out patient

        % KNN imputation with patient-ID-aware distance (prevents using
        % the held-out patient as a neighbor reference)
        id_tr_loo = id_list_impute(train_mask);
        id_te_loo = id_list_impute(loo_i);
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr_fold, X_te_fold, 5, id_tr_loo, id_te_loo);

        % Collinearity filtering on training data only (per-fold to avoid
        % using test-fold information for feature selection)
        keep_fold = filter_collinear_features(X_tr_imp, y_tr_fold);
        X_tr_kept = X_tr_imp(:, keep_fold);
        X_te_kept = X_te_imp(:, keep_fold);

        % Suppress expected warnings from elastic net on small folds
        w_state_loo = warning('off', 'stats:lassoGlm:PerfectSeparation');
        warning('off', 'stats:lassoGlm:IterationLimit');
        try
            % --- Inner loop: 5-fold CV on N-1 patients to select lambda ---
            % This nested CV ensures lambda selection is independent of the
            % held-out patient, preventing optimistic bias in risk scores.
            inn_fold_id = make_grouped_folds(id_tr_loo, y_tr_fold, 5);
            n_inn_folds = max(inn_fold_id);
            inn_Lambda = []; inn_deviance = [];
            for inn_i = 1:n_inn_folds
                inn_tr = (inn_fold_id ~= inn_i);
                inn_te = (inn_fold_id == inn_i);
                if inn_i == 1
                    % First inner fold: auto-generate lambda grid
                    [B_inn, FI_inn] = lassoglm(X_tr_kept(inn_tr,:), y_tr_fold(inn_tr), 'binomial', ...
                        'Alpha', 0.5, 'NumLambda', 10, 'Standardize', true, 'MaxIter', 1e7);
                    inn_Lambda = FI_inn.Lambda;
                    inn_deviance = zeros(numel(inn_Lambda), n_inn_folds);
                else
                    % Subsequent inner folds: reuse lambda grid from fold 1
                    [B_inn, FI_inn] = lassoglm(X_tr_kept(inn_tr,:), y_tr_fold(inn_tr), 'binomial', ...
                        'Alpha', 0.5, 'Lambda', inn_Lambda, 'Standardize', true, 'MaxIter', 1e7);
                end
                % Compute predicted probabilities via logistic sigmoid
                p_inn = 1 ./ (1 + exp(-(X_tr_kept(inn_te,:) * B_inn + FI_inn.Intercept)));
                p_inn = max(min(p_inn, 1 - 1e-10), 1e-10);  % clamp to avoid log(0)
                % Binomial deviance for lambda selection
                y_te_m = repmat(y_tr_fold(inn_te), 1, numel(inn_Lambda));
                inn_deviance(:, inn_i) = (-2 * sum(y_te_m .* log(p_inn) + (1 - y_te_m) .* log(1 - p_inn), 1))';
            end
            if isempty(inn_Lambda)
                error('grouped_cv:noLambda', 'No lambda path computed in inner CV.');
            end
            % Select lambda with minimum mean deviance across inner folds
            [~, best_idx] = min(mean(inn_deviance, 2));
            opt_lam_loo = inn_Lambda(best_idx);

            % Fit final model on all N-1 training patients with optimal lambda
            [B_loo, FitInfo_loo] = lassoglm(X_tr_kept, y_tr_fold, 'binomial', ...
                'Alpha', 0.5, 'Lambda', opt_lam_loo, 'Standardize', true, 'MaxIter', 1e7);
            coefs_loo = B_loo;
            intercept_loo = FitInfo_loo.Intercept;

            % --- Firth refit in LOOCV ---
            % Refit selected features with Firth to produce finite
            % coefficients when elastic net encounters separation.
            if use_firth
                nz_loo = find(coefs_loo ~= 0);
                if ~isempty(nz_loo)
                    try
                        mdl_firth_loo = fitglm(X_tr_kept(:, nz_loo), y_tr_fold, ...
                            'Distribution', 'binomial', ...
                            'LikelihoodPenalty', 'jeffreys-prior');
                        % Replace coefficients with Firth estimates
                        firth_coefs = zeros(size(coefs_loo));
                        firth_coefs(nz_loo) = mdl_firth_loo.Coefficients.Estimate(2:end);
                        coefs_loo = firth_coefs;
                        intercept_loo = mdl_firth_loo.Coefficients.Estimate(1);
                    catch
                        % Firth failed — keep elastic net coefficients
                    end
                end
            end
        catch
            % Mark this patient's risk score as NaN rather than producing a
            % meaningless zero-coefficient prediction that could corrupt
            % downstream median-based stratification and ROC analysis.
            coefs_loo = nan(size(X_tr_kept, 2), 1);
            intercept_loo = NaN;
        end
        warning(w_state_loo);

        % Compute the risk score for the held-out patient: linear predictor
        % (log-odds) from the elastic net model. Higher values = higher
        % predicted probability of local failure.
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
end
