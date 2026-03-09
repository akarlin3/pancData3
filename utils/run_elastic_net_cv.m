function [selected_indices, opt_lambda, common_Lambda, cv_failed, keep_fold_counts, coefs_en, final_feature_indices] = run_elastic_net_cv( ...
    X_impute, y_clean, id_list_impute, n_folds_en, use_firth, ...
    original_feature_indices, feat_names_lasso_full, fx_label)
% RUN_ELASTIC_NET_CV  5-fold elastic net CV + final model fitting.
%
%   Performs grouped 5-fold cross-validation with elastic net (alpha=0.5)
%   to select the optimal regularisation lambda, then fits a final model
%   on the full imputed dataset.  Optionally refits selected features
%   with Firth penalized logistic regression for bias-corrected coefficients.
%
% Inputs:
%   X_impute               - Feature matrix (patients x features), may contain NaN
%   y_clean                - Binary outcome vector (0/1)
%   id_list_impute         - Cell array of patient IDs for grouped folds
%   n_folds_en             - Number of CV folds (from make_grouped_folds)
%   use_firth              - Whether to refit with Firth penalty
%   original_feature_indices - Mapping from current columns to original 22 features
%   feat_names_lasso_full  - Full (unfiltered) feature names for display
%   fx_label               - Fraction label for console output
%
% Outputs:
%   selected_indices       - Indices (into original 22) of nonzero features
%   opt_lambda             - Optimal lambda from CV
%   common_Lambda          - Lambda grid used across folds
%   cv_failed              - True if elastic net failed to converge
%   keep_fold_counts       - Per-feature count of folds where feature was retained
%   coefs_en               - Final elastic net coefficient vector
%   final_feature_indices  - Original indices after consensus collinearity filtering

    n_lambdas = 10;  % lambda grid size for elastic net regularisation path

    common_Lambda = [];
    cv_failed = false;
    n_features_impute = size(X_impute, 2);
    keep_fold_counts = zeros(1, n_features_impute);  % track feature retention across folds

    % Fix random seed for reproducibility of fold assignments.
    % make_grouped_folds ensures patient-stratified folds: all scans from
    % one patient are in the same fold to prevent intra-patient data leakage
    % (where the model could memorise patient-specific patterns from training
    % data and exploit them when the same patient appears in the test fold).
    rng(42);
    fold_id_en = make_grouped_folds(id_list_impute, y_clean, 5);
    n_folds_en = max(fold_id_en);

    selected_indices = [];
    opt_lambda = [];
    coefs_en = [];
    final_feature_indices = [];

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
            p = 1 ./ (1 + exp(-eta));       % logistic transform: log-odds -> probability
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

            % --- Firth refit on elastic-net-selected features ---
            % Elastic net selects features via L1 sparsity; Firth produces
            % finite, bias-corrected coefficients on the selected subset.
            % This is especially important when perfect separation occurs
            % (common in small pancreatic cancer cohorts), where standard
            % MLE coefficients diverge to infinity.
            if use_firth && ~isempty(selected_indices)
                try
                    X_firth = X_clean_kept(:, nz_in_kept);
                    mdl_firth = fitglm(X_firth, y_clean, ...
                        'Distribution', 'binomial', ...
                        'LikelihoodPenalty', 'jeffreys-prior');
                    fprintf('  Firth refit successful for %s (%d features).\n', ...
                        fx_label, numel(selected_indices));
                catch ME_firth
                    fprintf('  ⚠️  Firth refit failed for %s: %s. Using elastic net coefficients.\n', ...
                        fx_label, ME_firth.message);
                end
            end
        catch
            cv_failed = true;
        end
        warning(w_state_final);
    end

    if cv_failed || isempty(common_Lambda)
        fprintf('Elastic Net failed to converge for %s. Fallback to empty selection.\n', fx_label);
        selected_indices = [];
    end
end
