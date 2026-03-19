function results = imputation_sensitivity(td_panel_raw, patient_ids, feature_names, ...
    y_clean, id_list_impute, dl_provenance, dtype, dtype_label, use_firth, output_folder)
% IMPUTATION_SENSITIVITY  Compare KNN imputation against alternatives.
%
%   Compares four imputation strategies on the time-dependent panel:
%     1. KNN (existing pipeline default)
%     2. LOCF (Last Observation Carried Forward)
%     3. Mean imputation (per-feature mean from non-missing values)
%     4. Linear interpolation between adjacent observed timepoints
%
%   For each method, imputes missing values in td_panel_raw, runs the
%   same elastic net + LOOCV pipeline, and records the AUC.  Returns a
%   summary struct with method names, AUC values, imputed counts, and a
%   concordance matrix (pairwise Spearman correlation of risk scores).
%
% Inputs:
%   td_panel_raw    - [n_obs x n_feat] raw panel with NaN for missing values
%   patient_ids     - [n_obs x 1] cell array of patient IDs per row
%   feature_names   - {1 x n_feat} cell array of feature names
%   y_clean         - [n_pts x 1] binary outcome for elastic net
%   id_list_impute  - {n_pts x 1} cell array of patient IDs for LOOCV
%   dl_provenance   - DL provenance struct for leakage checks
%   dtype           - DWI type index (1=Standard, 2=dnCNN, 3=IVIMnet)
%   dtype_label     - DWI type label string
%   use_firth       - Whether to use Firth penalty in logistic regression
%   output_folder   - Where to save figures
%
% Outputs:
%   results         - Struct with fields:
%                       method_names    - {1x4} cell array
%                       auc_values      - [1x4] AUC per method
%                       n_imputed_per_method - [1x4] count of imputed values
%                       concordance     - [4x4] Spearman correlation of risk scores
%                       selected_features - {1x4} cell of selected feature indices per method

    fprintf('\n  --- Imputation Sensitivity Analysis ---\n');

    % Define available imputation strategies
    strategies = get_imputation_strategies();
    method_names = {strategies.name};
    n_methods = length(strategies);
    
    nan_mask = isnan(td_panel_raw);
    total_missing = sum(nan_mask(:));

    auc_values = nan(1, n_methods);
    n_imputed = zeros(1, n_methods);
    risk_scores_all = nan(numel(y_clean), n_methods);
    selected_features = cell(1, n_methods);

    % Execute each imputation strategy
    for i = 1:n_methods
        strategy = strategies(i);
        fprintf('  [%d/%d] %s imputation...\n', i, n_methods, strategy.name);
        
        % Apply imputation strategy
        X_imputed = strategy.impute_func(td_panel_raw, patient_ids);
        
        % Count imputed values
        if strcmp(strategy.name, 'KNN')
            % For KNN, all missing values are imputed
            n_imputed(i) = total_missing;
        else
            n_imputed(i) = sum(~isnan(X_imputed(:)) & nan_mask(:));
        end
        
        % Evaluate imputed data
        [auc_values(i), risk_scores_all(:,i), selected_features{i}] = evaluate_imputed( ...
            X_imputed, y_clean, id_list_impute, dl_provenance, dtype, dtype_label, ...
            use_firth, feature_names, strategy.name);
    end

    % --- Concordance matrix (pairwise Spearman) ---
    concordance = eye(n_methods);
    for i = 1:n_methods
        for j = (i+1):n_methods
            valid = ~isnan(risk_scores_all(:,i)) & ~isnan(risk_scores_all(:,j));
            if sum(valid) >= 3
                concordance(i,j) = corr(risk_scores_all(valid,i), ...
                    risk_scores_all(valid,j), 'Type', 'Spearman');
                concordance(j,i) = concordance(i,j);
            end
        end
    end

    % --- Print comparison table ---
    fprintf('\n  Imputation Sensitivity Results:\n');
    fprintf('  %-20s  %8s  %12s\n', 'Method', 'AUC', 'N_Imputed');
    fprintf('  %s\n', repmat('-', 1, 44));
    for m = 1:n_methods
        fprintf('  %-20s  %8.3f  %12d\n', method_names{m}, auc_values(m), n_imputed(m));
    end

    fprintf('\n  Concordance (Spearman rho between risk scores):\n');
    fprintf('  %-20s', '');
    for m = 1:n_methods
        fprintf('  %8s', method_names{m});
    end
    fprintf('\n');
    for i = 1:n_methods
        fprintf('  %-20s', method_names{i});
        for j = 1:n_methods
            fprintf('  %8.3f', concordance(i,j));
        end
        fprintf('\n');
    end

    % --- Bar chart ---
    if ~isempty(output_folder)
        try
            fig = figure('Visible', 'off', 'Position', [100 100 600 400]);
            bar(auc_values);
            set(gca, 'XTickLabel', method_names, 'XTick', 1:n_methods);
            ylabel('AUC');
            title(sprintf('Imputation Sensitivity (%s)', dtype_label));
            ylim([0 1]);
            grid on;
            for m = 1:n_methods
                text(m, auc_values(m) + 0.02, sprintf('%.3f', auc_values(m)), ...
                    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            end
            fig_path = fullfile(output_folder, sprintf('imputation_sensitivity_%s.png', dtype_label));
            saveas(fig, fig_path);
            close(fig);
            fprintf('  📁 Imputation sensitivity plot saved: %s\n', fig_path);
        catch ME_fig
            fprintf('  ⚠️  Imputation sensitivity plot failed: %s\n', ME_fig.message);
        end
    end

    % --- Output struct ---
    results = struct();
    results.method_names = method_names;
    results.auc_values = auc_values;
    results.n_imputed_per_method = n_imputed;
    results.concordance = concordance;
    results.selected_features = selected_features;
end


%% ===== Imputation Strategy Registry =====

function strategies = get_imputation_strategies()
%GET_IMPUTATION_STRATEGIES  Returns array of available imputation strategies.
%
%   Each strategy is a struct with fields:
%     name        - String identifier for the strategy
%     impute_func - Function handle that takes (X_raw, patient_ids) and returns X_imputed
%
%   This centralized registry makes it easy to add new strategies or modify existing ones.

    strategies = [
        struct('name', 'KNN', 'impute_func', @impute_knn_strategy),
        struct('name', 'LOCF', 'impute_func', @impute_locf_strategy),
        struct('name', 'Mean', 'impute_func', @impute_mean_strategy),
        struct('name', 'Linear_Interp', 'impute_func', @impute_linear_strategy)
    ];
end


%% ===== Individual Imputation Strategies =====

function X_out = impute_knn_strategy(X_raw, patient_ids)
%IMPUTE_KNN_STRATEGY  KNN imputation using existing pipeline function.
    % Use all rows as both train and test for consistency
    [X_out, ~] = knn_impute_train_test(X_raw, zeros(0, size(X_raw, 2)), 5, patient_ids, {});
    % Fill any remaining NaN with column mean
    X_out = fill_remaining_nan_with_mean(X_out);
end

function X_out = impute_locf_strategy(X_raw, patient_ids)
%IMPUTE_LOCF_STRATEGY  Last observation carried forward per patient.
    X_out = X_raw;
    unique_pats = unique(patient_ids, 'stable');
    
    for p = 1:numel(unique_pats)
        rows = strcmp(patient_ids, unique_pats{p});
        row_idx = find(rows);
        
        for c = 1:size(X_out, 2)
            last_val = NaN;
            for r = 1:numel(row_idx)
                ri = row_idx(r);
                if ~isnan(X_out(ri, c))
                    last_val = X_out(ri, c);
                elseif ~isnan(last_val)
                    X_out(ri, c) = last_val;
                end
            end
        end
    end
    
    % Fill remaining NaN with column mean
    X_out = fill_remaining_nan_with_mean(X_out);
end

function X_out = impute_mean_strategy(X_raw, ~)
%IMPUTE_MEAN_STRATEGY  Per-feature mean from non-missing values.
    X_out = X_raw;
    for c = 1:size(X_out, 2)
        col = X_out(:, c);
        col_mean = nanmean_safe(col);
        col(isnan(col)) = col_mean;
        X_out(:, c) = col;
    end
end

function X_out = impute_linear_strategy(X_raw, patient_ids)
%IMPUTE_LINEAR_STRATEGY  Linear interpolation between adjacent observed timepoints.
    X_out = X_raw;
    unique_pats = unique(patient_ids, 'stable');
    
    for p = 1:numel(unique_pats)
        rows = strcmp(patient_ids, unique_pats{p});
        row_idx = find(rows);
        
        for c = 1:size(X_out, 2)
            vals = X_out(row_idx, c);
            observed = find(~isnan(vals));
            
            if numel(observed) >= 2
                % Interpolate between observed values
                interp_vals = interp1(observed, vals(observed), 1:numel(vals), 'linear', NaN);
                X_out(row_idx, c) = interp_vals(:);
            elseif numel(observed) == 1
                % Only one observation: fill with that value
                X_out(row_idx, c) = vals(observed);
            end
        end
    end
    
    % Fill remaining NaN with column mean
    X_out = fill_remaining_nan_with_mean(X_out);
end


%% ===== Utility Functions =====

function X_out = fill_remaining_nan_with_mean(X_in)
%FILL_REMAINING_NAN_WITH_MEAN  Fill any remaining NaN values with column means.
    X_out = X_in;
    for c = 1:size(X_out, 2)
        col = X_out(:, c);
        if any(isnan(col))
            col(isnan(col)) = nanmean_safe(col);
            X_out(:, c) = col;
        end
    end
end

function [auc_val, risk_scores, selected_indices] = evaluate_imputed( ...
    X_imputed, y_clean, id_list, dl_provenance, dtype, dtype_label, ...
    use_firth, feature_names, method_label)
%EVALUATE_IMPUTED  Run elastic net CV + LOOCV on imputed data and return AUC.
%
%   Delegates to run_elastic_net_cv for feature selection and
%   run_loocv_risk_scores for unbiased out-of-fold risk scores, matching
%   the main pipeline flow in metrics_stats_predictive.m.  This ensures
%   improvements (Firth refit, collinearity filtering, DL provenance
%   checks) automatically apply to all imputation sensitivity variants.

    n_pts = numel(y_clean);
    risk_scores = nan(n_pts, 1);
    auc_val = NaN;
    selected_indices = [];

    try
        % Feature selection via elastic net CV (mirrors main pipeline)
        original_feature_indices = 1:size(X_imputed, 2);
        [selected_indices] = run_elastic_net_cv( ...
            X_imputed, y_clean, id_list, 5, use_firth, ...
            original_feature_indices, feature_names, method_label);

        % Unbiased risk scores via nested LOOCV
        [risk_scores_oof, ~] = run_loocv_risk_scores( ...
            X_imputed, y_clean, id_list, dl_provenance, dtype, dtype_label, use_firth);
        risk_scores = risk_scores_oof;

        % Compute AUC
        valid = ~isnan(risk_scores) & ~isnan(y_clean);
        if sum(valid) >= 5 && numel(unique(y_clean(valid))) == 2
            [~, ~, ~, auc_val] = perfcurve(y_clean(valid), risk_scores(valid), 1);
        end
    catch ME
        fprintf('    ⚠️  Evaluation failed for %s: %s\n', method_label, ME.message);
    end
end