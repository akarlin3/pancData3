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

    fprintf('\n  --- Imputation Sensitivity Analysis ---\n');

    method_names = {'KNN', 'LOCF', 'Mean', 'Linear_Interp'};
    n_methods = 4;
    n_obs = size(td_panel_raw, 1);
    n_feat = size(td_panel_raw, 2);

    nan_mask = isnan(td_panel_raw);
    total_missing = sum(nan_mask(:));

    auc_values = nan(1, n_methods);
    n_imputed = zeros(1, n_methods);
    risk_scores_all = nan(numel(y_clean), n_methods);

    % --- Method 1: KNN (use existing pipeline function) ---
    fprintf('  [1/4] KNN imputation...\n');
    X_knn = impute_knn_wrapper(td_panel_raw, patient_ids);
    n_imputed(1) = total_missing;
    [auc_values(1), risk_scores_all(:,1)] = evaluate_imputed( ...
        X_knn, y_clean, id_list_impute, dl_provenance, dtype, dtype_label, use_firth);

    % --- Method 2: LOCF ---
    fprintf('  [2/4] LOCF imputation...\n');
    X_locf = impute_locf(td_panel_raw, patient_ids);
    n_imputed(2) = sum(~isnan(X_locf(:)) & nan_mask(:));
    [auc_values(2), risk_scores_all(:,2)] = evaluate_imputed( ...
        X_locf, y_clean, id_list_impute, dl_provenance, dtype, dtype_label, use_firth);

    % --- Method 3: Mean ---
    fprintf('  [3/4] Mean imputation...\n');
    X_mean = impute_mean(td_panel_raw);
    n_imputed(3) = sum(~isnan(X_mean(:)) & nan_mask(:));
    [auc_values(3), risk_scores_all(:,3)] = evaluate_imputed( ...
        X_mean, y_clean, id_list_impute, dl_provenance, dtype, dtype_label, use_firth);

    % --- Method 4: Linear Interpolation ---
    fprintf('  [4/4] Linear interpolation...\n');
    X_interp = impute_linear(td_panel_raw, patient_ids);
    n_imputed(4) = sum(~isnan(X_interp(:)) & nan_mask(:));
    [auc_values(4), risk_scores_all(:,4)] = evaluate_imputed( ...
        X_interp, y_clean, id_list_impute, dl_provenance, dtype, dtype_label, use_firth);

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
end


%% ===== Local helper functions =====

function X_out = impute_knn_wrapper(X_raw, patient_ids)
%IMPUTE_KNN_WRAPPER  KNN imputation using existing pipeline function.
    n = size(X_raw, 1);
    % Use all rows as both train and test for consistency
    [X_out, ~] = knn_impute_train_test(X_raw, zeros(0, size(X_raw, 2)), 5, patient_ids, {});
    % Fill any remaining NaN with column mean
    for c = 1:size(X_out, 2)
        col = X_out(:, c);
        if any(isnan(col))
            col(isnan(col)) = nanmean_safe(col);
            X_out(:, c) = col;
        end
    end
end

function X_out = impute_locf(X_raw, patient_ids)
%IMPUTE_LOCF  Last observation carried forward per patient.
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
    for c = 1:size(X_out, 2)
        col = X_out(:, c);
        if any(isnan(col))
            col(isnan(col)) = nanmean_safe(col);
            X_out(:, c) = col;
        end
    end
end

function X_out = impute_mean(X_raw)
%IMPUTE_MEAN  Per-feature mean from non-missing values.
    X_out = X_raw;
    for c = 1:size(X_out, 2)
        col = X_out(:, c);
        col_mean = nanmean_safe(col);
        col(isnan(col)) = col_mean;
        X_out(:, c) = col;
    end
end

function X_out = impute_linear(X_raw, patient_ids)
%IMPUTE_LINEAR  Linear interpolation between adjacent observed timepoints.
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
    for c = 1:size(X_out, 2)
        col = X_out(:, c);
        if any(isnan(col))
            col(isnan(col)) = nanmean_safe(col);
            X_out(:, c) = col;
        end
    end
end

function [auc_val, risk_scores] = evaluate_imputed(X_imputed, y_clean, id_list, dl_provenance, dtype, dtype_label, use_firth)
%EVALUATE_IMPUTED  Run elastic net + LOOCV on imputed data and return AUC.
    % Fill any remaining NaN before LOOCV
    for c = 1:size(X_imputed, 2)
        col = X_imputed(:, c);
        if any(isnan(col))
            col(isnan(col)) = nanmean_safe(col);
            X_imputed(:, c) = col;
        end
    end

    n_pts = numel(y_clean);
    risk_scores = nan(n_pts, 1);
    auc_val = NaN;

    try
        [risk_scores_oof, ~] = run_loocv_risk_scores( ...
            X_imputed, y_clean, id_list, dl_provenance, dtype, dtype_label, use_firth);
        risk_scores = risk_scores_oof;

        % Compute AUC
        valid = ~isnan(risk_scores) & ~isnan(y_clean);
        if sum(valid) >= 5 && numel(unique(y_clean(valid))) == 2
            [~, ~, ~, auc_val] = perfcurve(y_clean(valid), risk_scores(valid), 1);
        end
    catch ME
        fprintf('    ⚠️  LOOCV failed for this method: %s\n', ME.message);
    end
end
