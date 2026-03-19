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

    % Common validation
    validation_result = validate_inputs(td_panel_raw, patient_ids, y_clean, id_list_impute);
    if ~validation_result.valid
        error('Input validation failed: %s', validation_result.message);
    end

    % Missing data pattern analysis
    missing_pattern = analyze_missing_pattern(td_panel_raw, patient_ids);
    
    % Get imputation strategies
    strategies = ImputationRegistry.get_all_strategies();
    method_names = {strategies.name};
    n_methods = length(strategies);

    % Initialize results
    auc_values = nan(1, n_methods);
    n_imputed = zeros(1, n_methods);
    risk_scores_all = nan(numel(y_clean), n_methods);
    selected_features = cell(1, n_methods);

    % Execute each imputation strategy
    for i = 1:n_methods
        strategy = strategies(i);
        fprintf('  [%d/%d] %s imputation...\n', i, n_methods, strategy.name);
        
        % Apply imputation strategy
        X_imputed = strategy.impute(td_panel_raw, patient_ids);
        
        % Count imputed values
        n_imputed(i) = count_imputed_values(td_panel_raw, X_imputed, strategy.name, missing_pattern);
        
        % Evaluate imputed data
        [auc_values(i), risk_scores_all(:,i), selected_features{i}] = evaluate_imputed( ...
            X_imputed, y_clean, id_list_impute, dl_provenance, dtype, dtype_label, ...
            use_firth, feature_names, strategy.name);
    end

    % Compute concordance matrix
    concordance = compute_concordance_matrix(risk_scores_all, method_names);

    % Print results
    print_results_table(method_names, auc_values, n_imputed, concordance);

    % Generate visualization
    if ~isempty(output_folder)
        create_sensitivity_plot(auc_values, method_names, dtype_label, output_folder);
    end

    % Build output struct
    results = struct();
    results.method_names = method_names;
    results.auc_values = auc_values;
    results.n_imputed_per_method = n_imputed;
    results.concordance = concordance;
    results.selected_features = selected_features;
end


%% ===== Validation Functions =====

function result = validate_inputs(td_panel_raw, patient_ids, y_clean, id_list_impute)
%VALIDATE_INPUTS  Common input validation logic.
    result = struct('valid', true, 'message', '');
    
    if isempty(td_panel_raw) || isempty(patient_ids)
        result.valid = false;
        result.message = 'Input data or patient IDs are empty';
        return;
    end
    
    if size(td_panel_raw, 1) ~= numel(patient_ids)
        result.valid = false;
        result.message = 'Data rows must match patient ID count';
        return;
    end
    
    if numel(y_clean) ~= numel(id_list_impute)
        result.valid = false;
        result.message = 'Outcome vector must match patient list length';
        return;
    end
    
    if sum(~isnan(y_clean)) < 5
        result.valid = false;
        result.message = 'Insufficient non-missing outcomes for analysis';
        return;
    end
end

function pattern = analyze_missing_pattern(td_panel_raw, patient_ids)
%ANALYZE_MISSING_PATTERN  Analyze missing data patterns.
    nan_mask = isnan(td_panel_raw);
    
    pattern = struct();
    pattern.total_missing = sum(nan_mask(:));
    pattern.missing_by_feature = sum(nan_mask, 1);
    pattern.missing_by_patient = containers.Map();
    
    unique_pats = unique(patient_ids, 'stable');
    for p = 1:numel(unique_pats)
        rows = strcmp(patient_ids, unique_pats{p});
        pattern.missing_by_patient(unique_pats{p}) = sum(nan_mask(rows, :), 'all');
    end
end

function n_imputed = count_imputed_values(td_panel_raw, X_imputed, method_name, missing_pattern)
%COUNT_IMPUTED_VALUES  Count number of values imputed by strategy.
    if strcmp(method_name, 'KNN')
        % For KNN, all missing values are imputed
        n_imputed = missing_pattern.total_missing;
    else
        nan_mask_original = isnan(td_panel_raw);
        n_imputed = sum(~isnan(X_imputed(:)) & nan_mask_original(:));
    end
end


%% ===== Imputation Strategy Registry =====

classdef ImputationRegistry < handle
%IMPUTATIONREGISTRY  Central registry for imputation strategies.
    
    methods (Static)
        function strategies = get_all_strategies()
        %GET_ALL_STRATEGIES  Returns array of available imputation strategies.
            strategies = [
                KNNImputation(),
                LOCFImputation(),
                MeanImputation(),
                LinearImputation()
            ];
        end
    end
end


%% ===== Abstract Imputation Interface =====

classdef (Abstract) AbstractImputation < handle
%ABSTRACTIMPUTATION  Abstract base class for imputation strategies.
    
    properties (Abstract, Constant)
        name  % String identifier for the strategy
    end
    
    methods (Abstract)
        X_out = impute(obj, X_raw, patient_ids)
        % IMPUTE  Apply imputation strategy to raw data.
        %
        % Inputs:
        %   X_raw       - [n_obs x n_feat] raw data with NaN for missing values
        %   patient_ids - [n_obs x 1] cell array of patient IDs per row
        %
        % Outputs:
        %   X_out       - [n_obs x n_feat] imputed data
    end
    
    methods (Access = protected)
        function X_out = fill_remaining_nan_with_mean(~, X_in)
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
    end
end


%% ===== Concrete Imputation Implementations =====

classdef KNNImputation < AbstractImputation
%KNNIMPUTATION  KNN imputation using existing pipeline function.
    
    properties (Constant)
        name = 'KNN'
    end
    
    methods
        function X_out = impute(obj, X_raw, patient_ids)
            % Use all rows as both train and test for consistency
            [X_out, ~] = knn_impute_train_test(X_raw, zeros(0, size(X_raw, 2)), 5, patient_ids, {});
            % Fill any remaining NaN with column mean
            X_out = obj.fill_remaining_nan_with_mean(X_out);
        end
    end
end

classdef LOCFImputation < AbstractImputation
%LOCFIMPUTATION  Last observation carried forward per patient.
    
    properties (Constant)
        name = 'LOCF'
    end
    
    methods
        function X_out = impute(obj, X_raw, patient_ids)
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
            X_out = obj.fill_remaining_nan_with_mean(X_out);
        end
    end
end

classdef MeanImputation < AbstractImputation
%MEANIMPUTATION  Per-feature mean from non-missing values.
    
    properties (Constant)
        name = 'Mean'
    end
    
    methods
        function X_out = impute(~, X_raw, ~)
            X_out = X_raw;
            for c = 1:size(X_out, 2)
                col = X_out(:, c);
                col_mean = nanmean_safe(col);
                col(isnan(col)) = col_mean;
                X_out(:, c) = col;
            end
        end
    end
end

classdef LinearImputation < AbstractImputation
%LINEARIMPUTATION  Linear interpolation between adjacent observed timepoints.
    
    properties (Constant)
        name = 'Linear_Interp'
    end
    
    methods
        function X_out = impute(obj, X_raw, patient_ids)
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
            X_out = obj.fill_remaining_nan_with_mean(X_out);
        end
    end
end


%% ===== Analysis Functions =====

function concordance = compute_concordance_matrix(risk_scores_all, method_names)
%COMPUTE_CONCORDANCE_MATRIX  Compute pairwise Spearman correlation of risk scores.
    n_methods = size(risk_scores_all, 2);
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
end

function print_results_table(method_names, auc_values, n_imputed, concordance)
%PRINT_RESULTS_TABLE  Print comparison table and concordance matrix.
    n_methods = length(method_names);
    
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
end

function create_sensitivity_plot(auc_values, method_names, dtype_label, output_folder)
%CREATE_SENSITIVITY_PLOT  Generate bar chart visualization.
    try
        fig = figure('Visible', 'off', 'Position', [100 100 600 400]);
        bar(auc_values);
        set(gca, 'XTickLabel', method_names, 'XTick', 1:length(method_names));
        ylabel('AUC');
        title(sprintf('Imputation Sensitivity (%s)', dtype_label));
        ylim([0 1]);
        grid on;
        
        for m = 1:length(method_names)
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