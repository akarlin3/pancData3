function [risk_scores_all, is_high_risk, times_km, events_km] = metrics_stats(valid_pts, lf_group, metric_sets, set_names, time_labels, dtype_label, output_folder, dataloc, nTp, m_gtv_vol, adc_sd, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_pct, Dstar_pct, m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub, id_list, dtype, dl_provenance, x_labels, m_id_list, m_lf, m_total_time, m_total_follow_up_time)
% METRICS_STATS — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 4/5 of the metrics step.

fprintf('  --- SECTION 7: Univariate Analysis ---\n');

figure_titles = {
    '1. Absolute DWI/IVIM Metrics vs Local Failure (Wilcoxon Rank-Sum)', ...
    '2. Percent Change Metrics vs Local Failure (Wilcoxon Rank-Sum)', ...
    '3. Target Coverage (D95): Whole GTV vs Resistant Sub-volumes (Wilcoxon Rank-Sum)', ...
    '4. Target Coverage (V50): Whole GTV vs Resistant Sub-volumes (Wilcoxon Rank-Sum)'
};

p_val_store = struct('p_vals', {});

for s = 1:length(metric_sets)
    current_metrics = metric_sets{s};
    current_names = set_names{s};
    
    fig = figure('Name', [figure_titles{s} ' — ' dtype_label], 'Position', [50, 50, 1600, 1000]);
    sgtitle([figure_titles{s} ' (' dtype_label ')'], 'FontSize', 16, 'FontWeight', 'bold');
    
    num_rows = length(current_metrics);
    max_cols = 0;
    for m=1:num_rows
        max_cols = max(max_cols, size(current_metrics{m}, 2));
    end
    cols_to_plot = min(max_cols, length(time_labels));
    
    p_val_store(s).p_vals = nan(num_rows, cols_to_plot);

    plot_idx = 1;
    
    for m = 1:num_rows
        metric_data = current_metrics{m};
        
        for tp = 1:cols_to_plot
            subplot(num_rows, cols_to_plot, plot_idx);
            
            if tp > size(metric_data, 2)
                 axis off;
                 plot_idx = plot_idx + 1;
                 continue;
            end

            y_raw = metric_data(valid_pts, tp);
            has_data = ~isnan(y_raw);
            y = y_raw(has_data);
            g = lf_group(has_data);
            
            p = perform_statistical_test(y, g, 'ranksum');

            if ~isnan(p)
                p_val_store(s).p_vals(m, tp) = p;

                boxplot(y, g, 'Labels', {'LC (0)', 'LF (1)'});
                title_str = sprintf('%s - %s\np = %.3f', current_names{m}, time_labels{tp}, p);
                if p < 0.05
                    title(title_str, 'Color', 'r', 'FontWeight', 'bold');
                else
                    title(title_str, 'Color', 'k', 'FontWeight', 'normal');
                end
            else
                title(sprintf('%s - %s\n(Insufficient Data)', current_names{m}, time_labels{tp}), 'FontSize', 8);
                axis off; 
            end
            
            if m == num_rows
                xlabel('Outcome');
            end
            if tp == 1
                ylabel(current_names{m});
            end
            
            grid on;
            plot_idx = plot_idx + 1;
        end
    end
    subplot_scale = 0.92;
    allAx = findall(fig, 'Type', 'Axes');
    for k = 1:numel(allAx)
        pos = allAx(k).Position;
        allAx(k).Position = [pos(1), pos(2) * subplot_scale, pos(3), pos(4) * subplot_scale];
    end
    saveas(gcf, fullfile(output_folder, sprintf('Metric_Set_%d_%s.png', s, dtype_label)));
    close(gcf);
end

fprintf('  --- SECTION 8: Compile and Export Significant Results ---\n');
total_checks = 0;
for s = 1:length(metric_sets)
    total_checks = total_checks + length(metric_sets{s}) * length(time_labels);
end

sig_metric = cell(total_checks, 1);
sig_fraction = cell(total_checks, 1);
sig_pval = nan(total_checks, 1);
sig_mean_LC = nan(total_checks, 1);
sig_mean_LF = nan(total_checks, 1);
sig_count = 0;

for s = 1:length(metric_sets)
    current_metrics = metric_sets{s};
    current_names = set_names{s};
    num_metrics = length(current_metrics);
    
    for m = 1:num_metrics
        metric_data = current_metrics{m};
        cols_to_plot = min(size(metric_data, 2), length(time_labels));
        
        for tp = 1:cols_to_plot
            if tp <= size(p_val_store(s).p_vals, 2)
                p = p_val_store(s).p_vals(m, tp);
            else
                p = nan;
            end
            
            if ~isnan(p) && p < 0.05
                y_raw = metric_data(valid_pts, tp);
                has_data = ~isnan(y_raw);
                y = y_raw(has_data);
                g = lf_group(has_data);

                mean_LC = mean(y(g == 0), 'omitnan');
                mean_LF = mean(y(g == 1), 'omitnan');
                
                sig_count = sig_count + 1;
                sig_metric{sig_count, 1} = current_names{m};
                sig_fraction{sig_count, 1} = time_labels{tp};
                sig_pval(sig_count, 1) = p;
                sig_mean_LC(sig_count, 1) = mean_LC;
                sig_mean_LF(sig_count, 1) = mean_LF;
            end
        end
    end
end

if sig_count < total_checks
    sig_metric = sig_metric(1:sig_count, :);
    sig_fraction = sig_fraction(1:sig_count, :);
    sig_pval = sig_pval(1:sig_count, :);
    sig_mean_LC = sig_mean_LC(1:sig_count, :);
    sig_mean_LF = sig_mean_LF(1:sig_count, :);
end

if ~isempty(sig_pval)
    sig_results_table = table(sig_metric, sig_fraction, sig_pval, sig_mean_LC, sig_mean_LF, ...
        'VariableNames', {'Metric', 'Timepoint', 'P_Value', 'Mean_LC', 'Mean_LF'});
    sig_results_table = sortrows(sig_results_table, 'P_Value');
    
    disp('----- Significant Findings (p < 0.05) -----');
    disp(sig_results_table);
    
    export_filename = fullfile(dataloc, 'Significant_LF_Metrics.csv');
    writetable(sig_results_table, export_filename);
    fprintf('Saved significant results to: %s\n', export_filename);
else
    disp('No significant differences (p < 0.05) found between LC and LF groups for these metrics.');
end

fprintf('  --- SECTION 9: FDR Correction ---\n');
if ~isempty(sig_pval)
    fprintf('\n----- PER-TIMEPOINT FDR (Benjamini-Hochberg, Q < 0.05) -----\n');
    for tp = 1:length(time_labels)
        tp_pvals  = [];
        tp_labels = {};
        for s = 1:length(metric_sets)
            current_metrics = metric_sets{s};
            current_names = set_names{s};
            for mi = 1:length(current_metrics)
                if tp <= size(p_val_store(s).p_vals, 2)
                    p = p_val_store(s).p_vals(mi, tp);
                else
                    p = nan;
                end

                if ~isnan(p)
                    tp_pvals(end+1, 1)  = p;
                    tp_labels{end+1, 1} = current_names{mi};
                end
            end
        end
        
        if isempty(tp_pvals), continue; end
        
        n_tp = length(tp_pvals);
        [p_sort, sort_id] = sort(tp_pvals);
        q_tp = zeros(n_tp, 1);
        q_tp(n_tp) = p_sort(n_tp);
        for ii = n_tp-1:-1:1
            q_tp(ii) = min(q_tp(ii+1), p_sort(ii) * (n_tp / ii));
        end
        q_tp = min(q_tp, 1);
        q_unsorted = zeros(n_tp, 1);
        q_unsorted(sort_id) = q_tp;
        
        tp_table = table(tp_labels, tp_pvals, q_unsorted, ...
            'VariableNames', {'Metric', 'Raw_P', 'FDR_Q'});
        sig_tp = tp_table(tp_table.FDR_Q < 0.05, :);
        
        fprintf('\n  Timepoint: %s (family size = %d)\n', time_labels{tp}, n_tp);
        if isempty(sig_tp)
            fprintf('    None survived FDR correction.\n');
        else
            disp(sig_tp);
            writetable(sig_tp, fullfile(dataloc, sprintf('FDR_Sig_%s.csv', strrep(time_labels{tp},' ','_'))));
        end
    end
end

fprintf('  --- SECTION 10: Per-Timepoint Analysis Loop ---\n');

% Initialize risk outputs to be returned and used by survival
risk_scores_all = [];
is_high_risk = [];
times_km = [];
events_km = [];

for target_fx = 2:nTp
    fx_label = x_labels{target_fx};
    fprintf('\n=== Analyzing %s ===\n', fx_label);

    %% --- Elastic Net Feature Selection ---
    X_lasso_all = [ADC_abs(valid_pts, target_fx), D_abs(valid_pts, target_fx), ...
                   f_abs(valid_pts, target_fx),   Dstar_abs(valid_pts, target_fx), ...
                   ADC_pct(valid_pts, target_fx), D_pct(valid_pts, target_fx), ...
                   f_pct(valid_pts, target_fx),   Dstar_pct(valid_pts, target_fx), ...
                   m_d95_gtvp(valid_pts, target_fx), m_v50gy_gtvp(valid_pts, target_fx), ...
                   d95_adc_sub(valid_pts, target_fx), v50_adc_sub(valid_pts, target_fx), ...
                   d95_d_sub(valid_pts, target_fx),   v50_d_sub(valid_pts, target_fx), ...
                   d95_f_sub(valid_pts, target_fx),   v50_f_sub(valid_pts, target_fx), ...
                   d95_dstar_sub(valid_pts, target_fx), v50_dstar_sub(valid_pts, target_fx)];
               
    feat_names_lasso = {'ADC_Abs', 'D_Abs', 'f_Abs', 'Dstar_Abs', ...
                        'ADC_Pct', 'D_Pct', 'f_Pct', 'Dstar_Pct', ...
                        'D95_GTVp', 'V50_GTVp', ...
                        'D95_Sub_ADC', 'V50_Sub_ADC', ...
                        'D95_Sub_D', 'V50_Sub_D', ...
                        'D95_Sub_f', 'V50_Sub_f', ...
                        'D95_Sub_Dstar', 'V50_Sub_Dstar'};

    original_feature_indices = 1:18;

    if target_fx == 6
        X_lasso_all = X_lasso_all(:, 1:8);
        feat_names_lasso = feat_names_lasso(1:8);
        original_feature_indices = original_feature_indices(1:8);
    end
    
    valid_cols = ~all(isnan(X_lasso_all), 1);
    X_lasso_all = X_lasso_all(:, valid_cols);
    feat_names_lasso = feat_names_lasso(valid_cols);
    original_feature_indices = original_feature_indices(valid_cols);

    y_lasso_all = lf_group;

    base_cols = min(8, size(X_lasso_all, 2));
    has_any_imaging = any(~isnan(X_lasso_all(:, 1:base_cols)), 2);
    impute_mask = has_any_imaging & ~isnan(y_lasso_all);
    X_impute = X_lasso_all(impute_mask, :);
    y_clean  = y_lasso_all(impute_mask);

    id_list_valid = id_list(valid_pts);
    id_list_impute = id_list_valid(impute_mask);

    rng(42); 
    fold_id_en = make_grouped_folds(id_list_impute, y_clean, 5);
    n_folds_en = max(fold_id_en);
    n_lambdas = 10;
    
    common_Lambda = [];
    cv_failed = false;
    
    warning('off', 'all');
    for cv_i = 1:n_folds_en
        tr_idx = (fold_id_en ~= cv_i);
        te_idx = (fold_id_en == cv_i);
        X_tr = X_impute(tr_idx, :); y_tr = y_clean(tr_idx);
        X_te = X_impute(te_idx, :); y_te = y_clean(te_idx);
        
        id_tr = id_list_impute(tr_idx);
        id_te = id_list_impute(te_idx);
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, 5, id_tr, id_te);
        
        keep_fold = filter_collinear_features(X_tr_imp, y_tr);
        X_tr_kept = X_tr_imp(:, keep_fold);
        X_te_kept = X_te_imp(:, keep_fold);
        
        try
            if cv_i == 1
                [B_fold, FitInfo_fold] = lassoglm(X_tr_kept, y_tr, 'binomial', ...
                    'Alpha', 0.5, 'NumLambda', n_lambdas, 'Standardize', true, 'MaxIter', 10000);
                common_Lambda = FitInfo_fold.Lambda;
                all_deviance = zeros(length(common_Lambda), n_folds_en);
            else
                [B_fold, FitInfo_fold] = lassoglm(X_tr_kept, y_tr, 'binomial', ...
                    'Alpha', 0.5, 'Lambda', common_Lambda, 'Standardize', true, 'MaxIter', 10000);
            end
            
            eta = X_te_kept * B_fold + FitInfo_fold.Intercept;
            p = 1 ./ (1 + exp(-eta));
            p = max(min(p, 1 - 1e-10), 1e-10); 
            y_te_mat = repmat(y_te, 1, length(common_Lambda));
            
            dev_val = -2 * sum(y_te_mat .* log(p) + (1 - y_te_mat) .* log(1 - p), 1);
            all_deviance(:, cv_i) = dev_val';
        catch
            cv_failed = true;
            break;
        end
    end
    warning('on', 'all');
    
    if ~cv_failed && ~isempty(common_Lambda)
        mean_deviance = mean(all_deviance, 2);
        [~, idx_min] = min(mean_deviance);
        opt_lambda = common_Lambda(idx_min);
        
        X_clean_all = knn_impute_train_test(X_impute, [], 5, id_list_impute);
        
        try
            [B_final, FitInfo_final] = lassoglm(X_clean_all, y_clean, 'binomial', ...
                'Alpha', 0.5, 'Lambda', opt_lambda, 'Standardize', true, 'MaxIter', 10000);
            
            coefs_en = B_final;
            
            selected_indices = find(coefs_en ~= 0);
            selected_indices = original_feature_indices(selected_indices);
            
            fprintf('Elastic Net Selected Features for %s (Opt Lambda=%.4f): %s\n', ...
                fx_label, opt_lambda, strjoin(feat_names_lasso(selected_indices), ', '));
        catch
            cv_failed = true;
        end
    end
    
    if cv_failed || isempty(common_Lambda)
        fprintf('Elastic Net failed to converge for %s. Fallback to empty selection.\n', fx_label);
        selected_indices = [];
    end

    %% --- LOOCV For Risk Scores ---
    n_pts_impute = size(X_impute, 1);
    risk_scores_oof = nan(n_pts_impute, 1);
    is_high_risk_oof = false(n_pts_impute, 1);

    fprintf('  Generating unbiased out-of-fold risk scores via nested LOOCV...\n');
    for loo_i = 1:n_pts_impute
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
        
        keep_fold = filter_collinear_features(X_tr_imp, y_tr_fold);
        X_tr_kept = X_tr_imp(:, keep_fold);
        X_te_kept = X_te_imp(:, keep_fold);
        
        warning('off', 'all');
        try
            inn_fold_id = make_grouped_folds(id_tr_loo, y_tr_fold, 5);
            n_inn_folds = max(inn_fold_id);
            inn_Lambda = []; inn_deviance = [];
            for inn_i = 1:n_inn_folds
                inn_tr = (inn_fold_id ~= inn_i);
                inn_te = (inn_fold_id == inn_i);
                if inn_i == 1
                    [B_inn, FI_inn] = lassoglm(X_tr_kept(inn_tr,:), y_tr_fold(inn_tr), 'binomial', ...
                        'Alpha', 0.5, 'NumLambda', 10, 'Standardize', true, 'MaxIter', 10000);
                    inn_Lambda = FI_inn.Lambda;
                    inn_deviance = zeros(numel(inn_Lambda), n_inn_folds);
                else
                    [B_inn, FI_inn] = lassoglm(X_tr_kept(inn_tr,:), y_tr_fold(inn_tr), 'binomial', ...
                        'Alpha', 0.5, 'Lambda', inn_Lambda, 'Standardize', true, 'MaxIter', 10000);
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
                'Alpha', 0.5, 'Lambda', opt_lam_loo, 'Standardize', true, 'MaxIter', 10000);
            coefs_loo = B_loo;
            intercept_loo = FitInfo_loo.Intercept;
        catch
            coefs_loo = zeros(size(X_tr_kept, 2), 1);
            intercept_loo = 0;
        end
        warning('on', 'all');
        
        risk_scores_oof(loo_i) = X_te_kept * coefs_loo + intercept_loo;
        
        train_median = median(X_tr_kept * coefs_loo + intercept_loo);
        is_high_risk_oof(loo_i) = risk_scores_oof(loo_i) > train_median;
    end
    
    risk_scores_all_target = nan(sum(valid_pts), 1);
    risk_scores_all_target(impute_mask) = risk_scores_oof;
    
    is_high_risk_target = false(sum(valid_pts), 1);
    is_high_risk_target(impute_mask) = is_high_risk_oof;

    all_feat_data  = {ADC_abs,       D_abs,       f_abs,       Dstar_abs, ...
                      ADC_pct,       D_pct,       f_pct,       Dstar_pct, ...
                      m_d95_gtvp,    m_v50gy_gtvp, ...
                      d95_adc_sub,   v50_adc_sub, ...
                      d95_d_sub,     v50_d_sub, ...
                      d95_f_sub,     v50_f_sub, ...
                      d95_dstar_sub, v50_dstar_sub};
    
    all_feat_names = {'ADC',         'D',         'f',         'D*', ...
                      'ADC',         'D',         'f',         'D*', ...
                      'D95 GTVp',    'V50 GTVp', ...
                      'D95 Sub(ADC)','V50 Sub(ADC)', ...
                      'D95 Sub(D)',  'V50 Sub(D)', ...
                      'D95 Sub(f)',  'V50 Sub(f)', ...
                      'D95 Sub(D*)', 'V50 Sub(D*)'};
    
    all_feat_is_abs = [true          true         true         true  ...
                       false         false        false        false ...
                       true          false        true         false ...
                       true          false        true         false ...
                       true          false];
    
    all_feat_disp  = {'Abs ADC',     'Abs D',     'Abs f',     'Abs D*', ...
                      '\Delta ADC',  '\Delta D',  '\Delta f',  '\Delta D*', ...
                      'D95 GTVp',    'V50 GTVp', ...
                      'D95 Sub(ADC)','V50 Sub(ADC)', ...
                      'D95 Sub(D)',  'V50 Sub(D)', ...
                      'D95 Sub(f)',  'V50 Sub(f)', ...
                      'D95 Sub(D*)', 'V50 Sub(D*)'};
    
    all_feat_units = {'mm^2/s',      'mm^2/s',    'Fraction',  'mm^2/s', ...
                      '%',           '%',         '%',         '%', ...
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
    
    for si = 1:n_sig
        fi = selected_indices(si);              
        sig_data_selected{si} = all_feat_data{fi};
        sig_names{si}         = all_feat_names{fi};
        sig_is_abs(si)        = all_feat_is_abs(fi);
        sig_is_pct_imaging(si) = (fi >= 5 && fi <= 8);
        sig_disp_names{si}    = all_feat_disp{fi};
        sig_units{si}         = all_feat_units{fi};
        
        if fi <= 4
            sig_abs_data{si} = all_feat_data{fi};
            sig_pct_data{si} = all_feat_data{fi + 4};
        elseif fi >= 5 && fi <= 8
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
    times_km(m_lf==0) = m_total_follow_up_time(m_lf==0);
    events_km = m_lf;
    
    times_km = times_km(valid_pts);
    events_km = events_km(valid_pts);
    
    labels = lf_group; % 0 = LC, 1 = LF
    valid_roc = ~isnan(risk_scores_all_target) & ~isnan(labels);
    
    if sum(valid_roc) > 0
        [roc_X, roc_Y, roc_T, roc_AUC] = perfcurve(labels(valid_roc), risk_scores_all_target(valid_roc), 1);
        
        [~, roc_opt_idx] = max(roc_Y - roc_X);
        roc_opt_thresh = roc_T(roc_opt_idx);
        
        figure('Name', ['ROC Analysis - ' fx_label ' — ' dtype_label], 'Position', [200, 200, 700, 600]);
        hold on;
        
        plot(roc_X, roc_Y, 'b-', 'LineWidth', 2.5);
        leg_entries = {sprintf('LOOCV OOF Risk Score (AUC = %.3f)', roc_AUC)};
        
        plot(roc_X(roc_opt_idx), roc_Y(roc_opt_idx), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        leg_entries{end+1} = sprintf('Youden Cutoff (score = %.3f)', roc_opt_thresh);
        
        plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);
        leg_entries{end+1} = 'Random Guess';
        
        xlabel('False Positive Rate (1 - Specificity)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('True Positive Rate (Sensitivity)', 'FontSize', 12, 'FontWeight', 'bold');
        title(['PRIMARY ROC Curve: LOOCV Out-of-Fold Risk Score (' fx_label ', ' dtype_label ')'], 'FontSize', 14);
        
        legend(leg_entries, 'Location', 'SouthEast', 'FontSize', 11);
        grid on; box on;
        hold off;
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
        
        boxplot(vol_pct, lf_group, 'Labels', {'LC (0)', 'LF (1)'});
        ylabel(['% Change in GTV Volume (' fx_label ')']);
        title('Confounder Check: Volume', 'FontSize', 12, 'FontWeight', 'bold');
        p_vol = perform_statistical_test(vol_pct, lf_group, 'ranksum');
        
        y_lim = ylim;
        if numel(y_lim) >= 2 && all(isfinite(y_lim)) && y_lim(2) > y_lim(1)
            text(1.5, y_lim(1) + 0.9*(y_lim(2)-y_lim(1)), sprintf('p = %.3f', p_vol), ...
                'HorizontalAlignment', 'center', 'FontSize', 11);
        end
        grid on;
        if p_vol > 0.05, xlabel('Conclusion: No Volumetric Bias');
        else, xlabel('Warning: Volume is a Confounder'); end
        
        subplot(1, 3, 2);
        sd_fx1  = adc_sd(valid_pts, 1, 1);
        sd_fxN  = adc_sd(valid_pts, target_fx, 1);
        sd_delta = sd_fxN - sd_fx1;
        
        boxplot(sd_delta, lf_group, 'Labels', {'LC (0)', 'LF (1)'});
        ylabel(['\Delta ADC SD (' fx_label ') [mm\u00b2/s]']);
        title('Heterogeneity: ADC SD Change', 'FontSize', 12, 'FontWeight', 'bold');
        p_sd = perform_statistical_test(sd_delta, lf_group, 'ranksum');
        
        y_lim = ylim;
        if numel(y_lim) >= 2 && all(isfinite(y_lim)) && y_lim(2) > y_lim(1)
            text(1.5, y_lim(1) + 0.9*(y_lim(2)-y_lim(1)), sprintf('p = %.3f', p_sd), ...
                'HorizontalAlignment', 'center', 'FontSize', 11);
        end
        grid on;
        
        subplot(1, 3, 3);
        hold on;
        wcv_est = 2.8; 
        cor_est = 1.96 * sqrt(2) * wcv_est;
        
        x_scatter = ones(size(lf_group));
        x_scatter(lf_group==1) = 2;
        x_scatter = x_scatter + (rand(size(x_scatter))-0.5)*0.2;
        scatter(x_scatter, curr_sig_pct_full(valid_pts, target_fx), 50, 'filled', 'MarkerEdgeColor', 'k');
        
        base_idx = mod(selected_indices(vi)-1, 4) + 1;
        
        if sig_is_pct_imaging(vi) && base_idx == 1
            yfill = [-cor_est cor_est cor_est -cor_est];
            xfill = [0.5 0.5 2.5 2.5];
            fill(xfill, yfill, [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            yline(0, 'k-');
            yline(cor_est, 'k--', 'CoR (+7.8%)');
            yline(-cor_est, 'k--', 'CoR (-7.8%)');
        elseif ~sig_is_abs(vi) && contains(sig_units{vi}, '%')
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
            pos = allAx(k).Position;
            allAx(k).Position = [pos(1), pos(2) * 0.92, pos(3), pos(4) * 0.92];
        end
        safe_name = strrep(curr_sig_file, '*', 'star');
        saveas(gcf, fullfile(output_folder, ['Sanity_Checks_' safe_name '_' fx_label '_' dtype_label '.png']));
        close(gcf);
    end

    if n_sig >= 2
        for fi = 1:(n_sig-1)
            for fj = (fi+1):n_sig
                figure('Name', sprintf('2D Feature Space %s vs %s %s — %s', sig_names{fi}, sig_names{fj}, fx_label, dtype_label), 'Position', [100, 100, 800, 600]);
                hold on;
                
                x_val = sig_data_selected{fi}(valid_pts, target_fx);
                y_val = sig_data_selected{fj}(valid_pts, target_fx);
                group = lf_group;

                scatter(x_val(group==0), y_val(group==0), 80, 'b', 'filled', 'MarkerEdgeColor', 'k');
                scatter(x_val(group==1), y_val(group==1), 80, 'r', 'filled', 'MarkerEdgeColor', 'k');

                mdl = fitglm([x_val, y_val], group, 'Distribution', 'binomial', 'Options', statset('MaxIter', 10000));
                coefs = mdl.Coefficients.Estimate;
                x_range = linspace(min(x_val), max(x_val), 100);
                y_boundary = -(coefs(1) + coefs(2)*x_range) / coefs(3);
                plot(x_range, y_boundary, 'k--', 'LineWidth', 2);

                if sig_is_abs(fi), xl = sprintf('%s at %s (%s)', sig_names{fi}, fx_label, sig_units{fi}); else, xl = sprintf('\\Delta %s at %s (%s)', sig_names{fi}, fx_label, sig_units{fi}); end
                if sig_is_abs(fj), yl = sprintf('%s at %s (%s)', sig_names{fj}, fx_label, sig_units{fj}); else, yl = sprintf('\\Delta %s at %s (%s)', sig_names{fj}, fx_label, sig_units{fj}); end
                xlabel(xl, 'FontSize', 12, 'FontWeight', 'bold');
                ylabel(yl, 'FontSize', 12, 'FontWeight', 'bold');
                title(sprintf('Biomarker Interaction: Separation of LC vs LF (%s, %s)', fx_label, dtype_label), 'FontSize', 14);
                legend({'Local Control', 'Local Failure', 'Logistic Decision Boundary'}, 'Location', 'NorthWest');
                
                grid on; box on;
                xline(0, 'k-', 'Alpha', 0.2); yline(0, 'k-', 'Alpha', 0.2);
                
                safe_name1 = strrep(sig_names{fi}, '*', 'star');
                safe_name2 = strrep(sig_names{fj}, '*', 'star');
                saveas(gcf, fullfile(output_folder, sprintf('2D_Space_%s_vs_%s_%s_%s.png', safe_name1, safe_name2, fx_label, dtype_label)));
                close(gcf);
            end
        end
    else
        fprintf('Skipping 2D scatter plots: requires at least 2 significant variables.\n');
    end

    % Track values to pass to metrics_survival
    risk_scores_all = risk_scores_all_target;
    is_high_risk = is_high_risk_target;
end

% Let's also run the GLME Mixed-effects model from the end of the script
fprintf('\n--- LONGITUDINAL MIXED-EFFECTS MODEL (GLME) ---\n');
long_PatientID = [];
long_Timepoint = [];
long_ADC = [];
long_D = [];
long_f = [];
long_Dstar = [];
long_LF = [];

patient_indices = find(valid_pts);
for i = 1:length(patient_indices)
    p_idx = patient_indices(i);
    for t = 1:nTp
        if ~isnan(ADC_abs(p_idx, t)) || ~isnan(D_abs(p_idx, t)) || ~isnan(f_abs(p_idx, t)) || ~isnan(Dstar_abs(p_idx, t))
            long_PatientID = [long_PatientID; i]; 
            long_Timepoint = [long_Timepoint; t];
            long_ADC = [long_ADC; ADC_abs(p_idx, t)];
            long_D = [long_D; D_abs(p_idx, t)];
            long_f = [long_f; f_abs(p_idx, t)];
            long_Dstar = [long_Dstar; Dstar_abs(p_idx, t)];
            long_LF = [long_LF; lf_group(i)];
        end
    end
end

glme_table = table(categorical(long_PatientID), long_Timepoint, ...
    long_ADC, long_D, long_f, long_Dstar, long_LF, ...
    'VariableNames', {'PatientID', 'Timepoint', 'ADC', 'D', 'f', 'Dstar', 'LF'});

clean_idx = ~isnan(glme_table.ADC) & ~isnan(glme_table.D) & ~isnan(glme_table.f) & ~isnan(glme_table.Dstar);
glme_table_clean = glme_table(clean_idx, :);
baseline_idx = glme_table_clean.Timepoint == 1;

mean_ADC_base = mean(glme_table_clean.ADC(baseline_idx));
std_ADC_base = std(glme_table_clean.ADC(baseline_idx));
mean_D_base = mean(glme_table_clean.D(baseline_idx));
std_D_base = std(glme_table_clean.D(baseline_idx));
mean_f_base = mean(glme_table_clean.f(baseline_idx));
std_f_base = std(glme_table_clean.f(baseline_idx));
mean_Dstar_base = mean(glme_table_clean.Dstar(baseline_idx));
std_Dstar_base = std(glme_table_clean.Dstar(baseline_idx));

glme_table_clean.ADC_z = (glme_table_clean.ADC - mean_ADC_base) / std_ADC_base;
glme_table_clean.D_z = (glme_table_clean.D - mean_D_base) / std_D_base;
glme_table_clean.f_z = (glme_table_clean.f - mean_f_base) / std_f_base;
glme_table_clean.Dstar_z = (glme_table_clean.Dstar - mean_Dstar_base) / std_Dstar_base;

glme_table_clean.LF = categorical(glme_table_clean.LF);
glme_table_clean.Timepoint = categorical(glme_table_clean.Timepoint);

biomarkers = {'ADC_z', 'D_z', 'f_z', 'Dstar_z'};
warning('off', 'all');
for b = 1:length(biomarkers)
    bm = biomarkers{b};
    formula = sprintf('%s ~ 1 + LF * Timepoint + (1|PatientID)', bm);
    try
        glme = fitglme(glme_table_clean, formula, 'OptimizerOptions', statset('MaxIter', 10000));
        fprintf('\n--- %s ---\n', formula);
        
        anova_res = anova(glme);
        row_idx = find(strcmp(anova_res.Term, 'LF:Timepoint'));
        if ~isempty(row_idx)
            pval = anova_res.pValue(row_idx);
            fprintf('Interaction P-Value (LF * Timepoint): %.4f\n', pval);
            if pval < 0.05
                fprintf('  -> SIGNIFICANT DIFFERENCE in trajectory between LC and LF groups.\n');
            else
                fprintf('  -> No significant difference in trajectory between LC and LF groups.\n');
            end
        else
            disp(anova_res);
        end
    catch ME
        fprintf('GLME model for %s failed to converge: %s\n', bm, ME.message);
    end
end
warning('on', 'all');

end
