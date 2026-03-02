function [risk_scores_all, is_high_risk, times_km, events_km] = metrics_stats_predictive(valid_pts, lf_group, dtype_label, output_folder, dataloc, nTp, m_gtv_vol, adc_sd, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_pct, Dstar_pct, m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub, id_list, dtype, dl_provenance, x_labels, m_lf, m_total_time, m_total_follow_up_time)
% METRICS_STATS_PREDICTIVE — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 4b/5 of the metrics step: Predictive Modeling (Elastic Net & Cox prep).
% Performs multivariate feature selection via Elastic Net logistic regression
% and generates out-of-fold risk scores via nested Leave-One-Out Cross-Validation.
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
    
    if iscell(X_lasso_all)
        vars = {'ADC_abs', 'D_abs', 'f_abs', 'Dstar_abs', 'ADC_pct', 'D_pct', 'f_pct', 'Dstar_pct', 'm_d95_gtvp', 'm_v50gy_gtvp', 'd95_adc_sub', 'v50_adc_sub', 'd95_d_sub', 'v50_d_sub', 'd95_f_sub', 'v50_f_sub', 'd95_dstar_sub', 'v50_dstar_sub'};
        fp = fopen('debug_concat_error.txt', 'a');
        fprintf(fp, '\n--- X_lasso_all is a cell array at fx_label=%s ---\n', fx_label);
        for i_v = 1:length(vars)
            tmp_v = eval(vars{i_v});
            fprintf(fp, '%s -> Size: %s, Class: %s\n', vars{i_v}, mat2str(size(tmp_v)), class(tmp_v));
        end
        fclose(fp);
        error('Invalid data type: cell array detected in X_lasso_all concatenation. See debug_concat_error.txt');
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
        ylabel(['\Delta ADC SD (' fx_label ') [mm^2/s]']);
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

end
