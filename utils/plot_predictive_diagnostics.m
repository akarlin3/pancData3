function plot_predictive_diagnostics( ...
    selected_indices, n_sig, sig_data_selected, sig_names, sig_is_abs, ...
    sig_is_pct_imaging, sig_disp_names, sig_units, sig_col_idx, ...
    sig_abs_data, sig_pct_data, ...
    risk_scores_all_target, lf_group, valid_pts, ...
    m_gtv_vol, adc_sd, ADC_abs, ...
    target_fx, fx_label, dtype_label, dtype, output_folder, use_firth)
% PLOT_PREDICTIVE_DIAGNOSTICS  ROC curve, sanity checks, and 2D scatter plots.
%
%   Generates all diagnostic figures for the predictive modeling section:
%     1. ROC curve from LOOCV out-of-fold risk scores with AUC and Youden cutoff
%     2. Per-feature sanity check panels (volume confounder, ADC SD heterogeneity,
%        signal vs noise floor with CoR band)
%     3. Pairwise 2D feature space scatter plots with logistic decision boundary
%
%   These plots verify that selected biomarkers are not confounded and that
%   multivariate combinations provide meaningful separation between LC and LF.
%
% Inputs:
%   selected_indices       - Indices of selected features in the original 22
%   n_sig                  - Number of significant features
%   sig_data_selected      - Cell array of data matrices for selected features
%   sig_names              - Cell array of feature names
%   sig_is_abs             - Logical: true if feature is absolute (not change)
%   sig_is_pct_imaging     - Logical: true if feature is percent-change imaging
%   sig_disp_names         - Cell array of display names for features
%   sig_units              - Cell array of unit strings
%   sig_col_idx            - Column index to plot for each feature
%   sig_abs_data           - Cell array of absolute data for each feature
%   sig_pct_data           - Cell array of percent-change data for each feature
%   risk_scores_all_target - Out-of-fold risk scores for valid patients
%   lf_group               - Outcome labels (0=LC, 1=LF, 2=competing)
%   valid_pts              - Logical mask of valid patients
%   m_gtv_vol              - GTV volume array
%   adc_sd                 - ADC standard deviation array (3D: patients x timepoints x dtypes)
%   ADC_abs                - Absolute ADC array
%   target_fx              - Target fraction index
%   fx_label               - Fraction label string
%   dtype_label            - DWI type label string
%   dtype                  - DWI type index
%   output_folder          - Output folder path
%   use_firth              - Whether to use Firth penalty for decision boundary

    %% --- ROC Analysis ---
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

    %% --- 2D Feature Space Scatter Plots ---
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
                scatter_mask_2d = (group <= 1);
                x_val = x_val(scatter_mask_2d);
                y_val = y_val(scatter_mask_2d);
                group = group(scatter_mask_2d);

                scatter(x_val(group==0), y_val(group==0), 80, [0 0.4470 0.7410], 'filled', 'MarkerEdgeColor', 'k');
                scatter(x_val(group==1), y_val(group==1), 80, [0.8500 0.3250 0.0980], 'filled', 'MarkerEdgeColor', 'k');

                % NOTE: Decision boundary is fitted on the full displayed
                % dataset (not cross-validated) for visualization only.
                % Firth penalty produces stable boundaries under separation.
                w_state = warning('off', 'all');
                if use_firth
                    mdl = fitglm([x_val, y_val], group, 'Distribution', 'binomial', ...
                        'LikelihoodPenalty', 'jeffreys-prior', 'Options', statset('MaxIter', 1e7));
                else
                    mdl = fitglm([x_val, y_val], group, 'Distribution', 'binomial', 'Options', statset('MaxIter', 1e7));
                end
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
                if use_firth
                    boundary_label = 'Firth Boundary (illustrative, not CV)';
                else
                    boundary_label = 'Logistic Boundary (illustrative, not CV)';
                end
                if numel(coefs) >= 3 && coefs(3) ~= 0
                    legend({'Local Control', 'Local Failure', boundary_label}, 'Location', 'NorthWest');
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
end
