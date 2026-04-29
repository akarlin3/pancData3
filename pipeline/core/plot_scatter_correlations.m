function plot_scatter_correlations(dtype_label, dmean_gtvp, d95_gtvp, adc_mean, d_mean, f_mean, valid_pts, lf_group, dtype, output_folder)
% PLOT_SCATTER_CORRELATIONS
%
% ANALYTICAL OVERVIEW:
%   Explores the relationship between baseline tumour diffusion properties
%   and the radiation dose delivered to the GTV.  This analysis tests a
%   key hypothesis in adaptive radiotherapy: do tumours with restricted
%   diffusion (low ADC/D = dense cellularity) receive adequate radiation
%   dose, or are there systematic dose-response patterns that could inform
%   treatment planning?
%
%   Two dose endpoints are examined:
%     Mean GTV Dose — Average dose across the entire tumour volume.
%       A surrogate for overall treatment intensity.
%     D95 — Minimum dose to 95% of the GTV volume.
%       Clinically more relevant than mean dose because it reflects the
%       WORST-COVERED region of the tumour.  A low D95 means part of the
%       GTV is underdosed, which is a known risk factor for local failure.
%
%   Spearman (rank) correlation is used instead of Pearson because:
%     1. Dose-diffusion relationships may be non-linear
%     2. Spearman is robust to outliers (important with noisy DWI data)
%     3. Rank correlation captures monotonic trends without assuming
%        a specific functional form
%
%   Correlations and the trend line are computed on the full cohort
%   (LC + LF + CR pooled) so the reported r_s reflects the dose-diffusion
%   relationship across every patient with valid data, not within an
%   outcome subset. Markers remain colour-coded by outcome for visual
%   stratification (blue = LC, red = LF, gray = CR).

% Extract Fx1 (baseline) dose metrics for the valid patient subset.
% Baseline dose metrics reflect the planned dose distribution before any
% adaptive replanning.  In standard practice, the RT plan is designed
% once and delivered across all fractions, so Fx1 dose approximates the
% cumulative planned dose.
dose_mean_vec = dmean_gtvp(valid_pts, 1);   % mean dose inside GTV (Gy)
dose_d95_vec  = d95_gtvp(valid_pts, 1);     % D95 = dose to 95 % of GTV (Gy)

% Diffusion biomarkers to correlate with dose.
% D* is excluded from this analysis because its high physiological noise
% would produce meaningless correlations with dose metrics.
% ADC and D reflect cellularity (lower = denser tumour = potentially
% needs higher dose for local control).
% f reflects perfusion/oxygenation (lower = more hypoxic = more
% radioresistant, may need dose escalation).
diff_metrics = {adc_mean(valid_pts,1,dtype), d_mean(valid_pts,1,dtype), f_mean(valid_pts,1,dtype)};
diff_names   = {'Mean ADC', 'Mean D', 'Mean f'};
diff_units   = {'mm^2/s',   'mm^2/s',  ''};

% 2 rows (mean dose, D95) x 3 cols (ADC, D, f) = 6 scatter subplots
fig_scatter = figure('Visible', 'off', 'Name', ['Dose vs Diffusion Metrics — ' dtype_label], ...
       'Position', [150, 150, 1400, 500]);

plot_idx = 1;  % sequential subplot index across the 2x3 grid
n_diff_metrics = numel(diff_metrics);
for di = 1:n_diff_metrics
    text_progress_bar(di, n_diff_metrics, 'Generating scatter plots');
    y_vals = diff_metrics{di};

    % Loop over the two dose endpoints (mean dose, D95)
    for dose_type = 1:2
        if dose_type == 1
            x_vals = dose_mean_vec;
            x_label = 'Mean GTV Dose (Gy)';
        else
            x_vals = dose_d95_vec;
            x_label = 'D95 (Gy)';
        end

        subplot(2, numel(diff_metrics), plot_idx);

        % Ensure lf_group is a column vector to match clean, x_vals, and y_vals.
        % lf_group encoding: 0=local control (LC), 1=local failure (LF),
        % 2=competing risk (non-cancer death without prior LF).
        lf_group_col = lf_group(:);

        % Full cohort = every patient with a known outcome and valid x/y.
        % CR patients are full cohort members; they are kept in both the
        % pooled correlation and the trend line.
        eligible_all = ~isnan(lf_group_col);
        clean_all = eligible_all & ~isnan(x_vals) & ~isnan(y_vals);
        if sum(clean_all) < 3
            title([diff_names{di} ' — insufficient data']);
            plot_idx = plot_idx + 1;
            continue;
        end

        % Plot competing-risk patients (lf==2) first as gray squares so
        % they appear behind the LC/LF markers.
        cr_mask = clean_all & lf_group_col==2;
        if sum(cr_mask) > 0
            scatter(x_vals(cr_mask), y_vals(cr_mask), ...
                40, [0.5 0.5 0.5], 's', 'filled', 'MarkerEdgeColor', [0.3 0.3 0.3], ...
                'DisplayName', 'CR');
            hold on;
        end

        % Plot LC (blue, MATLAB default blue) and LF (red/orange, MATLAB default orange)
        % with black edge for visibility. Marker size 50 provides good contrast
        % for cohorts of ~20-40 patients without excessive overlap.
        lc_mask = clean_all & lf_group_col==0;
        lf_mask = clean_all & lf_group_col==1;
        scatter(x_vals(lc_mask), y_vals(lc_mask), ...
            50, [0 0.4470 0.7410], 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'LC'); hold on;
        scatter(x_vals(lf_mask), y_vals(lf_mask), ...
            50, [0.8500 0.3250 0.0980], 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'LF');

        % Pooled (full-cohort) linear trend line over LC + LF + CR.
        warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
        if sum(clean_all) >= 2
            x_line = linspace(min(x_vals(clean_all)), max(x_vals(clean_all)), 50);
            p_fit = polyfit(x_vals(clean_all), y_vals(clean_all), 1);
            plot(x_line, polyval(p_fit, x_line), '-', 'Color', [0.2 0.2 0.2], ...
                'LineWidth', 2, 'DisplayName', 'Cohort trend');
        end
        warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');
        hold off;

        % Full-cohort Spearman correlation across LC + LF + CR. A minimum of
        % 3 points is required for a meaningful estimate (fewer points
        % produce unreliable p-values).
        r_all = NaN; p_all = NaN;
        if sum(clean_all) >= 3
            if exist('OCTAVE_VERSION', 'builtin')
                % Octave uses spearman() and requires manual t-test for p-value.
                % t = r * sqrt((n-2) / (1-r^2)) follows t-distribution with n-2 df.
                % eps prevents division by zero when |r| = 1 (perfect correlation).
                r_all = spearman(x_vals(clean_all), y_vals(clean_all));
                n_all = sum(clean_all);
                t_all = r_all * sqrt((n_all - 2) / (1 - r_all^2 + eps));
                p_all = 2 * (1 - tcdf(abs(t_all), n_all - 2));
            else
                [r_all, p_all] = corr(x_vals(clean_all), y_vals(clean_all), 'Type', 'Spearman');
            end
        end

        xlabel(x_label);
        if isempty(diff_units{di})
            ylabel(diff_names{di});
        else
            ylabel([diff_names{di} ' (' diff_units{di} ')']);
        end
        % Per-outcome counts shown alongside the pooled correlation.
        n_lc = sum(lc_mask);
        n_lf = sum(lf_mask);
        n_cr = sum(cr_mask);
        n_total = sum(clean_all);

        if isnan(r_all)
            cohort_corr_str = sprintf('Cohort (n=%d): n<3', n_total);
        else
            cohort_corr_str = sprintf('Cohort (n=%d) r_s=%.2f %s', n_total, r_all, format_p_value(p_all));
        end
        if n_cr > 0
            breakdown_str = sprintf('LC: n=%d | LF: n=%d | CR: n=%d', n_lc, n_lf, n_cr);
        else
            breakdown_str = sprintf('LC: n=%d | LF: n=%d', n_lc, n_lf);
        end
        title(sprintf('%s vs Dose\n%s | %s', diff_names{di}, cohort_corr_str, breakdown_str), ...
            'FontSize', 9);
        if exist('OCTAVE_VERSION', 'builtin')
            if n_cr > 0
                legend('CR', 'LC', 'LF', 'Cohort trend', 'location', 'best');
            else
                legend('LC', 'LF', 'Cohort trend', 'location', 'best');
            end
        else
            if n_cr > 0
                lg = legend('CR', 'LC', 'LF', 'Location', 'best', 'FontSize', 8);
            else
                lg = legend('LC', 'LF', 'Location', 'best', 'FontSize', 8);
            end
            set(lg, 'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
        end
        grid on;

        plot_idx = plot_idx + 1;
    end
end
sgtitle(['RT Dose vs Diffusion Metrics (Fx1) (' dtype_label ')'], ...
        'FontSize', 14, 'FontWeight', 'bold');
set(findall(fig_scatter, 'Type', 'Axes'), 'Toolbar', []);
print(fig_scatter, fullfile(output_folder, ['Dose_vs_Diffusion_' dtype_label '.png']), '-dpng', '-r300');
close(fig_scatter);

fprintf('  Scatter plots generated (%s).\n', dtype_label);

end
