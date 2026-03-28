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
%   Per-group (LC/LF) trend lines and correlations avoid Simpson's paradox:
%   the pooled correlation across all patients could be misleading if LC
%   and LF groups have different dose distributions (e.g., LF patients
%   received lower doses AND had lower ADC, creating a spurious positive
%   correlation in the pooled data).
%
%  Each scatter panel is coloured by clinical outcome (blue = LC, red = LF),
%  with per-group linear trend-lines and Spearman correlation coefficients.

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
figure('Name', ['Dose vs Diffusion Metrics — ' dtype_label], ...
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

        % Exclude competing-risk patients (lf==2) from scatter and statistics
        % to prevent Simpson's paradox from masking group-specific trends.
        eligible = (lf_group_col <= 1);
        clean = eligible & ~isnan(x_vals) & ~isnan(y_vals);
        if sum(clean) < 3
            title([diff_names{di} ' — insufficient data']);
            plot_idx = plot_idx + 1;
            continue;
        end

        % Plot LC (blue, MATLAB default blue) and LF (red/orange, MATLAB default orange)
        % with black edge for visibility. Marker size 50 provides good contrast
        % for cohorts of ~20-40 patients without excessive overlap.
        scatter(x_vals(clean & lf_group_col==0), y_vals(clean & lf_group_col==0), ...
            50, [0 0.4470 0.7410], 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'LC'); hold on;
        scatter(x_vals(clean & lf_group_col==1), y_vals(clean & lf_group_col==1), ...
            50, [0.8500 0.3250 0.0980], 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'LF');

        % Overlay per-group linear trend lines to avoid pooled Simpson's paradox.
        % If LC and LF groups show opposite trends (e.g., LC: higher dose →
        % higher ADC; LF: higher dose → lower ADC), the pooled trend would
        % be meaningless and could mask the group-specific biology.
        warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
        lc_mask = clean & lf_group_col==0;
        lf_mask = clean & lf_group_col==1;
        if sum(lc_mask) >= 2
            x_line_lc = linspace(min(x_vals(lc_mask)), max(x_vals(lc_mask)), 50);
            p_fit_lc = polyfit(x_vals(lc_mask), y_vals(lc_mask), 1);
            plot(x_line_lc, polyval(p_fit_lc, x_line_lc), '-', 'Color', [0 0.4470 0.7410], ...
                'LineWidth', 2, 'DisplayName', 'LC trend');
        end
        if sum(lf_mask) >= 2
            x_line_lf = linspace(min(x_vals(lf_mask)), max(x_vals(lf_mask)), 50);
            p_fit_lf = polyfit(x_vals(lf_mask), y_vals(lf_mask), 1);
            plot(x_line_lf, polyval(p_fit_lf, x_line_lf), '--', 'Color', [0.8500 0.3250 0.0980], ...
                'LineWidth', 2, 'DisplayName', 'LF trend');
        end
        warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');
        hold off;

        % Compute per-group Spearman correlations to match the per-group
        % trend lines and avoid Simpson's paradox inflating pooled r_s.
        % A minimum of 3 points per group is required for a meaningful
        % correlation estimate (fewer points produce unreliable p-values).
        % Initialize Spearman rho and p-values to NaN (shown when too few data points)
        r_lc = NaN; p_lc = NaN; r_lf = NaN; p_lf = NaN;
        if sum(lc_mask) >= 3
            if exist('OCTAVE_VERSION', 'builtin')
                % Octave uses spearman() and requires manual t-test for p-value.
                % t = r * sqrt((n-2) / (1-r^2)) follows t-distribution with n-2 df.
                % eps prevents division by zero when |r| = 1 (perfect correlation).
                r_lc = spearman(x_vals(lc_mask), y_vals(lc_mask));
                n_lc = sum(lc_mask);
                t_lc = r_lc * sqrt((n_lc - 2) / (1 - r_lc^2 + eps));
                p_lc = 2 * (1 - tcdf(abs(t_lc), n_lc - 2));
            else
                [r_lc, p_lc] = corr(x_vals(lc_mask), y_vals(lc_mask), 'Type', 'Spearman');
            end
        end
        if sum(lf_mask) >= 3
            if exist('OCTAVE_VERSION', 'builtin')
                r_lf = spearman(x_vals(lf_mask), y_vals(lf_mask));
                n_lf = sum(lf_mask);
                t_lf = r_lf * sqrt((n_lf - 2) / (1 - r_lf^2 + eps));
                p_lf = 2 * (1 - tcdf(abs(t_lf), n_lf - 2));
            else
                [r_lf, p_lf] = corr(x_vals(lf_mask), y_vals(lf_mask), 'Type', 'Spearman');
            end
        end

        xlabel(x_label);
        if isempty(diff_units{di})
            ylabel(diff_names{di});
        else
            ylabel([diff_names{di} ' (' diff_units{di} ')']);
        end
        % Format correlation strings; show "n<3" when a group has insufficient data.
        if isnan(r_lc)
            lc_corr_str = 'LC: n<3';
        else
            lc_corr_str = sprintf('LC r_s=%.2f %s', r_lc, format_p_value(p_lc));
        end
        if isnan(r_lf)
            lf_corr_str = 'LF: n<3';
        else
            lf_corr_str = sprintf('LF r_s=%.2f %s', r_lf, format_p_value(p_lf));
        end
        title(sprintf('%s vs Dose\n%s | %s', diff_names{di}, lc_corr_str, lf_corr_str), ...
            'FontSize', 10);
        if exist('OCTAVE_VERSION', 'builtin')
            legend('LC', 'LF', 'LC trend', 'LF trend', 'location', 'best');
        else
            lg = legend('Location', 'best', 'FontSize', 8);
            set(lg, 'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
        end
        grid on;

        plot_idx = plot_idx + 1;
    end
end
sgtitle(['RT Dose vs Diffusion Metrics (Fx1) (' dtype_label ')'], ...
        'FontSize', 14, 'FontWeight', 'bold');
set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);
saveas(gcf, fullfile(output_folder, ['Dose_vs_Diffusion_' dtype_label '.png']));
close(gcf);

fprintf('  Scatter plots generated (%s).\n', dtype_label);

end
