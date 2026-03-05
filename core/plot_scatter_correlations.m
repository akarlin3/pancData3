function plot_scatter_correlations(dtype_label, dmean_gtvp, d95_gtvp, adc_mean, d_mean, f_mean, valid_pts, lf_group, dtype, output_folder)
% PLOT_SCATTER_CORRELATIONS
%
%  For each diffusion metric (ADC, D, f) plot it against two dose
%  endpoints — Mean GTV Dose and D95 — to explore potential dose–response
%  relationships.  Each scatter panel is coloured by clinical outcome
%  (blue = LC, red = LF), with a linear trend-line and Spearman
%  correlation coefficient annotated.

% Extract Fx1 dose metrics for the valid patient subset
dose_mean_vec = dmean_gtvp(valid_pts, 1);   % mean dose inside GTV (Gy)
dose_d95_vec  = d95_gtvp(valid_pts, 1);     % D95 = dose to 95 % of GTV (Gy)

% Diffusion biomarkers to correlate with dose
diff_metrics = {adc_mean(valid_pts,1,dtype), d_mean(valid_pts,1,dtype), f_mean(valid_pts,1,dtype)};
diff_names   = {'Mean ADC', 'Mean D', 'Mean f'};
diff_units   = {'mm^2/s',   'mm^2/s',  ''};

figure('Name', ['Dose vs Diffusion Metrics — ' dtype_label], ...
       'Position', [150, 150, 1400, 500]);

plot_idx = 1;
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

        % Ensure lf_group is a column vector to match clean, x_vals, and y_vals
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

        % Plot LC (blue) and LF (red) points with black edge
        scatter(x_vals(clean & lf_group_col==0), y_vals(clean & lf_group_col==0), ...
            50, [0 0.4470 0.7410], 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'LC'); hold on;
        scatter(x_vals(clean & lf_group_col==1), y_vals(clean & lf_group_col==1), ...
            50, [0.8500 0.3250 0.0980], 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'LF');

        % Overlay per-group linear trend lines to avoid pooled Simpson's paradox
        warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
        lc_mask = clean & lf_group_col==0;
        lf_mask = clean & lf_group_col==1;
        if sum(lc_mask) >= 2
            x_line_lc = linspace(min(x_vals(lc_mask)), max(x_vals(lc_mask)), 50);
            p_fit_lc = polyfit(x_vals(lc_mask), y_vals(lc_mask), 1);
            plot(x_line_lc, polyval(p_fit_lc, x_line_lc), '--', 'Color', [0 0.4470 0.7410], ...
                'LineWidth', 1.5, 'DisplayName', 'LC trend');
        end
        if sum(lf_mask) >= 2
            x_line_lf = linspace(min(x_vals(lf_mask)), max(x_vals(lf_mask)), 50);
            p_fit_lf = polyfit(x_vals(lf_mask), y_vals(lf_mask), 1);
            plot(x_line_lf, polyval(p_fit_lf, x_line_lf), '--', 'Color', [0.8500 0.3250 0.0980], ...
                'LineWidth', 1.5, 'DisplayName', 'LF trend');
        end
        warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');
        hold off;

        % Compute per-group Spearman correlations to match the per-group
        % trend lines and avoid Simpson's paradox inflating pooled r_s.
        r_lc = NaN; p_lc = NaN; r_lf = NaN; p_lf = NaN;
        if sum(lc_mask) >= 3
            if exist('OCTAVE_VERSION', 'builtin')
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
        title(sprintf('%s vs Dose\nLC r_s=%.2f %s | LF r_s=%.2f %s', ...
            diff_names{di}, r_lc, format_p_value(p_lc), r_lf, format_p_value(p_lf)), ...
            'FontSize', 10);
        if exist('OCTAVE_VERSION', 'builtin')
            legend('LC', 'LF', 'Linear fit', 'location', 'best');
        else
            legend('Location', 'best', 'FontSize', 7);
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
