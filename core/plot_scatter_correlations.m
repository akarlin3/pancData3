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
for di = 1:numel(diff_metrics)
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

        % Identify patients with both dose and diffusion data available
        clean = ~isnan(x_vals) & ~isnan(y_vals);
        if sum(clean) < 3
            title([diff_names{di} ' — insufficient data']);
            plot_idx = plot_idx + 1;
            continue;
        end

        % Plot LC (blue) and LF (red) points with black edge
        scatter(x_vals(clean & lf_group==0), y_vals(clean & lf_group==0), ...
            50, 'b', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'LC'); hold on;
        scatter(x_vals(clean & lf_group==1), y_vals(clean & lf_group==1), ...
            50, 'r', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'LF');

        % Overlay a first-order (linear) polynomial fit
        p_fit = polyfit(x_vals(clean), y_vals(clean), 1);
        x_line = linspace(min(x_vals(clean)), max(x_vals(clean)), 50);
        plot(x_line, polyval(p_fit, x_line), 'k--', 'LineWidth', 1.5, ...
            'DisplayName', 'Linear fit');
        hold off;

        % Compute Spearman rank correlation and annotate the title
        [r_sp, p_sp] = corr(x_vals(clean), y_vals(clean), 'Type', 'Spearman');

        xlabel(x_label);
        ylabel([diff_names{di} ' (' diff_units{di} ')']);
        title(sprintf('%s vs Dose\nr_s=%.2f, p=%.3f', diff_names{di}, r_sp, p_sp), ...
            'FontSize', 10);
        legend('Location', 'best', 'FontSize', 7);
        grid on;

        plot_idx = plot_idx + 1;
    end
end
sgtitle(['RT Dose vs Diffusion Metrics (Fx1) (' dtype_label ')'], ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(output_folder, ['Dose_vs_Diffusion_' dtype_label '.png']));
close(gcf);

fprintf('  Scatter plots generated (%s).\n', dtype_label);

end
