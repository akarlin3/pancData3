function metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_pct, Dstar_pct, nTp, dtype_label, output_folder)
% METRICS_LONGITUDINAL — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 2/5 of the metrics step.

fprintf('  --- SECTION 5: Longitudinal Metric Plotting ---\n');

% Group data for easy iteration
metrics_abs = {ADC_abs, D_abs, f_abs, Dstar_abs};
metrics_pct = {ADC_pct, D_pct, f_pct, Dstar_pct};
metric_names = {'ADC', 'D', 'f', 'D*'};
metric_units = {'mm^2/s', 'mm^2/s', 'Fraction', 'mm^2/s'};
x_labels = {'Fx1', 'Fx2', 'Fx3', 'Fx4', 'Fx5', 'Post'};
x_vals = 1:nTp;

% Create figure
figure('Name', ['Longitudinal Mean Metrics — ' dtype_label], 'Position', [100, 100, 1400, 700]);

for i = 1:4
    % -------------------------------------------------------------------
    % TOP ROW: Absolute Mean Values
    % -------------------------------------------------------------------
    subplot(2, 4, i);
    hold on;
    
    dat = metrics_abs{i};
    % Calculate population mean and standard error of the mean (SEM)
    pop_mean = mean(dat, 1, 'omitnan');
    pop_se   = std(dat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(dat), 1));
    
    % Plot individual patient trajectories (spaghetti plot)
    plot(x_vals, dat', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    
    % Plot population average with error bars
    errorbar(x_vals, pop_mean, pop_se, '-k', 'LineWidth', 2, ...
        'Marker', 'o', 'MarkerFaceColor', 'k');
    
    % Formatting
    set(gca, 'XTick', x_vals, 'XTickLabel', x_labels, 'FontSize', 10);
    title(['Mean ', metric_names{i}], 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(metric_units{i});
    xlim([0.5, nTp+0.5]);
    grid on; box on;
    
    % -------------------------------------------------------------------
    % BOTTOM ROW: Percent Change from Fx1
    % -------------------------------------------------------------------
    subplot(2, 4, i+4);
    hold on;
    
    dat_pct = metrics_pct{i};
    % Calculate population mean percent change and SEM
    pop_mean_pct = mean(dat_pct, 1, 'omitnan');
    pop_se_pct   = std(dat_pct, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(dat_pct), 1));
    
    % Plot individual patient trajectories
    plot(x_vals, dat_pct', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    
    % Plot population average with error bars
    errorbar(x_vals, pop_mean_pct, pop_se_pct, '-r', 'LineWidth', 2, ...
        'Marker', 's', 'MarkerFaceColor', 'r');
    
    % Add a baseline reference line at 0% change
    yline(0, 'k--', 'LineWidth', 1.5);
    
    % Formatting
    set(gca, 'XTick', x_vals, 'XTickLabel', x_labels, 'FontSize', 10);
    if strcmp(metric_names{i}, 'f')
        title(['\Delta ', metric_names{i}, ' (abs)'], 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Absolute Change');
    else
        title(['\Delta ', metric_names{i}, ' (%)'], 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('% Change from Fx1');
    end
    xlim([0.5, nTp+0.5]);
    grid on; box on;
end

sgtitle(['Longitudinal Evolution of DWI and IVIM Metrics (' dtype_label ')'], 'FontSize', 16, 'FontWeight', 'bold');
% Shift subplot axes down to create visual separation between super-title and subplot titles
subplot_scale = 0.92;
allAx = findall(gcf, 'Type', 'Axes');
for k = 1:numel(allAx)
    pos = allAx(k).Position;
    allAx(k).Position = [pos(1), pos(2) * subplot_scale, pos(3), pos(4) * subplot_scale];
end
saveas(gcf, fullfile(output_folder, ['Longitudinal_Mean_Metrics_' dtype_label '.png']));
close(gcf);

end
