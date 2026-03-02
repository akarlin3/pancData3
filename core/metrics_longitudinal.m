function metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_pct, Dstar_pct, nTp, dtype_label, output_folder)
% METRICS_LONGITUDINAL — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 2/5 of the metrics step. Visualizes longitudinal evolution of absolute
% parameters and their relative percent changes.
%
% Inputs:
%   *_abs             - Matrices of absolute values for ADC, D, f, D*
%   *_pct             - Matrices of percent change values for ADC, D, f, D*
%   nTp               - Number of timepoints expected
%   dtype_label       - String label indicating pipeline variant (Standard, etc)
%   output_folder     - Directory where generated figures will be saved
%
% Outputs:
%   None. Saves generated longitudinal plots to output_folder.
%

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
    title_str = ['Mean ', metric_names{i}];
    plot_metric_subplot(i, metrics_abs{i}, x_vals, x_labels, nTp, 'k', 'o', ...
        title_str, metric_units{i}, false);
    
    % -------------------------------------------------------------------
    % BOTTOM ROW: Percent Change from Fx1
    % -------------------------------------------------------------------
    if strcmp(metric_names{i}, 'f')
        title_str_pct = ['\Delta ', metric_names{i}, ' (abs)'];
        ylabel_pct = 'Absolute Change';
    else
        title_str_pct = ['\Delta ', metric_names{i}, ' (%)'];
        ylabel_pct = '% Change from Fx1';
    end
    plot_metric_subplot(i+4, metrics_pct{i}, x_vals, x_labels, nTp, 'r', 's', ...
        title_str_pct, ylabel_pct, true);
end

if exist('OCTAVE_VERSION', 'builtin')
    % Fallback for sgtitle in Octave
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.98, ['Longitudinal Evolution of DWI and IVIM Metrics (' dtype_label ')'], 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
else
    sgtitle(['Longitudinal Evolution of DWI and IVIM Metrics (' dtype_label ')'], 'FontSize', 16, 'FontWeight', 'bold');
end
% Shift subplot axes down to create visual separation between super-title and subplot titles
subplot_scale = 0.92;
allAx = findall(gcf, 'Type', 'Axes');
if ~exist('OCTAVE_VERSION', 'builtin')
    for k = 1:numel(allAx)
        pos = allAx(k).Position;
        allAx(k).Position = [pos(1), pos(2) * subplot_scale, pos(3), pos(4) * subplot_scale];
    end
end
saveas(gcf, fullfile(output_folder, ['Longitudinal_Mean_Metrics_' dtype_label '.png']));
close(gcf);

end

function plot_metric_subplot(idx, dat, x_vals, x_labels, nTp, color_spec, marker_style, title_str, y_label, add_zero_line)
% PLOT_METRIC_SUBPLOT — Helper to generate standardized longitudinal subplot
%
% Inputs:
%   idx           - Subplot index (1 to 8)
%   dat           - Data matrix (patients x timepoints)
%   x_vals        - X-axis values
%   x_labels      - X-axis tick labels
%   nTp           - Number of timepoints
%   color_spec    - Color specification string (e.g., 'k', 'r')
%   marker_style  - Marker string (e.g., 'o', 's')
%   title_str     - Subplot title string
%   y_label       - Y-axis label string
%   add_zero_line - Boolean indicating whether to add a dashed y=0 line

subplot(2, 4, idx);
hold on;

% Calculate population mean and standard error of the mean (SEM)
if exist('OCTAVE_VERSION', 'builtin')
    valid_mask = ~isnan(dat);
    N = sum(valid_mask, 1);
    dat_zero = dat;
    dat_zero(~valid_mask) = 0;
    pop_mean = sum(dat_zero, 1) ./ N;

    devs = dat - repmat(pop_mean, size(dat, 1), 1);
    devs_sq = devs.^2;
    devs_sq(~valid_mask) = 0;
    pop_std = sqrt(sum(devs_sq, 1) ./ max(1, N - 1));
    pop_se = pop_std ./ sqrt(N);
else
    pop_mean = mean(dat, 1, 'omitnan');
    pop_se   = std(dat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(dat), 1));
end

% Plot individual patient trajectories (spaghetti plot)
plot(x_vals, dat', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);

% Plot population average with error bars
if exist('OCTAVE_VERSION', 'builtin')
    % Octave's errorbar doesn't support the line spec string correctly with the Marker property
    % the same way MATLAB does in this specific overload, so we simplify for the mock
    errorbar(x_vals, pop_mean, pop_se);
else
    errorbar(x_vals, pop_mean, pop_se, ['-' color_spec], 'LineWidth', 2, ...
        'Marker', marker_style, 'MarkerFaceColor', color_spec);
end

if add_zero_line
    % Add a baseline reference line at 0% change
    if exist('OCTAVE_VERSION', 'builtin')
        % Fallback for yline in Octave
        xl = xlim;
        plot(xl, [0 0], 'k--', 'LineWidth', 1.5);
    else
        yline(0, 'k--', 'LineWidth', 1.5);
    end
end

% Formatting
set(gca, 'XTick', x_vals, 'XTickLabel', x_labels, 'FontSize', 10);
title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
ylabel(y_label);
xlim([0.5, nTp+0.5]);
grid on; box on;

end
