function plot_feature_distributions(dtype_label, adc_mean, d_mean, f_mean, dstar_mean, valid_pts, lf_group, dtype, output_folder)
% PLOT_FEATURE_DISTRIBUTIONS
%
%  For each DWI processing pipeline (Standard, dnCNN, IVIMnet):
%    2a. Histograms — overlay Local-Control vs Local-Failure distributions
%        for each of the four baseline (Fx1) diffusion biomarkers.
%    2b. Box plots  — side-by-side LC vs LF boxes with one-way ANOVA
%        p-values annotated on each panel.
%  The resulting figures are saved as PNG files in saved_figures/.

% Collect baseline (Fx1) biomarker values for the valid patient subset.
% Each cell element is an [nValid x 1] vector for one metric.
metric_data  = {adc_mean(valid_pts,1,dtype), ...
                d_mean(valid_pts,1,dtype), ...
                f_mean(valid_pts,1,dtype), ...
                dstar_mean(valid_pts,1,dtype)};
metric_names = {'Mean ADC', 'Mean D', 'Mean f', 'Mean D*'};
metric_units = {'mm^2/s',   'mm^2/s',  '',       'mm^2/s'};

% --- 2a. Histograms ---
% Use a taller figure (450 px) to prevent the sgtitle from overlapping
% the subplot titles, and shift subplots down for extra clearance.
figure('Name', ['Feature Distributions — Histograms — ' dtype_label], ...
       'Position', [100, 100, 1200, 450]);
for mi = 1:4
    subplot(1, 4, mi);
    vals = metric_data{mi};

    plot_feature_distribution(vals, lf_group, metric_names{mi}, metric_units{mi}, 'histogram');
end
sgtitle(['Baseline Feature Distributions by Outcome (' dtype_label ')'], ...
        'FontSize', 14, 'FontWeight', 'bold');
% Shift all subplot axes down slightly so the supertitle does not overlap
% the individual subplot titles.
allAx = findall(gcf, 'Type', 'Axes');
for k = 1:numel(allAx)
    pos = allAx(k).Position;
    allAx(k).Position = [pos(1), pos(2) * 0.92, pos(3), pos(4) * 0.92];
end
saveas(gcf, fullfile(output_folder, ['Feature_Histograms_' dtype_label '.png']));
close(gcf);

% --- 2b. Box Plots ---
% Same height increase and subplot shift applied here.
figure('Name', ['Feature Distributions — Box Plots — ' dtype_label], ...
       'Position', [100, 500, 1200, 450]);
for mi = 1:4
    subplot(1, 4, mi);
    vals = metric_data{mi};

    plot_feature_distribution(vals, lf_group, metric_names{mi}, metric_units{mi}, 'boxplot');
end
sgtitle(['Baseline Feature Box Plots by Outcome (' dtype_label ')'], ...
        'FontSize', 14, 'FontWeight', 'bold');
% Shift subplot axes down to avoid supertitle/subplot-title overlap
allAx = findall(gcf, 'Type', 'Axes');
for k = 1:numel(allAx)
    pos = allAx(k).Position;
    allAx(k).Position = [pos(1), pos(2) * 0.92, pos(3), pos(4) * 0.92];
end
saveas(gcf, fullfile(output_folder, ['Feature_BoxPlots_' dtype_label '.png']));
close(gcf);

fprintf('  Histograms and box plots generated for ADC, D, f, D* (%s).\n', dtype_label);

end
