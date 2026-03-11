function plot_feature_distributions(dtype_label, adc_mean, d_mean, f_mean, dstar_mean, valid_pts, lf_group, dtype, output_folder)
% PLOT_FEATURE_DISTRIBUTIONS
%
% ANALYTICAL OVERVIEW:
%   Visualises the baseline (pre-treatment, Fx1) distribution of each
%   diffusion biomarker, stratified by clinical outcome (LC vs LF).
%   This is the most fundamental exploratory analysis: if a biomarker's
%   baseline distribution does not differ between outcome groups, it is
%   unlikely to be useful as a pre-treatment predictor of local failure.
%
%   Two complementary visualisations are generated:
%     2a. Histograms — Show the full distributional shape (skewness,
%         bimodality, overlap between groups).  Overlapping distributions
%         indicate limited discriminative power.  Bimodal distributions
%         may suggest distinct tumour subtypes with different prognoses.
%     2b. Box plots  — Summarise location (median), spread (IQR), and
%         outliers.  Annotated with one-way ANOVA p-values to quantify
%         whether group means differ.  Note: ANOVA assumes normality;
%         the non-parametric Wilcoxon test in metrics_stats_comparisons
%         is the definitive test.
%
%   Expected biological patterns at baseline:
%     - ADC/D: LF tumours may have LOWER values (denser cellularity,
%       more treatment-resistant tissue).
%     - f: LF tumours may have LOWER perfusion (hypoxic, less radiosensitive).
%     - D*: Often too noisy for reliable pre-treatment prediction.
%
%  The resulting figures are saved as PNG files in the output folder.

% Collect baseline (Fx1) biomarker values for the valid patient subset.
% Each cell element is an [nValid x 1] vector for one metric.
% Only Fx1 (column 1) is used because this analysis characterises
% pre-treatment tumour properties — before radiation has induced any changes.
% The third index selects the DWI processing pipeline (1=Standard, etc.).
metric_data  = {adc_mean(valid_pts,1,dtype), ...
                d_mean(valid_pts,1,dtype), ...
                f_mean(valid_pts,1,dtype), ...
                dstar_mean(valid_pts,1,dtype)};
metric_names = {'Mean ADC', 'Mean D', 'Mean f', 'Mean D*'};
% f is dimensionless (volume fraction [0,1]), so no unit is displayed.
% ADC, D, and D* are all in mm^2/s (diffusion coefficient units).
metric_units = {'mm^2/s',   'mm^2/s',  '',       'mm^2/s'};

% --- 2a. Histograms ---
% Use a taller figure (450 px) to prevent the sgtitle from overlapping
% the subplot titles, and shift subplots down for extra clearance.
figure('Name', ['Feature Distributions — Histograms — ' dtype_label], ...
       'Position', [100, 100, 1200, 450]);
n_metrics_dist = 4;  % ADC, D, f, D*
for mi = 1:n_metrics_dist
    text_progress_bar(mi, n_metrics_dist, 'Generating histograms');
    subplot(1, 4, mi);
    vals = metric_data{mi};

    % plot_feature_distribution handles NaN removal, group stratification
    % (LC vs LF), bin selection, and ANOVA p-value annotation internally.
    plot_feature_distribution(vals, lf_group, metric_names{mi}, metric_units{mi}, 'histogram');
end
sgtitle(['Baseline Feature Distributions by Outcome (' dtype_label ')'], ...
        'FontSize', 14, 'FontWeight', 'bold');
% Shift all subplot axes down slightly so the supertitle does not overlap
% the individual subplot titles. The 0.92 multiplier reduces both the
% vertical position and height by 8%, creating clearance above each subplot.
allAx = findall(gcf, 'Type', 'Axes');
for k = 1:numel(allAx)
    pos = get(allAx(k), 'Position');
    set(allAx(k), 'Position', [pos(1), pos(2) * 0.92, pos(3), pos(4) * 0.92]);
end
% Disable axes toolbar widgets (zoom, pan, data tips) for clean figure export
set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);
saveas(gcf, fullfile(output_folder, ['Feature_Histograms_' dtype_label '.png']));
close(gcf);

% --- 2b. Box Plots ---
% Same height increase and subplot shift applied here.
figure('Name', ['Feature Distributions — Box Plots — ' dtype_label], ...
       'Position', [100, 500, 1200, 450]);
for mi = 1:n_metrics_dist
    text_progress_bar(mi, n_metrics_dist, 'Generating box plots');
    subplot(1, 4, mi);
    vals = metric_data{mi};

    plot_feature_distribution(vals, lf_group, metric_names{mi}, metric_units{mi}, 'boxplot');
end
sgtitle(['Baseline Feature Box Plots by Outcome (' dtype_label ')'], ...
        'FontSize', 14, 'FontWeight', 'bold');
% Shift subplot axes down to avoid supertitle/subplot-title overlap
allAx = findall(gcf, 'Type', 'Axes');
for k = 1:numel(allAx)
    pos = get(allAx(k), 'Position');
    set(allAx(k), 'Position', [pos(1), pos(2) * 0.92, pos(3), pos(4) * 0.92]);
end
set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);
saveas(gcf, fullfile(output_folder, ['Feature_BoxPlots_' dtype_label '.png']));
close(gcf);

fprintf('  Histograms and box plots generated for ADC, D, f, D* (%s).\n', dtype_label);

end
