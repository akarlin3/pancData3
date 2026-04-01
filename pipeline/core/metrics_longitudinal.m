function metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, nTp, dtype_label, output_folder, m_lf)
% METRICS_LONGITUDINAL — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 2/5 of the metrics step. Visualizes longitudinal evolution of absolute
% parameters and their relative percent changes.
%
% ANALYTICAL OVERVIEW:
%   Generates "spaghetti plots" showing individual patient trajectories
%   overlaid with population mean +/- SEM for each diffusion biomarker
%   across treatment fractions.  This visualization serves two purposes:
%
%   1. TREATMENT RESPONSE CHARACTERIZATION — Expected biological patterns:
%      - ADC/D should increase during RT as radiation kills tumour cells,
%        reducing cellularity and increasing extracellular water diffusion.
%      - f may decrease if radiation damages tumour microvasculature,
%        reducing the perfusion fraction.
%      - D* is physiologically noisy and may not show consistent trends.
%      Deviations from expected patterns may indicate treatment resistance
%      or inadequate dose delivery.
%
%   2. OUTLIER AND DATA QUALITY ASSESSMENT — Individual patient trajectories
%      (grey lines) reveal measurement artifacts (sudden spikes/drops from
%      poor image quality), protocol deviations, or genuinely unusual
%      biological responses that warrant clinical review.
%
%   The top row shows absolute parameter values (physiological scale),
%   while the bottom row shows changes from baseline (normalised scale),
%   making inter-patient response patterns comparable despite different
%   starting values.
%
%   When m_lf is provided, additional outcome-stratified figures are
%   generated: one per outcome group (LC, LF, Competing Risk) and one
%   combined overlay comparing group mean trajectories.
%
% Inputs:
%   *_abs             - Matrices of absolute values for ADC, D, f, D*
%   *_pct             - Matrices of percent change values for ADC, D, f, D*
%   nTp               - Number of timepoints expected
%   dtype_label       - String label indicating pipeline variant (Standard, etc)
%   output_folder     - Directory where generated figures will be saved
%   m_lf              - (Optional) Outcome vector (0=LC, 1=LF, 2=Competing Risk)
%
% Outputs:
%   None. Saves generated longitudinal plots to output_folder.
%

% Handle optional m_lf argument for backward compatibility
if nargin < 12, m_lf = []; end

fprintf('  --- SECTION 5: Longitudinal Metric Plotting ---\n');

% Diary: capture console output to output_folder
diary_file = fullfile(output_folder, ['metrics_longitudinal_output_' dtype_label '.txt']);
if exist(diary_file, 'file'), delete(diary_file); end
diary(diary_file);

% Group data for easy iteration.
% Order: ADC first (most clinically established biomarker), then IVIM
% parameters in order of decreasing reliability (D > f > D*).
% D* is the least reliable due to the ill-conditioned bi-exponential fit
% and high sensitivity to motion artifacts in the pancreatic bed.
metrics_abs = {ADC_abs, D_abs, f_abs, Dstar_abs};
metrics_pct = {ADC_pct, D_pct, f_delta, Dstar_pct};
metric_names = {'ADC', 'D', 'f', 'D*'};
metric_units = {'mm^2/s', 'mm^2/s', 'Fraction', 'mm^2/s'};
x_labels = [arrayfun(@(x) sprintf('Fx%d', x), 1:(nTp-1), 'UniformOutput', false), {'Post'}];
x_vals = 1:nTp;

% Create figure
figure('Name', ['Longitudinal Mean Metrics — ' dtype_label], 'Position', [100, 100, 1400, 700]);

n_metrics_long = 4;
for i = 1:n_metrics_long
    text_progress_bar(i, n_metrics_long, 'Plotting longitudinal metrics');
    % -------------------------------------------------------------------
    % TOP ROW: Absolute Mean Values
    % -------------------------------------------------------------------
    title_str = ['Mean ', metric_names{i}];
    plot_metric_subplot(i, metrics_abs{i}, x_vals, x_labels, nTp, 'k', 'o', ...
        title_str, metric_units{i}, false);

    % -------------------------------------------------------------------
    % BOTTOM ROW: Percent Change from Fx1 (or absolute change for f)
    % -------------------------------------------------------------------
    % f uses absolute change because its near-zero baseline values
    % (typical range 0.05-0.15) make percent change numerically unstable
    % and clinically uninterpretable.  All other parameters use percent
    % change to normalise for inter-patient baseline variation.
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
% Shift subplot axes down to create visual separation between super-title
% and subplot titles.  The 0.92 scale factor compresses the plot area
% vertically by 8% to make room for sgtitle without overlap.
subplot_scale = 0.92;
allAx = findall(gcf, 'Type', 'Axes');
if ~exist('OCTAVE_VERSION', 'builtin')
    for k = 1:numel(allAx)
        pos = allAx(k).Position;
        allAx(k).Position = [pos(1), pos(2) * subplot_scale, pos(3), pos(4) * subplot_scale];
    end
end
% Remove interactive toolbar from axes for cleaner PNG output
set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);
saveas(gcf, fullfile(output_folder, ['Longitudinal_Mean_Metrics_' dtype_label '.png']));
close(gcf);

% =====================================================================
% OUTCOME-STRATIFIED LONGITUDINAL PLOTS
% =====================================================================
% When clinical outcome data (m_lf) is available, generate:
%   1. Per-outcome spaghetti figures (same layout as all-patients figure)
%   2. Combined overlay figure comparing group mean trajectories
if ~isempty(m_lf) && numel(m_lf) == size(ADC_abs, 1)
    outcome_codes  = [0,                1,    2];
    outcome_labels = {'LC',             'LF', 'CompetingRisk'};
    outcome_titles = {'Local Control',  'Local Failure', 'Competing Risk'};

    % --- Per-outcome figures ---
    for g = 1:numel(outcome_codes)
        idx_g = find(m_lf == outcome_codes(g));
        if isempty(idx_g), continue; end

        fprintf('  Plotting longitudinal metrics for %s (n=%d)...\n', outcome_titles{g}, numel(idx_g));
        figure('Name', ['Longitudinal — ' outcome_titles{g} ' — ' dtype_label], ...
               'Position', [100, 100, 1400, 700]);

        for i = 1:n_metrics_long
            title_str = ['Mean ', metric_names{i}];
            plot_metric_subplot(i, metrics_abs{i}(idx_g, :), x_vals, x_labels, nTp, 'k', 'o', ...
                title_str, metric_units{i}, false);

            if strcmp(metric_names{i}, 'f')
                title_str_pct = ['\Delta ', metric_names{i}, ' (abs)'];
                ylabel_pct = 'Absolute Change';
            else
                title_str_pct = ['\Delta ', metric_names{i}, ' (%)'];
                ylabel_pct = '% Change from Fx1';
            end
            plot_metric_subplot(i+4, metrics_pct{i}(idx_g, :), x_vals, x_labels, nTp, 'r', 's', ...
                title_str_pct, ylabel_pct, true);
        end

        fig_title = sprintf('Longitudinal Metrics — %s (n=%d) [%s]', ...
            outcome_titles{g}, numel(idx_g), dtype_label);
        if exist('OCTAVE_VERSION', 'builtin')
            axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
            text(0.5, 0.98, fig_title, 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
        else
            sgtitle(fig_title, 'FontSize', 16, 'FontWeight', 'bold');
        end
        if ~exist('OCTAVE_VERSION', 'builtin')
            allAx = findall(gcf, 'Type', 'Axes');
            for k = 1:numel(allAx)
                pos = allAx(k).Position;
                allAx(k).Position = [pos(1), pos(2) * subplot_scale, pos(3), pos(4) * subplot_scale];
            end
        end
        set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);
        saveas(gcf, fullfile(output_folder, ['Longitudinal_Mean_Metrics_' dtype_label '_' outcome_labels{g} '.png']));
        close(gcf);
    end

    % --- Combined stratified overlay figure ---
    % Shows mean+SEM for each outcome group on the same axes for direct
    % comparison of treatment response trajectories between groups.
    % This is the key clinical figure: diverging trajectories between
    % LC (blue) and LF (red) groups indicate that diffusion biomarkers
    % can detect treatment resistance during the RT course.
    % Competing risk patients (grey) are shown for completeness but are
    % not the focus of the dose-response hypothesis.
    group_colors = {[0.0 0.4 0.8], [0.9 0.1 0.1], [0.5 0.5 0.5]};  % blue=LC, red=LF, grey=CR
    group_markers = {'o', 's', 'd'};
    % Small x-offsets to prevent error bar overlap between groups
    group_offsets = [-0.08, 0.0, 0.08];

    figure('Name', ['Longitudinal by Outcome — ' dtype_label], ...
           'Position', [100, 100, 1400, 700]);

    for i = 1:n_metrics_long
        % Top row: absolute values
        subplot(2, 4, i); hold on;
        legend_entries = {};
        h_lines = [];
        for g = 1:numel(outcome_codes)
            idx_g = find(m_lf == outcome_codes(g));
            if isempty(idx_g), continue; end
            dat_g = metrics_abs{i}(idx_g, :);
            h_g = plot_group_mean_sem(x_vals + group_offsets(g), dat_g, group_colors{g}, group_markers{g});
            h_lines(end+1) = h_g; %#ok<AGROW>
            legend_entries{end+1} = sprintf('%s (n=%d)', outcome_titles{g}, numel(idx_g)); %#ok<AGROW>
        end
        set(gca, 'XTick', x_vals, 'XTickLabel', x_labels, 'FontSize', 10);
        title(['Mean ', metric_names{i}], 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(metric_units{i});
        xlim([0.5, nTp+0.5]);
        grid on; box on;
        if i == 1 && ~isempty(legend_entries)
            legend(h_lines, legend_entries, 'Location', 'best', 'FontSize', 7, ...
                'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
        end

        % Bottom row: percent change (or absolute change for f)
        subplot(2, 4, i+4); hold on;
        for g = 1:numel(outcome_codes)
            idx_g = find(m_lf == outcome_codes(g));
            if isempty(idx_g), continue; end
            dat_g = metrics_pct{i}(idx_g, :);
            plot_group_mean_sem(x_vals + group_offsets(g), dat_g, group_colors{g}, group_markers{g});
        end
        if strcmp(metric_names{i}, 'f')
            title_str_pct = ['\Delta ', metric_names{i}, ' (abs)'];
            ylabel_pct = 'Absolute Change';
        else
            title_str_pct = ['\Delta ', metric_names{i}, ' (%)'];
            ylabel_pct = '% Change from Fx1';
        end
        if exist('OCTAVE_VERSION', 'builtin')
            xl = xlim;
            plot(xl, [0 0], 'k--', 'LineWidth', 1.5);
        else
            yline(0, 'k--', 'LineWidth', 1.5);
        end
        set(gca, 'XTick', x_vals, 'XTickLabel', x_labels, 'FontSize', 10);
        title(title_str_pct, 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(ylabel_pct);
        xlim([0.5, nTp+0.5]);
        grid on; box on;
    end

    fig_title = ['Longitudinal Metrics by Outcome (' dtype_label ')'];
    if exist('OCTAVE_VERSION', 'builtin')
        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        text(0.5, 0.98, fig_title, 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
    else
        sgtitle(fig_title, 'FontSize', 16, 'FontWeight', 'bold');
    end
    if ~exist('OCTAVE_VERSION', 'builtin')
        allAx = findall(gcf, 'Type', 'Axes');
        for k = 1:numel(allAx)
            pos = allAx(k).Position;
            allAx(k).Position = [pos(1), pos(2) * subplot_scale, pos(3), pos(4) * subplot_scale];
        end
    end
    set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);
    saveas(gcf, fullfile(output_folder, ['Longitudinal_Mean_Metrics_' dtype_label '_ByOutcome.png']));
    close(gcf);

    fprintf('  Outcome-stratified longitudinal plots saved.\n');
end

diary off;
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

% Calculate population mean and standard error of the mean (SEM).
% SEM (not SD) is used for error bars because we want to characterise
% uncertainty in the population mean trajectory, not individual patient
% variability.  In RT response assessment, the question is "does the
% average tumour respond?" (requires SEM), not "how variable are
% individual responses?" (requires SD).
if exist('OCTAVE_VERSION', 'builtin')
    valid_mask = ~isnan(dat);
    N = sum(valid_mask, 1);
    dat_zero = dat;
    dat_zero(~valid_mask) = 0;
    pop_mean = sum(dat_zero, 1) ./ N;

    devs = dat - repmat(pop_mean, size(dat, 1), 1);
    devs_sq = devs.^2;
    devs_sq(~valid_mask) = 0;
    pop_std = sqrt(sum(devs_sq, 1) ./ max(N - 1, 1));
    pop_se = pop_std ./ sqrt(max(N, 1));
    pop_se(N < 2) = NaN;
else
    pop_mean = mean(dat, 1, 'omitnan');
    N_pop = sum(~isnan(dat), 1);
    pop_se   = std(dat, 0, 1, 'omitnan') ./ sqrt(max(N_pop, 1));
    pop_se(N_pop < 2) = NaN;
end

% Plot individual patient trajectories (spaghetti plot).
% Light grey lines show each patient's evolution over treatment fractions.
% This reveals heterogeneous response patterns: some patients may show
% early ADC increases (responding), while others show stable or
% decreasing ADC (potentially resistant).  Divergent trajectories between
% responders and non-responders are the biological signal this pipeline
% aims to quantify.
h_indiv = plot(x_vals, dat', 'Color', [0.8 0.8 0.8 0.4], 'LineWidth', 0.5);

% Plot population average with error bars
if exist('OCTAVE_VERSION', 'builtin')
    % Octave's errorbar doesn't support the line spec string correctly with the Marker property
    % the same way MATLAB does in this specific overload, so we simplify for the mock
    h_mean = errorbar(x_vals, pop_mean, pop_se);
else
    h_mean = errorbar(x_vals, pop_mean, pop_se, ['-' color_spec], 'LineWidth', 2, ...
        'Marker', marker_style, 'MarkerFaceColor', color_spec, 'MarkerSize', 7, ...
        'CapSize', 8);
end

% Add legend on first subplot only to avoid clutter
if idx == 1 || idx == 5
    n_pts = size(dat, 1);
    legend([h_indiv(1), h_mean], {sprintf('Individual (n=%d)', n_pts), 'Mean \pm SEM'}, ...
        'Location', 'best', 'FontSize', 7, 'Box', 'on', 'EdgeColor', [0.5 0.5 0.5]);
end

if add_zero_line
    % Add a baseline reference line at 0% change.  Points above this line
    % indicate an increase from baseline (e.g., rising ADC = cell death),
    % while points below indicate a decrease (e.g., falling f = vascular
    % disruption).  The line helps quickly identify the fraction at which
    % treatment effects become detectable above measurement noise.
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

function h = plot_group_mean_sem(x_vals, dat, color_rgb, marker_style)
% PLOT_GROUP_MEAN_SEM — Plot mean+SEM line for a single outcome group
%
% Used by the combined stratified overlay figure.  Plots into the current
% axes (caller must have called subplot + hold on).
%
% Inputs:
%   x_vals       - X-axis values (1:nTp)
%   dat          - Data matrix (patients x timepoints) for this group
%   color_rgb    - 1x3 RGB color vector
%   marker_style - Marker string (e.g., 'o', 's', 'd')

if exist('OCTAVE_VERSION', 'builtin')
    valid_mask = ~isnan(dat);
    N = sum(valid_mask, 1);
    dat_zero = dat;
    dat_zero(~valid_mask) = 0;
    grp_mean = sum(dat_zero, 1) ./ N;

    devs = dat - repmat(grp_mean, size(dat, 1), 1);
    devs_sq = devs.^2;
    devs_sq(~valid_mask) = 0;
    grp_std = sqrt(sum(devs_sq, 1) ./ (N - 1));
    grp_se = grp_std ./ sqrt(N);
    grp_se(N < 2) = NaN;
else
    grp_mean = mean(dat, 1, 'omitnan');
    N_grp = sum(~isnan(dat), 1);
    grp_se   = std(dat, 0, 1, 'omitnan') ./ sqrt(max(N_grp, 1));
    grp_se(N_grp < 2) = NaN;
end

if exist('OCTAVE_VERSION', 'builtin')
    h = errorbar(x_vals, grp_mean, grp_se);
else
    h = errorbar(x_vals, grp_mean, grp_se, '-', 'Color', color_rgb, 'LineWidth', 2, ...
        'Marker', marker_style, 'MarkerFaceColor', color_rgb, 'MarkerSize', 7, ...
        'CapSize', 8);
end

end
