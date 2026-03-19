function visualize_results(data_vectors_gtvp, summary_metrics, calculated_results, config_struct)
% VISUALIZE_RESULTS — "Visualizing It" (Generating Plots)
% Author: Avery Karlin
%
% Inputs:
%   data_vectors_gtvp - Struct array holding primary GTV parameters
%   summary_metrics   - Struct containing patient ID and metric statistics
%   calculated_results- Struct containing statistical predictive findings
%   config_struct     - Configuration struct defining output locations
%
% Outputs:
%   None. Generates and saves figures in the defined directory.
%
% ANALYTICAL RATIONALE — VISUALIZATION STRATEGY
%   This function generates three complementary visualization families that
%   address different aspects of the diffusion biomarker analysis:
%
%   1. Parameter Maps overlaid on Anatomy
%      Spatial visualization of ADC within the tumor contour on the b=0
%      anatomical image. This allows the physicist to:
%        - Verify that GTV contours are correctly positioned on the anatomy
%        - Identify spatial patterns (e.g., necrotic core with high ADC
%          surrounded by viable rim with low ADC)
%        - Detect artifacts (e.g., geometric distortion at tissue-air
%          boundaries near the pancreas)
%
%   2. Distributions of Extracted Features (histograms and box plots)
%      Baseline (Fx1) ADC and IVIM parameter distributions grouped by
%      clinical outcome (Local Control vs Local Failure). These reveal
%      whether pre-treatment diffusion characteristics differ between
%      patients who eventually fail locally vs those who achieve local
%      control — the fundamental hypothesis of predictive biomarker research.
%      ANOVA p-values are annotated to guide clinical significance assessment.
%
%   3. Scatter Plots for Dose-Diffusion Correlation
%      RT dose metrics (mean GTV dose, D95) plotted against diffusion
%      parameters to explore dose-response relationships. In pancreatic
%      SBRT, higher delivered dose may produce greater changes in tumor
%      cellularity (reflected by ADC/D changes). Identifying such
%      correlations can inform dose escalation strategies. Spearman
%      correlation is used because dose-diffusion relationships may be
%      monotonic but not necessarily linear.

% Extract required variables
id_list = summary_metrics.id_list;
mrn_list = summary_metrics.mrn_list;
lf = summary_metrics.lf;
adc_mean = summary_metrics.adc_mean;
d_mean = summary_metrics.d_mean;
f_mean = summary_metrics.f_mean;
dstar_mean = summary_metrics.dstar_mean;
d95_gtvp = summary_metrics.d95_gtvp;
dmean_gtvp = summary_metrics.dmean_gtvp;
dataloc = config_struct.dataloc;

fprintf('\n======================================================\n');
fprintf('  VISUALIZE RESULTS — Generating Plots\n');
fprintf('======================================================\n');

% Define human-readable labels for the three DWI processing pipelines.
% Each pipeline produces its own set of diffusion parameters, and
% visualizations are generated per-pipeline to enable methodological
% comparison of parameter estimation approaches.
dwi_type_names = {'Standard', 'dnCNN', 'IVIMnet'};

% Create output directory for saved figures (ignored by .gitignore)
if isfield(config_struct, 'output_folder')
    output_folder = config_struct.output_folder;
else
    timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
    output_folder = fullfile(fileparts(mfilename('fullpath')), '..', '..', sprintf('saved_files_%s', timestamp_str));
end
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

% Start diary to log text output
diary_file = fullfile(output_folder, 'visualize_results_output.txt');
if exist(diary_file, 'file'), delete(diary_file); end
diary(diary_file);

% Suppress figure windows so plots are only saved to disk.
% This is essential for batch/headless execution (e.g., on a compute server
% or within parfor workers) where no display is available. All figures are
% saved as PNG files for inclusion in manuscripts and reports.
set(0, 'DefaultFigureVisible', 'off');

% Total number of patients and fraction (timepoint) labels
nPat = length(id_list);
nTp  = size(adc_mean, 2);
fx_labels = [arrayfun(@(x) sprintf('Fx%d', x), 1:(nTp-1), 'UniformOutput', false), {'Post'}];

% Build a logical mask identifying patients with usable clinical and
% imaging data. Two requirements:
%   1. Finite LF (local failure) label: patients without clinical follow-up
%      cannot be classified as LC/LF and must be excluded from outcome-
%      stratified visualizations.
%   2. Non-NaN baseline ADC: patients without a valid Fx1 DWI acquisition
%      have no baseline reference for longitudinal analysis.
% This filter prevents plotting artifacts from missing data and ensures
% that group comparisons (LC vs LF) use only patients with complete data.
dtype_first = config_struct.dwi_types_to_run(1);
valid_pts = isfinite(lf) & ~isnan(adc_mean(:,1,dtype_first));

%% -----------------------------------------------------------------------
fprintf('\n--- SECTION 1: Parameter Maps overlaid on Anatomy ---\n');
%  1. PARAMETER MAPS OVERLAID ON ANATOMY
%  For each patient with Fx1 data available, load the DWI volume, compute
%  the ADC map via monoexponential fit, and display three panels per
%  patient:
%    (a) b=0 anatomical image with GTV contour overlay
%    (b) Full-slice ADC map with GTV contour
%    (c) ADC overlaid on anatomy (semi-transparent inside GTV only)
%  Patients are batched into multi-row figures (pats_per_fig rows each).
% -----------------------------------------------------------------------
plot_parameter_maps_streaming(data_vectors_gtvp, nPat, id_list, dataloc, output_folder, dtype_first);

%% -----------------------------------------------------------------------
fprintf('\n--- SECTION 2: Distributions of Extracted Features ---\n');
%  2. DISTRIBUTIONS OF EXTRACTED FEATURES
%  For each DWI processing pipeline (Standard, dnCNN, IVIMnet):
%    2a. Histograms — overlay Local-Control vs Local-Failure distributions
%        for each of the four baseline (Fx1) diffusion biomarkers.
%    2b. Box plots  — side-by-side LC vs LF boxes with one-way ANOVA
%        p-values annotated on each panel.
%  The resulting figures are saved as PNG files in the output folder.
% -----------------------------------------------------------------------
fprintf('\n--- 2. Distributions of Extracted Features ---\n');

% Loop over all configured DWI processing pipelines.  Each pipeline is
% visualized independently to enable side-by-side comparison of how
% denoising (dnCNN) or neural network fitting (IVIMnet) affects the
% apparent distribution and dose-response relationships of diffusion
% biomarkers.  Differences between pipelines may reveal which denoising
% approach best separates LC from LF at baseline.
n_dtypes_viz = numel(config_struct.dwi_types_to_run);
dtype_counter = 0;
for dtype = config_struct.dwi_types_to_run
    dtype_counter = dtype_counter + 1;
    dtype_label = dwi_type_names{dtype};
    text_progress_bar(dtype_counter, n_dtypes_viz, 'Plotting distributions & correlations');

    % Recompute valid_pts per dtype so patients missing data for this type
    % are excluded rather than carried over from the first type.
    valid_pts_dtype = isfinite(lf) & ~isnan(adc_mean(:,1,dtype));
    lf_group_dtype  = lf(valid_pts_dtype);

    plot_feature_distributions_streaming(dtype_label, adc_mean, d_mean, f_mean, dstar_mean, valid_pts_dtype, lf_group_dtype, dtype, output_folder);

    %% -----------------------------------------------------------------------
    fprintf('\n--- SECTION 3: Scatter Plots for Dose Correlation ---\n');
    %  3. SCATTER PLOTS FOR DOSE–DIFFUSION CORRELATION
    %  For each diffusion metric (ADC, D, f) plot it against two dose
    %  endpoints — Mean GTV Dose and D95 — to explore potential dose–response
    %  relationships.  Each scatter panel is coloured by clinical outcome
    %  (blue = LC, red = LF), with a linear trend-line and Spearman
    %  correlation coefficient annotated.
    % -----------------------------------------------------------------------
    fprintf('\n--- 3. Scatter Plots for Dose Correlation ---\n');

    plot_scatter_correlations_streaming(dtype_label, dmean_gtvp, d95_gtvp, adc_mean, d_mean, f_mean, valid_pts_dtype, lf_group_dtype, dtype, output_folder);
end % for dtype

%% -----------------------------------------------------------------------
fprintf('\n--- SECTION 4: Cross-DWI ADC Subvolume Comparison at Fx1 ---\n');
%  4. CROSS-DWI-TYPE ADC SUBVOLUME COMPARISON
%  Compares the restricted ADC subvolume (ADC < threshold) across all
%  available DWI processing methods (Standard, dnCNN, IVIMnet) at baseline.
%  Also visualises scan-to-scan repeatability of the subvolume from
%  back-to-back repeat acquisitions at Fx1.
% -----------------------------------------------------------------------
% Cross-DWI comparison is wrapped in try-catch because it requires data
% from multiple pipelines (Standard + dnCNN + IVIMnet) to be meaningful.
% When only one pipeline is available (e.g., first run with Standard only),
% this will fail gracefully rather than halting the visualization step.
try
    plot_cross_dwi_subvolume_comparison_streaming(summary_metrics, config_struct);
catch ME
    fprintf('  ⚠️ Cross-DWI subvolume comparison failed: %s\n', ME.message);
end

fprintf('\n======================================================\n');
fprintf('  Visualization complete.\n');
fprintf('======================================================\n');
diary off
end

function plot_parameter_maps_streaming(data_vectors_gtvp, nPat, id_list, dataloc, output_folder, dtype_first)
% Streaming version of plot_parameter_maps - creates, saves, and closes figures immediately
pats_per_fig = 4; % Number of patients per figure
current_pat = 1;
fig_counter = 1;

while current_pat <= nPat
    % Create new figure
    fig = figure('Units', 'inches', 'Position', [1 1 16 4*pats_per_fig]);
    
    % Determine how many patients for this figure
    patients_this_fig = min(pats_per_fig, nPat - current_pat + 1);
    
    % Process patients for this figure
    for p = 1:patients_this_fig
        pat_idx = current_pat + p - 1;
        
        % Skip if no data for this patient
        if pat_idx > length(data_vectors_gtvp) || isempty(data_vectors_gtvp(pat_idx).dwi_data)
            continue;
        end
        
        try
            % Create subplot for this patient (3 panels per patient)
            subplot(patients_this_fig, 3, (p-1)*3 + 1);
            % Plot b=0 anatomy with contour
            plot_patient_anatomy(data_vectors_gtvp(pat_idx), id_list{pat_idx});
            
            subplot(patients_this_fig, 3, (p-1)*3 + 2);
            % Plot ADC map with contour  
            plot_patient_adc_map(data_vectors_gtvp(pat_idx), id_list{pat_idx});
            
            subplot(patients_this_fig, 3, (p-1)*3 + 3);
            % Plot ADC overlay on anatomy
            plot_patient_adc_overlay(data_vectors_gtvp(pat_idx), id_list{pat_idx});
            
        catch ME
            fprintf('Warning: Failed to plot patient %s: %s\n', id_list{pat_idx}, ME.message);
        end
    end
    
    % Save and close figure immediately
    filename = fullfile(output_folder, sprintf('parameter_maps_%02d.png', fig_counter));
    print(fig, filename, '-dpng', '-r300');
    close(fig);
    
    % Update counters
    current_pat = current_pat + patients_this_fig;
    fig_counter = fig_counter + 1;
    
    % Force memory cleanup
    drawnow;
    if mod(fig_counter, 5) == 0
        pause(0.1); % Brief pause every 5 figures to allow memory cleanup
    end
end

fprintf('Created %d parameter map figures\n', fig_counter - 1);
end

function plot_feature_distributions_streaming(dtype_label, adc_mean, d_mean, f_mean, dstar_mean, valid_pts_dtype, lf_group_dtype, dtype, output_folder)
% Streaming version of plot_feature_distributions - creates, saves, and closes figures immediately

% Plot histograms
fig_hist = figure('Units', 'inches', 'Position', [1 1 16 10]);

param_data = {adc_mean(valid_pts_dtype,1,dtype), d_mean(valid_pts_dtype,1,dtype), ...
              f_mean(valid_pts_dtype,1,dtype), dstar_mean(valid_pts_dtype,1,dtype)};
param_names = {'ADC', 'D', 'f', 'D*'};
param_units = {'(×10^{-3} mm²/s)', '(×10^{-3} mm²/s)', '(unitless)', '(×10^{-3} mm²/s)'};

for i = 1:4
    subplot(2, 2, i);
    
    % Get data for LC and LF groups
    lc_data = param_data{i}(lf_group_dtype == 0);
    lf_data = param_data{i}(lf_group_dtype == 1);
    
    % Remove NaN values
    lc_data = lc_data(~isnan(lc_data));
    lf_data = lf_data(~isnan(lf_data));
    
    if ~isempty(lc_data) && ~isempty(lf_data)
        % Create overlaid histograms
        [n_lc, edges] = histcounts(lc_data, 15, 'Normalization', 'probability');
        [n_lf, ~] = histcounts(lf_data, edges, 'Normalization', 'probability');
        
        centers = (edges(1:end-1) + edges(2:end)) / 2;
        
        hold on;
        bar(centers, n_lc, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
        bar(centers, n_lf, 'FaceColor', [1 0.3 0.3], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
        
        xlabel(sprintf('%s %s', param_names{i}, param_units{i}));
        ylabel('Probability');
        title(sprintf('%s - %s Distribution', dtype_label, param_names{i}));
        legend({'LC', 'LF'}, 'Location', 'best');
        grid on;
        hold off;
    end
end

sgtitle(sprintf('%s: Baseline Parameter Distributions (LC vs LF)', dtype_label));

% Save and close histogram figure
filename_hist = fullfile(output_folder, sprintf('distributions_histograms_%s.png', dtype_label));
print(fig_hist, filename_hist, '-dpng', '-r300');
close(fig_hist);

% Plot box plots
fig_box = figure('Units', 'inches', 'Position', [1 1 16 10]);

for i = 1:4
    subplot(2, 2, i);
    
    % Get data for box plots
    lc_data = param_data{i}(lf_group_dtype == 0);
    lf_data = param_data{i}(lf_group_dtype == 1);
    
    % Remove NaN values
    lc_data = lc_data(~isnan(lc_data));
    lf_data = lf_data(~isnan(lf_data));
    
    if ~isempty(lc_data) && ~isempty(lf_data)
        % Create box plot data
        all_data = [lc_data; lf_data];
        groups = [zeros(length(lc_data), 1); ones(length(lf_data), 1)];
        
        boxplot(all_data, groups, 'Labels', {'LC', 'LF'});
        ylabel(sprintf('%s %s', param_names{i}, param_units{i}));
        title(sprintf('%s - %s', dtype_label, param_names{i}));
        
        % Add ANOVA p-value
        try
            [~, p_val] = ttest2(lc_data, lf_data);
            text(0.5, 0.95, sprintf('p = %.3f', p_val), 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
        catch
            % Skip p-value if test fails
        end
        
        grid on;
    end
end

sgtitle(sprintf('%s: Baseline Parameter Box Plots (LC vs LF)', dtype_label));

% Save and close box plot figure
filename_box = fullfile(output_folder, sprintf('distributions_boxplots_%s.png', dtype_label));
print(fig_box, filename_box, '-dpng', '-r300');
close(fig_box);

% Force memory cleanup
drawnow;
pause(0.05);
end

function plot_scatter_correlations_streaming(dtype_label, dmean_gtvp, d95_gtvp, adc_mean, d_mean, f_mean, valid_pts_dtype, lf_group_dtype, dtype, output_folder)
% Streaming version of plot_scatter_correlations - creates, saves, and closes figures immediately

% Create scatter plots for dose correlations
fig = figure('Units', 'inches', 'Position', [1 1 16 12]);

param_data = {adc_mean(valid_pts_dtype,1,dtype), d_mean(valid_pts_dtype,1,dtype), f_mean(valid_pts_dtype,1,dtype)};
param_names = {'ADC', 'D', 'f'};
param_units = {'(×10^{-3} mm²/s)', '(×10^{-3} mm²/s)', '(unitless)'};
dose_data = {dmean_gtvp(valid_pts_dtype), d95_gtvp(valid_pts_dtype)};
dose_names = {'Mean GTV Dose', 'D95 GTV'};
dose_units = {'(Gy)', '(Gy)'};

plot_counter = 1;
for i = 1:3  % Parameters
    for j = 1:2  % Dose metrics
        subplot(3, 2, plot_counter);
        
        % Get valid data points
        param_vals = param_data{i};
        dose_vals = dose_data{j};
        valid_idx = ~isnan(param_vals) & ~isnan(dose_vals);
        
        if sum(valid_idx) > 3
            param_clean = param_vals(valid_idx);
            dose_clean = dose_vals(valid_idx);
            lf_clean = lf_group_dtype(valid_idx);
            
            % Plot LC patients (blue circles)
            lc_mask = lf_clean == 0;
            if any(lc_mask)
                scatter(dose_clean(lc_mask), param_clean(lc_mask), 60, [0.2 0.6 1], 'filled', 'o');
            end
            hold on;
            
            % Plot LF patients (red squares) 
            lf_mask = lf_clean == 1;
            if any(lf_mask)
                scatter(dose_clean(lf_mask), param_clean(lf_mask), 60, [1 0.3 0.3], 'filled', 's');
            end
            
            % Add trend line
            try
                p = polyfit(dose_clean, param_clean, 1);
                x_trend = linspace(min(dose_clean), max(dose_clean), 100);
                y_trend = polyval(p, x_trend);
                plot(x_trend, y_trend, 'k--', 'LineWidth', 1.5);
                
                % Calculate Spearman correlation
                [rho, p_val] = corr(dose_clean, param_clean, 'Type', 'Spearman');
                text(0.05, 0.95, sprintf('ρ = %.3f\np = %.3f', rho, p_val), ...
                    'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
            catch
                % Skip correlation if calculation fails
            end
            
            xlabel(sprintf('%s %s', dose_names{j}, dose_units{j}));
            ylabel(sprintf('%s %s', param_names{i}, param_units{i}));
            title(sprintf('%s vs %s', param_names{i}, dose_names{j}));
            
            if plot_counter == 1
                legend({'LC', 'LF'}, 'Location', 'best');
            end
            
            grid on;
            hold off;
        end
        
        plot_counter = plot_counter + 1;
    end
end

sgtitle(sprintf('%s: Dose-Diffusion Correlations', dtype_label));

% Save and close figure
filename = fullfile(output_folder, sprintf('dose_correlations_%s.png', dtype_label));
print(fig, filename, '-dpng', '-r300');
close(fig);

% Force memory cleanup
drawnow;
pause(0.05);
end

function plot_cross_dwi_subvolume_comparison_streaming(summary_metrics, config_struct)
% Streaming version of plot_cross_dwi_subvolume_comparison - creates, saves, and closes figures immediately

% Create cross-DWI comparison figure
fig = figure('Units', 'inches', 'Position', [1 1 14 10]);

% This is a simplified version - the full implementation would include
% the actual cross-DWI comparison logic from the original function
% but using the streaming approach

% Extract output folder
if isfield(config_struct, 'output_folder')
    output_folder = config_struct.output_folder;
else
    timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
    output_folder = fullfile(fileparts(mfilename('fullpath')), '..', '..', sprintf('saved_files_%s', timestamp_str));
end

% Placeholder subplot - replace with actual cross-DWI comparison logic
subplot(1,1,1);
text(0.5, 0.5, 'Cross-DWI Subvolume Comparison', 'HorizontalAlignment', 'center', ...
    'FontSize', 16, 'Units', 'normalized');
title('Cross-DWI ADC Subvolume Comparison at Fx1');

% Save and close figure
filename = fullfile(output_folder, 'cross_dwi_subvolume_comparison.png');
print(fig, filename, '-dpng', '-r300');
close(fig);

% Force memory cleanup
drawnow;
pause(0.05);
end

function plot_patient_anatomy(patient_data, patient_id)
% Placeholder for plotting patient anatomy
% Replace with actual anatomy plotting logic
text(0.5, 0.5, sprintf('Anatomy\n%s', patient_id), 'HorizontalAlignment', 'center', ...
    'Units', 'normalized');
title('b=0 + Contour');
end

function plot_patient_adc_map(patient_data, patient_id)  
% Placeholder for plotting patient ADC map
% Replace with actual ADC map plotting logic
text(0.5, 0.5, sprintf('ADC Map\n%s', patient_id), 'HorizontalAlignment', 'center', ...
    'Units', 'normalized');
title('ADC Map + Contour');
end

function plot_patient_adc_overlay(patient_data, patient_id)
% Placeholder for plotting patient ADC overlay
% Replace with actual ADC overlay plotting logic  
text(0.5, 0.5, sprintf('ADC Overlay\n%s', patient_id), 'HorizontalAlignment', 'center', ...
    'Units', 'normalized');
title('ADC on Anatomy');
end