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
plot_parameter_maps(data_vectors_gtvp, nPat, id_list, dataloc, output_folder, dtype_first);

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

    plot_feature_distributions(dtype_label, adc_mean, d_mean, f_mean, dstar_mean, valid_pts_dtype, lf_group_dtype, dtype, output_folder);

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

    plot_scatter_correlations(dtype_label, dmean_gtvp, d95_gtvp, adc_mean, d_mean, f_mean, valid_pts_dtype, lf_group_dtype, dtype, output_folder);
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
    plot_cross_dwi_subvolume_comparison(summary_metrics, config_struct);
catch ME
    fprintf('  ⚠️ Cross-DWI subvolume comparison failed: %s\n', ME.message);
end

%% -----------------------------------------------------------------------
%  5. TRAJECTORY PLOTS (Waterfall, Swimmer, Spider)
% -----------------------------------------------------------------------
if isfield(config_struct, 'run_trajectory_plots') && config_struct.run_trajectory_plots
    fprintf('\n--- SECTION 5: Trajectory Plots (Waterfall / Swimmer / Spider) ---\n');
    try
        plot_waterfall_chart(calculated_results, config_struct);
    catch ME
        fprintf('  ⚠️ Waterfall chart failed: %s\n', ME.message);
    end
    try
        plot_swimmer_chart(calculated_results, config_struct);
    catch ME
        fprintf('  ⚠️ Swimmer chart failed: %s\n', ME.message);
    end
    try
        plot_spider_chart(calculated_results, config_struct);
    catch ME
        fprintf('  ⚠️ Spider chart failed: %s\n', ME.message);
    end
end

fprintf('\n======================================================\n');
fprintf('  Visualization complete.\n');
fprintf('======================================================\n');
diary off
end

% plot_parameter_maps_streaming removed — was using placeholder helpers
% that produced blank images. Now calls the full plot_parameter_maps.m
% which loads NIfTI volumes and renders actual ADC maps.

% plot_feature_distributions_streaming removed — was a simplified
% reimplementation that produced blank box plots. Now calls the full
% plot_feature_distributions.m which handles single-group cases, Octave
% compatibility, and proper statistical annotation.

% plot_scatter_correlations_streaming removed — was a simplified
% reimplementation that produced blank scatter plots. Now calls the full
% plot_scatter_correlations.m which handles per-group trend lines,
% competing risk exclusion, and Octave compatibility.

% plot_cross_dwi_subvolume_comparison_streaming removed — was a placeholder
% that produced blank graphs. Now calls the real utility directly (line 175).

% Placeholder helper functions (plot_patient_anatomy, plot_patient_adc_map,
% plot_patient_adc_overlay) removed — they produced blank images with only
% text labels. The full plot_parameter_maps.m handles all rendering.

function plot_waterfall_chart(calculated_results, config_struct)
% PLOT_WATERFALL_CHART  Sorted bar chart of best ADC response per patient.
%
%   Requires calculated_results.percent_deltas with fields:
%     patient_ids  — cell array of patient IDs
%     adc_mean     — vector of percent-change values (best response)

output_folder = config_struct.output_folder;
if ~isfield(calculated_results, 'percent_deltas')
    error('visualize_results:missingField', 'percent_deltas not found in calculated_results.');
end
pd = calculated_results.percent_deltas;
if ~isfield(pd, 'patient_ids') || ~isfield(pd, 'adc_mean')
    error('visualize_results:missingField', 'percent_deltas must contain patient_ids and adc_mean.');
end

pct = pd.adc_mean(:);
ids = pd.patient_ids(:);
valid = isfinite(pct);
pct = pct(valid);
ids = ids(valid);
if isempty(pct), return; end

[pct_sorted, idx] = sort(pct, 'ascend');
ids_sorted = ids(idx);

% Colour by response category
colors = zeros(numel(pct_sorted), 3);
for k = 1:numel(pct_sorted)
    if pct_sorted(k) <= -20
        colors(k, :) = [0.2 0.6 0.2];   % responder (green)
    elseif pct_sorted(k) >= 20
        colors(k, :) = [0.8 0.2 0.2];   % progressor (red)
    else
        colors(k, :) = [0.5 0.5 0.5];   % stable (grey)
    end
end

fig = figure('Visible', 'off', 'Position', [100 100 900 500]);
b = bar(pct_sorted, 'FaceColor', 'flat');
b.CData = colors;
hold on;
yline(-20, '--k', 'LineWidth', 1);
yline(20, '--k', 'LineWidth', 1);
hold off;
set(gca, 'XTick', 1:numel(ids_sorted), 'XTickLabel', ids_sorted, ...
    'XTickLabelRotation', 90, 'FontSize', 7);
ylabel('Best ADC Response (% change)');
title(sprintf('Waterfall Plot — ADC Response (%s)', config_struct.dwi_type));
grid on;

fname = fullfile(output_folder, 'waterfall_adc_response');
print(fig, fname, '-dpng', '-r300');
print(fig, fname, '-depsc');
close(fig);
drawnow; pause(0.05);
fprintf('  📁 Saved %s.png\n', fname);
end

function plot_swimmer_chart(calculated_results, config_struct)
% PLOT_SWIMMER_CHART  Horizontal timeline per patient with scan markers.
%
%   Requires calculated_results.scan_timeline with fields:
%     patient_ids  — cell array of patient IDs
%     scan_days    — cell array of numeric vectors (days from baseline)
%     follow_up    — vector of total follow-up days

output_folder = config_struct.output_folder;
if ~isfield(calculated_results, 'scan_timeline')
    error('visualize_results:missingField', 'scan_timeline not found in calculated_results.');
end
st = calculated_results.scan_timeline;
if ~isfield(st, 'patient_ids') || ~isfield(st, 'scan_days') || ~isfield(st, 'follow_up')
    error('visualize_results:missingField', 'scan_timeline must contain patient_ids, scan_days, and follow_up.');
end

ids = st.patient_ids(:);
follow_up = st.follow_up(:);
scan_days = st.scan_days(:);
nPat = numel(ids);
if nPat == 0, return; end

[~, order] = sort(follow_up, 'descend');

fig = figure('Visible', 'off', 'Position', [100 100 900 max(400, nPat * 18)]);
hold on;
for i = 1:nPat
    row = order(i);
    barh(i, follow_up(row), 0.6, 'FaceColor', [0.7 0.85 1.0], 'EdgeColor', 'none');
    days_vec = scan_days{row};
    if ~isempty(days_vec)
        plot(days_vec, repmat(i, size(days_vec)), 'k^', 'MarkerSize', 5, 'MarkerFaceColor', [0.2 0.4 0.8]);
    end
end
hold off;
set(gca, 'YTick', 1:nPat, 'YTickLabel', ids(order), 'FontSize', 7, 'YDir', 'reverse');
xlabel('Days from Baseline');
title(sprintf('Swimmer Plot — Patient Timelines (%s)', config_struct.dwi_type));

fname = fullfile(output_folder, 'swimmer_timeline');
print(fig, fname, '-dpng', '-r300');
print(fig, fname, '-depsc');
close(fig);
drawnow; pause(0.05);
fprintf('  📁 Saved %s.png\n', fname);
end

function plot_spider_chart(calculated_results, config_struct)
% PLOT_SPIDER_CHART  Per-patient ADC trajectories normalised to cohort median baseline.
%
%   Requires calculated_results.longitudinal_trajectories with fields:
%     patient_ids  — cell array of patient IDs
%     timepoints   — numeric vector of common timepoints (e.g. scan days)
%     adc_values   — nPatients x nTimepoints matrix of ADC values

output_folder = config_struct.output_folder;
if ~isfield(calculated_results, 'longitudinal_trajectories')
    error('visualize_results:missingField', 'longitudinal_trajectories not found in calculated_results.');
end
lt = calculated_results.longitudinal_trajectories;
if ~isfield(lt, 'patient_ids') || ~isfield(lt, 'timepoints') || ~isfield(lt, 'adc_values')
    error('visualize_results:missingField', 'longitudinal_trajectories must contain patient_ids, timepoints, and adc_values.');
end

tp = lt.timepoints(:)';
vals = lt.adc_values;  % nPat x nTP
nPat = size(vals, 1);
if nPat == 0 || isempty(tp), return; end

% Normalise: percent change from cohort median at baseline
baseline_median = nanmedian(vals(:, 1));
if baseline_median == 0 || isnan(baseline_median)
    baseline_median = 1;  % avoid division by zero
end
pct_change = ((vals - baseline_median) ./ baseline_median) * 100;

fig = figure('Visible', 'off', 'Position', [100 100 800 500]);
hold on;
cmap = lines(min(nPat, 64));
for k = 1:nPat
    cidx = mod(k - 1, size(cmap, 1)) + 1;
    plot(tp, pct_change(k, :), '-o', 'Color', cmap(cidx, :), ...
        'MarkerSize', 4, 'LineWidth', 1);
end
yline(0, '--k', 'LineWidth', 1);
hold off;
xlabel('Timepoint (days)');
ylabel('ADC Change from Cohort Median (%)');
title(sprintf('Spider Plot — Longitudinal ADC Trajectories (%s)', config_struct.dwi_type));
grid on;

fname = fullfile(output_folder, 'spider_adc_trajectories');
print(fig, fname, '-dpng', '-r300');
print(fig, fname, '-depsc');
close(fig);
drawnow; pause(0.05);
fprintf('  📁 Saved %s.png\n', fname);
end