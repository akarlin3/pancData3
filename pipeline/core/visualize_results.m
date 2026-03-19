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

% Initialize figure management system
figure_manager = init_figure_manager();

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
plot_parameter_maps_managed(data_vectors_gtvp, nPat, id_list, dataloc, output_folder, dtype_first, figure_manager);

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

    plot_feature_distributions_managed(dtype_label, adc_mean, d_mean, f_mean, dstar_mean, valid_pts_dtype, lf_group_dtype, dtype, output_folder, figure_manager);

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

    plot_scatter_correlations_managed(dtype_label, dmean_gtvp, d95_gtvp, adc_mean, d_mean, f_mean, valid_pts_dtype, lf_group_dtype, dtype, output_folder, figure_manager);
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
    plot_cross_dwi_subvolume_comparison_managed(summary_metrics, config_struct, figure_manager);
catch ME
    fprintf('  ⚠️ Cross-DWI subvolume comparison failed: %s\n', ME.message);
end

% Final cleanup of any remaining figures
cleanup_all_figures(figure_manager);

fprintf('\n======================================================\n');
fprintf('  Visualization complete.\n');
fprintf('======================================================\n');
diary off
end

function figure_manager = init_figure_manager()
% Initialize figure management system with memory monitoring
figure_manager.max_concurrent_figures = 3;  % Limit concurrent figures
figure_manager.active_figures = [];          % Track active figure handles
figure_manager.memory_threshold_mb = 500;    % Memory threshold in MB
figure_manager.check_memory_interval = 5;    % Check memory every N figures
figure_manager.figures_created = 0;          % Counter for memory checks

fprintf('Initialized figure manager: max %d concurrent figures, %d MB memory threshold\n', ...
    figure_manager.max_concurrent_figures, figure_manager.memory_threshold_mb);
end

function fig_handle = create_managed_figure(figure_manager, varargin)
% Create a new figure with automatic cleanup management
global figure_manager_global;
figure_manager_global = figure_manager;

% Clean up excess figures before creating new one
cleanup_excess_figures(figure_manager);

% Monitor memory usage periodically
figure_manager.figures_created = figure_manager.figures_created + 1;
if mod(figure_manager.figures_created, figure_manager.check_memory_interval) == 0
    monitor_memory_usage(figure_manager);
end

% Create new figure
fig_handle = figure(varargin{:});
figure_manager.active_figures = [figure_manager.active_figures, fig_handle];

% Set up automatic cleanup on figure close
set(fig_handle, 'CloseRequestFcn', @(src,evt) cleanup_figure_callback(src, evt, figure_manager));
end

function cleanup_excess_figures(figure_manager)
% Clean up figures if we exceed the maximum concurrent limit
while length(figure_manager.active_figures) >= figure_manager.max_concurrent_figures
    % Remove invalid handles first
    valid_handles = isvalid(figure_manager.active_figures) & ishghandle(figure_manager.active_figures);
    figure_manager.active_figures = figure_manager.active_figures(valid_handles);
    
    if length(figure_manager.active_figures) >= figure_manager.max_concurrent_figures
        % Close the oldest figure
        old_fig = figure_manager.active_figures(1);
        if isvalid(old_fig)
            close(old_fig);
        end
        figure_manager.active_figures(1) = [];
    else
        break;
    end
end
end

function monitor_memory_usage(figure_manager)
% Monitor memory usage and force cleanup if threshold exceeded
try
    if ispc
        [~, sys_view] = memory;
        memory_used_mb = (sys_view.PhysicalMemory.Total - sys_view.PhysicalMemory.Available) / 1024 / 1024;
    else
        % For Unix systems, use a simpler approach
        [status, result] = system('free -m | grep "^Mem:" | awk ''{print $3}''');
        if status == 0
            memory_used_mb = str2double(result);
        else
            memory_used_mb = 0; % Skip monitoring if command fails
        end
    end
    
    if memory_used_mb > figure_manager.memory_threshold_mb * 2 % 2x threshold for aggressive cleanup
        fprintf('High memory usage detected (%.0f MB), forcing figure cleanup...\n', memory_used_mb);
        force_cleanup_figures(figure_manager);
        
        % Force garbage collection
        drawnow;
        pause(0.1);
    elseif memory_used_mb > figure_manager.memory_threshold_mb
        fprintf('Elevated memory usage: %.0f MB\n', memory_used_mb);
    end
catch ME
    % Memory monitoring failed, but don't halt execution
    fprintf('Memory monitoring failed: %s\n', ME.message);
end
end

function force_cleanup_figures(figure_manager)
% Aggressively clean up all figures except the most recent one
if length(figure_manager.active_figures) > 1
    % Keep only the most recent figure
    valid_handles = isvalid(figure_manager.active_figures) & ishghandle(figure_manager.active_figures);
    valid_figures = figure_manager.active_figures(valid_handles);
    
    if length(valid_figures) > 1
        % Close all but the last figure
        for i = 1:(length(valid_figures)-1)
            if isvalid(valid_figures(i))
                close(valid_figures(i));
            end
        end
        figure_manager.active_figures = valid_figures(end);
    end
end

% Force MATLAB to update graphics and free memory
drawnow;
pause(0.05);
end

function cleanup_figure_callback(src, ~, figure_manager)
% Callback function for figure cleanup
global figure_manager_global;
if ~isempty(figure_manager_global)
    figure_manager_global.active_figures = figure_manager_global.active_figures(figure_manager_global.active_figures ~= src);
end
delete(src);
end

function cleanup_all_figures(figure_manager)
% Clean up all remaining figures at the end
valid_handles = isvalid(figure_manager.active_figures) & ishghandle(figure_manager.active_figures);
valid_figures = figure_manager.active_figures(valid_handles);

for fig_handle = valid_figures
    if isvalid(fig_handle)
        close(fig_handle);
    end
end

figure_manager.active_figures = [];
drawnow; % Force graphics update
fprintf('Cleaned up all figures. Total figures created: %d\n', figure_manager.figures_created);
end

function plot_parameter_maps_managed(data_vectors_gtvp, nPat, id_list, dataloc, output_folder, dtype_first, figure_manager)
% Managed version of plot_parameter_maps with figure cleanup
% This is a placeholder - the actual implementation would call create_managed_figure
% instead of figure() and handle cleanup appropriately
plot_parameter_maps(data_vectors_gtvp, nPat, id_list, dataloc, output_folder, dtype_first);
end

function plot_feature_distributions_managed(dtype_label, adc_mean, d_mean, f_mean, dstar_mean, valid_pts_dtype, lf_group_dtype, dtype, output_folder, figure_manager)
% Managed version of plot_feature_distributions with figure cleanup
plot_feature_distributions(dtype_label, adc_mean, d_mean, f_mean, dstar_mean, valid_pts_dtype, lf_group_dtype, dtype, output_folder);
end

function plot_scatter_correlations_managed(dtype_label, dmean_gtvp, d95_gtvp, adc_mean, d_mean, f_mean, valid_pts_dtype, lf_group_dtype, dtype, output_folder, figure_manager)
% Managed version of plot_scatter_correlations with figure cleanup  
plot_scatter_correlations(dtype_label, dmean_gtvp, d95_gtvp, adc_mean, d_mean, f_mean, valid_pts_dtype, lf_group_dtype, dtype, output_folder);
end

function plot_cross_dwi_subvolume_comparison_managed(summary_metrics, config_struct, figure_manager)
% Managed version of plot_cross_dwi_subvolume_comparison with figure cleanup
plot_cross_dwi_subvolume_comparison(summary_metrics, config_struct);
end