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
    output_folder = fullfile(fileparts(mfilename('fullpath')), '..', sprintf('saved_files_%s', timestamp_str));
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

fprintf('\n======================================================\n');
fprintf('  Visualization complete.\n');
fprintf('======================================================\n');
diary off
end
