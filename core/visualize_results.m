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

% This function generates three families of visualizations:
%   1. Parameter Maps overlaid on Anatomy
%      – ADC maps computed from the Fx1 DWI volume, overlaid on the b=0
%        anatomical image with the GTV contour
%   2. Distributions of Extracted Features (histograms & box plots)
%      – Baseline (Fx1) ADC and IVIM parameters grouped by clinical
%        outcome (Local Control vs Local Failure)
%   3. Scatter Plots for Dose–Diffusion Correlation
%      – RT Dose (mean GTV dose and D95) plotted against diffusion
%        metrics to identify clusters, thresholds, or linear trends

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

% Define human-readable labels for the three DWI processing pipelines
dwi_type_names = {'Standard', 'dnCNN', 'IVIMnet'};

% Create output directory for saved figures (ignored by .gitignore)
if isfield(config_struct, 'output_folder')
    output_folder = config_struct.output_folder;
else
    output_folder = fullfile(pwd, 'saved_figures');
end
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

% Start diary to log text output
diary_file = fullfile(output_folder, 'visualize_results_output.txt');
if exist(diary_file, 'file'), delete(diary_file); end
diary(diary_file);

% Suppress figure windows so plots are only saved to disk
set(0, 'DefaultFigureVisible', 'off');

% Total number of patients and fraction (timepoint) labels
nPat = length(id_list);
fx_labels = {'Fx1','Fx2','Fx3','Fx4','Fx5','Post'};

% Build a logical mask identifying patients with usable clinical and
% imaging data.  Require a finite LF label and a
% non-NaN baseline ADC value (Standard DWI, Fx1).
valid_pts = isfinite(lf) & ~isnan(adc_mean(:,1,1));
% Subset the outcome labels to the valid patients for later grouping
lf_group = lf(valid_pts);

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
plot_parameter_maps(data_vectors_gtvp, nPat, id_list, dataloc, output_folder);

%% -----------------------------------------------------------------------
fprintf('\n--- SECTION 2: Distributions of Extracted Features ---\n');
%  2. DISTRIBUTIONS OF EXTRACTED FEATURES
%  For each DWI processing pipeline (Standard, dnCNN, IVIMnet):
%    2a. Histograms — overlay Local-Control vs Local-Failure distributions
%        for each of the four baseline (Fx1) diffusion biomarkers.
%    2b. Box plots  — side-by-side LC vs LF boxes with one-way ANOVA
%        p-values annotated on each panel.
%  The resulting figures are saved as PNG files in saved_figures/.
% -----------------------------------------------------------------------
fprintf('\n--- 2. Distributions of Extracted Features ---\n');

for dtype = config_struct.dwi_types_to_run
    dtype_label = dwi_type_names{dtype};

    plot_feature_distributions(dtype_label, adc_mean, d_mean, f_mean, dstar_mean, valid_pts, lf_group, dtype, output_folder);

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

    plot_scatter_correlations(dtype_label, dmean_gtvp, d95_gtvp, adc_mean, d_mean, f_mean, valid_pts, lf_group, dtype, output_folder);
end % for dtype

fprintf('\n======================================================\n');
fprintf('  Visualization complete.\n');
fprintf('======================================================\n');
diary off
end
