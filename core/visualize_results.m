function visualize_results(data_vectors_gtvp, summary_metrics, calculated_results, config_struct)
% VISUALIZE_RESULTS — "Visualizing It" (Generating Plots)
% Author: Avery Karlin
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
fprintf('\n--- 1. Parameter Maps overlaid on Anatomy ---\n');

% Track how many patients have been plotted so far
patients_plotted = 0;
% Maximum number of patient rows per figure panel
pats_per_fig     = 5;

% Pre-count eligible patients (those with Fx1 DWI, bval, and GTV files)
% so that each figure can be allocated the correct number of subplot rows.
n_eligible = 0;
for j = 1:nPat
    if j > size(data_vectors_gtvp, 1)
        continue;
    end
    s = data_vectors_gtvp(j, 1, 1);
    if isempty(s.adc_vector), continue; end
    basefolder = fullfile(dataloc, id_list{j});
    nii_path   = fullfile(basefolder, 'nii');
    if exist(fullfile(nii_path, 'fx1_dwi1.nii.gz'), 'file') && ...
       exist(fullfile(nii_path, 'fx1_dwi1.bval'),   'file') && ...
       exist(fullfile(nii_path, 'fx1_gtv1.nii.gz'), 'file')
        n_eligible = n_eligible + 1;
    end
end

n_rows_cur_fig = 0;  % will be set when each new figure is created

for j = 1:nPat
    if j > size(data_vectors_gtvp, 1)
        continue;
    end
    % Skip patients without extracted ADC data at Fx1
    s = data_vectors_gtvp(j, 1, 1);
    if isempty(s.adc_vector), continue; end

    % Build file paths for the Fx1 DWI volume, b-value table, and GTV mask
    basefolder = fullfile(dataloc, id_list{j});
    nii_path   = fullfile(basefolder, 'nii');
    dwi_file   = fullfile(nii_path, 'fx1_dwi1.nii.gz');
    bval_file  = fullfile(nii_path, 'fx1_dwi1.bval');
    gtv_file   = fullfile(nii_path, 'fx1_gtv1.nii.gz');

    % All three files must exist on disk to proceed
    if ~exist(dwi_file, 'file') || ~exist(bval_file, 'file') || ~exist(gtv_file, 'file')
        continue;
    end

    % Load the 4-D DWI NIfTI volume and the 3-D GTV binary mask.
    % rot90 reorients from NIfTI radiological convention to display coords.
    dwi_info = niftiinfo(dwi_file);
    dwi_img = double(rot90(niftiread(dwi_info)));
    gtv_info = niftiinfo(gtv_file);
    gtv_img = double(rot90(niftiread(gtv_info)));

    % Read b-values from the accompanying text file (one line, space-delimited)
    fid = fopen(bval_file); bvals = sscanf(fgetl(fid), '%f')'; fclose(fid);
    bvals = bvals(:);

    % Validate b-values against expected study protocol
    expected_bvals = [0; 30; 150; 550];
    if size(dwi_img, 4) ~= length(bvals) || ~isequal(sort(bvals), expected_bvals)
        fprintf('  Protocol deviation: Pt %d has non-standard b-values %s — excluding from comparative mapping.\n', ...
            j, mat2str(bvals'));
        continue;
    end

    % Compute a voxel-wise ADC map using a monoexponential fit:
    %   S(b) = S0 * exp(-b * ADC)
    % Clamp non-physical values to the range [0, 3e-3] mm^2/s.
    adc_map = fit_adc_mono(dwi_img, bvals);
    adc_map(adc_map < 0)    = 0;
    adc_map(adc_map > 3e-3) = 3e-3;

    % Pick the axial slice containing the largest cross-sectional GTV area
    gtv_areas = squeeze(sum(sum(gtv_img, 1), 2));
    [~, z_slice] = max(gtv_areas);

    % Extract 2-D slices for display
    b0_slice  = squeeze(dwi_img(:,:,z_slice,1));   % b=0 anatomical image
    adc_slice = squeeze(adc_map(:,:,z_slice));      % ADC map
    gtv_slice = squeeze(gtv_img(:,:,z_slice));      % GTV mask

    patients_plotted = patients_plotted + 1;
    row_in_fig = mod(patients_plotted - 1, pats_per_fig);
    fig_num    = ceil(patients_plotted / pats_per_fig);

    % --- Start a new figure every pats_per_fig patients ---
    if row_in_fig == 0
        % Finalise and save the previous figure (if any)
        if patients_plotted > 1
            sgtitle('Parameter Maps Overlaid on Anatomy (Fx1)', 'FontSize', 14, 'FontWeight', 'bold');
            cb = colorbar('Position', [0.93 0.11 0.015 0.8]);
            ylabel(cb, 'ADC (mm^2/s)');
            saveas(gcf, fullfile(output_folder, sprintf('Parameter_Maps_%d.png', fig_num - 1)));
            close(gcf);
        end
        % Determine how many rows this new figure needs
        n_rows_cur_fig = min(pats_per_fig, n_eligible - (fig_num - 1) * pats_per_fig);
        fig_height = max(300, n_rows_cur_fig * 150);
        figure('Name', sprintf('ADC Maps on Anatomy (%d)', fig_num), ...
               'Position', [50, 50, 1400, fig_height]);
        colormap(jet);
    end

    % --- Column 1: b=0 anatomy with GTV contour (red) ---
    subplot(n_rows_cur_fig, 3, row_in_fig*3 + 1);
    imagesc(b0_slice); axis image; axis off; colormap(gca, gray);
    hold on;
    contour(gtv_slice, [0.5 0.5], 'r', 'LineWidth', 1.5);
    hold off;
    title(sprintf('%s — b0 (GTV contour)', id_list{j}), 'Interpreter', 'none', 'FontSize', 9);

    % --- Column 2: ADC map with GTV contour (white) ---
    subplot(n_rows_cur_fig, 3, row_in_fig*3 + 2);
    imagesc(adc_slice, [0 2.5e-3]); axis image; axis off;
    hold on;
    contour(gtv_slice, [0.5 0.5], 'w', 'LineWidth', 1.5);
    hold off;
    title('ADC map', 'FontSize', 9);

    % --- Column 3: ADC overlaid on anatomy (inside GTV only) ---
    subplot(n_rows_cur_fig, 3, row_in_fig*3 + 3);
    % Normalise b=0 image to [0,1] so it serves as a greyscale underlay
    b0_norm = mat2gray(b0_slice);
    imshow(b0_norm, []); hold on;
    % Mask ADC outside the GTV to NaN so only tumour voxels are coloured
    adc_overlay = adc_slice;
    adc_overlay(gtv_slice < 0.5) = NaN;
    h_ov = imagesc(adc_overlay, [0 2.5e-3]);
    % Use 60 % opacity for the colour overlay
    set(h_ov, 'AlphaData', ~isnan(adc_overlay) * 0.6);
    contour(gtv_slice, [0.5 0.5], 'w', 'LineWidth', 1.5);
    hold off;
    title('ADC on Anatomy', 'FontSize', 9);
end

% Save and close the final parameter-maps figure
if patients_plotted > 0
    sgtitle('Parameter Maps Overlaid on Anatomy (Fx1)', 'FontSize', 14, 'FontWeight', 'bold');
    cb = colorbar('Position', [0.93 0.11 0.015 0.8]);
    ylabel(cb, 'ADC (mm^2/s)');
    saveas(gcf, fullfile(output_folder, sprintf('Parameter_Maps_%d.png', fig_num)));
    close(gcf);
    fprintf('  Plotted %d patients.\n', patients_plotted);
else
    fprintf('  No patients with complete NIfTI data found on disk. Skipping.\n');
end

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

    % Split values by clinical outcome
    vals_lc = vals(lf_group == 0);   % Local Control patients
    vals_lf = vals(lf_group == 1);   % Local Failure patients

    % Create 15 equally-spaced bins spanning the combined value range
    edges = linspace(min(vals, [], 'omitnan'), max(vals, [], 'omitnan'), 16);

    % Overlay semi-transparent histograms for each outcome group
    histogram(vals_lc, edges, 'FaceColor', [0.2 0.4 0.8], 'FaceAlpha', 0.6, ...
        'EdgeColor', 'none', 'DisplayName', 'Local Control'); hold on;
    histogram(vals_lf, edges, 'FaceColor', [0.8 0.2 0.2], 'FaceAlpha', 0.6, ...
        'EdgeColor', 'none', 'DisplayName', 'Local Failure');
    hold off;

    xlabel(metric_units{mi}); ylabel('Count');
    title(metric_names{mi}, 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
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

    % Remove NaN entries so boxplot renders correctly
    has_data = ~isnan(vals);
    boxplot(vals(has_data), lf_group(has_data), 'Labels', {'LC (0)', 'LF (1)'});

    ylabel(metric_units{mi});
    title(metric_names{mi}, 'FontSize', 11, 'FontWeight', 'bold');
    grid on;

    % Annotate with a one-way ANOVA p-value comparing LC vs LF
    if sum(has_data) > 2 && numel(unique(lf_group(has_data))) > 1
        p = anova1(vals(has_data), lf_group(has_data), 'off');
        yl = ylim;
        text(1.5, yl(2)*0.95, sprintf('p = %.3f', p), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
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

% Extract Fx1 dose metrics for the valid patient subset
dose_mean_vec = dmean_gtvp(valid_pts, 1);   % mean dose inside GTV (Gy)
dose_d95_vec  = d95_gtvp(valid_pts, 1);     % D95 = dose to 95 % of GTV (Gy)

% Diffusion biomarkers to correlate with dose
diff_metrics = {adc_mean(valid_pts,1,dtype), d_mean(valid_pts,1,dtype), f_mean(valid_pts,1,dtype)};
diff_names   = {'Mean ADC', 'Mean D', 'Mean f'};
diff_units   = {'mm^2/s',   'mm^2/s',  ''};

figure('Name', ['Dose vs Diffusion Metrics — ' dtype_label], ...
       'Position', [150, 150, 1400, 500]);

plot_idx = 1;
for di = 1:numel(diff_metrics)
    y_vals = diff_metrics{di};

    % Loop over the two dose endpoints (mean dose, D95)
    for dose_type = 1:2
        if dose_type == 1
            x_vals = dose_mean_vec;
            x_label = 'Mean GTV Dose (Gy)';
        else
            x_vals = dose_d95_vec;
            x_label = 'D95 (Gy)';
        end

        subplot(2, numel(diff_metrics), plot_idx);

        % Identify patients with both dose and diffusion data available
        clean = ~isnan(x_vals) & ~isnan(y_vals);
        if sum(clean) < 3
            title([diff_names{di} ' — insufficient data']);
            plot_idx = plot_idx + 1;
            continue;
        end

        % Plot LC (blue) and LF (red) points with black edge
        scatter(x_vals(clean & lf_group==0), y_vals(clean & lf_group==0), ...
            50, 'b', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'LC'); hold on;
        scatter(x_vals(clean & lf_group==1), y_vals(clean & lf_group==1), ...
            50, 'r', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'LF');

        % Overlay a first-order (linear) polynomial fit
        p_fit = polyfit(x_vals(clean), y_vals(clean), 1);
        x_line = linspace(min(x_vals(clean)), max(x_vals(clean)), 50);
        plot(x_line, polyval(p_fit, x_line), 'k--', 'LineWidth', 1.5, ...
            'DisplayName', 'Linear fit');
        hold off;

        % Compute Spearman rank correlation and annotate the title
        [r_sp, p_sp] = corr(x_vals(clean), y_vals(clean), 'Type', 'Spearman');

        xlabel(x_label);
        ylabel([diff_names{di} ' (' diff_units{di} ')']);
        title(sprintf('%s vs Dose\nr_s=%.2f, p=%.3f', diff_names{di}, r_sp, p_sp), ...
            'FontSize', 10);
        legend('Location', 'best', 'FontSize', 7);
        grid on;

        plot_idx = plot_idx + 1;
    end
end
sgtitle(['RT Dose vs Diffusion Metrics (Fx1) (' dtype_label ')'], ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(output_folder, ['Dose_vs_Diffusion_' dtype_label '.png']));
close(gcf);

fprintf('  Scatter plots generated (%s).\n', dtype_label);

end % for dtype

fprintf('\n======================================================\n');
fprintf('  Visualization complete.\n');
fprintf('======================================================\n');
diary off
end