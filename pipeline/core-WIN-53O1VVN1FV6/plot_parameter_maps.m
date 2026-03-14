function plot_parameter_maps(data_vectors_gtvp, nPat, id_list, dataloc, output_folder, dtype)
% PLOT_PARAMETER_MAPS — "Parameter Maps overlaid on Anatomy"
%
% ANALYTICAL OVERVIEW:
%   Generates spatial visualisations of the ADC map overlaid on anatomical
%   DWI images for quality assurance and clinical interpretation.  These
%   maps show WHERE within the tumour diffusion is restricted (low ADC =
%   dark blue) vs. unrestricted (high ADC = warm colours), enabling:
%
%   1. QUALITY ASSURANCE — Verify that the GTV contour correctly delineates
%      the tumour on the b=0 image.  Misregistration between the contour
%      and anatomy indicates a segmentation or co-registration error that
%      would corrupt all downstream voxel-level analyses.
%
%   2. SPATIAL HETEROGENEITY ASSESSMENT — Tumours with spatially uniform
%      ADC may respond differently to RT than those with focal low-ADC
%      regions (potential treatment-resistant subclones).  The overlay
%      (panel c) reveals this heterogeneity within anatomical context.
%
%   3. PROTOCOL COMPLIANCE — B-value validation ensures each patient was
%      scanned with the expected protocol (b = 0, 30, 150, 550 s/mm^2).
%      Deviating protocols invalidate the monoexponential ADC fit because
%      the b-value range determines the sensitivity to diffusion vs.
%      perfusion contributions (low b < 200 captures perfusion effects).
%
%   The ADC map is computed via weighted least squares (WLS) monoexponential
%   fit: S(b) = S0 * exp(-b * ADC), consistent with fit_models.m.  WLS
%   weights by S^2 to account for the heteroscedastic noise in the log-
%   transformed signal (Rician noise becomes heteroscedastic after log
%   transformation).
%
%  For each patient with Fx1 data available, load the DWI volume, compute
%  the ADC map via monoexponential fit, and display three panels per
%  patient:
%    (a) b=0 anatomical image with GTV contour overlay
%    (b) Full-slice ADC map with GTV contour
%    (c) ADC overlaid on anatomy (semi-transparent inside GTV only)
%  Patients are batched into multi-row figures (pats_per_fig rows each).

if nargin < 6 || isempty(dtype)
    dtype = 1;  % default to Standard (1=Standard, 2=DnCNN, 3=IVIMnet)
end

fprintf('\n--- 1. Parameter Maps overlaid on Anatomy ---\n');

% Track how many patients have been plotted so far (across all figures)
patients_plotted = 0;
% Maximum number of patient rows per figure panel. Each row contains 3
% subplots (b0, ADC map, ADC overlay), so a 5-row figure has 15 subplots.
pats_per_fig     = 5;

% Expected b-value protocol for validation.
% This 4-point protocol is designed for IVIM analysis:
%   b=0    — reference image (no diffusion weighting), provides S0
%   b=30   — low b-value, captures perfusion (pseudo-diffusion) contribution
%   b=150  — intermediate b-value, transition region between perfusion and diffusion
%   b=550  — high b-value, dominated by true tissue diffusion (f contribution negligible)
% The ratio of high-to-low b-values determines the sensitivity of the
% IVIM fit to separate D (true diffusion) from f*D* (perfusion).
expected_bvals = [0; 30; 150; 550];

% Diagnostic counters for skip reasons
diag_out_of_range = 0;
diag_empty_adc    = 0;
diag_missing_file = 0;
diag_bad_bval     = 0;

% Pre-count eligible patients (those with Fx1 DWI, bval, and GTV files)
% so that each figure can be allocated the correct number of subplot rows.
n_eligible = 0;
for j = 1:nPat
    if j > size(data_vectors_gtvp, 1)
        continue;
    end
    if dtype > size(data_vectors_gtvp, 3), continue; end
    s = data_vectors_gtvp(j, 1, dtype);
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
    text_progress_bar(j, nPat, 'Generating parameter maps');
    if j > size(data_vectors_gtvp, 1)
        diag_out_of_range = diag_out_of_range + 1;
        fprintf('  💡 Pt %d: index exceeds data_vectors_gtvp rows (%d)\n', j, size(data_vectors_gtvp, 1));
        continue;
    end
    % Skip patients whose data lacks the requested DWI type
    if dtype > size(data_vectors_gtvp, 3)
        diag_empty_adc = diag_empty_adc + 1;
        fprintf('  💡 Pt %d (%s): dtype %d unavailable (max %d) — skipping\n', j, id_list{j}, dtype, size(data_vectors_gtvp, 3));
        continue;
    end
    % Skip patients without extracted ADC data at Fx1
    s = data_vectors_gtvp(j, 1, dtype);
    if isempty(s.adc_vector)
        diag_empty_adc = diag_empty_adc + 1;
        fprintf('  💡 Pt %d (%s): empty adc_vector at Fx1 — skipping\n', j, id_list{j});
        continue;
    end

    % Build file paths for the Fx1 DWI volume, b-value table, and GTV mask
    basefolder = fullfile(dataloc, id_list{j});
    nii_path   = fullfile(basefolder, 'nii');
    dwi_file   = fullfile(nii_path, 'fx1_dwi1.nii.gz');
    bval_file  = fullfile(nii_path, 'fx1_dwi1.bval');
    gtv_file   = fullfile(nii_path, 'fx1_gtv1.nii.gz');

    % All three files must exist on disk to proceed
    if ~exist(dwi_file, 'file') || ~exist(bval_file, 'file') || ~exist(gtv_file, 'file')
        diag_missing_file = diag_missing_file + 1;
        fprintf('  💡 Pt %d (%s): missing NIfTI — dwi=%d bval=%d gtv=%d\n', ...
            j, id_list{j}, exist(dwi_file,'file')~=0, exist(bval_file,'file')~=0, exist(gtv_file,'file')~=0);
        continue;
    end

    % Load the 4-D DWI NIfTI volume (x, y, z, b-value) and the 3-D GTV binary mask.
    % rot90 reorients from NIfTI radiological convention to MATLAB display coords.
    % double() ensures consistent arithmetic (NIfTI may store as int16/uint16).
    dwi_info = niftiinfo(dwi_file);
    dwi_img = double(rot90(niftiread(dwi_info)));
    gtv_info = niftiinfo(gtv_file);
    gtv_img = double(rot90(niftiread(gtv_info)));

    % Read b-values from the accompanying text file (one line, space-delimited).
    % Force to column vector for consistent indexing below.
    fid = fopen(bval_file); bvals = sscanf(fgetl(fid), '%f')'; fclose(fid);
    bvals = bvals(:);

    % B-value validation: skip patients whose b-values do not match
    % the expected acquisition protocol.
    if ~isequal(sort(bvals), sort(expected_bvals))
        diag_bad_bval = diag_bad_bval + 1;
        fprintf('  ⚠️  Pt %d (%s): Protocol Deviation — b-values %s do not match expected %s — skipping\n', ...
            j, id_list{j}, mat2str(bvals'), mat2str(expected_bvals'));
        continue;
    end

    % Validate that the 4-D volume has one frame per b-value and that
    % a b=0 image is present (required for monoexponential ADC fitting).
    if size(dwi_img, 4) ~= length(bvals)
        diag_bad_bval = diag_bad_bval + 1;
        fprintf('  💡 Pt %d (%s): DWI volume has %d frames but %d b-values — skipping\n', ...
            j, id_list{j}, size(dwi_img, 4), length(bvals));
        continue;
    end
    if ~any(bvals == 0)
        diag_bad_bval = diag_bad_bval + 1;
        fprintf('  💡 Pt %d (%s): no b=0 image found (b-values: %s) — skipping\n', ...
            j, id_list{j}, mat2str(bvals'));
        continue;
    end

    % Compute a voxel-wise ADC map using the same WLS monoexponential fit
    % as the main pipeline (fit_models.m) for consistency.
    %   Model: S(b) = S0 * exp(-b * ADC)
    %   Log-linearised: ln(S/S0) = -b * ADC
    %   WLS weights: W = S^2 (accounts for Rician noise heteroscedasticity)
    % Clamp non-physical values to the range [0, 3e-3] mm^2/s:
    %   0 mm^2/s      = no diffusion (solid tissue / bone)
    %   3e-3 mm^2/s   = free water at body temperature (upper physiological limit)
    %   Values outside this range indicate fitting artifacts (noise, motion).
    % Find the b=0 index (minimum b-value = no diffusion weighting = S0 reference)
    [~, b0_idx] = min(bvals);
    sz_dwi = size(dwi_img);
    % Flatten spatial dims to (voxels x b-values) for vectorized WLS fit
    dwi_flat = reshape(dwi_img, [prod(sz_dwi(1:3)), sz_dwi(4)]);
    % Only fit voxels with strictly positive signal at ALL b-values (log requires >0)
    all_pos = all(dwi_flat > 0, 2);
    adc_map = nan(sz_dwi(1:3));
    if any(all_pos)
        S_a = dwi_flat(all_pos, :);             % signal intensities at all b-values
        non_b0 = true(1, length(bvals)); non_b0(b0_idx) = false;  % exclude b=0 from regression
        A_b = -bvals(non_b0);                    % design matrix: negative b-values
        Y = log(S_a(:, non_b0) ./ S_a(:, b0_idx));  % log(S/S0) for linearized model
        W = S_a(:, non_b0).^2;                  % WLS weights correct log-domain heteroscedasticity
        numer = sum(W .* Y .* A_b', 2);         % weighted numerator per voxel
        denom = sum(W .* (A_b'.^2), 2);         % weighted denominator per voxel
        adc_vals = numer ./ denom;               % closed-form WLS estimate of ADC
        % Clamp to physiological range [0, 3e-3] mm^2/s
        adc_vals(adc_vals < 0 | ~isfinite(adc_vals)) = 0;
        adc_vals(adc_vals > 3e-3) = 3e-3;
        adc_map_flat = zeros(prod(sz_dwi(1:3)), 1);
        adc_map_flat(all_pos) = adc_vals;
        adc_map = reshape(adc_map_flat, sz_dwi(1:3));
    end

    % Pick the axial slice containing the largest cross-sectional GTV area.
    % This is the most representative slice for visual inspection because
    % it shows the widest extent of the tumour, maximising the number of
    % intra-tumoural voxels visible and reducing partial-volume effects
    % from tumour edges.
    gtv_areas = squeeze(sum(sum(gtv_img, 1), 2));
    [~, z_slice] = max(gtv_areas);

    % Extract 2-D slices for display
    b0_slice  = squeeze(dwi_img(:,:,z_slice,1));   % b=0 anatomical image
    adc_slice = squeeze(adc_map(:,:,z_slice));      % ADC map
    gtv_slice = squeeze(gtv_img(:,:,z_slice));      % GTV mask

    patients_plotted = patients_plotted + 1;
    row_in_fig = mod(patients_plotted - 1, pats_per_fig);  % 0-indexed row within current figure
    fig_num    = ceil(patients_plotted / pats_per_fig);     % 1-indexed figure number

    % --- Start a new figure every pats_per_fig patients ---
    if row_in_fig == 0
        % Finalise and save the previous figure (if any)
        if patients_plotted > 1
            sgtitle('Parameter Maps Overlaid on Anatomy (Fx1)', 'FontSize', 16, 'FontWeight', 'bold');
            cb = colorbar('Position', [0.93 0.11 0.015 0.8]);
            ylabel(cb, 'ADC (\times10^{-3} mm^2/s)', 'FontSize', 11);
            set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);
            print(gcf, fullfile(output_folder, sprintf('Parameter_Maps_%d.png', fig_num - 1)), '-dpng', '-r150');
            close(gcf);
        end
        % Determine how many rows this new figure needs
        n_rows_cur_fig = min(pats_per_fig, n_eligible - (fig_num - 1) * pats_per_fig);
        fig_height = max(400, n_rows_cur_fig * 200);
        figure('Name', sprintf('ADC Maps on Anatomy (%d)', fig_num), ...
               'Position', [50, 50, 1800, fig_height]);
        % Use turbo colormap (perceptually uniform, better than jet) in MATLAB;
        % fall back to jet in Octave where turbo is not available.
        if exist('OCTAVE_VERSION', 'builtin')
            colormap(jet);
        else
            colormap(turbo);
        end
    end

    % --- Column 1: b=0 anatomy with GTV contour (red) ---
    % The b=0 image provides the best soft-tissue contrast for anatomical
    % localisation (no diffusion weighting = T2-weighted contrast).
    % Red GTV contour enables verification of tumour delineation accuracy.
    subplot(n_rows_cur_fig, 3, row_in_fig*3 + 1);
    imagesc(b0_slice); axis image; axis off; colormap(gca, gray);
    hold on;
    contour(gtv_slice, [0.5 0.5], 'r', 'LineWidth', 1.5);
    hold off;
    title(sprintf('%s — b0 (GTV contour)', id_list{j}), 'Interpreter', 'none', 'FontSize', 11);

    % --- Column 2: ADC map with GTV contour (white) ---
    % Full-slice ADC map (fixed colour scale [0, 3] in ×10^-3 mm^2/s) shows
    % diffusion properties across all tissues.  White GTV contour used
    % because red would be invisible against warm-coloured high-ADC regions.
    % Scale to ×10^-3 so the colorbar label is self-contained (no separate exponent).
    subplot(n_rows_cur_fig, 3, row_in_fig*3 + 2);
    imagesc(adc_slice * 1e3, [0 3]); axis image; axis off;
    hold on;
    contour(gtv_slice, [0.5 0.5], 'w', 'LineWidth', 1.5);
    hold off;
    title('ADC map', 'FontSize', 11);

    % --- Column 3: ADC overlaid on anatomy (inside GTV only) ---
    % This fusion view is the most clinically useful: it shows the spatial
    % distribution of ADC within the tumour in anatomical context.  Regions
    % of low ADC (cool colours) within the GTV may represent dense,
    % treatment-resistant tumour tissue that could benefit from dose
    % escalation in an adaptive RT strategy.
    subplot(n_rows_cur_fig, 3, row_in_fig*3 + 3);
    % Normalise b=0 image to [0,1] and make it an RGB TrueColor image
    % so that the figure's colormap does not apply to it
    b0_range = max(b0_slice(:)) - min(b0_slice(:));
    if b0_range == 0, b0_range = 1; end
    b0_norm = (b0_slice - min(b0_slice(:))) ./ b0_range;
    b0_rgb = cat(3, b0_norm, b0_norm, b0_norm);
    imagesc(b0_rgb); axis image; axis off; hold on;
    % Mask ADC outside the GTV to NaN so only tumour voxels are coloured
    adc_overlay = adc_slice * 1e3;
    adc_overlay(gtv_slice < 0.5) = NaN;
    h_ov = imagesc(adc_overlay, [0 3]);
    % Use 60 % opacity for the colour overlay
    set(h_ov, 'AlphaData', ~isnan(adc_overlay) * 0.6);
    contour(gtv_slice, [0.5 0.5], 'w', 'LineWidth', 1.5);
    hold off;
    title('ADC on Anatomy', 'FontSize', 11);
end

% Save and close the final parameter-maps figure
if patients_plotted > 0
    sgtitle('Parameter Maps Overlaid on Anatomy (Fx1)', 'FontSize', 16, 'FontWeight', 'bold');
    cb = colorbar('Position', [0.93 0.11 0.015 0.8]);
    ylabel(cb, 'ADC (\times10^{-3} mm^2/s)', 'FontSize', 11);
    set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);
    print(gcf, fullfile(output_folder, sprintf('Parameter_Maps_%d.png', fig_num)), '-dpng', '-r150');
    close(gcf);
    fprintf('  Plotted %d patients.\n', patients_plotted);
else
    fprintf('  No patients with complete NIfTI data found on disk. Skipping.\n');
end

% Diagnostic summary of skip reasons
if diag_out_of_range + diag_empty_adc + diag_missing_file + diag_bad_bval > 0
    fprintf('  --- Parameter Maps skip summary ---\n');
    fprintf('  nPat=%d, data_vectors_gtvp rows=%d\n', nPat, size(data_vectors_gtvp, 1));
    if diag_out_of_range > 0, fprintf('    Index out of range: %d\n', diag_out_of_range); end
    if diag_empty_adc    > 0, fprintf('    Empty adc_vector:   %d\n', diag_empty_adc); end
    if diag_missing_file > 0, fprintf('    Missing NIfTI file: %d\n', diag_missing_file); end
    if diag_bad_bval     > 0, fprintf('    b-value mismatch:   %d\n', diag_bad_bval); end
    fprintf('    Plotted:            %d\n', patients_plotted);
end

end
