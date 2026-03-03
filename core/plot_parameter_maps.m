function plot_parameter_maps(data_vectors_gtvp, nPat, id_list, dataloc, output_folder)
% PLOT_PARAMETER_MAPS — "Parameter Maps overlaid on Anatomy"
%
%  For each patient with Fx1 data available, load the DWI volume, compute
%  the ADC map via monoexponential fit, and display three panels per
%  patient:
%    (a) b=0 anatomical image with GTV contour overlay
%    (b) Full-slice ADC map with GTV contour
%    (c) ADC overlaid on anatomy (semi-transparent inside GTV only)
%  Patients are batched into multi-row figures (pats_per_fig rows each).

fprintf('\n--- 1. Parameter Maps overlaid on Anatomy ---\n');

% Track how many patients have been plotted so far
patients_plotted = 0;
% Maximum number of patient rows per figure panel
pats_per_fig     = 5;

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
        diag_out_of_range = diag_out_of_range + 1;
        fprintf('  💡 Pt %d: index exceeds data_vectors_gtvp rows (%d)\n', j, size(data_vectors_gtvp, 1));
        continue;
    end
    % Skip patients without extracted ADC data at Fx1
    s = data_vectors_gtvp(j, 1, 1);
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
        diag_bad_bval = diag_bad_bval + 1;
        fprintf('  💡 Pt %d (%s): protocol deviation — b-values %s (expected %s)\n', ...
            j, id_list{j}, mat2str(bvals'), mat2str(expected_bvals'));
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
        if exist('OCTAVE_VERSION', 'builtin')
            colormap(jet);
        else
            colormap(turbo);
        end
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
    % Normalise b=0 image to [0,1] and make it an RGB TrueColor image
    % so that the figure's colormap does not apply to it
    b0_norm = (b0_slice - min(b0_slice(:))) ./ (max(b0_slice(:)) - min(b0_slice(:)));
    b0_rgb = cat(3, b0_norm, b0_norm, b0_norm);
    imagesc(b0_rgb); axis image; axis off; hold on;
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
