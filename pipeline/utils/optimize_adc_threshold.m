function opt_results = optimize_adc_threshold(data_vectors_gtvp, config_struct, id_list, gtv_locations)
%OPTIMIZE_ADC_THRESHOLD Sweep ADC threshold and find value maximizing Fx1 repeat Dice.
%
%   opt_results = optimize_adc_threshold(data_vectors_gtvp, config_struct, ...
%                                        id_list, gtv_locations)
%
%   Sweeps ADC threshold from 0.8e-3 to 2.0e-3 in steps of 0.1e-3 (13 values).
%   For each threshold, for each patient with >=2 Fx1 repeat scans, computes
%   pairwise Dice coefficients between threshold-defined binary sub-volumes
%   (after morphological cleanup) and averages across pairs per patient.
%   The optimal threshold is the one that maximizes the median Dice across
%   patients.
%
%   Outputs are saved to opt_results struct and a figure is generated with
%   dual y-axis (median Dice left, median volume fraction right). Vertical
%   dashed lines mark the current default (0.001), proposed (0.0016), and
%   optimal thresholds.
%
%   Inputs:
%       data_vectors_gtvp - struct array [nPatients x nTp x nRepeats]
%       config_struct     - pipeline configuration (must contain output_folder,
%                           dwi_types_to_run, dwi_type_name). morph_se and
%                           morph_min_cc are honoured if present; otherwise
%                           sensible defaults are used.
%       id_list           - cell array of patient identifiers
%       gtv_locations     - cell array {nPatients x nTp x nRepeats} of
%                           GTV mask file paths
%
%   Output struct fields:
%       thresholds        - [1x13] threshold values in mm^2/s
%       median_dice       - [1x13] median Dice across patients per threshold
%       mean_dice         - [1x13] mean Dice across patients per threshold
%       std_dice          - [1x13] std Dice across patients per threshold
%       median_vol_frac   - [1x13] median sub-volume fraction per threshold
%       mean_vol_frac     - [1x13] mean sub-volume fraction per threshold
%       n_patients        - scalar: patient count with >=2 repeats
%       optimal_thresh    - scalar: threshold maximizing median Dice
%       optimal_dice      - scalar
%       optimal_vol_frac  - scalar
%       per_patient_dice  - [nPatients x 13] per-patient mean Dice per threshold

    % --- Threshold sweep ---
    thresholds = 0.8e-3 : 0.1e-3 : 2.0e-3;  % 13 values
    n_thresh = numel(thresholds);
    n_patients = numel(id_list);

    if isfield(config_struct, 'dwi_types_to_run') && ~isempty(config_struct.dwi_types_to_run)
        dwi_type = config_struct.dwi_types_to_run(1);
    else
        dwi_type = 1;
    end

    if isfield(config_struct, 'dwi_type_name') && ~isempty(config_struct.dwi_type_name)
        dwi_label = config_struct.dwi_type_name;
    else
        dwi_label = 'Standard';
    end

    % Morphological cleanup parameters
    if isfield(config_struct, 'morph_se') && ~isempty(config_struct.morph_se)
        morph_se = config_struct.morph_se;
    else
        try
            morph_se = strel('sphere', 1);
        catch
            morph_se = strel('square', 3);
        end
    end
    if isfield(config_struct, 'morph_min_cc') && ~isempty(config_struct.morph_min_cc)
        morph_min_cc = config_struct.morph_min_cc;
    else
        morph_min_cc = 10;
    end

    % --- Pre-allocate ---
    per_patient_dice = nan(n_patients, n_thresh);
    per_patient_vol_frac = nan(n_patients, n_thresh);
    n_valid_patients = 0;

    fprintf('\n\xf0\x9f\x94\xac Optimizing ADC threshold over %d values for %s...\n', ...
        n_thresh, dwi_label);

    % --- Loop over patients first, then thresholds (so we load GTV masks once) ---
    for j = 1:n_patients
        if exist('text_progress_bar', 'file')
            text_progress_bar(j, n_patients, 'ADC threshold sweep');
        end

        % Collect valid Fx1 repeats and their ADC vectors.
        valid_rpis = [];
        rpt_adc = {};
        for rpi = 1:size(data_vectors_gtvp, 3)
            [adc_vec, ~, ~, ~] = select_dwi_vectors(data_vectors_gtvp, j, 1, rpi, dwi_type);
            if ~isempty(adc_vec)
                valid_rpis(end+1) = rpi; %#ok<AGROW>
                rpt_adc{end+1} = adc_vec; %#ok<AGROW>
            end
        end

        if numel(valid_rpis) < 2
            continue;
        end

        % Load 3D GTV masks for each valid repeat.  When only one shared
        % Fx1 contour exists (gtv_locations{j,1,1} populated but {j,1,2..N}
        % empty), fall back to that shared path so the sweep is not skipped
        % for every patient in the cohort.
        shared_fx1_path = '';
        if ~isempty(gtv_locations) && size(gtv_locations, 1) >= j
            for ri_shared = 1:size(gtv_locations, 3)
                candidate = gtv_locations{j, 1, ri_shared};
                if ~isempty(candidate); shared_fx1_path = candidate; break; end
            end
        end
        rpt_masks_3d = cell(numel(valid_rpis), 1);
        has_all_3d = true;
        for ri = 1:numel(valid_rpis)
            rpi_idx = valid_rpis(ri);
            if isempty(gtv_locations) || size(gtv_locations, 1) < j || ...
                    size(gtv_locations, 3) < rpi_idx
                has_all_3d = false; break;
            end
            gtv_path = gtv_locations{j, 1, rpi_idx};
            if isempty(gtv_path); gtv_path = shared_fx1_path; end
            if isempty(gtv_path) || ~exist(gtv_path, 'file')
                has_all_3d = false; break;
            end
            try
                rpt_masks_3d{ri} = safe_load_mask(gtv_path, 'Stvol3d');
            catch
                rpt_masks_3d{ri} = [];
            end
            if isempty(rpt_masks_3d{ri})
                has_all_3d = false; break;
            end
        end

        if ~has_all_3d
            continue;
        end

        % Determine voxel dimensions for Dice (size only matters for HD, not Dice).
        vox_dims = data_vectors_gtvp(j, 1, 1).vox_dims;
        if isempty(vox_dims) || ~isnumeric(vox_dims) || numel(vox_dims) ~= 3
            vox_dims = [1, 1, 1];
        end

        n_valid_patients = n_valid_patients + 1;

        % --- Loop over thresholds ---
        for ti = 1:n_thresh
            thresh = thresholds(ti);

            % Pre-compute threshold masks per repeat.
            subvols_3d = cell(numel(valid_rpis), 1);
            subvol_fracs = nan(numel(valid_rpis), 1);
            for ri = 1:numel(valid_rpis)
                mask_3d = rpt_masks_3d{ri};
                adc_vec = rpt_adc{ri};
                n_gtv_vox = sum(mask_3d(:) == 1);
                if numel(adc_vec) ~= n_gtv_vox
                    subvols_3d{ri} = [];
                    continue;
                end

                subvol = false(size(mask_3d));
                subvol(mask_3d == 1) = adc_vec < thresh;
                % Morphological cleanup (same as compute_spatial_repeatability)
                try
                    subvol = imclose(imopen(subvol, morph_se), morph_se);
                    subvol = bwareaopen(subvol, morph_min_cc);
                catch
                    % Skip cleanup if image processing functions fail
                end
                subvols_3d{ri} = subvol;
                if n_gtv_vox > 0
                    subvol_fracs(ri) = sum(subvol(:)) / n_gtv_vox;
                end
            end

            % Pairwise Dice
            pair_dice = [];
            for ri1 = 1:numel(valid_rpis)-1
                for ri2 = ri1+1:numel(valid_rpis)
                    a = subvols_3d{ri1};
                    b = subvols_3d{ri2};
                    if isempty(a) || isempty(b) || ~isequal(size(a), size(b))
                        continue;
                    end
                    [d_val, ~, ~] = compute_dice_hausdorff(a, b, vox_dims);
                    if ~isnan(d_val)
                        pair_dice(end+1) = d_val; %#ok<AGROW>
                    end
                end
            end

            if ~isempty(pair_dice)
                per_patient_dice(j, ti) = mean(pair_dice);
            end
            if any(~isnan(subvol_fracs))
                per_patient_vol_frac(j, ti) = mean(subvol_fracs, 'omitnan');
            end
        end
    end

    % --- Aggregate statistics across patients ---
    median_dice = nan(1, n_thresh);
    mean_dice = nan(1, n_thresh);
    std_dice = nan(1, n_thresh);
    median_vol_frac = nan(1, n_thresh);
    mean_vol_frac = nan(1, n_thresh);
    for ti = 1:n_thresh
        col_d = per_patient_dice(:, ti);
        col_d = col_d(~isnan(col_d));
        if ~isempty(col_d)
            median_dice(ti) = median(col_d);
            mean_dice(ti) = mean(col_d);
            std_dice(ti) = std(col_d);
        end
        col_v = per_patient_vol_frac(:, ti);
        col_v = col_v(~isnan(col_v));
        if ~isempty(col_v)
            median_vol_frac(ti) = median(col_v);
            mean_vol_frac(ti) = mean(col_v);
        end
    end

    % --- Optimal threshold ---
    if all(isnan(median_dice))
        optimal_idx = 1;
    else
        [~, optimal_idx] = max(median_dice);
    end
    optimal_thresh = thresholds(optimal_idx);
    optimal_dice = median_dice(optimal_idx);
    optimal_vol_frac = median_vol_frac(optimal_idx);

    % --- Console output ---
    fprintf('\n  Optimal ADC threshold: %.4f mm\xc2\xb2/s (Dice=%.2f, vol_frac=%.1f%%)\n', ...
        optimal_thresh, optimal_dice, 100 * optimal_vol_frac);

    [~, curr_idx] = min(abs(thresholds - 0.001));
    [~, prop_idx] = min(abs(thresholds - 0.0016));
    fprintf('  Current default: 0.0010 (Dice=%.2f, vol_frac=%.1f%%)\n', ...
        safe_val(median_dice, curr_idx), 100 * safe_val(median_vol_frac, curr_idx));
    fprintf('  Proposed 0.0016: (Dice=%.2f, vol_frac=%.1f%%)\n', ...
        safe_val(median_dice, prop_idx), 100 * safe_val(median_vol_frac, prop_idx));

    % --- Figure ---
    try
        fig = figure('Visible', 'off', 'Position', [100, 100, 900, 500]);

        yyaxis left;
        err_mask = ~isnan(median_dice);
        errorbar(thresholds(err_mask), median_dice(err_mask), std_dice(err_mask), ...
            '-o', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, ...
            'MarkerFaceColor', [0 0.4470 0.7410]);
        ylabel('Median Dice (Fx1 repeat)');
        ylim([0, 1]);

        yyaxis right;
        vf_mask = ~isnan(median_vol_frac);
        plot(thresholds(vf_mask), 100 * median_vol_frac(vf_mask), '-s', ...
            'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, ...
            'MarkerFaceColor', [0.8500 0.3250 0.0980]);
        ylabel('Median Sub-Volume Fraction (%)');

        xlabel('ADC Threshold (mm^2/s)');
        title(sprintf('ADC Threshold Optimization (%s, N=%d patients)', ...
            dwi_label, n_valid_patients));
        grid on;

        % Reference lines
        hold on;
        yl = ylim;
        yyaxis left;
        yl_left = ylim;
        plot([0.001, 0.001], yl_left, 'k--', 'LineWidth', 1.2);
        text(0.001, yl_left(2)*0.95, ' Current', 'Color', 'k', 'FontSize', 9);

        plot([0.0016, 0.0016], yl_left, '--', 'Color', [0 0.6 0], 'LineWidth', 1.2);
        text(0.0016, yl_left(2)*0.88, ' Proposed', 'Color', [0 0.6 0], 'FontSize', 9);

        plot([optimal_thresh, optimal_thresh], yl_left, 'm--', 'LineWidth', 1.5);
        text(optimal_thresh, yl_left(2)*0.80, sprintf(' Optimal (%.4f)', optimal_thresh), ...
            'Color', 'm', 'FontSize', 9);
        hold off;

        if isfield(config_struct, 'output_folder') && ~isempty(config_struct.output_folder) && ...
                exist(config_struct.output_folder, 'dir')
            png_path = fullfile(config_struct.output_folder, ...
                sprintf('adc_threshold_optimization_%s.png', dwi_label));
            print(fig, png_path, '-dpng', '-r300');
        end
        close(fig);
    catch ME_fig
        fprintf('  \xe2\x9a\xa0\xef\xb8\x8f Figure generation failed: %s\n', ME_fig.message);
    end

    % --- Pack output struct ---
    opt_results = struct();
    opt_results.thresholds = thresholds;
    opt_results.median_dice = median_dice;
    opt_results.mean_dice = mean_dice;
    opt_results.std_dice = std_dice;
    opt_results.median_vol_frac = median_vol_frac;
    opt_results.mean_vol_frac = mean_vol_frac;
    opt_results.n_patients = n_valid_patients;
    opt_results.optimal_thresh = optimal_thresh;
    opt_results.optimal_dice = optimal_dice;
    opt_results.optimal_vol_frac = optimal_vol_frac;
    opt_results.per_patient_dice = per_patient_dice;
end


function v = safe_val(vec, idx)
    if idx < 1 || idx > numel(vec) || isnan(vec(idx))
        v = NaN;
    else
        v = vec(idx);
    end
end
