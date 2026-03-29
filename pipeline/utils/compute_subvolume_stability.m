function stability = compute_subvolume_stability(data_vectors_gtvp, config_struct, id_list, gtv_locations)
%COMPUTE_SUBVOLUME_STABILITY Dice between each fraction's core mask and Fx1 baseline.
%
%   For each active core method x DWI pipeline x patient, computes
%   extract_tumor_core at every timepoint and compares to the Fx1 mask via Dice.
%
%   Inputs:
%       data_vectors_gtvp  - struct array [nPatients x nTp x nRepeats]
%       config_struct      - pipeline configuration struct
%       id_list            - cell array of patient identifiers
%       gtv_locations      - cell array of GTV mask file paths
%
%   Outputs:
%       stability.dice_vs_baseline - [nMethods x nTp x nPatients] Dice values
%                                    (Fx1 column is always 1.0 by definition)
%       stability.method_names     - cell array of method names
%       stability.n_patients       - scalar
%       stability.n_timepoints     - scalar

    % --- Diary ---
    output_folder = config_struct.output_folder;
    diary_file = fullfile(output_folder, ...
        sprintf('subvolume_stability_output_%s.txt', config_struct.dwi_type_name));
    diary(diary_file);

    fprintf('\n🔬 Sub-Volume Stability Over Fractions (%s)\n', config_struct.dwi_type_name);

    % --- Constants ---
    ALL_METHODS = {'adc_threshold', 'd_threshold', 'df_intersection', ...
        'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
        'percentile', 'spectral', 'fdm'};

    if isfield(config_struct, 'active_core_methods') && ~isempty(config_struct.active_core_methods)
        active_methods = config_struct.active_core_methods;
    else
        active_methods = ALL_METHODS;
    end

    n_methods = numel(active_methods);
    n_patients = numel(id_list);
    nTp = size(data_vectors_gtvp, 2);
    dwi_type = config_struct.dwi_types_to_run(1);

    % --- Pre-allocate ---
    dice_vs_baseline = nan(n_methods, nTp, n_patients);

    % Suppress warnings during batch extraction
    prev_warn_state = suppress_core_warnings();
    cleanup_obj = onCleanup(@() warning(prev_warn_state));

    rng(42);

    % --- Main loop ---
    for j = 1:n_patients
        text_progress_bar(j, n_patients, 'Sub-volume stability');

        % Load 3D GTV mask for Fx1
        has_3d_fx1 = false;
        gtv_mask_3d_fx1 = [];
        if ~isempty(gtv_locations) && ...
                size(gtv_locations, 1) >= j && size(gtv_locations, 2) >= 1
            gtv_mat = gtv_locations{j, 1, 1};
            if ~isempty(gtv_mat) && exist(gtv_mat, 'file')
                gtv_mask_3d_fx1 = safe_load_mask(gtv_mat, 'Stvol3d');
                [adc_ref, ~, ~, ~] = select_dwi_vectors(data_vectors_gtvp, j, 1, 1, dwi_type);
                if ~isempty(gtv_mask_3d_fx1) && ~isempty(adc_ref) && ...
                        sum(gtv_mask_3d_fx1(:) == 1) == numel(adc_ref)
                    has_3d_fx1 = true;
                end
            end
        end

        for m = 1:n_methods
            temp_config = config_struct;
            temp_config.core_method = active_methods{m};

            % --- Compute baseline (Fx1) mask ---
            [adc_fx1, d_fx1, f_fx1, dstar_fx1] = select_dwi_vectors( ...
                data_vectors_gtvp, j, 1, 1, dwi_type);
            if isempty(adc_fx1)
                continue;
            end

            core_opts_fx1 = struct('timepoint_index', 1);
            rng(42);
            [baseline_mask, ~] = extract_tumor_core(temp_config, ...
                adc_fx1, d_fx1, f_fx1, dstar_fx1, has_3d_fx1, gtv_mask_3d_fx1, core_opts_fx1);

            % Fx1 Dice vs itself = 1.0 (by definition)
            if any(baseline_mask)
                dice_vs_baseline(m, 1, j) = 1.0;
            else
                dice_vs_baseline(m, 1, j) = NaN;
            end

            % --- Compute mask at each subsequent timepoint ---
            for k = 2:nTp
                [adc_k, d_k, f_k, dstar_k] = select_dwi_vectors( ...
                    data_vectors_gtvp, j, k, 1, dwi_type);
                if isempty(adc_k)
                    continue;
                end

                % Ensure same number of voxels for Dice comparison
                if numel(adc_k) ~= numel(adc_fx1)
                    continue;
                end

                core_opts_k = struct('timepoint_index', k);
                if k > 1
                    core_opts_k.baseline_adc_vec = adc_fx1;
                    core_opts_k.baseline_d_vec = d_fx1;
                end

                % Load timepoint-specific 3D mask if available
                has_3d_k = false;
                gtv_mask_3d_k = [];
                if ~isempty(gtv_locations) && ...
                        size(gtv_locations, 1) >= j && size(gtv_locations, 2) >= k
                    gtv_mat_k = gtv_locations{j, k, 1};
                    if ~isempty(gtv_mat_k) && exist(gtv_mat_k, 'file')
                        gtv_mask_3d_k = safe_load_mask(gtv_mat_k, 'Stvol3d');
                        if ~isempty(gtv_mask_3d_k) && ...
                                sum(gtv_mask_3d_k(:) == 1) == numel(adc_k)
                            has_3d_k = true;
                        end
                    end
                end
                if ~has_3d_k
                    has_3d_k = has_3d_fx1;
                    gtv_mask_3d_k = gtv_mask_3d_fx1;
                end

                rng(42);
                [mask_k, ~] = extract_tumor_core(temp_config, ...
                    adc_k, d_k, f_k, dstar_k, has_3d_k, gtv_mask_3d_k, core_opts_k);

                % Compute Dice
                if ~any(baseline_mask) || ~any(mask_k)
                    dice_vs_baseline(m, k, j) = NaN;
                else
                    intersection = sum(baseline_mask & mask_k);
                    dice_vs_baseline(m, k, j) = 2 * intersection / ...
                        (sum(baseline_mask) + sum(mask_k));
                end
            end
        end
    end

    % --- Package results ---
    stability.dice_vs_baseline = dice_vs_baseline;
    stability.method_names = active_methods;
    stability.n_patients = n_patients;
    stability.n_timepoints = nTp;

    % --- Print summary ---
    fprintf('\n\nMean Dice vs Fx1 Baseline (%s):\n', config_struct.dwi_type_name);
    tp_labels = cell(1, nTp);
    for k = 1:nTp
        if k <= 5
            tp_labels{k} = sprintf('Fx%d', k);
        else
            tp_labels{k} = 'Post';
        end
    end
    fprintf('%-25s', 'Method');
    for k = 1:nTp
        fprintf(' %8s', tp_labels{k});
    end
    fprintf('\n%s\n', repmat('-', 1, 25 + 9*nTp));
    for m = 1:n_methods
        fprintf('%-25s', active_methods{m});
        for k = 1:nTp
            mean_dice = nanmean_safe(squeeze(dice_vs_baseline(m, k, :)));
            fprintf(' %8.3f', mean_dice);
        end
        fprintf('\n');
    end

    % --- Generate figure: line plot of mean Dice across fractions ---
    try
        fig = figure('Visible', 'off', 'Position', [100 100 800 500]);
        colors = lines(n_methods);
        hold on;
        for m = 1:n_methods
            mean_dice_per_tp = nan(1, nTp);
            for k = 1:nTp
                mean_dice_per_tp(k) = nanmean_safe(squeeze(dice_vs_baseline(m, k, :)));
            end
            plot(1:nTp, mean_dice_per_tp, '-o', 'Color', colors(m,:), ...
                'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', active_methods{m});
        end
        hold off;
        set(gca, 'XTick', 1:nTp, 'XTickLabel', tp_labels, 'FontSize', 9);
        xlabel('Fraction');
        ylabel('Mean Dice vs Fx1 Baseline');
        ylim([0 1.05]);
        title(sprintf('Sub-Volume Stability Over Fractions (%s)', config_struct.dwi_type_name));
        legend('Location', 'southwestoutside', 'FontSize', 7);
        saveas(fig, fullfile(output_folder, ...
            sprintf('subvolume_stability_%s.png', config_struct.dwi_type_name)));
        close(fig);
    catch ME_fig
        fprintf('⚠️ Could not generate stability figure: %s\n', ME_fig.message);
    end

    fprintf('\n✅ Sub-volume stability analysis complete: %d patients, %d timepoints.\n', ...
        n_patients, nTp);
    diary off;
end
