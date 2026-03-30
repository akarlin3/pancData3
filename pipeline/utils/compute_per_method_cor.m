function cor_results = compute_per_method_cor(data_vectors_gtvp, config_struct, id_list, gtv_locations)
%COMPUTE_PER_METHOD_COR Coefficient of Reproducibility for each core method's sub-volume.
%
%   For patients with >=2 Fx1 repeat scans, computes the sub-volume fraction
%   (% of GTV) for each repeat using each core method, then derives wCV and CoR.
%
%   Inputs:
%       data_vectors_gtvp  - struct array [nPatients x nTp x nRepeats]
%       config_struct      - pipeline configuration struct
%       id_list            - cell array of patient identifiers
%       gtv_locations      - cell array of GTV mask file paths
%
%   Outputs:
%       cor_results.method_names            - {nMethods x 1}
%       cor_results.median_wcv              - [nMethods x 1] median within-subject CV
%       cor_results.cor                     - [nMethods x 1] Coefficient of Reproducibility (%)
%       cor_results.n_patients_with_repeats - scalar
%       cor_results.per_patient_wcv         - [nMethods x nPatients] wCV per patient (NaN if <2 repeats)
%       cor_results.subvol_per_repeat       - [nMethods x nPatients x nRepeats] sub-volume fraction

    % --- Diary ---
    output_folder = config_struct.output_folder;
    diary_file = fullfile(output_folder, ...
        sprintf('per_method_cor_output_%s.txt', config_struct.dwi_type_name));
    diary(diary_file);

    fprintf('\n🔬 Per-Method Coefficient of Reproducibility (%s)\n', config_struct.dwi_type_name);

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
    n_repeats = size(data_vectors_gtvp, 3);
    dwi_type = config_struct.dwi_types_to_run(1);

    % --- Pre-allocate ---
    subvol_per_repeat = nan(n_methods, n_patients, max(n_repeats, 1));
    per_patient_wcv = nan(n_methods, n_patients);

    % Suppress warnings during batch extraction
    prev_warn_state = suppress_core_warnings();
    cleanup_obj = onCleanup(@() warning(prev_warn_state));

    rng(42);

    % --- Identify patients with repeats ---
    has_repeat = false(n_patients, 1);
    for j = 1:n_patients
        if n_repeats >= 2
            [adc_r2, ~, ~, ~] = select_dwi_vectors(data_vectors_gtvp, j, 1, 2, dwi_type);
            if ~isempty(adc_r2)
                has_repeat(j) = true;
            end
        end
    end
    n_with_repeats = sum(has_repeat);
    fprintf('   Found %d/%d patients with Fx1 repeat scans.\n', n_with_repeats, n_patients);

    if n_with_repeats == 0
        fprintf('   ⚠️ No repeat scans found. CoR cannot be computed.\n');
        cor_results.method_names = active_methods;
        cor_results.median_wcv = nan(n_methods, 1);
        cor_results.cor = nan(n_methods, 1);
        cor_results.n_patients_with_repeats = 0;
        cor_results.per_patient_wcv = per_patient_wcv;
        cor_results.subvol_per_repeat = subvol_per_repeat;
        diary off;
        return;
    end

    % --- Main loop ---
    for j = 1:n_patients
        text_progress_bar(j, n_patients, 'Per-method CoR');

        if ~has_repeat(j)
            continue;
        end

        % Build core_opts (Fx1, timepoint_index=1)
        core_opts = struct('timepoint_index', 1);

        % Load 3D GTV mask
        has_3d = false;
        gtv_mask_3d = [];
        if ~isempty(gtv_locations) && ...
                size(gtv_locations, 1) >= j && size(gtv_locations, 2) >= 1
            gtv_mat = gtv_locations{j, 1, 1};
            if ~isempty(gtv_mat) && exist(gtv_mat, 'file')
                gtv_mask_3d = safe_load_mask(gtv_mat, 'Stvol3d');
                [adc_ref, ~, ~, ~] = select_dwi_vectors(data_vectors_gtvp, j, 1, 1, dwi_type);
                if ~isempty(gtv_mask_3d) && ~isempty(adc_ref) && ...
                        sum(gtv_mask_3d(:) == 1) == numel(adc_ref)
                    has_3d = true;
                end
            end
        end

        for m = 1:n_methods
            fracs = [];
            temp_config = config_struct;
            temp_config.core_method = active_methods{m};

            for rpi = 1:n_repeats
                [adc_vec, d_vec, f_vec, dstar_vec] = select_dwi_vectors( ...
                    data_vectors_gtvp, j, 1, rpi, dwi_type);

                if isempty(adc_vec)
                    continue;
                end

                rng(42);
                [core_mask, ~] = extract_tumor_core(temp_config, ...
                    adc_vec, d_vec, f_vec, dstar_vec, has_3d, gtv_mask_3d, core_opts);

                n_valid = sum(~isnan(adc_vec));
                if n_valid > 0
                    frac = sum(core_mask) / n_valid;
                else
                    frac = NaN;
                end

                subvol_per_repeat(m, j, rpi) = frac;
                fracs = [fracs; frac]; %#ok<AGROW>
            end

            % Compute within-subject CV
            if numel(fracs) >= 2 && mean(fracs) > 0
                per_patient_wcv(m, j) = std(fracs) / mean(fracs);
            end
        end
    end

    % --- Aggregate: median wCV and CoR ---
    median_wcv = nan(n_methods, 1);
    cor_val = nan(n_methods, 1);
    for m = 1:n_methods
        valid_wcv = per_patient_wcv(m, :);
        valid_wcv = valid_wcv(~isnan(valid_wcv));
        if ~isempty(valid_wcv)
            median_wcv(m) = median(valid_wcv);
            cor_val(m) = 1.96 * sqrt(2) * median_wcv(m) * 100;
        end
    end

    % --- Package results ---
    cor_results.method_names = active_methods;
    cor_results.median_wcv = median_wcv;
    cor_results.cor = cor_val;
    cor_results.n_patients_with_repeats = n_with_repeats;
    cor_results.per_patient_wcv = per_patient_wcv;
    cor_results.subvol_per_repeat = subvol_per_repeat;

    % --- Print summary ---
    fprintf('\n\nPer-Method CoR Summary (%s pipeline):\n', config_struct.dwi_type_name);
    fprintf('%-25s %12s %10s\n', 'Method', 'Median wCV', 'CoR (%)');
    fprintf('%s\n', repmat('-', 1, 50));
    for m = 1:n_methods
        fprintf('%-25s %11.3f %9.1f%%\n', active_methods{m}, median_wcv(m), cor_val(m));
    end

    % --- Generate bar chart ---
    try
        fig = figure('Visible', 'off', 'Position', [100 100 800 450]);
        bar(cor_val);
        set(gca, 'XTick', 1:n_methods, 'XTickLabel', active_methods, ...
            'FontSize', 7, 'XTickLabelRotation', 45);
        ylabel('CoR (%)');
        title(sprintf('Per-Method Coefficient of Reproducibility (%s)', config_struct.dwi_type_name));
        hold on;
        yline(15, '--g', 'Reproducible (<15%)');
        yline(30, '--r', 'Caution (>30%)');
        hold off;
        saveas(fig, fullfile(output_folder, ...
            sprintf('per_method_cor_%s.png', config_struct.dwi_type_name)));
        close(fig);
    catch ME_fig
        fprintf('⚠️ Could not generate CoR figure: %s\n', ME_fig.message);
    end

    fprintf('\n✅ Per-method CoR analysis complete: %d patients with repeats.\n', n_with_repeats);
    diary off;
end
