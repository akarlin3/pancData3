function failure_table = compute_core_failure_rates(data_vectors_gtvp, config_struct, id_list, gtv_locations)
%COMPUTE_CORE_FAILURE_RATES Aggregate failure rates for all 11 core methods across the cohort.
%
%   Runs extract_tumor_core (with fit_info output) for all 11 methods x 3 DWI
%   pipelines x all patients x all timepoints. Aggregates four failure modes:
%     1. Fallback triggered (method fell back to ADC threshold)
%     2. Empty mask (core mask is all-false)
%     3. Insufficient voxels (core has < min_vox_hist voxels)
%     4. All-NaN input (parameter vectors were all-NaN, fit never ran)
%
%   This is a NON-FATAL pipeline step: errors are caught by the orchestrator.
%
% Inputs:
%   data_vectors_gtvp  - struct array [nPatients x nTp x 1]
%   config_struct      - pipeline configuration struct
%   id_list            - cell array of patient identifiers
%   gtv_locations      - cell array of GTV mask file paths
%
% Outputs:
%   failure_table      - struct with per-method, per-pipeline failure rates

    % --- Diary ---
    output_folder = config_struct.output_folder;
    diary_file = fullfile(output_folder, ...
        sprintf('core_failure_rates_output_%s.txt', config_struct.dwi_type_name));
    diary(diary_file);

    fprintf('\n🔬 Core Method Failure Rates (%s)\n', config_struct.dwi_type_name);
    fprintf('   Aggregating failure modes for all 11 methods x 3 pipelines.\n\n');

    % --- Constants ---
    ALL_METHODS = {'adc_threshold', 'd_threshold', 'df_intersection', ...
        'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
        'percentile', 'spectral', 'fdm'};
    PIPELINE_NAMES = {'Standard', 'DnCNN', 'IVIMNet'};
    n_methods = 11;
    n_pipelines = 3;
    n_patients = numel(id_list);
    nTp = size(data_vectors_gtvp, 2);

    % --- Pre-allocate ---
    was_fallback = false(n_methods, n_pipelines, n_patients, nTp);
    was_empty = false(n_methods, n_pipelines, n_patients, nTp);
    was_insufficient = false(n_methods, n_pipelines, n_patients, nTp);
    was_all_nan = false(n_methods, n_pipelines, n_patients, nTp);
    core_voxel_counts = nan(n_methods, n_pipelines, n_patients, nTp);
    pipe_skipped = false(n_pipelines, n_patients, nTp);

    % --- Main loop ---
    for j = 1:n_patients
        text_progress_bar(j, n_patients, 'Core failure rates');
        for k = 1:nTp
            entry = data_vectors_gtvp(j, k, 1);

            % Extract parameter vectors for each pipeline
            pipe_vecs = cell(n_pipelines, 4);

            % Standard
            pipe_vecs(1,:) = {entry.adc_vector, entry.d_vector, entry.f_vector, entry.dstar_vector};

            % DnCNN
            dncnn_adc = []; dncnn_d = []; dncnn_f = []; dncnn_dstar = [];
            if isfield(entry, 'adc_vector_dncnn'), dncnn_adc = entry.adc_vector_dncnn; end
            if isfield(entry, 'd_vector_dncnn'), dncnn_d = entry.d_vector_dncnn; end
            if isfield(entry, 'f_vector_dncnn'), dncnn_f = entry.f_vector_dncnn; end
            if isfield(entry, 'dstar_vector_dncnn'), dncnn_dstar = entry.dstar_vector_dncnn; end
            pipe_vecs(2,:) = {dncnn_adc, dncnn_d, dncnn_f, dncnn_dstar};

            % IVIMNet (ADC shared with Standard)
            ivimnet_d = []; ivimnet_f = []; ivimnet_dstar = [];
            if isfield(entry, 'd_vector_ivimnet'), ivimnet_d = entry.d_vector_ivimnet; end
            if isfield(entry, 'f_vector_ivimnet'), ivimnet_f = entry.f_vector_ivimnet; end
            if isfield(entry, 'dstar_vector_ivimnet'), ivimnet_dstar = entry.dstar_vector_ivimnet; end
            pipe_vecs(3,:) = {entry.adc_vector, ivimnet_d, ivimnet_f, ivimnet_dstar};

            % Check which pipelines have data
            for p = 1:n_pipelines
                if isempty(pipe_vecs{p,1}) || all(isnan(pipe_vecs{p,1}))
                    pipe_skipped(p, j, k) = true;
                end
            end

            % Load 3D GTV mask
            has_3d = false;
            gtv_mask_3d = [];
            if ~isempty(gtv_locations) && ...
                    size(gtv_locations, 1) >= j && size(gtv_locations, 2) >= k
                gtv_mat = gtv_locations{j, k, 1};
                if ~isempty(gtv_mat) && exist(gtv_mat, 'file')
                    gtv_mask_3d = safe_load_mask(gtv_mat, 'Stvol3d');
                    adc_ref = pipe_vecs{1,1};
                    if ~isempty(gtv_mask_3d) && ~isempty(adc_ref) && ...
                            sum(gtv_mask_3d(:) == 1) == numel(adc_ref)
                        has_3d = true;
                    end
                end
                % Fx1 mask fallback for later timepoints
                if ~has_3d && k > 1
                    fx1_mat = gtv_locations{j, 1, 1};
                    if ~isempty(fx1_mat) && exist(fx1_mat, 'file')
                        fx1_mask_3d = safe_load_mask(fx1_mat, 'Stvol3d');
                        adc_ref = pipe_vecs{1,1};
                        if ~isempty(fx1_mask_3d) && ~isempty(adc_ref) && ...
                                sum(fx1_mask_3d(:) == 1) == numel(adc_ref)
                            gtv_mask_3d = fx1_mask_3d;
                            has_3d = true;
                        end
                    end
                end
            end

            % Build core_opts for fDM
            core_opts = struct('timepoint_index', k);
            if k > 1
                % Baseline vectors for fDM (use Standard pipeline baseline)
                core_opts.baseline_adc_vec = data_vectors_gtvp(j, 1, 1).adc_vector;
                core_opts.baseline_d_vec = data_vectors_gtvp(j, 1, 1).d_vector;
            end

            % Suppress warnings during method loop
            prev_warn_state = suppress_core_warnings();

            for p = 1:n_pipelines
                if pipe_skipped(p, j, k)
                    continue;
                end

                adc_p = pipe_vecs{p,1};
                d_p = pipe_vecs{p,2};
                f_p = pipe_vecs{p,3};
                dstar_p = pipe_vecs{p,4};

                for m = 1:n_methods
                    rng(42);
                    temp_config = config_struct;
                    temp_config.core_method = ALL_METHODS{m};

                    [~, fit_info] = extract_tumor_core(temp_config, ...
                        adc_p, d_p, f_p, dstar_p, has_3d, gtv_mask_3d, core_opts);

                    was_fallback(m, p, j, k) = fit_info.fallback;
                    was_empty(m, p, j, k) = fit_info.empty_mask;
                    was_insufficient(m, p, j, k) = fit_info.insufficient_voxels;
                    was_all_nan(m, p, j, k) = fit_info.all_nan_input;
                    core_voxel_counts(m, p, j, k) = fit_info.n_core_voxels;
                end
            end

            % Restore warnings
            warning(prev_warn_state);
        end
    end

    % --- Compute aggregate rates ---
    fallback_rate = nan(n_methods, n_pipelines);
    empty_rate = nan(n_methods, n_pipelines);
    insufficient_rate = nan(n_methods, n_pipelines);
    all_nan_rate = nan(n_methods, n_pipelines);
    any_failure_rate = nan(n_methods, n_pipelines);
    median_core_voxels = nan(n_methods, n_pipelines);

    for m = 1:n_methods
        for p = 1:n_pipelines
            % Count non-skipped entries
            skipped = squeeze(pipe_skipped(p, :, :));
            n_total = sum(~skipped(:));
            if n_total == 0
                continue;
            end

            fb = squeeze(was_fallback(m, p, :, :));
            em = squeeze(was_empty(m, p, :, :));
            ins = squeeze(was_insufficient(m, p, :, :));
            an = squeeze(was_all_nan(m, p, :, :));
            vc = squeeze(core_voxel_counts(m, p, :, :));

            % Only count non-skipped entries
            valid = ~skipped(:);
            fallback_rate(m, p) = sum(fb(valid)) / n_total;
            empty_rate(m, p) = sum(em(valid)) / n_total;
            insufficient_rate(m, p) = sum(ins(valid)) / n_total;
            all_nan_rate(m, p) = sum(an(valid)) / n_total;
            any_failure_rate(m, p) = sum(fb(valid) | em(valid) | ins(valid) | an(valid)) / n_total;

            % Median core voxels for successful runs
            successful = valid & ~fb(:) & ~em(:) & ~an(:);
            if any(successful)
                median_core_voxels(m, p) = nanmedian(vc(successful));
            end
        end
    end

    % --- Package results ---
    failure_table.method_names = ALL_METHODS;
    failure_table.pipeline_names = PIPELINE_NAMES;
    failure_table.n_patients = n_patients;
    failure_table.n_timepoints = nTp;
    failure_table.fallback_rate = fallback_rate;
    failure_table.empty_rate = empty_rate;
    failure_table.insufficient_rate = insufficient_rate;
    failure_table.all_nan_rate = all_nan_rate;
    failure_table.any_failure_rate = any_failure_rate;
    failure_table.median_core_voxels = median_core_voxels;
    failure_table.core_voxel_counts = core_voxel_counts;
    failure_table.was_fallback = was_fallback;
    failure_table.was_empty = was_empty;
    failure_table.was_insufficient = was_insufficient;
    failure_table.was_all_nan = was_all_nan;

    % --- Print summary table ---
    for p = 1:n_pipelines
        fprintf('\nCore Method Failure Rates (%s pipeline):\n', PIPELINE_NAMES{p});
        fprintf('%-25s %8s %8s %8s %8s %10s\n', ...
            'Method', 'Fallback', 'Empty', 'Insuff', 'NaN', 'Total Fail');
        fprintf('%s\n', repmat('-', 1, 75));
        for m = 1:n_methods
            fprintf('%-25s %7.1f%% %7.1f%% %7.1f%% %7.1f%% %9.1f%%\n', ...
                ALL_METHODS{m}, ...
                fallback_rate(m, p) * 100, ...
                empty_rate(m, p) * 100, ...
                insufficient_rate(m, p) * 100, ...
                all_nan_rate(m, p) * 100, ...
                any_failure_rate(m, p) * 100);
        end
    end

    % --- Generate stacked bar chart ---
    try
        fig = figure('Visible', 'off', 'Position', [100 100 900 500]);
        % Use Standard pipeline (p=1) for the figure
        p_idx = 1;
        bar_data = [fallback_rate(:, p_idx), empty_rate(:, p_idx), ...
                    insufficient_rate(:, p_idx), all_nan_rate(:, p_idx)] * 100;
        b = bar(bar_data, 'stacked');
        b(1).FaceColor = [1.0 0.6 0.2];  % orange - fallback
        b(2).FaceColor = [0.9 0.2 0.2];  % red - empty
        b(3).FaceColor = [1.0 0.9 0.3];  % yellow - insufficient
        b(4).FaceColor = [0.6 0.6 0.6];  % gray - all_nan
        set(gca, 'XTick', 1:n_methods, 'XTickLabel', ALL_METHODS, ...
            'FontSize', 7, 'XTickLabelRotation', 45);
        ylabel('Failure Rate (%)');
        title(sprintf('Core Method Failure Rates (%s)', config_struct.dwi_type_name));
        legend({'Fallback', 'Empty mask', 'Insufficient voxels', 'All-NaN'}, ...
            'Location', 'northeastoutside');
        ylim([0 max(100, max(sum(bar_data, 2)) * 1.1)]);
        saveas(fig, fullfile(output_folder, ...
            sprintf('core_failure_rates_%s.png', config_struct.dwi_type_name)));
        close(fig);
    catch ME_fig
        fprintf('⚠️ Could not generate failure rates figure: %s\n', ME_fig.message);
    end

    fprintf('\n✅ Core failure rate analysis complete: %d patients, %d timepoints.\n', ...
        n_patients, nTp);

    diary off;
end
