function dice_results = compute_cross_pipeline_dice(data_vectors_gtvp, config_struct, id_list, gtv_locations)
% COMPUTE_CROSS_PIPELINE_DICE  Cross-pipeline Dice coefficients at fraction 1.
%
% Computes pairwise Dice coefficients between the 3 DWI processing
% pipelines (Standard, DnCNN, IVIMNet) for each of the 11 tumor core
% delineation methods at fraction 1 (baseline).
%
% This measures whether the choice of DWI processing pipeline — raw
% (Standard), DnCNN-denoised, or IVIMNet neural-network-estimated —
% changes the spatial location of the treatment-resistant sub-volume.
% High Dice (>0.7) indicates that the core identification is robust to
% the denoising/fitting method; low Dice (<0.5) indicates that
% pipeline choice materially affects which voxels are classified as
% resistant, with downstream implications for dose escalation planning.
%
% This is a NON-FATAL pipeline step: errors are caught by the orchestrator.
%
% NOTE: This does NOT modify or replace compare_core_methods.m, which
% handles cross-METHOD Dice (same pipeline, different core algorithms).
% This function handles cross-PIPELINE Dice (same core algorithm,
% different DWI processing pipelines).
%
% Inputs:
%   data_vectors_gtvp  - struct array [nPatients x nTp x 1] with voxel-
%                         level parameter vectors for all 3 pipelines
%   config_struct      - pipeline configuration struct
%   id_list            - cell array of patient identifiers
%   gtv_locations      - cell array of GTV mask file paths
%
% Outputs:
%   dice_results       - struct with fields:
%     .dice               - [11 x 3 x nPatients] Dice per (method, pair, patient)
%     .method_names       - cell array of 11 method names
%     .pipeline_names     - {'Standard', 'DnCNN', 'IVIMNet'}
%     .pipeline_pair_labels - {'Std_vs_DnCNN', 'Std_vs_IVIMNet', 'DnCNN_vs_IVIMNet'}
%     .patient_ids        - cell array of patient identifiers
%     .n_voxels           - [11 x 3 x nPatients] core mask voxel counts
%     .fallback_flags     - [11 x 3 x nPatients] logical fallback indicators
%     .n_patients         - scalar

    % --- Diary ---
    output_folder = config_struct.output_folder;
    diary_file = fullfile(output_folder, ...
        sprintf('cross_pipeline_dice_output_%s.txt', config_struct.dwi_type_name));
    diary(diary_file);

    fprintf('\n🔬 Cross-Pipeline Dice Coefficients at Fx1 (%s)\n', config_struct.dwi_type_name);
    fprintf('   Comparing Standard vs DnCNN vs IVIMNet for all 11 core methods.\n\n');

    % --- Constants ---
    ALL_METHODS = {'adc_threshold', 'd_threshold', 'df_intersection', ...
        'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
        'percentile', 'spectral', 'fdm'};
    PIPELINE_NAMES = {'Standard', 'DnCNN', 'IVIMNet'};
    PIPELINE_PAIR_LABELS = {'Std_vs_DnCNN', 'Std_vs_IVIMNet', 'DnCNN_vs_IVIMNet'};
    n_methods = 11;
    n_pipelines = 3;
    n_pairs = 3;
    n_patients = numel(id_list);

    % --- Pre-allocate ---
    dice = nan(n_methods, n_pairs, n_patients);
    n_voxels = nan(n_methods, n_pipelines, n_patients);
    fallback_flags = false(n_methods, n_pipelines, n_patients);

    % --- Main loop: patients ---
    for j = 1:n_patients
        text_progress_bar(j, n_patients, 'Cross-pipeline Dice');

        % Fix random seed per patient for reproducibility of stochastic methods
        rng(42);

        % --- Extract parameter vectors for each pipeline at Fx1 (k=1) ---
        entry = data_vectors_gtvp(j, 1, 1);

        % Standard pipeline
        std_adc = entry.adc_vector;
        std_d = entry.d_vector;
        std_f = entry.f_vector;
        std_dstar = entry.dstar_vector;

        % DnCNN pipeline
        dncnn_adc = [];
        dncnn_d = [];
        dncnn_f = [];
        dncnn_dstar = [];
        if isfield(entry, 'adc_vector_dncnn')
            dncnn_adc = entry.adc_vector_dncnn;
        end
        if isfield(entry, 'd_vector_dncnn')
            dncnn_d = entry.d_vector_dncnn;
        end
        if isfield(entry, 'f_vector_dncnn')
            dncnn_f = entry.f_vector_dncnn;
        end
        if isfield(entry, 'dstar_vector_dncnn')
            dncnn_dstar = entry.dstar_vector_dncnn;
        end

        % IVIMNet pipeline (ADC shared with Standard)
        ivimnet_adc = entry.adc_vector;  % shared
        ivimnet_d = [];
        ivimnet_f = [];
        ivimnet_dstar = [];
        if isfield(entry, 'd_vector_ivimnet')
            ivimnet_d = entry.d_vector_ivimnet;
        end
        if isfield(entry, 'f_vector_ivimnet')
            ivimnet_f = entry.f_vector_ivimnet;
        end
        if isfield(entry, 'dstar_vector_ivimnet')
            ivimnet_dstar = entry.dstar_vector_ivimnet;
        end

        % Package pipeline vectors: {adc, d, f, dstar} per pipeline
        pipe_vecs = cell(n_pipelines, 4);
        pipe_vecs(1,:) = {std_adc, std_d, std_f, std_dstar};
        pipe_vecs(2,:) = {dncnn_adc, dncnn_d, dncnn_f, dncnn_dstar};
        pipe_vecs(3,:) = {ivimnet_adc, ivimnet_d, ivimnet_f, ivimnet_dstar};

        % Check which pipelines have data
        pipe_has_data = false(n_pipelines, 1);
        for p = 1:n_pipelines
            pipe_has_data(p) = ~isempty(pipe_vecs{p,1}) && ~all(isnan(pipe_vecs{p,1}));
        end

        % Need at least 2 pipelines to compute any pairwise Dice
        if sum(pipe_has_data) < 2
            continue;
        end

        % --- Load 3D GTV mask (follow compare_core_methods pattern) ---
        has_3d = false;
        gtv_mask_3d = [];
        if ~isempty(gtv_locations) && ...
                size(gtv_locations, 1) >= j && size(gtv_locations, 2) >= 1
            gtv_mat = gtv_locations{j, 1, 1};
            if ~isempty(gtv_mat) && exist(gtv_mat, 'file')
                gtv_mask_3d = safe_load_mask(gtv_mat, 'Stvol3d');
                if ~isempty(gtv_mask_3d) && ~isempty(std_adc) && ...
                        sum(gtv_mask_3d(:) == 1) == numel(std_adc)
                    has_3d = true;
                end
            end
        end

        % Build core_opts for fDM (at Fx1, no baseline available)
        core_opts = struct('timepoint_index', 1);

        % Suppress expected warnings during method loop
        prev_warn_state = suppress_core_warnings();

        % --- Compute adc_threshold masks first (reference for fallback detection) ---
        masks_adc_thresh = cell(n_pipelines, 1);
        temp_config = config_struct;
        temp_config.core_method = 'adc_threshold';
        for p = 1:n_pipelines
            if ~pipe_has_data(p)
                masks_adc_thresh{p} = [];
                continue;
            end
            try
                masks_adc_thresh{p} = extract_tumor_core(temp_config, ...
                    pipe_vecs{p,1}, pipe_vecs{p,2}, pipe_vecs{p,3}, pipe_vecs{p,4}, ...
                    has_3d, gtv_mask_3d, core_opts);
            catch
                masks_adc_thresh{p} = false(size(pipe_vecs{p,1}));
            end
            masks_adc_thresh{p} = masks_adc_thresh{p}(:);
        end

        % --- Run all 11 methods across all 3 pipelines ---
        for m = 1:n_methods
            temp_config = config_struct;
            temp_config.core_method = ALL_METHODS{m};

            masks = cell(n_pipelines, 1);
            for p = 1:n_pipelines
                if ~pipe_has_data(p)
                    masks{p} = [];
                    continue;
                end

                % Reuse adc_threshold result for method 1
                if m == 1
                    masks{p} = masks_adc_thresh{p};
                else
                    try
                        masks{p} = extract_tumor_core(temp_config, ...
                            pipe_vecs{p,1}, pipe_vecs{p,2}, pipe_vecs{p,3}, pipe_vecs{p,4}, ...
                            has_3d, gtv_mask_3d, core_opts);
                    catch
                        masks{p} = false(size(pipe_vecs{p,1}));
                    end
                    masks{p} = masks{p}(:);
                end

                n_voxels(m, p, j) = sum(masks{p});

                % Detect fallback to adc_threshold
                if m > 1 && ~isempty(masks_adc_thresh{p}) && isequal(masks{p}, masks_adc_thresh{p})
                    fallback_flags(m, p, j) = true;
                end
            end

            % Record n_voxels for method 1 (adc_threshold)
            if m == 1
                for p = 1:n_pipelines
                    if ~isempty(masks{p})
                        n_voxels(m, p, j) = sum(masks{p});
                    end
                end
            end

            % --- Compute pairwise Dice between pipelines ---
            % Pair 1: Standard vs DnCNN
            if ~isempty(masks{1}) && ~isempty(masks{2})
                dice(m, 1, j) = dice_coeff(masks{1}, masks{2});
            end
            % Pair 2: Standard vs IVIMNet
            if ~isempty(masks{1}) && ~isempty(masks{3})
                dice(m, 2, j) = dice_coeff(masks{1}, masks{3});
            end
            % Pair 3: DnCNN vs IVIMNet
            if ~isempty(masks{2}) && ~isempty(masks{3})
                dice(m, 3, j) = dice_coeff(masks{2}, masks{3});
            end
        end

        % Restore warnings
        warning(prev_warn_state);
    end

    % --- Package results ---
    dice_results.dice = dice;
    dice_results.method_names = ALL_METHODS;
    dice_results.pipeline_names = PIPELINE_NAMES;
    dice_results.pipeline_pair_labels = PIPELINE_PAIR_LABELS;
    dice_results.patient_ids = id_list;
    dice_results.n_voxels = n_voxels;
    dice_results.fallback_flags = fallback_flags;
    dice_results.n_patients = n_patients;

    % --- Save results ---
    save(fullfile(output_folder, ...
        sprintf('cross_pipeline_dice_%s.mat', config_struct.dwi_type_name)), 'dice_results');

    fprintf('\n✅ Cross-pipeline Dice complete: %d patients, %d methods, %d pipeline pairs.\n', ...
        n_patients, n_methods, n_pairs);

    diary off;
end


function d = dice_coeff(mask_a, mask_b)
% DICE_COEFF  Compute Dice coefficient between two binary masks.
%   Returns NaN if both masks are empty (Dice is undefined).
    sum_ab = sum(mask_a(:)) + sum(mask_b(:));
    if sum_ab == 0
        d = NaN;
    else
        d = 2 * sum(mask_a(:) & mask_b(:)) / sum_ab;
    end
end
