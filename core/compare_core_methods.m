function compare_results = compare_core_methods(data_vectors_gtvp, summary_metrics, config_struct)
% COMPARE_CORE_METHODS  Pairwise comparison of all 11 tumor core methods.
%
% Runs extract_tumor_core once per method for each patient/timepoint,
% computes pairwise Dice coefficients and (where 3D masks are available)
% Hausdorff distances, and generates summary heatmaps and volume bar charts.
%
% This is a NON-FATAL pipeline step: errors are caught by the orchestrator.
%
% Inputs:
%   data_vectors_gtvp  - struct array with voxel-level parameter vectors
%   summary_metrics    - struct with id_list, gtv_locations, etc.
%   config_struct      - pipeline configuration struct
%
% Outputs:
%   compare_results    - struct with method_names, mean_dice_matrix,
%                        mean_hd95_matrix, volume_fractions, fallback_flags,
%                        all_dice, all_hd95

    % --- Diary ---
    output_folder = config_struct.output_folder;
    diary_file = fullfile(output_folder, ...
        sprintf('compare_core_methods_output_%s.txt', config_struct.dwi_type_name));
    diary(diary_file);

    fprintf('\n🔬 Core Method Comparison (%s)\n', config_struct.dwi_type_name);
    fprintf('   Comparing all 11 tumor core delineation methods.\n\n');

    % --- Constants ---
    % All 11 tumor core delineation methods implemented in extract_tumor_core.m.
    % Each method identifies the "treatment-resistant core" — the subregion
    % of the GTV with the most restricted diffusion (low ADC/D), highest
    % cellularity, or most distinct tissue characteristics. Different methods
    % use different assumptions about what constitutes the core:
    %   - Threshold methods (adc/d/df): hard cutoff on parameter values
    %   - Statistical methods (otsu/gmm/kmeans): unsupervised clustering
    %   - Spatial methods (region_growing/active_contours): 3D connectivity
    %   - Temporal methods (fdm): change from baseline (functional diffusion map)
    ALL_METHODS = {'adc_threshold', 'd_threshold', 'df_intersection', ...
        'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
        'percentile', 'spectral', 'fdm'};
    n_methods = numel(ALL_METHODS);

    % Fix random seed for reproducible clustering results (GMM, k-means, spectral).
    % Without this, different runs could produce different cluster assignments,
    % making the Dice/Hausdorff comparison non-reproducible.
    rng(42);

    id_list = summary_metrics.id_list;         % patient identifiers for logging
    gtv_locations = summary_metrics.gtv_locations;  % cell array of GTV .mat file paths
    n_patients = numel(id_list);
    nTp = size(data_vectors_gtvp, 2);          % number of timepoints (fractions)
    dwi_type = config_struct.dwi_types_to_run; % 1=Standard, 2=DnCNN, 3=IVIMnet

    % --- Pre-allocate ---
    % Running sums and counts for incremental mean computation:
    % mean = sum / count, computed after the main loop completes.
    dice_sum = zeros(n_methods, n_methods);    % running sum of pairwise Dice coefficients
    dice_count = zeros(n_methods, n_methods);  % number of valid Dice observations per pair
    hd95_sum = zeros(n_methods, n_methods);    % running sum of pairwise HD95 distances (mm)
    hd95_count = zeros(n_methods, n_methods);  % number of valid HD95 observations per pair

    % Per-patient-timepoint storage for detailed analysis
    volume_fractions = nan(n_patients, nTp, n_methods);   % core volume as fraction of GTV
    fallback_flags = false(n_patients, nTp, n_methods);   % true when method fell back to adc_threshold

    % Cell arrays store the full pairwise matrices for each patient x timepoint
    all_dice = cell(n_patients, nTp);   % each cell: n_methods x n_methods Dice matrix
    all_hd95 = cell(n_patients, nTp);   % each cell: n_methods x n_methods HD95 matrix (mm)

    % --- Main loop: patient x timepoint ---
    for j = 1:n_patients
        text_progress_bar(j, n_patients, 'Comparing core methods');
        for k = 1:nTp
            % Extract the appropriate voxel-level parameter vectors based on
            % the DWI processing pipeline. Each type stores vectors in
            % different struct fields within data_vectors_gtvp.
            switch dwi_type
                case 1  % Standard (conventional fitting on raw DWI)
                    adc_vec = data_vectors_gtvp(j,k,1).adc_vector;
                    d_vec = data_vectors_gtvp(j,k,1).d_vector;
                    f_vec = data_vectors_gtvp(j,k,1).f_vector;
                    dstar_vec = data_vectors_gtvp(j,k,1).dstar_vector;
                    adc_baseline = data_vectors_gtvp(j,1,1).adc_vector;
                    d_baseline = data_vectors_gtvp(j,1,1).d_vector;
                case 2  % dnCNN (fitting on DnCNN-denoised DWI)
                    adc_vec = data_vectors_gtvp(j,k,1).adc_vector_dncnn;
                    d_vec = data_vectors_gtvp(j,k,1).d_vector_dncnn;
                    f_vec = data_vectors_gtvp(j,k,1).f_vector_dncnn;
                    dstar_vec = data_vectors_gtvp(j,k,1).dstar_vector_dncnn;
                    adc_baseline = data_vectors_gtvp(j,1,1).adc_vector_dncnn;
                    d_baseline = data_vectors_gtvp(j,1,1).d_vector_dncnn;
                case 3  % IVIMnet (neural network IVIM estimation; ADC still from standard)
                    adc_vec = data_vectors_gtvp(j,k,1).adc_vector;
                    d_vec = data_vectors_gtvp(j,k,1).d_vector_ivimnet;
                    f_vec = data_vectors_gtvp(j,k,1).f_vector_ivimnet;
                    dstar_vec = data_vectors_gtvp(j,k,1).dstar_vector_ivimnet;
                    adc_baseline = data_vectors_gtvp(j,1,1).adc_vector;
                    d_baseline = data_vectors_gtvp(j,1,1).d_vector_ivimnet;
            end

            % Skip empty/all-NaN data
            if isempty(adc_vec) || all(isnan(adc_vec))
                continue;
            end

            % Load 3D GTV mask.  The mask must have exactly as many
            % true-voxels as the parameter vector length, because spatial
            % methods (active_contours, region_growing) reconstruct 3D
            % maps via: map_3d(mask == 1) = vec.  A size mismatch would
            % assign voxels to wrong spatial positions.
            has_3d = false;
            gtv_mask_3d = [];
            if ~isempty(gtv_locations) && ...
                    size(gtv_locations, 1) >= j && size(gtv_locations, 2) >= k
                gtv_mat = gtv_locations{j, k, 1};
                if ~isempty(gtv_mat) && exist(gtv_mat, 'file')
                    gtv_mask_3d = safe_load_mask(gtv_mat, 'Stvol3d');
                    if ~isempty(gtv_mask_3d) && sum(gtv_mask_3d(:) == 1) == numel(adc_vec)
                        has_3d = true;
                    end
                end
                % Fx1 mask fallback: at Fx2+, parameter maps are warped
                % to Fx1 geometry via DIR (imregdemons), so vectors were
                % extracted using the Fx1 mask.  gtv_locations{j,k,1}
                % points to the native fraction mask (different voxel
                % count).  The Fx1 mask has the matching geometry.
                if ~has_3d && k > 1
                    fx1_mat = gtv_locations{j, 1, 1};
                    if ~isempty(fx1_mat) && exist(fx1_mat, 'file')
                        fx1_mask_3d = safe_load_mask(fx1_mat, 'Stvol3d');
                        if ~isempty(fx1_mask_3d) && sum(fx1_mask_3d(:) == 1) == numel(adc_vec)
                            gtv_mask_3d = fx1_mask_3d;
                            has_3d = true;
                        end
                    end
                end
                % Diagnostic: log when 3D mask is unavailable so the
                % source of active_contours/region_growing fallbacks
                % can be traced to specific patients/timepoints.
                if ~has_3d && ~isempty(gtv_mask_3d)
                    fprintf('  [3D mask] %s Fx%d: mask has %d voxels, vector has %d (mismatch)\n', ...
                        id_list{j}, k, sum(gtv_mask_3d(:) == 1), numel(adc_vec));
                elseif ~has_3d
                    fprintf('  [3D mask] %s Fx%d: no mask loaded\n', id_list{j}, k);
                end
            end

            % Voxel dimensions for Hausdorff
            vox_dims = [1 1 1];
            if isfield(data_vectors_gtvp(j,k,1), 'vox_dims') && ...
                    isnumeric(data_vectors_gtvp(j,k,1).vox_dims) && ...
                    numel(data_vectors_gtvp(j,k,1).vox_dims) == 3
                vox_dims = data_vectors_gtvp(j,k,1).vox_dims;
            end

            % Build opts for fDM (functional Diffusion Map) method.
            % fDM identifies core voxels by comparing current-fraction
            % parameter values to baseline (Fx1). At Fx1 (k=1), no baseline
            % comparison is possible, so fDM falls back to adc_threshold.
            core_opts = struct('timepoint_index', k);
            if k > 1
                core_opts.baseline_adc_vec = adc_baseline;
                core_opts.baseline_d_vec = d_baseline;
            end

            % --- Run all 11 methods ---
            % Reuse pre-computed masks from compute_summary_metrics when
            % available (run_all_core_methods + store_core_masks were true).
            has_precomputed = isfield(summary_metrics, 'all_core_metrics') && ...
                isstruct(summary_metrics.all_core_metrics) && ...
                isfield(summary_metrics.all_core_metrics, 'adc_threshold') && ...
                isfield(summary_metrics.all_core_metrics.adc_threshold, 'core_masks');

            masks_1d = cell(n_methods, 1);
            temp_config = config_struct;

            % Suppress expected fallback warnings during comparison
            prev_warn_state = warning('query');
            warning('off', 'extract_tumor_core:tooFewForSpectral');
            warning('off', 'extract_tumor_core:no3DForActiveContours');
            warning('off', 'extract_tumor_core:no3DForRegionGrowing');
            warning('off', 'extract_tumor_core:fdmBaseline');
            warning('off', 'extract_tumor_core:fdmNoBaseline');
            warning('off', 'extract_tumor_core:noDValues');
            warning('off', 'extract_tumor_core:noIVIMValues');
            warning('off', 'extract_tumor_core:noSpectralCluster');

            for m = 1:n_methods
                % Try pre-computed mask first
                if has_precomputed && isfield(summary_metrics.all_core_metrics, ALL_METHODS{m}) && ...
                        j <= size(summary_metrics.all_core_metrics.(ALL_METHODS{m}).core_masks, 1) && ...
                        k <= size(summary_metrics.all_core_metrics.(ALL_METHODS{m}).core_masks, 2) && ...
                        ~isempty(summary_metrics.all_core_metrics.(ALL_METHODS{m}).core_masks{j,k})
                    masks_1d{m} = summary_metrics.all_core_metrics.(ALL_METHODS{m}).core_masks{j,k};
                else
                    temp_config.core_method = ALL_METHODS{m};
                    try
                        masks_1d{m} = extract_tumor_core(temp_config, adc_vec, ...
                            d_vec, f_vec, dstar_vec, has_3d, gtv_mask_3d, core_opts);
                    catch
                        masks_1d{m} = false(size(adc_vec));
                    end
                end
                % Force column vector to prevent implicit expansion in
                % pairwise Dice/Hausdorff computations (row & column would
                % produce a matrix via broadcasting).
                masks_1d{m} = masks_1d{m}(:);

                % Volume fraction
                n_finite = sum(~isnan(adc_vec));
                if n_finite > 0
                    volume_fractions(j, k, m) = sum(masks_1d{m}) / n_finite;
                end

                % Detect fallback to adc_threshold
                if m > 1 && isequal(masks_1d{m}, masks_1d{1})
                    needs_3d = ismember(ALL_METHODS{m}, {'region_growing', 'active_contours'});
                    needs_baseline = strcmp(ALL_METHODS{m}, 'fdm') && k <= 1;
                    if (needs_3d && ~has_3d) || needs_baseline
                        fallback_flags(j, k, m) = true;
                    end
                end
            end

            % Restore warning state
            warning(prev_warn_state);

            % --- Pairwise Dice (1D) ---
            % Dice coefficient measures spatial overlap between two binary masks:
            %   Dice = 2 * |A AND B| / (|A| + |B|)
            % Dice = 1.0 → identical masks, Dice = 0 → no overlap.
            % Computed on 1D (voxel-vector) masks; does not require 3D geometry.
            dice_matrix = nan(n_methods, n_methods);
            for a = 1:n_methods
                for b = a:n_methods
                    sum_a = sum(masks_1d{a});   % number of core voxels in method A
                    sum_b = sum(masks_1d{b});   % number of core voxels in method B
                    if sum_a == 0 && sum_b == 0
                        dice_matrix(a, b) = NaN;  % both empty: Dice is undefined
                    elseif sum_a == 0 || sum_b == 0
                        dice_matrix(a, b) = 0;    % one empty: no overlap possible
                    else
                        dice_matrix(a, b) = 2 * sum(masks_1d{a} & masks_1d{b}) / (sum_a + sum_b);
                    end
                    dice_matrix(b, a) = dice_matrix(a, b);  % symmetric
                end
            end
            all_dice{j, k} = dice_matrix;

            valid = ~isnan(dice_matrix);
            dice_sum(valid) = dice_sum(valid) + dice_matrix(valid);
            dice_count(valid) = dice_count(valid) + 1;

            % --- Pairwise Hausdorff (3D, when available) ---
            % HD95 (95th percentile Hausdorff distance) measures the maximum
            % boundary disagreement between two masks in mm, robust to outliers.
            % Requires 3D mask reconstruction from 1D vectors, which is only
            % possible when we have the original 3D GTV mask with matching voxel count.
            hd95_matrix = nan(n_methods, n_methods);
            n_gtv_voxels = 0;
            if has_3d
                n_gtv_voxels = sum(gtv_mask_3d(:) == 1);
            end
            if has_3d && n_gtv_voxels == numel(masks_1d{1})
                % Reconstruct 1D core masks back into 3D volumes by placing each
                % voxel-vector value into its spatial position within the GTV.
                masks_3d = cell(n_methods, 1);
                for m = 1:n_methods
                    vol_3d = false(size(gtv_mask_3d));
                    vol_3d(gtv_mask_3d == 1) = masks_1d{m};  % map 1D→3D via GTV linear indices
                    masks_3d{m} = vol_3d;
                end

                for a = 1:n_methods
                    for b = (a+1):n_methods
                        [~, ~, hd95_val] = compute_dice_hausdorff(masks_3d{a}, masks_3d{b}, vox_dims);
                        hd95_matrix(a, b) = hd95_val;
                        hd95_matrix(b, a) = hd95_val;
                    end
                    hd95_matrix(a, a) = 0;
                end

                valid_hd = ~isnan(hd95_matrix) & ~isinf(hd95_matrix);
                hd95_sum(valid_hd) = hd95_sum(valid_hd) + hd95_matrix(valid_hd);
                hd95_count(valid_hd) = hd95_count(valid_hd) + 1;
            end
            all_hd95{j, k} = hd95_matrix;
        end
    end

    % --- Mean matrices ---
    % Compute element-wise mean from running sums. max(count, 1) prevents
    % division by zero; pairs with zero observations are then set to NaN.
    mean_dice_matrix = dice_sum ./ max(dice_count, 1);
    mean_dice_matrix(dice_count == 0) = NaN;

    mean_hd95_matrix = hd95_sum ./ max(hd95_count, 1);
    mean_hd95_matrix(hd95_count == 0) = NaN;

    % --- Package results ---
    compare_results = struct();
    compare_results.method_names = ALL_METHODS;
    compare_results.mean_dice_matrix = mean_dice_matrix;
    compare_results.mean_hd95_matrix = mean_hd95_matrix;
    compare_results.dice_count = dice_count;
    compare_results.hd95_count = hd95_count;
    compare_results.volume_fractions = volume_fractions;
    compare_results.fallback_flags = fallback_flags;
    compare_results.all_dice = all_dice;
    compare_results.all_hd95 = all_hd95;
    compare_results.n_patients = n_patients;
    compare_results.nTp = nTp;
    compare_results.dwi_type_name = config_struct.dwi_type_name;

    % --- Save MAT file ---
    results_mat = fullfile(output_folder, ...
        sprintf('compare_core_results_%s.mat', config_struct.dwi_type_name));
    save(results_mat, 'compare_results');
    fprintf('  📁 Saved comparison results to %s\n', results_mat);

    % --- Figures ---
    generate_comparison_figures(compare_results, output_folder, config_struct.dwi_type_name);

    fprintf('✅ Core method comparison complete.\n');
    diary off;
end


%% ---- Local helper functions ----

function generate_comparison_figures(results, output_folder, dwi_type_name)
% Generate heatmaps and bar charts for core method comparison.

    method_labels = strrep(results.method_names, '_', ' ');
    n_methods = numel(results.method_names);

    % --- Figure 1: Mean Dice Heatmap ---
    fig1 = figure('Visible', 'off', 'Position', [100 100 800 700]);
    imagesc(results.mean_dice_matrix);
    colorbar;
    colormap(parula);
    caxis([0 1]);
    set(gca, 'XTick', 1:n_methods, 'XTickLabel', method_labels, ...
        'YTick', 1:n_methods, 'YTickLabel', method_labels, ...
        'XTickLabelRotation', 45, 'FontSize', 8);
    title(sprintf('Mean Pairwise Dice Coefficient (%s)', dwi_type_name));

    % Overlay numeric values
    for i = 1:n_methods
        for j = 1:n_methods
            val = results.mean_dice_matrix(i, j);
            if ~isnan(val)
                if val > 0.5
                    txt_color = [0 0 0];
                else
                    txt_color = [1 1 1];
                end
                text(j, i, sprintf('%.2f', val), 'HorizontalAlignment', 'center', ...
                    'FontSize', 7, 'Color', txt_color);
            end
        end
    end

    saveas(fig1, fullfile(output_folder, sprintf('core_method_dice_heatmap_%s.png', dwi_type_name)));
    close(fig1);

    % --- Figure 2: Mean HD95 Heatmap (if data available) ---
    if any(results.hd95_count(:) > 0)
        fig2 = figure('Visible', 'off', 'Position', [100 100 800 700]);
        imagesc(results.mean_hd95_matrix);
        colorbar;
        colormap(hot);
        set(gca, 'XTick', 1:n_methods, 'XTickLabel', method_labels, ...
            'YTick', 1:n_methods, 'YTickLabel', method_labels, ...
            'XTickLabelRotation', 45, 'FontSize', 8);
        title(sprintf('Mean Pairwise HD95 in mm (%s)', dwi_type_name));

        for i = 1:n_methods
            for j = 1:n_methods
                val = results.mean_hd95_matrix(i, j);
                if ~isnan(val)
                    text(j, i, sprintf('%.1f', val), 'HorizontalAlignment', 'center', ...
                        'FontSize', 7);
                end
            end
        end

        saveas(fig2, fullfile(output_folder, sprintf('core_method_hd95_heatmap_%s.png', dwi_type_name)));
        close(fig2);
    end

    % --- Figure 3: Volume Fraction Bar Chart ---
    % Reshape to (observations x methods) for nanmean/nanstd
    vf_2d = reshape(results.volume_fractions, [], n_methods);
    mean_vol_frac = nanmean(vf_2d, 1);
    std_vol_frac = nanstd(vf_2d, 0, 1);

    fig3 = figure('Visible', 'off', 'Position', [100 100 900 500]);
    bar(mean_vol_frac(:) * 100);
    set(gca, 'XTick', 1:n_methods, 'XTickLabel', method_labels, ...
        'XTickLabelRotation', 45, 'FontSize', 8);
    ylabel('Mean Core Volume (% of GTV)');
    title(sprintf('Core Volume by Method (%s)', dwi_type_name));
    hold on;
    errorbar(1:n_methods, mean_vol_frac(:) * 100, std_vol_frac(:) * 100, '.k');
    hold off;

    saveas(fig3, fullfile(output_folder, sprintf('core_method_volume_comparison_%s.png', dwi_type_name)));
    close(fig3);

    % --- Figure 4: Fallback Summary (if any) ---
    n_fallbacks = squeeze(sum(sum(results.fallback_flags, 1), 2));
    if any(n_fallbacks > 0)
        n_total = results.n_patients * results.nTp;
        fig4 = figure('Visible', 'off', 'Position', [100 100 900 400]);
        bar(n_fallbacks(:));
        set(gca, 'XTick', 1:n_methods, 'XTickLabel', method_labels, ...
            'XTickLabelRotation', 45, 'FontSize', 8);
        ylabel('Number of Fallbacks');
        title(sprintf('Fallbacks to ADC Threshold (%s, N=%d patient-scans)', dwi_type_name, n_total));
        saveas(fig4, fullfile(output_folder, sprintf('core_method_fallbacks_%s.png', dwi_type_name)));
        close(fig4);
    end

    fprintf('  📁 Generated comparison figures in %s\n', output_folder);
end
