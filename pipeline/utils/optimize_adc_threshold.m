function opt_results = optimize_adc_threshold(data_vectors_gtvp, config_struct, id_list, gtv_locations)
%OPTIMIZE_ADC_THRESHOLD Sweep ADC threshold and pick optimum via three tactics.
%
%   opt_results = optimize_adc_threshold(data_vectors_gtvp, config_struct, ...
%                                        id_list, gtv_locations)
%
%   Sweeps ADC threshold from 0.8e-3 to 2.0e-3 in steps of 0.1e-3 (13 values).
%   For each candidate threshold the function computes, per patient, a
%   pairwise Fx1 repeat Dice and a sub-volume fraction, then aggregates
%   across the cohort.  Three independent optimisation tactics are then
%   reported; each has a different rationale, and reading them together
%   tells the analyst which threshold to favour for which goal.
%
%     Tactic 1 — Reproducibility (primary):
%       Threshold that maximises median Fx1-repeat Dice.  Picks the most
%       test-retest stable cut.
%
%     Tactic 2 — Volume inflection (secondary):
%       Threshold at the knee of the median sub-volume-fraction curve
%       (the most concave point of vol_frac vs threshold, found via
%       discrete second derivative on a 3-point smoothed curve).  Marks
%       the biological transition between dense tumour core and
%       surrounding tissue, where each additional 0.0001 in threshold
%       starts adding sharply less volume.
%
%     Tactic 3 — Outcome significance (tertiary):
%       Threshold that minimises the Wilcoxon rank-sum p-value for
%       sub-volume fraction between LC (LF==0) and LF (LF==1) patients.
%       Picks the cut that best discriminates outcomes.  Requires at
%       least 3 patients per group at that threshold.
%
%   Outputs are saved to opt_results and a figure is generated with
%   dual y-axis (median Dice left, median volume fraction right).
%   Vertical dashed lines mark the current default (0.001), proposed
%   (0.0016), and the three tactic-specific optima.
%
%   Inputs:
%       data_vectors_gtvp - struct array [nPatients x nTp x nRepeats].
%                           Each element should carry .LF (0 = local
%                           control, 1 = local failure) for tactic 3.
%       config_struct     - pipeline configuration (must contain output_folder,
%                           dwi_types_to_run, dwi_type_name). morph_se and
%                           morph_min_cc are honoured if present; otherwise
%                           sensible defaults are used.
%       id_list           - cell array of patient identifiers
%       gtv_locations     - cell array {nPatients x nTp x nRepeats} of
%                           GTV mask file paths
%
%   Output struct fields (Tactic 1 — Dice):
%       thresholds        - [1x13] threshold values in mm^2/s
%       median_dice       - [1x13] median Dice across patients per threshold
%       mean_dice         - [1x13] mean Dice across patients per threshold
%       std_dice          - [1x13] std Dice across patients per threshold
%       median_vol_frac   - [1x13] median sub-volume fraction per threshold
%       mean_vol_frac     - [1x13] mean sub-volume fraction per threshold
%       n_patients        - scalar: patient count with >=2 repeats
%       optimal_thresh    - scalar: threshold maximising median Dice
%       optimal_dice      - scalar
%       optimal_vol_frac  - scalar
%       per_patient_dice  - [nPatients x 13] per-patient mean Dice per threshold
%
%   Output struct fields (Tactic 2 — Inflection):
%       inflection_thresh    - scalar: threshold at the knee of vol_frac curve
%                              (NaN if curve has fewer than 3 finite points)
%       inflection_idx       - scalar: index of the knee in `thresholds`
%       inflection_curvature - scalar: discrete 2nd derivative at the knee
%                              (negative magnitude = concave-down saturation)
%       vol_frac_curvature   - [1x13] discrete 2nd derivative of smoothed
%                              median_vol_frac (endpoints are NaN)
%
%   Output struct fields (Tactic 3 — Significance):
%       significance_pvalues   - [1x13] Wilcoxon rank-sum p-values per threshold
%       significance_thresh    - scalar: threshold with smallest valid p-value
%                                (NaN if no threshold has both LC>=3 and LF>=3)
%       significance_pvalue    - scalar: minimum p-value
%       significance_n_lc      - scalar: LC group size used at that threshold
%       significance_n_lf      - scalar: LF group size used at that threshold
%       significance_metric    - string: 'wilcoxon ranksum on vol_frac'

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
    per_patient_lf = nan(n_patients, 1);  % 0 = LC, 1 = LF (Tactic 3)
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
            if isempty(gtv_path)
                has_all_3d = false; break;
            end
            gtv_path = normalize_path_preserving_roots(gtv_path);
            if ~exist(gtv_path, 'file')
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

        % Capture LF status for the significance tactic (Tactic 3).
        % data_vectors_gtvp(j,1,1).LF is the patient-level outcome.
        try
            lf_val = data_vectors_gtvp(j, 1, 1).LF;
            if isscalar(lf_val) && (lf_val == 0 || lf_val == 1)
                per_patient_lf(j) = lf_val;
            end
        catch
            % Field absent or unreadable — leave as NaN; this patient
            % will not contribute to the significance scan.
        end

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

    % --- Tactic 1: Optimal threshold (Dice) ---
    if all(isnan(median_dice))
        optimal_idx = 1;
    else
        [~, optimal_idx] = max(median_dice);
    end
    optimal_thresh = thresholds(optimal_idx);
    optimal_dice = median_dice(optimal_idx);
    optimal_vol_frac = median_vol_frac(optimal_idx);

    % --- Tactic 2: Volume-fraction inflection ---
    % The volume curve V(t) grows monotonically with threshold.  Its knee
    % (most concave point) marks the saturation transition between dense
    % tumour core and surrounding tissue.  We approximate the second
    % derivative on a 3-point smoothed copy of median_vol_frac and pick
    % the index that maximises -d^2V/dt^2 (most negative curvature).
    [inflection_thresh, inflection_idx, inflection_curvature, vol_frac_curvature] = ...
        local_inflection(thresholds, median_vol_frac);

    % --- Tactic 3: Outcome-significance ---
    % Per-threshold Wilcoxon rank-sum on per-patient sub-volume fraction
    % between LC (LF==0) and LF (LF==1) groups.  The optimal threshold
    % is the one with the smallest valid p-value (both group sizes >=3).
    [significance_thresh, significance_idx, significance_pvalue, ...
     significance_pvalues, significance_n_lc, significance_n_lf] = ...
        local_significance(thresholds, per_patient_vol_frac, per_patient_lf);
    significance_metric = 'wilcoxon ranksum on vol_frac';

    % --- Console output ---
    fprintf('\n  [Tactic 1] Reproducibility: %.4f mm\xc2\xb2/s (Dice=%.2f, vol_frac=%.1f%%)\n', ...
        optimal_thresh, optimal_dice, 100 * optimal_vol_frac);
    if ~isnan(inflection_thresh)
        fprintf('  [Tactic 2] Volume inflection: %.4f mm\xc2\xb2/s (curvature=%.4g, vol_frac=%.1f%%)\n', ...
            inflection_thresh, inflection_curvature, ...
            100 * safe_val(median_vol_frac, inflection_idx));
    else
        fprintf('  [Tactic 2] Volume inflection: insufficient data\n');
    end
    if ~isnan(significance_thresh)
        fprintf('  [Tactic 3] Outcome significance: %.4f mm\xc2\xb2/s (p=%.3g, n_LC=%d, n_LF=%d)\n', ...
            significance_thresh, significance_pvalue, significance_n_lc, significance_n_lf);
    else
        fprintf('  [Tactic 3] Outcome significance: insufficient LC/LF data\n');
    end

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
        text(optimal_thresh, yl_left(2)*0.80, sprintf(' T1: Dice (%.4f)', optimal_thresh), ...
            'Color', 'm', 'FontSize', 9);

        if ~isnan(inflection_thresh)
            plot([inflection_thresh, inflection_thresh], yl_left, '--', ...
                'Color', [0.4 0.2 0.8], 'LineWidth', 1.5);
            text(inflection_thresh, yl_left(2)*0.72, ...
                sprintf(' T2: Inflection (%.4f)', inflection_thresh), ...
                'Color', [0.4 0.2 0.8], 'FontSize', 9);
        end
        if ~isnan(significance_thresh)
            plot([significance_thresh, significance_thresh], yl_left, '--', ...
                'Color', [0.85 0.4 0.0], 'LineWidth', 1.5);
            text(significance_thresh, yl_left(2)*0.64, ...
                sprintf(' T3: Significance (%.4f)', significance_thresh), ...
                'Color', [0.85 0.4 0.0], 'FontSize', 9);
        end
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
    % Tactic 2 — Inflection
    opt_results.inflection_thresh    = inflection_thresh;
    opt_results.inflection_idx       = inflection_idx;
    opt_results.inflection_curvature = inflection_curvature;
    opt_results.vol_frac_curvature   = vol_frac_curvature;
    % Tactic 3 — Significance
    opt_results.significance_pvalues = significance_pvalues;
    opt_results.significance_thresh  = significance_thresh;
    opt_results.significance_pvalue  = significance_pvalue;
    opt_results.significance_n_lc    = significance_n_lc;
    opt_results.significance_n_lf    = significance_n_lf;
    opt_results.significance_metric  = significance_metric;
end


function [knee_thresh, knee_idx, knee_curvature, curvature] = local_inflection(thresholds, vol_frac)
% LOCAL_INFLECTION  Find the knee of vol_frac(threshold) via discrete d^2/dt^2.
%
% Returns NaNs everywhere when the curve has fewer than 3 finite samples
% (no second-derivative defined).  Curvature is computed on a 3-point
% moving-average smoothed copy to suppress noise from the small (13-pt)
% sweep.  The knee is the index that maximises -d^2V/dt^2 (most concave
% saturation point).  Endpoints of `curvature` are NaN.

    n = numel(vol_frac);
    curvature = nan(1, n);
    knee_thresh    = NaN;
    knee_idx       = NaN;
    knee_curvature = NaN;

    finite_mask = ~isnan(vol_frac);
    if sum(finite_mask) < 3
        return;
    end

    % 3-point moving-average smoothing (NaN-aware).  Boundary points
    % use a 2-point average so we don't lose them.
    v = vol_frac;
    sm = nan(1, n);
    for i = 1:n
        lo = max(1, i - 1);
        hi = min(n, i + 1);
        win = v(lo:hi);
        win = win(~isnan(win));
        if ~isempty(win)
            sm(i) = mean(win);
        end
    end

    % Discrete second derivative.  Uniform spacing is assumed (the sweep
    % is fixed step), so the divisor is constant and the argmax of the
    % unscaled finite difference is the same as that of the true d^2V/dt^2.
    for i = 2:n-1
        if ~isnan(sm(i-1)) && ~isnan(sm(i)) && ~isnan(sm(i+1))
            curvature(i) = sm(i+1) - 2*sm(i) + sm(i-1);
        end
    end

    if all(isnan(curvature))
        return;
    end
    % Pick the most concave (most negative curvature) interior point.
    [~, knee_idx] = min(curvature);
    knee_thresh    = thresholds(knee_idx);
    knee_curvature = curvature(knee_idx);
end


function [best_thresh, best_idx, best_p, pvalues, n_lc_at_best, n_lf_at_best] = ...
    local_significance(thresholds, per_patient_vol_frac, per_patient_lf)
% LOCAL_SIGNIFICANCE  Per-threshold Wilcoxon p-value of vol_frac LC vs LF.
%
% Returns NaNs when no threshold has both LC>=3 and LF>=3 finite values
% (Wilcoxon undefined / underpowered below this floor).

    n_thresh = numel(thresholds);
    pvalues  = nan(1, n_thresh);

    best_thresh   = NaN;
    best_idx      = NaN;
    best_p        = NaN;
    n_lc_at_best  = 0;
    n_lf_at_best  = 0;

    n_lc_arr = zeros(1, n_thresh);
    n_lf_arr = zeros(1, n_thresh);

    for ti = 1:n_thresh
        col = per_patient_vol_frac(:, ti);
        finite_mask = ~isnan(col) & ~isnan(per_patient_lf);
        lc_mask = finite_mask & per_patient_lf == 0;
        lf_mask = finite_mask & per_patient_lf == 1;
        n_lc = sum(lc_mask);
        n_lf = sum(lf_mask);
        n_lc_arr(ti) = n_lc;
        n_lf_arr(ti) = n_lf;
        if n_lc < 3 || n_lf < 3
            continue;
        end
        data = col(finite_mask);
        groups = per_patient_lf(finite_mask);
        try
            pvalues(ti) = perform_statistical_test(data, groups, 'ranksum');
        catch
            % Test failed — leave NaN at this threshold.
        end
    end

    valid = ~isnan(pvalues);
    if ~any(valid)
        return;
    end
    [best_p, best_idx] = min(pvalues);
    best_thresh   = thresholds(best_idx);
    n_lc_at_best  = n_lc_arr(best_idx);
    n_lf_at_best  = n_lf_arr(best_idx);
end


function v = safe_val(vec, idx)
    if idx < 1 || idx > numel(vec) || isnan(vec(idx))
        v = NaN;
    else
        v = vec(idx);
    end
end
