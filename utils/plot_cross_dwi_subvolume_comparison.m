function plot_cross_dwi_subvolume_comparison(summary_metrics, config_struct)
% PLOT_CROSS_DWI_SUBVOLUME_COMPARISON — Cross-DWI-type ADC subvolume comparison at Fx1
%
%   Loads summary_metrics checkpoint files for all available DWI processing
%   types and generates a comparison visualization of the restricted ADC
%   subvolume (ADC < adc_thresh) at baseline (Fx1).
%
%   ANALYTICAL RATIONALE:
%   Different DWI processing methods (Standard, DnCNN denoising, IVIMnet)
%   produce different ADC maps from the same raw DWI data.  This visualization
%   addresses two methodological questions:
%
%   1. CROSS-DWI-TYPE EFFECT: Does the choice of processing pipeline
%      systematically alter the measured restricted diffusion subvolume?
%      DnCNN denoising may sharpen or shift ADC distributions, changing how
%      many voxels fall below the threshold.  IVIMnet reuses Standard ADC
%      (mono-exponential fit is pipeline-independent), so Standard and
%      IVIMnet should be identical for ADC-based subvolume.
%
%   2. SCAN-TO-SCAN REPEATABILITY: How reproducible is the ADC subvolume
%      across back-to-back repeat acquisitions at Fx1?  This quantifies
%      the measurement noise floor for subvolume as a biomarker — only
%      longitudinal changes exceeding this variability can be attributed
%      to treatment effects.
%
%   Inputs:
%       summary_metrics  - Summary metrics struct from current DWI type run
%       config_struct    - Configuration struct with dataloc, output_folder, etc.
%
%   Outputs:
%       None.  Saves figure as PNG.

    dwi_type_names = {'Standard', 'dnCNN', 'IVIMnet'};

    % --- 1. Load summary_metrics for all available DWI types ---
    all_sm = cell(1, 3);
    types_available = false(1, 3);

    for t = 1:3
        sm = load_type_checkpoint(dwi_type_names{t}, config_struct.dataloc);
        if ~isempty(sm) && isfield(sm, 'adc_sub_vol_pc')
            col = sm.adc_sub_vol_pc(:, 1, t);
            if any(~isnan(col))
                all_sm{t} = sm;
                types_available(t) = true;
            end
        end
    end

    n_available = sum(types_available);
    if n_available < 2
        fprintf('  💡 Cross-DWI ADC subvolume comparison skipped: need >= 2 types, found %d.\n', n_available);
        return;
    end

    % --- 2. Verify patient alignment across checkpoints ---
    avail_idx = find(types_available);
    ref_ids = all_sm{avail_idx(1)}.id_list;
    for ai = 2:numel(avail_idx)
        other_ids = all_sm{avail_idx(ai)}.id_list;
        if ~isequal(ref_ids, other_ids)
            fprintf('  ⚠️ Patient ID mismatch between DWI type checkpoints. Using intersection.\n');
            [ref_ids, ia, ~] = intersect(ref_ids, other_ids, 'stable');
            for aj = 1:numel(avail_idx)
                [~, idx] = ismember(ref_ids, all_sm{avail_idx(aj)}.id_list);
                all_sm{avail_idx(aj)}.adc_sub_vol_pc = all_sm{avail_idx(aj)}.adc_sub_vol_pc(idx, :, :);
                if isfield(all_sm{avail_idx(aj)}, 'adc_sub_vol_pc_rpt')
                    all_sm{avail_idx(aj)}.adc_sub_vol_pc_rpt = all_sm{avail_idx(aj)}.adc_sub_vol_pc_rpt(idx, :, :);
                end
            end
            break;
        end
    end
    nPat = numel(ref_ids);

    % --- 3. Extract Fx1 subvolume data ---
    % adc_sub_vol_pc is stored as a fraction (0-1); convert to percentage
    fx1_data = nan(nPat, 3);
    for t = avail_idx
        fx1_data(:, t) = all_sm{t}.adc_sub_vol_pc(:, 1, t) * 100;
    end

    % Detect Standard-IVIMnet identity
    std_ivimnet_identical = false;
    if types_available(1) && types_available(3)
        d = fx1_data(:, 1) - fx1_data(:, 3);
        if all(d == 0 | (isnan(fx1_data(:,1)) & isnan(fx1_data(:,3))))
            std_ivimnet_identical = true;
        end
    end

    % Check for repeat scan data
    has_rpt_data = false;
    rpt_data = [];
    rpt_type_idx = avail_idx(1);  % Use first available type for repeat panel
    if isfield(all_sm{rpt_type_idx}, 'adc_sub_vol_pc_rpt')
        rpt_raw = all_sm{rpt_type_idx}.adc_sub_vol_pc_rpt(:, :, rpt_type_idx);
        % Count patients with >= 2 repeat scans
        n_valid_rpts = sum(~isnan(rpt_raw), 2);
        if any(n_valid_rpts >= 2)
            has_rpt_data = true;
            rpt_data = rpt_raw * 100;  % Convert to percentage
        end
    end

    % --- 4. Determine output folder ---
    if isfield(config_struct, 'master_output_folder') && exist(config_struct.master_output_folder, 'dir')
        save_folder = config_struct.master_output_folder;
    elseif isfield(config_struct, 'output_folder') && exist(config_struct.output_folder, 'dir')
        save_folder = config_struct.output_folder;
    else
        save_folder = pwd;
    end

    % --- 5. Create figure ---
    if has_rpt_data
        n_panels = 3;
        fig_width = 1800;
    else
        n_panels = 2;
        fig_width = 1200;
    end

    fig = figure('Name', 'Cross-DWI ADC Subvolume Comparison at Fx1', ...
                 'Position', [50, 50, fig_width, 500]);

    % ----- Panel 1: Cross-DWI-type box plots -----
    subplot(1, n_panels, 1);
    plot_cross_dwi_boxes(fx1_data, avail_idx, dwi_type_names, std_ivimnet_identical, types_available);

    % ----- Panel 2 (if repeat data): Repeat-scan variability -----
    if has_rpt_data
        subplot(1, n_panels, 2);
        plot_repeat_variability(rpt_data, dwi_type_names{rpt_type_idx});
    end

    % ----- Last panel: Paired difference Standard vs dnCNN -----
    subplot(1, n_panels, n_panels);
    plot_paired_difference(fx1_data, types_available);

    sgtitle('ADC Restricted Subvolume (ADC < threshold) — Fx1 Cross-DWI Comparison', ...
            'FontSize', 13, 'FontWeight', 'bold');
    % Shift subplots down to avoid supertitle overlap
    allAx = findall(gcf, 'Type', 'Axes');
    for k = 1:numel(allAx)
        pos = get(allAx(k), 'Position');
        set(allAx(k), 'Position', [pos(1), pos(2) * 0.92, pos(3), pos(4) * 0.92]);
    end
    set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);

    out_file = fullfile(save_folder, 'Cross_DWI_ADC_Subvolume_Fx1.png');
    saveas(gcf, out_file);
    close(gcf);
    fprintf('  📁 Saved cross-DWI subvolume comparison: %s\n', out_file);
end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function sm = load_type_checkpoint(type_name, dataloc)
% Load summary_metrics checkpoint file for a specific DWI type.
    sm = [];
    f = fullfile(dataloc, sprintf('summary_metrics_%s.mat', type_name));
    if exist(f, 'file')
        tmp = load(f, 'summary_metrics');
        if isfield(tmp, 'summary_metrics')
            sm = tmp.summary_metrics;
        end
    end
end

function plot_cross_dwi_boxes(fx1_data, avail_idx, dwi_type_names, std_ivimnet_identical, types_available)
% Left panel: box plots of subvolume % across DWI types with patient lines.
    hold on;

    n_types = numel(avail_idx);
    % Build labels
    labels = cell(1, n_types);
    for i = 1:n_types
        t = avail_idx(i);
        labels{i} = dwi_type_names{t};
        if std_ivimnet_identical && t == 3
            labels{i} = [labels{i} '*'];
        end
    end

    % Prepare data matrix for boxplot (nPat x n_types)
    box_data = fx1_data(:, avail_idx);
    positions = 1:n_types;

    boxplot(box_data, 'Labels', labels, 'Positions', positions, ...
            'Widths', 0.5, 'Colors', [0.3 0.3 0.3]);

    % Overlay individual patient points with jitter
    colors = lines(n_types);
    jitter_width = 0.15;
    for i = 1:n_types
        vals = box_data(:, i);
        valid = ~isnan(vals);
        x_jitter = positions(i) + (rand(sum(valid), 1) - 0.5) * jitter_width * 2;
        scatter(x_jitter, vals(valid), 25, colors(i, :), 'filled', 'MarkerFaceAlpha', 0.5);
    end

    % Connect same-patient points
    for p = 1:size(box_data, 1)
        pt_vals = box_data(p, :);
        valid = ~isnan(pt_vals);
        if sum(valid) >= 2
            plot(positions(valid), pt_vals(valid), '-', 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
        end
    end

    % Signed-rank test between Standard and dnCNN if both available
    if types_available(1) && types_available(2)
        std_vals = fx1_data(:, 1);
        dncnn_vals = fx1_data(:, 2);
        paired_valid = ~isnan(std_vals) & ~isnan(dncnn_vals);
        n_paired = sum(paired_valid);
        if n_paired >= 5
            p_val = signrank(std_vals(paired_valid), dncnn_vals(paired_valid));
            p_str = format_p_value(p_val);
            % Find positions of Standard and dnCNN in avail_idx
            pos_std = find(avail_idx == 1);
            pos_dncnn = find(avail_idx == 2);
            y_max = max([std_vals; dncnn_vals], [], 'omitnan');
            y_top = y_max * 1.08;
            plot([positions(pos_std), positions(pos_dncnn)], [y_top, y_top], 'k-', 'LineWidth', 1);
            text(mean([positions(pos_std), positions(pos_dncnn)]), y_top * 1.02, ...
                 sprintf('p = %s (n=%d)', p_str, n_paired), ...
                 'HorizontalAlignment', 'center', 'FontSize', 8);
        end
    end

    ylabel('Restricted Subvolume (% of GTV)');
    title('Across DWI Types');
    if std_ivimnet_identical
        % Add footnote
        yl = ylim;
        text(0.02, 0.02, '*IVIMnet uses Standard ADC', ...
             'Units', 'normalized', 'FontSize', 7, 'FontAngle', 'italic', 'Color', [0.5 0.5 0.5]);
    end
    hold off;
    box on;
end

function plot_repeat_variability(rpt_data, type_label)
% Middle panel: subvolume across Fx1 repeat scans.
    hold on;

    nRpt = size(rpt_data, 2);
    % Determine max repeat with data
    max_rpt = 0;
    for r = 1:nRpt
        if any(~isnan(rpt_data(:, r)))
            max_rpt = r;
        end
    end
    if max_rpt < 2
        text(0.5, 0.5, 'Insufficient repeat data', ...
             'Units', 'normalized', 'HorizontalAlignment', 'center');
        title(sprintf('Repeat Scan Variability (%s)', type_label));
        hold off;
        return;
    end

    rpt_data = rpt_data(:, 1:max_rpt);
    positions = 1:max_rpt;

    % Build labels
    rpt_labels = arrayfun(@(x) sprintf('Scan %d', x), positions, 'UniformOutput', false);

    boxplot(rpt_data, 'Labels', rpt_labels, 'Positions', positions, ...
            'Widths', 0.5, 'Colors', [0.3 0.3 0.3]);

    % Overlay individual patient points with jitter
    colors = lines(max_rpt);
    jitter_width = 0.12;
    for r = 1:max_rpt
        vals = rpt_data(:, r);
        valid = ~isnan(vals);
        x_jitter = positions(r) + (rand(sum(valid), 1) - 0.5) * jitter_width * 2;
        scatter(x_jitter, vals(valid), 25, colors(r, :), 'filled', 'MarkerFaceAlpha', 0.5);
    end

    % Connect same-patient points across repeats
    for p = 1:size(rpt_data, 1)
        pt_vals = rpt_data(p, :);
        valid = ~isnan(pt_vals);
        if sum(valid) >= 2
            plot(positions(valid), pt_vals(valid), '-', 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
        end
    end

    % Compute and display wCV
    pts_with_rpts = sum(~isnan(rpt_data), 2) >= 2;
    if sum(pts_with_rpts) >= 3
        rpt_subset = rpt_data(pts_with_rpts, :);
        if exist('OCTAVE_VERSION', 'builtin')
            pt_means = mean(rpt_subset, 2, 'omitnan');
            pt_sds = std(rpt_subset, 0, 2, 'omitnan');
        else
            pt_means = mean(rpt_subset, 2, 'omitnan');
            pt_sds = std(rpt_subset, 0, 2, 'omitnan');
        end
        pt_means(pt_means == 0) = NaN;
        wcv_vals = pt_sds ./ pt_means;
        median_wcv = median(wcv_vals, 'omitnan') * 100;
        text(0.98, 0.98, sprintf('Median wCV = %.1f%%', median_wcv), ...
             'Units', 'normalized', 'HorizontalAlignment', 'right', ...
             'VerticalAlignment', 'top', 'FontSize', 9, 'FontWeight', 'bold');
    end

    ylabel('Restricted Subvolume (% of GTV)');
    title(sprintf('Fx1 Repeat Scans (%s)', type_label));
    hold off;
    box on;
end

function plot_paired_difference(fx1_data, types_available)
% Right panel: Bland-Altman-style paired difference (Standard - dnCNN).
    hold on;

    if ~types_available(1) || ~types_available(2)
        text(0.5, 0.5, 'Standard and/or dnCNN not available', ...
             'Units', 'normalized', 'HorizontalAlignment', 'center');
        title('Paired Difference (Standard - dnCNN)');
        hold off;
        return;
    end

    std_vals = fx1_data(:, 1);
    dncnn_vals = fx1_data(:, 2);
    paired_valid = ~isnan(std_vals) & ~isnan(dncnn_vals);
    n_paired = sum(paired_valid);

    if n_paired < 3
        text(0.5, 0.5, sprintf('Insufficient paired data (n=%d)', n_paired), ...
             'Units', 'normalized', 'HorizontalAlignment', 'center');
        title('Paired Difference (Standard - dnCNN)');
        hold off;
        return;
    end

    diff_vals = std_vals(paired_valid) - dncnn_vals(paired_valid);
    mean_vals = (std_vals(paired_valid) + dncnn_vals(paired_valid)) / 2;

    scatter(mean_vals, diff_vals, 40, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.6);

    % Mean difference line
    mean_diff = mean(diff_vals, 'omitnan');
    xl = xlim;
    plot(xl, [mean_diff mean_diff], 'r-', 'LineWidth', 1.5);
    plot(xl, [0 0], 'k--', 'LineWidth', 0.8);

    % Limits of agreement (mean +/- 1.96*SD)
    sd_diff = std(diff_vals, 'omitnan');
    loa_upper = mean_diff + 1.96 * sd_diff;
    loa_lower = mean_diff - 1.96 * sd_diff;
    plot(xl, [loa_upper loa_upper], 'r--', 'LineWidth', 0.8);
    plot(xl, [loa_lower loa_lower], 'r--', 'LineWidth', 0.8);

    % Signed-rank p-value
    p_str = '';
    if n_paired >= 5
        p_val = signrank(std_vals(paired_valid), dncnn_vals(paired_valid));
        p_str = format_p_value(p_val);
    end

    % Annotation text
    ann_lines = {sprintf('Mean diff = %.2f%%', mean_diff), ...
                 sprintf('n = %d', n_paired)};
    if ~isempty(p_str)
        ann_lines{end+1} = sprintf('p = %s', p_str);
    end
    text(0.02, 0.98, strjoin(ann_lines, '\n'), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'FontSize', 8, 'BackgroundColor', [1 1 1 0.8]);

    xlabel('Mean Subvolume (% GTV)');
    ylabel('Difference: Standard - dnCNN (% GTV)');
    title('Paired Difference');
    hold off;
    box on;
end
