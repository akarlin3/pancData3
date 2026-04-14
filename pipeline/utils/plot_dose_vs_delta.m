function plot_dose_vs_delta(baseline_results, dosimetry_results, config_struct)
%PLOT_DOSE_VS_DELTA Scatter plots of GTV D95 vs percent change in D from baseline.
%
%   plot_dose_vs_delta(baseline_results, dosimetry_results, config_struct)
%
%   For each DWI type and each on-treatment fraction (Fx2, Fx3), generates a
%   1x2 figure:
%       Left  panel: Whole-GTV D95 vs Delta D (%)
%       Right panel: D-defined sub-volume D95 vs Delta D (%) (if available)
%
%   Outcomes are colour-coded:
%       LC = blue circles    [0 0.4470 0.7410]
%       LF = orange triangles [0.8500 0.3250 0.0980]
%       CR = grey squares    [0.5 0.5 0.5]
%   All markers are filled, black-edged, size 50.
%
%   Linear trend lines are fitted for LC and LF groups when n >= 3 per group
%   (solid for LC, dashed for LF). Spearman rho + p-value annotations are
%   computed per group; groups with n<3 are reported as "n<3".
%
%   PNG files are saved as
%       dose_vs_deltaD_Fx{fx}_{dtype_label}.png
%   at 300 DPI in config_struct.output_folder.
%
%   Inputs:
%       baseline_results - struct with fields D_pct, m_d95_gtvp, m_lf,
%                          m_id_list (fields may be missing; the function
%                          skips plots quietly when data is unavailable).
%       dosimetry_results- struct optionally containing d95_d_sub.
%       config_struct    - config struct with fields output_folder and
%                          dwi_types (cell array of DWI type labels).

    % --- Graceful early exit if required inputs are missing ---
    if ~isstruct(baseline_results) || ~isfield(baseline_results, 'D_pct') || ...
            ~isfield(baseline_results, 'm_d95_gtvp') || ~isfield(baseline_results, 'm_lf')
        fprintf('  \xe2\x9a\xa0\xef\xb8\x8f plot_dose_vs_delta: required baseline fields missing. Skipping.\n');
        return;
    end

    D_pct = baseline_results.D_pct;
    m_d95_gtvp = baseline_results.m_d95_gtvp;
    m_lf = baseline_results.m_lf(:);

    if isempty(D_pct) || isempty(m_d95_gtvp) || isempty(m_lf)
        fprintf('  \xe2\x9a\xa0\xef\xb8\x8f plot_dose_vs_delta: empty baseline data. Skipping.\n');
        return;
    end

    % Determine DWI type labels to iterate over
    if isfield(config_struct, 'dwi_types') && ~isempty(config_struct.dwi_types)
        dwi_types = config_struct.dwi_types;
        if ischar(dwi_types)
            dwi_types = {dwi_types};
        end
    elseif isfield(config_struct, 'dwi_type_name') && ~isempty(config_struct.dwi_type_name)
        dwi_types = {config_struct.dwi_type_name};
    else
        % Fall back to sensible defaults based on 3rd dimension size
        n_dtypes = size(D_pct, 3);
        default_names = {'Standard', 'dnCNN', 'IVIMnet'};
        dwi_types = default_names(1:min(n_dtypes, numel(default_names)));
    end

    % Optional sub-volume D95 matrix
    has_subvol = isstruct(dosimetry_results) && ...
        isfield(dosimetry_results, 'd95_d_sub') && ~isempty(dosimetry_results.d95_d_sub);
    if has_subvol
        d95_d_sub = dosimetry_results.d95_d_sub;
    else
        d95_d_sub = [];
    end

    output_folder = '';
    if isfield(config_struct, 'output_folder') && ~isempty(config_struct.output_folder)
        output_folder = config_struct.output_folder;
    end
    if isempty(output_folder) || ~exist(output_folder, 'dir')
        fprintf('  \xe2\x9a\xa0\xef\xb8\x8f plot_dose_vs_delta: output_folder missing. Skipping.\n');
        return;
    end

    n_dtypes = size(D_pct, 3);
    n_timepoints = size(D_pct, 2);

    fx_targets = [2, 3];

    for dt = 1:n_dtypes
        if dt > numel(dwi_types)
            dtype_label = sprintf('DWI%d', dt);
        else
            dtype_label = dwi_types{dt};
        end

        for fi = 1:numel(fx_targets)
            fx = fx_targets(fi);
            if fx > n_timepoints
                continue;
            end

            y_delta = D_pct(:, fx, dt);
            x_whole = m_d95_gtvp(:, fx);
            if all(isnan(y_delta)) || all(isnan(x_whole))
                % No usable data; skip quietly for this fx
                continue;
            end

            x_sub = [];
            if has_subvol && size(d95_d_sub, 2) >= fx
                x_sub = d95_d_sub(:, fx);
                if all(isnan(x_sub))
                    x_sub = [];
                end
            end

            fig = figure('Visible', 'off', 'Position', [150, 150, 1000, 400]);

            % --- Left panel: whole GTV D95 ---
            subplot(1, 2, 1);
            plot_single_panel(x_whole, y_delta, m_lf, ...
                'Whole GTV D95 (Gy)', '\Delta D (%)', 'Whole GTV D95 vs \Delta D');

            % --- Right panel: sub-volume D95 ---
            subplot(1, 2, 2);
            if ~isempty(x_sub)
                plot_single_panel(x_sub, y_delta, m_lf, ...
                    'Sub(D) D95 (Gy)', '\Delta D (%)', 'Sub(D) D95 vs \Delta D');
            else
                axis off;
                text(0.5, 0.5, 'Sub-volume D95 unavailable', ...
                    'HorizontalAlignment', 'center', 'FontSize', 11);
            end

            sgtitle(sprintf('Dose vs \\Delta D at Fx%d (%s)', fx, dtype_label), ...
                'FontSize', 13, 'FontWeight', 'bold');

            png_name = sprintf('dose_vs_deltaD_Fx%d_%s.png', fx, dtype_label);
            print(fig, fullfile(output_folder, png_name), '-dpng', '-r300');
            close(fig);
        end
    end

    fprintf('  \xe2\x9c\x85 plot_dose_vs_delta: scatter plots saved.\n');
end


function plot_single_panel(x_vals, y_vals, lf_group, x_label, y_label, panel_title)
%PLOT_SINGLE_PANEL Colour-coded scatter + per-group trends + Spearman annotation.

    x_vals = x_vals(:);
    y_vals = y_vals(:);
    lf_group_col = lf_group(:);

    % Align lengths defensively
    n_common = min([numel(x_vals), numel(y_vals), numel(lf_group_col)]);
    x_vals = x_vals(1:n_common);
    y_vals = y_vals(1:n_common);
    lf_group_col = lf_group_col(1:n_common);

    eligible_all = ~isnan(lf_group_col);
    clean_all = eligible_all & ~isnan(x_vals) & ~isnan(y_vals);
    clean = clean_all & (lf_group_col <= 1);

    hold on;

    % --- CR markers first (back layer) ---
    cr_mask = clean_all & lf_group_col == 2;
    if any(cr_mask)
        scatter(x_vals(cr_mask), y_vals(cr_mask), ...
            50, [0.5 0.5 0.5], 's', 'filled', ...
            'MarkerEdgeColor', 'k', 'DisplayName', 'CR');
    end

    % --- LC and LF markers ---
    lc_mask = clean & lf_group_col == 0;
    lf_mask = clean & lf_group_col == 1;

    if any(lc_mask)
        scatter(x_vals(lc_mask), y_vals(lc_mask), ...
            50, [0 0.4470 0.7410], 'o', 'filled', ...
            'MarkerEdgeColor', 'k', 'DisplayName', 'LC');
    end
    if any(lf_mask)
        scatter(x_vals(lf_mask), y_vals(lf_mask), ...
            50, [0.8500 0.3250 0.0980], '^', 'filled', ...
            'MarkerEdgeColor', 'k', 'DisplayName', 'LF');
    end

    % --- Trend lines (require n>=3 per group) ---
    warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    if sum(lc_mask) >= 3
        xl = linspace(min(x_vals(lc_mask)), max(x_vals(lc_mask)), 50);
        p_fit = polyfit(x_vals(lc_mask), y_vals(lc_mask), 1);
        plot(xl, polyval(p_fit, xl), '-', 'Color', [0 0.4470 0.7410], ...
            'LineWidth', 2, 'DisplayName', 'LC trend');
    end
    if sum(lf_mask) >= 3
        xl = linspace(min(x_vals(lf_mask)), max(x_vals(lf_mask)), 50);
        p_fit = polyfit(x_vals(lf_mask), y_vals(lf_mask), 1);
        plot(xl, polyval(p_fit, xl), '--', 'Color', [0.8500 0.3250 0.0980], ...
            'LineWidth', 2, 'DisplayName', 'LF trend');
    end
    warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    hold off;

    % --- Spearman per group ---
    r_lc = NaN; p_lc = NaN; r_lf = NaN; p_lf = NaN;
    if sum(lc_mask) >= 3
        [r_lc, p_lc] = safe_spearman(x_vals(lc_mask), y_vals(lc_mask));
    end
    if sum(lf_mask) >= 3
        [r_lf, p_lf] = safe_spearman(x_vals(lf_mask), y_vals(lf_mask));
    end

    n_lc = sum(lc_mask);
    n_lf = sum(lf_mask);
    n_cr = sum(clean_all & lf_group_col == 2);

    if isnan(r_lc)
        lc_str = sprintf('LC (n=%d): n<3', n_lc);
    else
        lc_str = sprintf('LC (n=%d) r_s=%.2f %s', n_lc, r_lc, format_p_value(p_lc));
    end
    if isnan(r_lf)
        lf_str = sprintf('LF (n=%d): n<3', n_lf);
    else
        lf_str = sprintf('LF (n=%d) r_s=%.2f %s', n_lf, r_lf, format_p_value(p_lf));
    end

    xlabel(x_label);
    ylabel(y_label);
    if n_cr > 0
        cap = sprintf('%s | %s | CR: n=%d', lc_str, lf_str, n_cr);
    else
        cap = sprintf('%s | %s', lc_str, lf_str);
    end
    if nargin < 6 || isempty(panel_title)
        title(cap, 'FontSize', 9);
    else
        title({panel_title, cap}, 'FontSize', 9);
    end
    grid on;
end


function [r, p] = safe_spearman(x, y)
%SAFE_SPEARMAN Octave-compatible Spearman rho with t-distribution p-value.
    r = NaN; p = NaN;
    if numel(x) < 3 || numel(y) < 3
        return;
    end
    try
        if exist('OCTAVE_VERSION', 'builtin')
            r = spearman(x, y);
            n = numel(x);
            t = r * sqrt((n - 2) / (1 - r^2 + eps));
            p = 2 * (1 - tcdf(abs(t), n - 2));
        else
            [r, p] = corr(x, y, 'Type', 'Spearman');
        end
    catch
        r = NaN; p = NaN;
    end
end
