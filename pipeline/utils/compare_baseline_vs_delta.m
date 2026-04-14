function comparison = compare_baseline_vs_delta(baseline_results, config_struct)
%COMPARE_BASELINE_VS_DELTA Univariable Cox PH: baseline vs delta diffusion predictors.
%
%   comparison = compare_baseline_vs_delta(baseline_results, config_struct)
%
%   For each diffusion parameter (ADC, D, f, D*) at each timepoint (Fx2, Fx3),
%   fits two univariable Cox PH models:
%       Model 1: baseline value (timepoint 1)
%       Model 2: percent / absolute change from baseline to the timepoint
%   Harrell's C-index is computed on the training data for each model and
%   the better predictor (higher C-index) is identified per row.
%
%   A grouped bar chart is produced and saved as
%       baseline_vs_delta_{dwi_label}.png
%   in config_struct.output_folder.
%
%   Inputs:
%       baseline_results - struct with fields ADC_abs, D_abs, f_abs, Dstar_abs
%                          (baseline absolute values) and ADC_pct, D_pct,
%                          f_delta, Dstar_pct (percent/absolute changes);
%                          plus m_lf (0=LC,1=LF,2=CR) and a survival time
%                          field such as td_time, m_total_time, or
%                          m_total_follow_up_time.
%       config_struct    - config struct with output_folder and
%                          dwi_types_to_run (scalar DWI index).
%
%   Output:
%       comparison.parameters  = {'ADC','D','f','Dstar'}
%       comparison.timepoints  = [2, 3]
%       comparison.results     = struct array with per-row fields
%                                .parameter, .timepoint,
%                                .baseline_hr, .baseline_p, .baseline_cindex,
%                                .delta_hr, .delta_p, .delta_cindex,
%                                .better_predictor  ('baseline' or 'delta'),
%                                .n_events, .n_patients
%
%   The function returns a populated struct even if data is insufficient; rows
%   where Cox fitting is infeasible are filled with NaN values.

    comparison = struct();
    comparison.parameters = {'ADC', 'D', 'f', 'Dstar'};
    comparison.timepoints = [2, 3];
    comparison.results = [];

    if ~isstruct(baseline_results) || ~isfield(baseline_results, 'D_pct')
        fprintf('  \xe2\x9a\xa0\xef\xb8\x8f compare_baseline_vs_delta: required field D_pct missing. Skipping.\n');
        return;
    end

    required_fields = {'m_lf', 'ADC_abs', 'D_abs', 'f_abs', 'Dstar_abs', ...
        'ADC_pct', 'D_pct', 'f_delta', 'Dstar_pct'};
    for i = 1:numel(required_fields)
        if ~isfield(baseline_results, required_fields{i})
            fprintf('  \xe2\x9a\xa0\xef\xb8\x8f compare_baseline_vs_delta: required field %s missing. Skipping.\n', ...
                required_fields{i});
            return;
        end
    end

    % DWI type index
    if isfield(config_struct, 'dwi_types_to_run') && ~isempty(config_struct.dwi_types_to_run)
        dtype = config_struct.dwi_types_to_run(1);
    else
        dtype = 1;
    end

    % DWI label
    dwi_names = {'Standard', 'dnCNN', 'IVIMnet'};
    if isfield(config_struct, 'dwi_type_name') && ~isempty(config_struct.dwi_type_name)
        dwi_label = config_struct.dwi_type_name;
    elseif dtype >= 1 && dtype <= 3
        dwi_label = dwi_names{dtype};
    else
        dwi_label = 'DWI';
    end

    % Resolve survival time vector
    if isfield(baseline_results, 'td_time') && ~isempty(baseline_results.td_time)
        T_all = baseline_results.td_time(:);
    elseif isfield(baseline_results, 'm_total_time') && ~isempty(baseline_results.m_total_time)
        T_all = baseline_results.m_total_time(:);
    elseif isfield(baseline_results, 'm_total_follow_up_time') && ...
            ~isempty(baseline_results.m_total_follow_up_time)
        T_all = baseline_results.m_total_follow_up_time(:);
    else
        fprintf('  \xe2\x9a\xa0\xef\xb8\x8f compare_baseline_vs_delta: no follow-up time field. Skipping.\n');
        return;
    end

    % Event: treat CR (==2) as censored
    m_lf = baseline_results.m_lf(:);
    event = m_lf;
    event(event == 2) = 0;
    event_binary = event == 1;  % true = LF event, false = censored

    baseline_fields = {'ADC_abs', 'D_abs', 'f_abs', 'Dstar_abs'};
    delta_fields    = {'ADC_pct', 'D_pct', 'f_delta', 'Dstar_pct'};
    param_labels    = {'ADC', 'D', 'f', 'Dstar'};
    timepoints = [2, 3];

    ws = warning('off', 'stats:coxphfit:FitWarning');
    cleanupWarn = onCleanup(@() warning(ws));

    results = struct('parameter', {}, 'timepoint', {}, ...
        'baseline_hr', {}, 'baseline_p', {}, 'baseline_cindex', {}, ...
        'delta_hr', {}, 'delta_p', {}, 'delta_cindex', {}, ...
        'better_predictor', {}, 'n_events', {}, 'n_patients', {});

    row = 0;
    for pi = 1:numel(param_labels)
        X_base_full = baseline_results.(baseline_fields{pi});
        X_delta_full = baseline_results.(delta_fields{pi});

        if isempty(X_base_full) || isempty(X_delta_full)
            continue;
        end

        % Ensure 3D: [pts x Tp x dtype]
        if ndims(X_base_full) < 3
            X_base_fx1 = X_base_full(:, 1);
        else
            X_base_fx1 = X_base_full(:, 1, min(dtype, size(X_base_full, 3)));
        end

        for ti = 1:numel(timepoints)
            fx = timepoints(ti);
            row = row + 1;

            if ndims(X_delta_full) < 3 || size(X_delta_full, 2) < fx
                if size(X_delta_full, 2) < fx
                    X_delta = nan(size(X_delta_full, 1), 1);
                else
                    X_delta = X_delta_full(:, fx);
                end
            else
                X_delta = X_delta_full(:, fx, min(dtype, size(X_delta_full, 3)));
            end

            % --- Baseline model ---
            [b_hr, b_p, b_c, b_n, b_ne] = fit_cox_and_cindex(X_base_fx1, T_all, event_binary);

            % --- Delta model ---
            [d_hr, d_p, d_c, d_n, d_ne] = fit_cox_and_cindex(X_delta, T_all, event_binary);

            if ~isnan(b_c) && ~isnan(d_c)
                if d_c > b_c
                    better = 'delta';
                else
                    better = 'baseline';
                end
            elseif ~isnan(d_c)
                better = 'delta';
            elseif ~isnan(b_c)
                better = 'baseline';
            else
                better = 'baseline';  % default when both undefined
            end

            n_patients_used = max(b_n, d_n);
            n_events_used = max(b_ne, d_ne);

            results(row).parameter = param_labels{pi};
            results(row).timepoint = fx;
            results(row).baseline_hr = b_hr;
            results(row).baseline_p = b_p;
            results(row).baseline_cindex = b_c;
            results(row).delta_hr = d_hr;
            results(row).delta_p = d_p;
            results(row).delta_cindex = d_c;
            results(row).better_predictor = better;
            results(row).n_events = n_events_used;
            results(row).n_patients = n_patients_used;
        end
    end

    comparison.results = results;

    % --- Figure: grouped bar chart of C-index ---
    try
        n_rows = numel(results);
        if n_rows > 0
            cindex_mat = nan(n_rows, 2);
            p_mat = nan(n_rows, 2);
            row_labels = cell(n_rows, 1);
            for r = 1:n_rows
                cindex_mat(r, 1) = results(r).baseline_cindex;
                cindex_mat(r, 2) = results(r).delta_cindex;
                p_mat(r, 1) = results(r).baseline_p;
                p_mat(r, 2) = results(r).delta_p;
                row_labels{r} = sprintf('%s Fx%d', results(r).parameter, results(r).timepoint);
            end

            fig = figure('Visible', 'off', 'Position', [100, 100, 1000, 500]);
            h = bar(cindex_mat, 'grouped');
            if numel(h) >= 2
                set(h(1), 'FaceColor', [0.2 0.4 0.7]);
                set(h(2), 'FaceColor', [0.85 0.4 0.2]);
            end

            set(gca, 'XTick', 1:n_rows, 'XTickLabel', row_labels, ...
                'XTickLabelRotation', 30);
            ylabel('Harrell''s C-Index');
            ylim([0, 1]);
            hold on;
            yline_val = 0.5;
            plot(xlim, [yline_val, yline_val], 'k--', 'LineWidth', 1.0);
            legend({'Baseline', 'Delta'}, 'Location', 'best');
            title(sprintf('Baseline vs Delta Cox PH C-Index (%s)', dwi_label));
            grid on;

            % Annotate with p-values above bars
            for r = 1:n_rows
                if ~isnan(cindex_mat(r, 1))
                    text(r - 0.15, cindex_mat(r, 1) + 0.03, ...
                        format_p_compact(p_mat(r, 1)), 'FontSize', 7, ...
                        'HorizontalAlignment', 'center');
                end
                if ~isnan(cindex_mat(r, 2))
                    text(r + 0.15, cindex_mat(r, 2) + 0.03, ...
                        format_p_compact(p_mat(r, 2)), 'FontSize', 7, ...
                        'HorizontalAlignment', 'center');
                end
            end
            hold off;

            if isfield(config_struct, 'output_folder') && ~isempty(config_struct.output_folder) && ...
                    exist(config_struct.output_folder, 'dir')
                png_path = fullfile(config_struct.output_folder, ...
                    sprintf('baseline_vs_delta_%s.png', dwi_label));
                print(fig, png_path, '-dpng', '-r300');
            end
            close(fig);
        end
    catch ME_fig
        fprintf('  \xe2\x9a\xa0\xef\xb8\x8f Figure generation failed: %s\n', ME_fig.message);
    end
end


function [hr, p, cindex, n_used, n_events_used] = fit_cox_and_cindex(X, T, event_binary)
%FIT_COX_AND_CINDEX Univariable Cox PH fit + C-index on training set.
    hr = NaN; p = NaN; cindex = NaN; n_used = 0; n_events_used = 0;

    X = X(:); T = T(:); event_binary = logical(event_binary(:));
    n_common = min([numel(X), numel(T), numel(event_binary)]);
    X = X(1:n_common);
    T = T(1:n_common);
    event_binary = event_binary(1:n_common);

    good = ~isnan(X) & ~isnan(T) & T > 0;
    X = X(good); T = T(good); event_binary = event_binary(good);

    if numel(X) < 5 || sum(event_binary) < 2 || std(X) == 0
        return;
    end

    n_used = numel(X);
    n_events_used = sum(event_binary);

    try
        censoring = ~event_binary;
        [b, ~, ~, stats] = coxphfit(X, T, 'Censoring', censoring);
        hr = exp(b);
        if isfield(stats, 'p')
            p = stats.p;
        end
        % Linear predictor = beta * X (for univariable, sign is the same as beta)
        risk_scores = b * X;
        cindex = harrell_cindex(risk_scores, T, event_binary);
    catch
        % Leave outputs as NaN
    end
end


function c = harrell_cindex(risk, T, event)
%HARRELL_CINDEX  Simple concordance index for univariable risk scores.
%   Higher risk should indicate shorter survival.
    risk = risk(:); T = T(:); event = logical(event(:));
    n = numel(risk);
    num = 0;
    den = 0;
    for i = 1:n
        if ~event(i)
            continue;
        end
        for j = 1:n
            if i == j
                continue;
            end
            if T(j) > T(i)
                den = den + 1;
                if risk(i) > risk(j)
                    num = num + 1;
                elseif risk(i) == risk(j)
                    num = num + 0.5;
                end
            end
        end
    end
    if den == 0
        c = NaN;
    else
        c = num / den;
    end
end


function s = format_p_compact(p)
    if isnan(p)
        s = 'p=NA';
    elseif p < 0.001
        s = 'p<.001';
    else
        s = sprintf('p=%.2f', p);
    end
end
