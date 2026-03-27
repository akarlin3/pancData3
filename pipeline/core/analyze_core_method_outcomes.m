function results = analyze_core_method_outcomes(per_method_dosimetry, baseline_results, active_methods, config_struct)
%ANALYZE_CORE_METHOD_OUTCOMES Test dose coverage vs local control per core method.
%
%   For each method in active_methods, extracts D95 and V50 dose metrics from
%   per_method_dosimetry, pairs with outcome data from baseline_results, and
%   runs univariable Cox PH + Kaplan-Meier analysis.
%
%   This is a NON-FATAL pipeline step: errors are caught by the orchestrator.
%
% Inputs:
%   per_method_dosimetry - struct from metrics_dosimetry (run_all_core_methods=true)
%   baseline_results     - struct from metrics_baseline
%   active_methods       - cell array of surviving method names
%   config_struct        - pipeline config
%
% Output:
%   results - struct with .method_results, .ranking, .summary_table, .active_methods

    % --- Diary ---
    output_folder = config_struct.output_folder;
    diary_file = fullfile(output_folder, ...
        sprintf('core_method_outcomes_output_%s.txt', config_struct.dwi_type_name));
    diary(diary_file);

    fprintf('\n🔬 Core Method Outcome Analysis (%s)\n', config_struct.dwi_type_name);
    fprintf('   Testing dose coverage vs local control for %d methods.\n\n', numel(active_methods));

    % --- Initialize output ---
    results = struct();
    results.active_methods = active_methods;
    results.method_results = struct([]);
    results.ranking = {};
    results.summary_table = {};

    % --- Handle empty dosimetry ---
    if isempty(per_method_dosimetry) || isempty(fieldnames(per_method_dosimetry))
        fprintf('⚠️ per_method_dosimetry is empty. Set run_all_core_methods=true.\n');
        diary off;
        return;
    end

    % --- Prepare outcome data ---
    valid = baseline_results.valid_pts;
    m_lf = baseline_results.m_lf;
    m_time = baseline_results.m_total_time;
    m_follow = baseline_results.m_total_follow_up_time;

    td_lf = m_lf(valid);
    td_time = m_time(valid);
    follow_valid = m_follow(valid);

    % Use follow-up time for censored/competing patients
    cens_mask = (td_lf == 0 | td_lf == 2) & ~isnan(follow_valid);
    td_time(cens_mask) = follow_valid(cens_mask);

    % Cause-specific hazard: competing risk -> censored
    event_csh = td_lf;
    event_csh(event_csh == 2) = 0;

    n_events = sum(event_csh == 1);
    n_patients = sum(valid);

    fprintf('   Patients: %d, LF events: %d\n', n_patients, n_events);

    % Check feasibility
    if n_events < 3
        fprintf('⚠️ Insufficient events (%d < 3). Skipping outcome analysis.\n', n_events);
        diary off;
        return;
    end

    % --- Select timepoint ---
    tp = 1;  % Default: Fx1
    if isfield(config_struct, 'best_risk_fx')
        tp = config_struct.best_risk_fx;
    end

    % --- Dose metric names ---
    METRIC_NAMES = {'d95_adc_sub', 'v50_adc_sub', 'd95_d_sub', 'v50_d_sub'};
    n_metrics = numel(METRIC_NAMES);

    % --- Loop over active methods ---
    m_idx = 0;
    best_p_values = [];
    best_metric_names = {};
    method_names_found = {};

    ws = warning('off', 'stats:coxphfit:FitWarning');
    cleanupWarn = onCleanup(@() warning(ws));

    for i = 1:numel(active_methods)
        mname = active_methods{i};

        % Check if method has dosimetry data
        if ~isfield(per_method_dosimetry, mname)
            fprintf('   ⏭️ %s: no dosimetry data, skipping.\n', mname);
            continue;
        end

        method_dos = per_method_dosimetry.(mname);
        m_idx = m_idx + 1;

        % --- Univariable Cox for each dose metric ---
        univar = struct([]);
        for mi = 1:n_metrics
            met_name = METRIC_NAMES{mi};
            if ~isfield(method_dos, met_name)
                continue;
            end

            dose_vals_all = method_dos.(met_name);
            if size(dose_vals_all, 2) < tp
                continue;
            end
            dose_vals = dose_vals_all(valid, tp);

            % Filter to complete cases
            finite_mask = ~isnan(dose_vals) & ~isnan(td_time) & ~isnan(event_csh);
            n_finite = sum(finite_mask);
            if n_finite < 3 || sum(event_csh(finite_mask) == 1) < 1
                continue;
            end

            x = dose_vals(finite_mask);
            t = td_time(finite_mask);
            e = event_csh(finite_mask);

            % Fit univariable Cox PH
            try
                [b, ~, ~, stats] = coxphfit(x, t, 'Censoring', e == 0);
                hr = exp(b);
                se = stats.se;
                hr_ci = exp([b - 1.96*se, b + 1.96*se]);
                p_val = stats.p;
            catch
                continue;
            end

            entry = struct();
            entry.metric_name = met_name;
            entry.hr = hr;
            entry.hr_ci = hr_ci;
            entry.p_value = p_val;
            entry.n = n_finite;
            entry.n_events = sum(e == 1);

            if isempty(univar)
                univar = entry;
            else
                univar(end+1) = entry; %#ok<AGROW>
            end
        end

        % --- Find best metric (lowest p) ---
        km_result = struct('best_metric', '', 'logrank_p', NaN, ...
            'logrank_chi2', NaN, 'median_high', NaN, 'median_low', NaN);

        if ~isempty(univar)
            p_vals = [univar.p_value];
            [~, best_idx] = min(p_vals);
            best_met = univar(best_idx).metric_name;

            % --- KM + log-rank for best metric ---
            dose_vals = method_dos.(best_met)(valid, tp);
            finite_mask = ~isnan(dose_vals) & ~isnan(td_time) & ~isnan(event_csh);
            x = dose_vals(finite_mask);
            t = td_time(finite_mask);
            e = event_csh(finite_mask);

            med_val = nanmedian(x);
            group = x >= med_val;  % high = true

            % Log-rank test
            try
                [lr_p, lr_chi2] = simple_logrank(t, e, group);
            catch
                lr_p = NaN;
                lr_chi2 = NaN;
            end

            % Median survival per group
            med_high = compute_median_survival(t(group), e(group));
            med_low = compute_median_survival(t(~group), e(~group));

            km_result.best_metric = best_met;
            km_result.logrank_p = lr_p;
            km_result.logrank_chi2 = lr_chi2;
            km_result.median_high = med_high;
            km_result.median_low = med_low;

            % --- KM plot ---
            try
                fig = figure('Visible', 'off', 'Position', [100 100 600 450]);
                [f_h, x_h] = ecdf(t(group), 'function', 'survivor', 'Censoring', e(group) == 0);
                [f_l, x_l] = ecdf(t(~group), 'function', 'survivor', 'Censoring', e(~group) == 0);
                stairs(x_h, f_h, 'b-', 'LineWidth', 2); hold on;
                stairs(x_l, f_l, 'r-', 'LineWidth', 2);
                xlabel('Time (days)'); ylabel('Survival Probability');
                title(sprintf('KM: %s (%s) — %s', mname, best_met, config_struct.dwi_type_name));
                legend({sprintf('High coverage (n=%d)', sum(group)), ...
                        sprintf('Low coverage (n=%d)', sum(~group))}, ...
                    'Location', 'southwest');
                text(0.95, 0.05, sprintf('Log-rank p = %.3f', lr_p), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 9);
                saveas(fig, fullfile(output_folder, ...
                    sprintf('km_%s_%s.png', mname, config_struct.dwi_type_name)));
                close(fig);
            catch
                % Non-fatal: skip figure on error
            end

            best_p_values(end+1) = min(p_vals); %#ok<AGROW>
        else
            best_p_values(end+1) = NaN; %#ok<AGROW>
        end

        % Package method result
        mr = struct();
        mr.method_name = mname;
        mr.univariable = univar;
        mr.km = km_result;
        mr.n_patients = n_patients;
        mr.n_events = n_events;

        if isempty(results.method_results)
            results.method_results = mr;
        else
            results.method_results(end+1) = mr;
        end

        method_names_found{end+1} = mname; %#ok<AGROW>
        best_metric_names{end+1} = km_result.best_metric; %#ok<AGROW>
    end

    % --- Rank methods by best p-value ---
    if ~isempty(best_p_values)
        [sorted_p, sort_idx] = sort(best_p_values, 'ascend');
        results.ranking = method_names_found(sort_idx);

        % Summary table
        summary = cell(numel(sort_idx), 5);
        for si = 1:numel(sort_idx)
            idx = sort_idx(si);
            mr = results.method_results(idx);
            summary{si, 1} = mr.method_name;
            summary{si, 2} = best_metric_names{idx};
            if ~isempty(mr.univariable)
                p_vals = [mr.univariable.p_value];
                [~, bi] = min(p_vals);
                summary{si, 3} = mr.univariable(bi).hr;
                summary{si, 4} = mr.univariable(bi).p_value;
                summary{si, 5} = mr.univariable(bi).n_events;
            else
                summary{si, 3} = NaN;
                summary{si, 4} = NaN;
                summary{si, 5} = 0;
            end
        end
        results.summary_table = summary;
    end

    % --- Print summary ---
    fprintf('\nCore Method Outcome Analysis (%s):\n', config_struct.dwi_type_name);
    fprintf('%-22s %-16s %-22s %-10s %s\n', 'Method', 'Best Metric', 'HR (95% CI)', 'p-value', 'Events');
    fprintf('%s\n', repmat('-', 1, 80));
    for si = 1:numel(results.ranking)
        idx = sort_idx(si);
        mr = results.method_results(idx);
        if ~isempty(mr.univariable)
            p_vals = [mr.univariable.p_value];
            [~, bi] = min(p_vals);
            uv = mr.univariable(bi);
            fprintf('%-22s %-16s %5.2f (%4.2f-%4.2f)     %-10.4f %d\n', ...
                mr.method_name, uv.metric_name, uv.hr, uv.hr_ci(1), uv.hr_ci(2), ...
                uv.p_value, uv.n_events);
        else
            fprintf('%-22s %-16s %-22s %-10s %s\n', mr.method_name, '—', '—', '—', '—');
        end
    end

    if ~isempty(results.ranking)
        fprintf('\nBest method: %s\n', results.ranking{1});
    end

    % --- Summary bar chart ---
    try
        valid_p = best_p_values(sort_idx);
        valid_names = method_names_found(sort_idx);
        finite = ~isnan(valid_p) & valid_p > 0;
        if any(finite)
            fig = figure('Visible', 'off', 'Position', [100 100 700 400]);
            barh(-log10(valid_p(finite)));
            set(gca, 'YTick', 1:sum(finite), 'YTickLabel', valid_names(finite), 'FontSize', 8);
            xlabel('-log_{10}(p-value)');
            title(sprintf('Core Method Outcomes (%s)', config_struct.dwi_type_name));
            hold on;
            xline(-log10(0.05), 'r--', 'p=0.05', 'LineWidth', 1.5);
            saveas(fig, fullfile(output_folder, ...
                sprintf('core_method_outcomes_summary_%s.png', config_struct.dwi_type_name)));
            close(fig);
        end
    catch
        % Non-fatal
    end

    fprintf('\n✅ Core method outcome analysis complete.\n');
    diary off;
end


%% ===== Local helper functions =====

function [p, chi2] = simple_logrank(time, event, group)
%SIMPLE_LOGRANK Two-sample log-rank (Mantel-Haenszel) test.
%   group is logical: true = group 1, false = group 0.

    event_times = unique(time(event == 1));
    event_times = sort(event_times);

    O1 = 0; E1 = 0;  % observed / expected in group 1

    for ti = 1:numel(event_times)
        t_i = event_times(ti);

        % At risk at time t_i (patients with time >= t_i)
        at_risk = time >= t_i;
        n = sum(at_risk);
        n1 = sum(at_risk & group);

        % Events at exactly t_i
        events_at_t = (time == t_i) & (event == 1);
        d = sum(events_at_t);
        d1 = sum(events_at_t & group);

        if n == 0, continue; end

        O1 = O1 + d1;
        E1 = E1 + d * n1 / n;
    end

    % Variance under null
    V = 0;
    for ti = 1:numel(event_times)
        t_i = event_times(ti);
        at_risk = time >= t_i;
        n = sum(at_risk);
        n1 = sum(at_risk & group);
        n0 = n - n1;
        events_at_t = (time == t_i) & (event == 1);
        d = sum(events_at_t);

        if n <= 1, continue; end
        V = V + d * n1 * n0 * (n - d) / (n^2 * (n - 1));
    end

    if V > 0
        chi2 = (O1 - E1)^2 / V;
        p = 1 - chi2cdf(chi2, 1);
    else
        chi2 = 0;
        p = 1;
    end
end


function med = compute_median_survival(time, event)
%COMPUTE_MEDIAN_SURVIVAL Compute median survival from KM curve.
    if isempty(time)
        med = NaN;
        return;
    end
    try
        [f, x] = ecdf(time, 'function', 'survivor', 'Censoring', event == 0);
        idx = find(f <= 0.5, 1, 'first');
        if ~isempty(idx)
            med = x(idx);
        else
            med = NaN;  % Median not reached
        end
    catch
        med = NaN;
    end
end
