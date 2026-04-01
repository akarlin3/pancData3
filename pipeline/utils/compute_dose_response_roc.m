function roc_results = compute_dose_response_roc(per_method_dosimetry, baseline_results, active_methods, config_struct)
%COMPUTE_DOSE_RESPONSE_ROC ROC analysis on sub-volume dose metrics to find optimal D95/V50 cutoff.
%
%   For each surviving core method, run ROC analysis on D95 (and V50) of the
%   sub-volume to find the optimal dose cutoff separating LC from LF.
%
%   Inputs:
%       per_method_dosimetry - struct from metrics_dosimetry (run_all_core_methods=true)
%       baseline_results     - struct from metrics_baseline
%       active_methods       - cell array of active method names
%       config_struct        - pipeline configuration struct
%
%   Outputs:
%       roc_results.method_results - struct array per method
%       roc_results.ranking        - methods sorted by best AUC descending

    % --- Diary ---
    output_folder = config_struct.output_folder;
    diary_file = fullfile(output_folder, ...
        sprintf('dose_response_roc_output_%s.txt', config_struct.dwi_type_name));
    diary(diary_file);

    fprintf('\n🔬 Dose-Response ROC Analysis (%s)\n', config_struct.dwi_type_name);

    n_methods = numel(active_methods);
    m_lf = baseline_results.m_lf;

    % Prepare outcome: event = LF (m_lf == 1), exclude competing risks (m_lf == 2)
    valid_idx = (m_lf == 0 | m_lf == 1);
    event = m_lf(valid_idx) == 1;

    if sum(event) < 2 || sum(~event) < 2
        fprintf('   ⚠️ Insufficient events for ROC analysis (need ≥2 in each group).\n');
        roc_results.method_results = struct([]);
        roc_results.ranking = {};
        diary off;
        return;
    end

    dose_metric_names = {'d95_adc_sub', 'v50_adc_sub', 'd95_d_sub', 'v50_d_sub'};
    dose_metric_labels = {'D95 ADC Sub', 'V50 ADC Sub', 'D95 D Sub', 'V50 D Sub'};

    method_results = struct([]);

    for mi = 1:n_methods
        method_name = active_methods{mi};
        fprintf('\n   Processing method: %s\n', method_name);

        if ~isfield(per_method_dosimetry, method_name)
            fprintf('      ⚠️ No dosimetry data for %s. Skipping.\n', method_name);
            continue;
        end
        method_dose = per_method_dosimetry.(method_name);

        result = struct();
        result.method_name = method_name;
        result.metrics = struct([]);
        best_auc = 0;
        best_metric = '';

        for di = 1:numel(dose_metric_names)
            dm_name = dose_metric_names{di};
            if ~isfield(method_dose, dm_name)
                continue;
            end

            dose_vals = method_dose.(dm_name);
            % Use Fx1 (column 1) if matrix
            if size(dose_vals, 2) > 1
                dose_vals = dose_vals(:, 1);
            end

            dose_valid = dose_vals(valid_idx);

            % Remove NaN entries
            nan_mask = isnan(dose_valid);
            dose_clean = dose_valid(~nan_mask);
            event_clean = event(~nan_mask);

            if numel(dose_clean) < 5
                continue;
            end

            % ROC via perfcurve
            try
                [fpr, tpr, thresholds, auc_val] = perfcurve(event_clean, dose_clean, true);
                % If AUC < 0.5, the predictor works in the opposite direction
                % (e.g., lower D95 predicts LF). Recompute with negated scores
                % so that Youden index and optimal threshold are correct.
                if auc_val < 0.5
                    [fpr, tpr, thresholds, auc_val] = perfcurve(event_clean, -dose_clean, true);
                    thresholds = -thresholds;
                end
            catch
                continue;
            end

            % Bootstrap AUC CI (1000 resamples)
            rng(42);
            n_boot = 1000;
            n_obs = numel(dose_clean);
            boot_aucs = nan(n_boot, 1);
            for bi = 1:n_boot
                idx = randi(n_obs, n_obs, 1);
                boot_event = event_clean(idx);
                boot_dose = dose_clean(idx);
                if numel(unique(boot_event)) < 2
                    continue;
                end
                try
                    [~, ~, ~, boot_auc] = perfcurve(boot_event, boot_dose, true);
                    boot_aucs(bi) = boot_auc;
                catch
                    % skip failed bootstrap
                end
            end
            boot_aucs = boot_aucs(~isnan(boot_aucs));
            if numel(boot_aucs) >= 100
                auc_ci = [prctile(boot_aucs, 2.5), prctile(boot_aucs, 97.5)];
            else
                auc_ci = [NaN, NaN];
            end

            % Optimal threshold via Youden index
            youden = tpr - fpr;
            [~, opt_idx] = max(youden);
            optimal_threshold = thresholds(opt_idx);
            opt_sens = tpr(opt_idx);
            opt_spec = 1 - fpr(opt_idx);

            metric_result = struct();
            metric_result.metric_name = dm_name;
            metric_result.auc = auc_val;
            metric_result.auc_ci = auc_ci;
            metric_result.optimal_threshold = optimal_threshold;
            metric_result.sensitivity = opt_sens;
            metric_result.specificity = opt_spec;
            metric_result.tpr = tpr;
            metric_result.fpr = fpr;

            if isempty(result.metrics)
                result.metrics = metric_result;
            else
                result.metrics(end+1) = metric_result;
            end

            if auc_val > best_auc
                best_auc = auc_val;
                best_metric = dm_name;
            end

            fprintf('      %s: AUC=%.3f [%.3f, %.3f], threshold=%.2f, sens=%.1f%%, spec=%.1f%%\n', ...
                dm_name, auc_val, auc_ci(1), auc_ci(2), optimal_threshold, ...
                opt_sens*100, opt_spec*100);
        end

        result.best_metric = best_metric;
        result.best_auc = best_auc;

        if isempty(method_results)
            method_results = result;
        else
            method_results(end+1) = result; %#ok<AGROW>
        end
    end

    % --- Ranking: sort by best AUC descending ---
    if ~isempty(method_results)
        best_aucs = [method_results.best_auc];
        [~, sort_idx] = sort(best_aucs, 'descend');
        ranking = {method_results(sort_idx).method_name};
    else
        ranking = {};
    end

    roc_results.method_results = method_results;
    roc_results.ranking = ranking;

    % --- Generate summary bar chart ---
    if ~isempty(method_results)
        try
            fig = figure('Visible', 'off', 'Position', [100 100 800 450]);
            method_labels = {method_results.method_name};
            aucs = [method_results.best_auc];
            bar(aucs);
            set(gca, 'XTick', 1:numel(method_labels), 'XTickLabel', method_labels, ...
                'FontSize', 7, 'XTickLabelRotation', 45);
            ylabel('Best AUC');
            title(sprintf('Dose-Response ROC: Best AUC per Method (%s)', config_struct.dwi_type_name));
            hold on;
            yline(0.7, '--g', 'AUC ≥ 0.7');
            hold off;
            ylim([0 1]);
            saveas(fig, fullfile(output_folder, ...
                sprintf('dose_roc_summary_%s.png', config_struct.dwi_type_name)));
            close(fig);
        catch ME_fig
            fprintf('⚠️ Could not generate ROC summary figure: %s\n', ME_fig.message);
        end
    end

    fprintf('\n✅ Dose-response ROC analysis complete.\n');
    diary off;
end
