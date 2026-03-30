function concordance = compute_risk_dose_concordance(predictive_results, per_method_dosimetry, baseline_results, active_methods, config_struct)
%COMPUTE_RISK_DOSE_CONCORDANCE Compare elastic net risk model vs dose coverage stratification.
%
%   Inputs:
%       predictive_results  - struct with .risk_scores_all, .is_high_risk
%       per_method_dosimetry - from metrics_dosimetry (run_all_core_methods=true)
%       baseline_results     - from metrics_baseline
%       active_methods       - cell array of surviving method names
%       config_struct        - pipeline configuration struct
%
%   Outputs:
%       concordance.method_results - struct array per method
%       concordance.summary        - overall interpretation text

    % --- Diary ---
    output_folder = config_struct.output_folder;
    diary_file = fullfile(output_folder, ...
        sprintf('risk_dose_concordance_output_%s.txt', config_struct.dwi_type_name));
    diary(diary_file);

    fprintf('\n🔬 Risk-Dose Concordance Analysis (%s)\n', config_struct.dwi_type_name);

    % Check for predictive results
    if ~isfield(predictive_results, 'is_high_risk') || all(isnan(predictive_results.is_high_risk))
        fprintf('   ⚠️ No predictive model results available. Skipping.\n');
        concordance.method_results = struct([]);
        concordance.summary = 'No predictive model results available.';
        diary off;
        return;
    end

    is_high_risk = logical(predictive_results.is_high_risk);
    risk_scores = predictive_results.risk_scores_all;
    m_lf = baseline_results.m_lf;
    n_methods = numel(active_methods);

    % Exclude competing risks for combined AUC
    valid_idx = (m_lf == 0 | m_lf == 1);
    event = m_lf(valid_idx) == 1;

    method_results = struct([]);

    for mi = 1:n_methods
        method_name = active_methods{mi};
        fprintf('\n   Processing method: %s\n', method_name);

        result = struct();
        result.method_name = method_name;
        result.best_dose_metric = '';
        result.cohen_kappa = NaN;
        result.concordance_pct = NaN;
        result.confusion_matrix = nan(2, 2);
        result.n_complementary = NaN;
        result.combined_auc = NaN;

        if ~isfield(per_method_dosimetry, method_name)
            if isempty(method_results)
                method_results = result;
            else
                method_results(end+1) = result; %#ok<AGROW>
            end
            continue;
        end
        method_dose = per_method_dosimetry.(method_name);

        % Use d95_adc_sub at Fx1 as default dose metric
        d95_field = 'd95_adc_sub';
        if ~isfield(method_dose, d95_field)
            if isempty(method_results)
                method_results = result;
            else
                method_results(end+1) = result; %#ok<AGROW>
            end
            continue;
        end
        dose_vals = method_dose.(d95_field);
        if size(dose_vals, 2) > 1
            dose_vals = dose_vals(:, 1);
        end
        result.best_dose_metric = d95_field;

        % Median-split dose into high/low coverage
        valid_dose = ~isnan(dose_vals) & ~isnan(is_high_risk);
        if sum(valid_dose) < 6
            if isempty(method_results)
                method_results = result;
            else
                method_results(end+1) = result; %#ok<AGROW>
            end
            continue;
        end

        dose_median = nanmedian(dose_vals(valid_dose));
        is_low_coverage = dose_vals < dose_median;

        % Build 2x2 confusion matrix: is_high_risk x is_low_coverage
        % [both_low_risk_high_dose, risk_only; dose_only, both_high_risk_low_dose]
        both_low = sum(valid_dose & ~is_high_risk & ~is_low_coverage);
        risk_only = sum(valid_dose & is_high_risk & ~is_low_coverage);
        dose_only = sum(valid_dose & ~is_high_risk & is_low_coverage);
        both_high = sum(valid_dose & is_high_risk & is_low_coverage);

        cm = [both_low, risk_only; dose_only, both_high];
        result.confusion_matrix = cm;
        n_total = sum(cm(:));

        % Cohen's kappa
        p_o = (both_low + both_high) / n_total;  % observed agreement
        p_risk = (risk_only + both_high) / n_total;
        p_dose = (dose_only + both_high) / n_total;
        p_no_risk = (both_low + dose_only) / n_total;
        p_no_dose = (both_low + risk_only) / n_total;
        p_e = p_risk * p_dose + p_no_risk * p_no_dose;  % expected by chance

        if p_e < 1
            result.cohen_kappa = (p_o - p_e) / (1 - p_e);
        end
        result.concordance_pct = p_o * 100;
        result.n_complementary = risk_only + dose_only;

        fprintf('      Concordance: %.1f%%, Kappa=%.3f, Complementary=%d patients\n', ...
            result.concordance_pct, result.cohen_kappa, result.n_complementary);

        % Combined AUC: logistic regression with risk_score + dose as predictors
        valid_combined = valid_idx & valid_dose;
        risk_clean = risk_scores(valid_combined);
        dose_clean = dose_vals(valid_combined);
        event_combined = m_lf(valid_combined) == 1;

        if sum(valid_combined) >= 10 && sum(event_combined) >= 3 && sum(~event_combined) >= 3
            try
                X_comb = [risk_clean(:), dose_clean(:)];
                % Simple combined score: standardize and sum
                z_risk = (risk_clean(:) - nanmean_safe(risk_clean(:))) / max(nanstd_safe(risk_clean(:)), eps);
                z_dose = (dose_clean(:) - nanmean_safe(dose_clean(:))) / max(nanstd_safe(dose_clean(:)), eps);
                combined_score = z_risk - z_dose;  % high risk + low dose = worst
                [~, ~, ~, combined_auc] = perfcurve(event_combined, combined_score, true);
                result.combined_auc = combined_auc;
                fprintf('      Combined AUC: %.3f\n', combined_auc);
            catch ME_auc
                fprintf('      ⚠️ Combined AUC failed: %s\n', ME_auc.message);
            end
        end

        if isempty(method_results)
            method_results = result;
        else
            method_results(end+1) = result; %#ok<AGROW>
        end
    end

    concordance.method_results = method_results;

    % Summarize
    if ~isempty(method_results)
        kappas = [method_results.cohen_kappa];
        mean_kappa = nanmean_safe(kappas(:));
        if mean_kappa > 0.6
            concordance.summary = sprintf('Models are largely concordant (mean kappa=%.2f). Risk model and dose coverage identify similar patients.', mean_kappa);
        elseif mean_kappa > 0.2
            concordance.summary = sprintf('Models show moderate agreement (mean kappa=%.2f). Combining both may improve prediction.', mean_kappa);
        else
            concordance.summary = sprintf('Models are largely complementary (mean kappa=%.2f). Each identifies different at-risk patients; combined use is recommended.', mean_kappa);
        end
    else
        concordance.summary = 'No method results computed.';
    end

    fprintf('\n   Summary: %s\n', concordance.summary);

    % Generate summary heatmap
    if ~isempty(method_results)
        try
            fig = figure('Visible', 'off', 'Position', [100 100 600 400]);
            labels = {method_results.method_name};
            kappas = [method_results.cohen_kappa];
            bar(kappas);
            set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, ...
                'FontSize', 7, 'XTickLabelRotation', 45);
            ylabel('Cohen''s Kappa');
            title(sprintf('Risk-Dose Concordance (%s)', config_struct.dwi_type_name));
            hold on;
            yline(0.6, '--g', 'Substantial');
            yline(0.2, '--r', 'Poor');
            hold off;
            saveas(fig, fullfile(output_folder, ...
                sprintf('risk_dose_concordance_%s.png', config_struct.dwi_type_name)));
            close(fig);
        catch ME_fig
            fprintf('⚠️ Could not generate concordance figure: %s\n', ME_fig.message);
        end
    end

    fprintf('\n✅ Risk-dose concordance analysis complete.\n');
    diary off;
end
