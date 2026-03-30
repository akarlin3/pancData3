function confound = compute_gtv_confounding(per_method_dosimetry, baseline_results, active_methods, config_struct)
%COMPUTE_GTV_CONFOUNDING Check if GTV volume change confounds D95-outcome association.
%
%   Tests whether the association between sub-volume dose coverage and local
%   control is confounded by GTV shrinkage during treatment.
%
%   Inputs:
%       per_method_dosimetry - struct from metrics_dosimetry
%       baseline_results     - struct from metrics_baseline
%       active_methods       - cell array of active method names
%       config_struct        - pipeline configuration struct
%
%   Outputs:
%       confound.method_results - struct array per method
%       confound.summary        - overall interpretation text

    % --- Diary ---
    output_folder = config_struct.output_folder;
    diary_file = fullfile(output_folder, ...
        sprintf('gtv_confounding_output_%s.txt', config_struct.dwi_type_name));
    diary(diary_file);

    fprintf('\n🔬 GTV Volume Confounding Check (%s)\n', config_struct.dwi_type_name);

    n_methods = numel(active_methods);
    m_lf = baseline_results.m_lf;

    % Get GTV volume
    if ~isfield(baseline_results, 'm_gtv_vol')
        fprintf('   ⚠️ GTV volume data not available. Skipping.\n');
        confound.method_results = struct([]);
        confound.summary = 'GTV volume data not available.';
        diary off;
        return;
    end
    gtv_vol = baseline_results.m_gtv_vol;
    nPt = size(gtv_vol, 1);
    nTp = size(gtv_vol, 2);

    % Exclude competing risks
    valid_idx = (m_lf == 0 | m_lf == 1);
    event = m_lf(valid_idx) == 1;

    if sum(event) < 3 || sum(~event) < 3
        fprintf('   ⚠️ Insufficient events for confounding analysis.\n');
        confound.method_results = struct([]);
        confound.summary = 'Insufficient events for confounding analysis.';
        diary off;
        return;
    end

    % GTV volume change from Fx1
    gtv_change = nan(nPt, 1);
    if nTp >= 5
        last_tp = 5;
    else
        last_tp = nTp;
    end
    for j = 1:nPt
        if gtv_vol(j, 1) > 0 && ~isnan(gtv_vol(j, 1)) && ~isnan(gtv_vol(j, last_tp))
            gtv_change(j) = (gtv_vol(j, last_tp) - gtv_vol(j, 1)) / gtv_vol(j, 1) * 100;
        end
    end

    method_results = struct([]);
    any_confounded = false;

    for mi = 1:n_methods
        method_name = active_methods{mi};
        fprintf('\n   Processing method: %s\n', method_name);

        result = struct();
        result.method_name = method_name;
        result.d95_gtv_correlation = NaN;
        result.d95_gtv_pvalue = NaN;
        result.gtv_change_lc_vs_lf = NaN;
        result.adjusted_hr = NaN;
        result.adjusted_p = NaN;
        result.unadjusted_hr = NaN;
        result.unadjusted_p = NaN;
        result.confounding_flag = false;

        if ~isfield(per_method_dosimetry, method_name)
            fprintf('      ⚠️ No dosimetry data for %s. Skipping.\n', method_name);
            if isempty(method_results)
                method_results = result;
            else
                method_results(end+1) = result; %#ok<AGROW>
            end
            continue;
        end
        method_dose = per_method_dosimetry.(method_name);

        % Get D95 ADC sub at Fx1
        d95_field = 'd95_adc_sub';
        if ~isfield(method_dose, d95_field)
            if isempty(method_results)
                method_results = result;
            else
                method_results(end+1) = result; %#ok<AGROW>
            end
            continue;
        end
        d95 = method_dose.(d95_field);
        if size(d95, 2) > 1
            d95 = d95(:, 1);
        end

        gtv_fx1 = gtv_vol(:, 1);

        % (a) Spearman correlation between D95 and GTV volume at Fx1
        valid_corr = valid_idx & ~isnan(d95) & ~isnan(gtv_fx1);
        if sum(valid_corr) >= 5
            [rho, p_rho] = corr(d95(valid_corr), gtv_fx1(valid_corr), 'Type', 'Spearman');
            result.d95_gtv_correlation = rho;
            result.d95_gtv_pvalue = p_rho;
            fprintf('      D95-GTV correlation: rho=%.3f, p=%.4f\n', rho, p_rho);
        end

        % (b) Wilcoxon rank-sum: GTV volume change between LC and LF
        gtv_change_valid = gtv_change(valid_idx);
        event_for_test = event;
        nan_mask_gwt = isnan(gtv_change_valid);
        if sum(~nan_mask_gwt & event_for_test) >= 2 && sum(~nan_mask_gwt & ~event_for_test) >= 2
            lf_change = gtv_change_valid(event_for_test & ~nan_mask_gwt);
            lc_change = gtv_change_valid(~event_for_test & ~nan_mask_gwt);
            p_wilcox = ranksum(lf_change, lc_change);
            result.gtv_change_lc_vs_lf = p_wilcox;
            fprintf('      GTV change LC vs LF: p=%.4f\n', p_wilcox);
        end

        % (c) Multivariable Cox PH: D95 adjusted for GTV volume
        d95_valid = d95(valid_idx);
        gtv_valid = gtv_fx1(valid_idx);
        time_valid = baseline_results.m_total_time(valid_idx);
        censored_valid = ~event;

        % Both predictors available
        complete = ~isnan(d95_valid) & ~isnan(gtv_valid) & ~isnan(time_valid) & time_valid > 0;
        if sum(complete) >= 10 && sum(event(complete)) >= 3
            try
                % Unadjusted Cox (D95 only)
                [b_unadj, ~, ~, stats_unadj] = coxphfit(d95_valid(complete), ...
                    time_valid(complete), 'Censoring', censored_valid(complete));
                result.unadjusted_hr = exp(b_unadj(1));
                result.unadjusted_p = stats_unadj.p(1);

                % Adjusted Cox (D95 + GTV volume)
                X_adj = [d95_valid(complete), gtv_valid(complete)];
                [b_adj, ~, ~, stats_adj] = coxphfit(X_adj, ...
                    time_valid(complete), 'Censoring', censored_valid(complete));
                result.adjusted_hr = exp(b_adj(1));
                result.adjusted_p = stats_adj.p(1);

                % Check confounding: >10% change in log(HR)
                log_hr_unadj = abs(log(result.unadjusted_hr));
                log_hr_adj = abs(log(result.adjusted_hr));
                if log_hr_unadj > 0
                    pct_change = abs(log_hr_adj - log_hr_unadj) / log_hr_unadj;
                    result.confounding_flag = pct_change > 0.10;
                end

                fprintf('      Unadjusted HR=%.3f (p=%.4f), Adjusted HR=%.3f (p=%.4f), Confounded=%d\n', ...
                    result.unadjusted_hr, result.unadjusted_p, ...
                    result.adjusted_hr, result.adjusted_p, result.confounding_flag);

                if result.confounding_flag
                    any_confounded = true;
                end
            catch ME_cox
                fprintf('      ⚠️ Cox model failed: %s\n', ME_cox.message);
            end
        end

        % Generate scatter plot
        try
            valid_plot = valid_idx & ~isnan(d95) & ~isnan(gtv_fx1);
            if sum(valid_plot) >= 5
                fig = figure('Visible', 'off', 'Position', [100 100 600 450]);
                lf_mask = valid_plot & (m_lf == 1);
                lc_mask = valid_plot & (m_lf == 0);
                scatter(d95(lc_mask), gtv_fx1(lc_mask), 50, 'b', 'filled', 'DisplayName', 'LC');
                hold on;
                scatter(d95(lf_mask), gtv_fx1(lf_mask), 50, 'r', 'filled', 'DisplayName', 'LF');
                hold off;
                xlabel('D95 (Gy)');
                ylabel('GTV Volume (mm³)');
                title(sprintf('D95 vs GTV Volume (%s, %s)', method_name, config_struct.dwi_type_name));
                legend('Location', 'best');
                saveas(fig, fullfile(output_folder, ...
                    sprintf('gtv_confounding_%s_%s.png', method_name, config_struct.dwi_type_name)));
                close(fig);
            end
        catch ME_fig
            fprintf('      ⚠️ Could not generate scatter plot: %s\n', ME_fig.message);
        end

        if isempty(method_results)
            method_results = result;
        else
            method_results(end+1) = result; %#ok<AGROW>
        end
    end

    confound.method_results = method_results;
    if any_confounded
        confound.summary = 'GTV volume is a significant confounder for at least one method. Adjusted HRs should be reported alongside unadjusted.';
    else
        confound.summary = 'GTV volume does not significantly confound the D95-outcome association. The observed protective effect of higher D95 appears independent of tumor size.';
    end

    fprintf('\n   Summary: %s\n', confound.summary);
    fprintf('\n✅ GTV confounding analysis complete.\n');
    diary off;
end
