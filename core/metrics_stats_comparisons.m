function metrics_stats_comparisons(valid_pts, lf_group, metric_sets, set_names, time_labels, dtype_label, output_folder, dataloc, nTp, ADC_abs, D_abs, f_abs, Dstar_abs)
% METRICS_STATS_COMPARISONS — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 4a/5 of the metrics step: Comparisons (Rank-sum & Mixed-Effects)
% Performs univariate hypothesis testing (Wilcoxon Rank-Sum) to compare parameter
% distributions between Local Control (LC) and Local Failure (LF) groups.
%
% Inputs:
%   valid_pts         - Logical mask of patients with valid survival/failure data
%   lf_group          - Local failure grouping variable (0 or 1)
%   metric_sets       - Cell array grouped by measurement type (e.g. Abs, Pct)
%   set_names         - String descriptors for the grouped metrics
%   time_labels       - String labels for timepoints (Fx1, Fx2, etc.)
%   dtype_label       - String indicating DWI pipeline
%   output_folder     - Output directory for generated figures
%   dataloc           - Input/output directory for data CSVs
%   nTp               - Total timepoints
%   *_abs             - Absolute values used in the GLME models
%
% Outputs:
%   None. Exports statistical tables (CSV) and figures.
%

fprintf('  --- SECTION 7: Univariate Analysis ---\n');

figure_titles = {
    '1. Absolute DWI/IVIM Metrics vs Local Failure (Wilcoxon Rank-Sum)', ...
    '2. Percent Change Metrics vs Local Failure (Wilcoxon Rank-Sum)', ...
    '3. Target Coverage (D95): Whole GTV vs Resistant Sub-volumes (Wilcoxon Rank-Sum)', ...
    '4. Target Coverage (V50): Whole GTV vs Resistant Sub-volumes (Wilcoxon Rank-Sum)'
};

p_val_store = struct('p_vals', {});

for s = 1:length(metric_sets)
    current_metrics = metric_sets{s};
    current_names = set_names{s};
    
    fig = figure('Name', [figure_titles{s} ' — ' dtype_label], 'Position', [50, 50, 1600, 1000]);
    sgtitle([figure_titles{s} ' (' dtype_label ')'], 'FontSize', 16, 'FontWeight', 'bold');
    
    num_rows = length(current_metrics);
    max_cols = 0;
    for m=1:num_rows
        max_cols = max(max_cols, size(current_metrics{m}, 2));
    end
    cols_to_plot = min(max_cols, length(time_labels));
    
    p_val_store(s).p_vals = nan(num_rows, cols_to_plot);

    plot_idx = 1;
    
    for m = 1:num_rows
        metric_data = current_metrics{m};
        
        for tp = 1:cols_to_plot
            subplot(num_rows, cols_to_plot, plot_idx);
            
            if tp > size(metric_data, 2)
                 axis off;
                 plot_idx = plot_idx + 1;
                 continue;
            end

            y_raw = metric_data(valid_pts, tp);
            has_data = ~isnan(y_raw);
            y = y_raw(has_data);
            g = lf_group(has_data);
            
            p = perform_statistical_test(y, g, 'ranksum');

            if ~isnan(p)
                p_val_store(s).p_vals(m, tp) = p;

                boxplot(y, g, 'Labels', {'LC (0)', 'LF (1)'});
                title_str = sprintf('%s - %s\np = %.3f', current_names{m}, time_labels{tp}, p);
                if p < 0.05
                    title(title_str, 'Color', 'r', 'FontWeight', 'bold');
                else
                    title(title_str, 'Color', 'k', 'FontWeight', 'normal');
                end
            else
                title(sprintf('%s - %s\n(Insufficient Data)', current_names{m}, time_labels{tp}), 'FontSize', 8);
                axis off; 
            end
            
            if m == num_rows
                xlabel('Outcome');
            end
            if tp == 1
                ylabel(current_names{m});
            end
            
            grid on;
            plot_idx = plot_idx + 1;
        end
    end
    subplot_scale = 0.92;
    allAx = findall(fig, 'Type', 'Axes');
    for k = 1:numel(allAx)
        pos = allAx(k).Position;
        allAx(k).Position = [pos(1), pos(2) * subplot_scale, pos(3), pos(4) * subplot_scale];
    end
    saveas(gcf, fullfile(output_folder, sprintf('Metric_Set_%d_%s.png', s, dtype_label)));
    close(gcf);
end

fprintf('  --- SECTION 8: Compile and Export Significant Results ---\n');
total_checks = 0;
for s = 1:length(metric_sets)
    total_checks = total_checks + length(metric_sets{s}) * length(time_labels);
end

sig_metric = cell(total_checks, 1);
sig_fraction = cell(total_checks, 1);
sig_pval = nan(total_checks, 1);
sig_mean_LC = nan(total_checks, 1);
sig_mean_LF = nan(total_checks, 1);
sig_count = 0;

for s = 1:length(metric_sets)
    current_metrics = metric_sets{s};
    current_names = set_names{s};
    num_metrics = length(current_metrics);
    
    for m = 1:num_metrics
        metric_data = current_metrics{m};
        cols_to_plot = min(size(metric_data, 2), length(time_labels));
        
        for tp = 1:cols_to_plot
            if tp <= size(p_val_store(s).p_vals, 2)
                p = p_val_store(s).p_vals(m, tp);
            else
                p = nan;
            end
            
            if ~isnan(p) && p < 0.05
                y_raw = metric_data(valid_pts, tp);
                has_data = ~isnan(y_raw);
                y = y_raw(has_data);
                g = lf_group(has_data);

                mean_LC = mean(y(g == 0), 'omitnan');
                mean_LF = mean(y(g == 1), 'omitnan');
                
                sig_count = sig_count + 1;
                sig_metric{sig_count, 1} = current_names{m};
                sig_fraction{sig_count, 1} = time_labels{tp};
                sig_pval(sig_count, 1) = p;
                sig_mean_LC(sig_count, 1) = mean_LC;
                sig_mean_LF(sig_count, 1) = mean_LF;
            end
        end
    end
end

if sig_count < total_checks
    sig_metric = sig_metric(1:sig_count, :);
    sig_fraction = sig_fraction(1:sig_count, :);
    sig_pval = sig_pval(1:sig_count, :);
    sig_mean_LC = sig_mean_LC(1:sig_count, :);
    sig_mean_LF = sig_mean_LF(1:sig_count, :);
end

if ~isempty(sig_pval)
    sig_results_table = table(sig_metric, sig_fraction, sig_pval, sig_mean_LC, sig_mean_LF, ...
        'VariableNames', {'Metric', 'Timepoint', 'P_Value', 'Mean_LC', 'Mean_LF'});
    sig_results_table = sortrows(sig_results_table, 'P_Value');
    
    disp('----- Significant Findings (p < 0.05) -----');
    disp(sig_results_table);
    
    export_filename = fullfile(dataloc, 'Significant_LF_Metrics.csv');
    writetable(sig_results_table, export_filename);
    fprintf('Saved significant results to: %s\n', export_filename);
else
    disp('No significant differences (p < 0.05) found between LC and LF groups for these metrics.');
end

fprintf('  --- SECTION 9: FDR Correction ---\n');
if ~isempty(sig_pval)
    fprintf('\n----- PER-TIMEPOINT FDR (Benjamini-Hochberg, Q < 0.05) -----\n');
    for tp = 1:length(time_labels)
        tp_pvals  = [];
        tp_labels = {};
        for s = 1:length(metric_sets)
            current_metrics = metric_sets{s};
            current_names = set_names{s};
            for mi = 1:length(current_metrics)
                if tp <= size(p_val_store(s).p_vals, 2)
                    p = p_val_store(s).p_vals(mi, tp);
                else
                    p = nan;
                end

                if ~isnan(p)
                    tp_pvals(end+1, 1)  = p;
                    tp_labels{end+1, 1} = current_names{mi};
                end
            end
        end
        
        if isempty(tp_pvals), continue; end
        
        n_tp = length(tp_pvals);
        [p_sort, sort_id] = sort(tp_pvals);
        q_tp = zeros(n_tp, 1);
        q_tp(n_tp) = p_sort(n_tp);
        for ii = n_tp-1:-1:1
            q_tp(ii) = min(q_tp(ii+1), p_sort(ii) * (n_tp / ii));
        end
        q_tp = min(q_tp, 1);
        q_unsorted = zeros(n_tp, 1);
        q_unsorted(sort_id) = q_tp;
        
        tp_table = table(tp_labels, tp_pvals, q_unsorted, ...
            'VariableNames', {'Metric', 'Raw_P', 'FDR_Q'});
        sig_tp = tp_table(tp_table.FDR_Q < 0.05, :);
        
        fprintf('\n  Timepoint: %s (family size = %d)\n', time_labels{tp}, n_tp);
        if isempty(sig_tp)
            fprintf('    None survived FDR correction.\n');
        else
            disp(sig_tp);
            writetable(sig_tp, fullfile(dataloc, sprintf('FDR_Sig_%s.csv', strrep(time_labels{tp},' ','_'))));
        end
    end
end

% Let's also run the GLME Mixed-effects model from the end of the script
fprintf('\n--- LONGITUDINAL MIXED-EFFECTS MODEL (GLME) ---\n');
long_PatientID = [];
long_Timepoint = [];
long_ADC = [];
long_D = [];
long_f = [];
long_Dstar = [];
long_LF = [];

patient_indices = find(valid_pts);
for i = 1:length(patient_indices)
    p_idx = patient_indices(i);
    for t = 1:nTp
        if ~isnan(ADC_abs(p_idx, t)) || ~isnan(D_abs(p_idx, t)) || ~isnan(f_abs(p_idx, t)) || ~isnan(Dstar_abs(p_idx, t))
            long_PatientID = [long_PatientID; i]; 
            long_Timepoint = [long_Timepoint; t];
            long_ADC = [long_ADC; ADC_abs(p_idx, t)];
            long_D = [long_D; D_abs(p_idx, t)];
            long_f = [long_f; f_abs(p_idx, t)];
            long_Dstar = [long_Dstar; Dstar_abs(p_idx, t)];
            long_LF = [long_LF; lf_group(i)];
        end
    end
end

glme_table = table(categorical(long_PatientID), long_Timepoint, ...
    long_ADC, long_D, long_f, long_Dstar, long_LF, ...
    'VariableNames', {'PatientID', 'Timepoint', 'ADC', 'D', 'f', 'Dstar', 'LF'});

clean_idx = ~isnan(glme_table.ADC) & ~isnan(glme_table.D) & ~isnan(glme_table.f) & ~isnan(glme_table.Dstar);
glme_table_clean = glme_table(clean_idx, :);
baseline_idx = glme_table_clean.Timepoint == 1;

mean_ADC_base = mean(glme_table_clean.ADC(baseline_idx));
std_ADC_base = std(glme_table_clean.ADC(baseline_idx));
mean_D_base = mean(glme_table_clean.D(baseline_idx));
std_D_base = std(glme_table_clean.D(baseline_idx));
mean_f_base = mean(glme_table_clean.f(baseline_idx));
std_f_base = std(glme_table_clean.f(baseline_idx));
mean_Dstar_base = mean(glme_table_clean.Dstar(baseline_idx));
std_Dstar_base = std(glme_table_clean.Dstar(baseline_idx));

glme_table_clean.ADC_z = (glme_table_clean.ADC - mean_ADC_base) / std_ADC_base;
glme_table_clean.D_z = (glme_table_clean.D - mean_D_base) / std_D_base;
glme_table_clean.f_z = (glme_table_clean.f - mean_f_base) / std_f_base;
glme_table_clean.Dstar_z = (glme_table_clean.Dstar - mean_Dstar_base) / std_Dstar_base;

glme_table_clean.LF = categorical(glme_table_clean.LF);
glme_table_clean.Timepoint = categorical(glme_table_clean.Timepoint);

biomarkers = {'ADC_z', 'D_z', 'f_z', 'Dstar_z'};
warning('off', 'all');
for b = 1:length(biomarkers)
    bm = biomarkers{b};
    formula = sprintf('%s ~ 1 + LF * Timepoint + (1|PatientID)', bm);
    try
        glme = fitglme(glme_table_clean, formula, 'OptimizerOptions', statset('MaxIter', 10000));
        fprintf('\n--- %s ---\n', formula);
        
        anova_res = anova(glme);
        row_idx = find(strcmp(anova_res.Term, 'LF:Timepoint'));
        if ~isempty(row_idx)
            pval = anova_res.pValue(row_idx);
            fprintf('Interaction P-Value (LF * Timepoint): %.4f\n', pval);
            if pval < 0.05
                fprintf('  -> SIGNIFICANT DIFFERENCE in trajectory between LC and LF groups.\n');
            else
                fprintf('  -> No significant difference in trajectory between LC and LF groups.\n');
            end
        else
            disp(anova_res);
        end
    catch ME
        fprintf('GLME model for %s failed to converge: %s\n', bm, ME.message);
    end
end
warning('on', 'all');

end
