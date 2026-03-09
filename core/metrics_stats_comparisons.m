function metrics_stats_comparisons(valid_pts, lf_group, metric_sets, set_names, time_labels, dtype_label, output_folder, dataloc, nTp, ADC_abs, D_abs, f_abs, Dstar_abs)
% METRICS_STATS_COMPARISONS — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 4a/5 of the metrics step: Comparisons (Rank-sum & Mixed-Effects)
% Performs univariate hypothesis testing (Wilcoxon Rank-Sum) to compare parameter
% distributions between Local Control (LC) and Local Failure (LF) groups.
%
% ANALYTICAL OVERVIEW:
%   This module performs two complementary statistical analyses:
%
%   1. WILCOXON RANK-SUM TESTS (univariate, per-timepoint)
%      Tests whether the distribution of each diffusion biomarker differs
%      between patients who achieved local control (LC) vs. those who
%      experienced local failure (LF).  Wilcoxon rank-sum (Mann-Whitney U)
%      is used instead of the t-test because:
%        - Diffusion parameters are often non-normally distributed
%          (right-skewed ADC, bimodal f in heterogeneous tumours)
%        - Small sample sizes make normality assumptions unreliable
%        - Rank-based tests are robust to outliers
%      The clinical hypothesis: patients whose tumours show smaller ADC/D
%      increases during RT (reflecting persistent cellularity) are more
%      likely to experience local failure.
%
%   2. GENERALISED LINEAR MIXED-EFFECTS MODEL (GLME, longitudinal)
%      Models the LF*Timepoint interaction to test whether LC and LF groups
%      have divergent biomarker trajectories over the course of treatment.
%      The random intercept/slope per patient accounts for within-patient
%      correlation across repeated measures.  A significant interaction
%      means the treatment response curves differ between groups — the key
%      biological signal for early response prediction.
%
%   Multiple comparison correction is critical: with ~48 rank-sum tests
%   (4 metrics x 2 sets x 6 timepoints), Benjamini-Hochberg FDR controls
%   the expected proportion of false discoveries.  The 4 GLME interaction
%   tests use the more conservative Holm-Bonferroni step-down procedure
%   to control family-wise error rate.
%
% Inputs:
%   valid_pts         - Logical mask of patients with valid survival/failure data
%   lf_group          - Local failure grouping variable (0=LC, 1=LF, 2=competing risk)
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

% Diary: capture console output to output_folder
diary_file = fullfile(output_folder, ['metrics_stats_comparisons_output_' dtype_label '.txt']);
if exist(diary_file, 'file'), delete(diary_file); end
diary(diary_file);

figure_titles = {
    '1. Absolute DWI/IVIM Metrics vs Local Failure (Wilcoxon Rank-Sum)', ...
    '2. Percent Change Metrics vs Local Failure (Wilcoxon Rank-Sum)', ...
    '3. Target Coverage (D95): Whole GTV vs Resistant Sub-volumes (Wilcoxon Rank-Sum)', ...
    '4. Target Coverage (V50): Whole GTV vs Resistant Sub-volumes (Wilcoxon Rank-Sum)'
};

% Storage for all raw p-values from Wilcoxon rank-sum tests, indexed by
% metric set.  These are collected across all timepoints and metrics,
% then passed to the global BH FDR correction in Section 8.
p_val_store = struct('p_vals', {});

% Iterate over metric sets (Set 1: absolute values, Set 2: percent changes,
% and optionally Sets 3-4 for dose coverage metrics).
n_metric_sets = length(metric_sets);
for s = 1:n_metric_sets
    text_progress_bar(s, n_metric_sets, 'Univariate analysis');
    current_metrics = metric_sets{s};
    current_names = set_names{s};
    
    fig = figure('Name', [figure_titles{s} ' — ' dtype_label], 'Position', [50, 50, 1600, 1000]);
    if exist('OCTAVE_VERSION', 'builtin')
        % sgtitle not supported in Octave
    else
        sgtitle([figure_titles{s} ' (' dtype_label ')'], 'FontSize', 16, 'FontWeight', 'bold');
    end
    
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

            % Exclude competing risk patients (lf==2) from the two-group
            % comparison.  These patients died of non-cancer causes before
            % local failure could be assessed — they are neither LC nor LF,
            % and including them in either group would bias the comparison.
            % Wilcoxon rank-sum requires exactly two groups; competing risk
            % patients are handled separately in the survival analysis
            % (metrics_survival.m) via the CSH framework.
            non_competing = (g <= 1);
            y = y(non_competing);
            g = g(non_competing);

            % Wilcoxon rank-sum test: non-parametric comparison of biomarker
            % distributions between LC and LF groups at this timepoint.
            p = perform_statistical_test(y, g, 'ranksum');

            if ~isnan(p)
                p_val_store(s).p_vals(m, tp) = p;

                if exist('OCTAVE_VERSION', 'builtin')
                    try
                        boxplot(y, g);
                    catch
                        plot(g + 1, y, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', [0 0.4470 0.7410]);
                        xlim([0.5, 2.5]);
                        set(gca, 'XTick', [1, 2]);
                        set(gca, 'XTickLabel', {'LC (0)', 'LF (1)'});
                    end
                else
                    boxplot(y, g, 'Labels', {'LC (0)', 'LF (1)'});
                end
                title_str = sprintf('%s - %s\np = %.3f', current_names{m}, time_labels{tp}, p);
                % Highlight nominally significant results in red for visual
                % screening.  Note: these are uncorrected p-values — FDR
                % correction is applied globally in Section 8 below to
                % determine which results survive multiple testing.
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
    if exist('OCTAVE_VERSION', 'builtin')
        % findall with 'Type' 'Axes' works differently in Octave and returns handles that might not have Position
        % skip subplot scaling in mock
    else
        allAx = findall(fig, 'Type', 'Axes');
        for k = 1:numel(allAx)
            pos = allAx(k).Position;
            allAx(k).Position = [pos(1), pos(2) * subplot_scale, pos(3), pos(4) * subplot_scale];
        end
    end
    set(findall(gcf, 'Type', 'Axes'), 'Toolbar', []);
    saveas(gcf, fullfile(output_folder, sprintf('Metric_Set_%d_%s.png', s, dtype_label)));
    close(gcf);
end

% ---- SECTION 8 & 9: FDR Correction FIRST, Then Report Significant Results ----
% With ~48 tests (4 metrics × 2 sets × 6 timepoints), reporting at nominal
% alpha=0.05 before correction yields ~2.4 expected false positives.  We
% therefore compute Benjamini-Hochberg FDR correction on the full family of
% comparisons BEFORE flagging any individual result as significant.
fprintf('  --- SECTION 8: Global FDR Correction (Benjamini-Hochberg) ---\n');

% Collect ALL p-values across every timepoint and metric set into one
% family so that the BH procedure controls the false discovery rate over
% the full set of comparisons (n_sets * n_metrics * n_timepoints).
max_metrics = sum(cellfun(@length, metric_sets));
max_total   = max_metrics * length(time_labels);
all_pvals   = nan(max_total, 1);
all_labels  = cell(max_total, 1);
all_tp_idx  = zeros(max_total, 1);
all_set_idx = zeros(max_total, 1);
all_met_idx = zeros(max_total, 1);
total_count = 0;

for tp = 1:length(time_labels)
    for s = 1:length(metric_sets)
        current_names = set_names{s};
        for mi = 1:length(metric_sets{s})
            if tp <= size(p_val_store(s).p_vals, 2)
                p = p_val_store(s).p_vals(mi, tp);
            else
                p = nan;
            end
            if ~isnan(p)
                total_count = total_count + 1;
                all_pvals(total_count)  = p;
                all_labels{total_count} = sprintf('%s @ %s', current_names{mi}, time_labels{tp});
                all_tp_idx(total_count) = tp;
                all_set_idx(total_count) = s;
                all_met_idx(total_count) = mi;
            end
        end
    end
end

% Compute FDR-corrected q-values for the entire family
q_unsorted = ones(total_count, 1);
if total_count > 0
    all_pvals   = all_pvals(1:total_count);
    all_labels  = all_labels(1:total_count);
    all_tp_idx  = all_tp_idx(1:total_count);
    all_set_idx = all_set_idx(1:total_count);
    all_met_idx = all_met_idx(1:total_count);

    % Benjamini-Hochberg FDR correction on the full family of comparisons.
    q_unsorted = benjamini_hochberg_fdr(all_pvals);
    n_all = length(all_pvals);

    fprintf('  Global family size = %d comparisons\n', n_all);
end

fprintf('  --- SECTION 9: Compile and Export FDR-Significant Results ---\n');

% Now report only results that survive FDR correction (q < 0.05)
sig_metric = cell(total_count, 1);
sig_fraction = cell(total_count, 1);
sig_pval = nan(total_count, 1);
sig_qval = nan(total_count, 1);
sig_mean_LC = nan(total_count, 1);
sig_mean_LF = nan(total_count, 1);
sig_count = 0;

for idx = 1:total_count
    if q_unsorted(idx) < 0.05
        s  = all_set_idx(idx);
        mi = all_met_idx(idx);
        tp = all_tp_idx(idx);
        current_metrics = metric_sets{s};
        current_names = set_names{s};
        metric_data = current_metrics{mi};

        y_raw = metric_data(valid_pts, tp);
        has_data = ~isnan(y_raw);
        y = y_raw(has_data);
        g = lf_group(has_data);

        % Exclude competing risk patients for group summary
        non_competing = (g <= 1);
        y = y(non_competing);
        g = g(non_competing);

        mean_LC = mean(y(g == 0), 'omitnan');
        mean_LF = mean(y(g == 1), 'omitnan');

        sig_count = sig_count + 1;
        sig_metric{sig_count, 1} = current_names{mi};
        sig_fraction{sig_count, 1} = time_labels{tp};
        sig_pval(sig_count, 1) = all_pvals(idx);
        sig_qval(sig_count, 1) = q_unsorted(idx);
        sig_mean_LC(sig_count, 1) = mean_LC;
        sig_mean_LF(sig_count, 1) = mean_LF;
    end
end

if sig_count > 0
    sig_metric = sig_metric(1:sig_count, :);
    sig_fraction = sig_fraction(1:sig_count, :);
    sig_pval = sig_pval(1:sig_count, :);
    sig_qval = sig_qval(1:sig_count, :);
    sig_mean_LC = sig_mean_LC(1:sig_count, :);
    sig_mean_LF = sig_mean_LF(1:sig_count, :);

    sig_results_table = table(sig_metric, sig_fraction, sig_pval, sig_qval, sig_mean_LC, sig_mean_LF, ...
        'VariableNames', {'Metric', 'Timepoint', 'Raw_P', 'FDR_Q', 'Mean_LC', 'Mean_LF'});
    sig_results_table = sortrows(sig_results_table, 'FDR_Q');

    disp('----- FDR-Significant Findings (BH q < 0.05) -----');
    disp(sig_results_table);

    export_filename = fullfile(output_folder, 'Significant_LF_Metrics.csv');
    writetable(sig_results_table, export_filename);
    fprintf('Saved FDR-significant results to: %s\n', export_filename);
else
    disp('No significant differences survived FDR correction (BH q < 0.05).');
end

% Full FDR table for reference — export all tests (significant and not)
% as a CSV for downstream analysis scripts (parse_csv_results.py) to
% compare significance patterns across DWI processing pipelines.
if total_count > 0
    fdr_table = table(all_labels, all_pvals, q_unsorted, ...
        'VariableNames', {'Metric_Timepoint', 'Raw_P', 'FDR_Q'});
    sig_global = fdr_table(fdr_table.FDR_Q < 0.05, :);

    if ~isempty(sig_global)
        writetable(sig_global, fullfile(output_folder, 'FDR_Sig_Global.csv'));
    end

    % Per-timepoint breakdown for readability
    for tp = 1:length(time_labels)
        tp_mask = (all_tp_idx == tp);
        sig_tp = fdr_table(tp_mask & fdr_table.FDR_Q < 0.05, :);
        fprintf('\n  Timepoint: %s — %d significant (global FDR)\n', time_labels{tp}, height(sig_tp));
        if ~isempty(sig_tp)
            disp(sig_tp);
        end
    end
end

% GLME Mixed-effects model requires categorical() and fitglme() which are
% not available in Octave. Skip this section in Octave.
if exist('OCTAVE_VERSION', 'builtin')
    fprintf('\n--- LONGITUDINAL MIXED-EFFECTS MODEL (GLME) ---\n');
    fprintf('  Skipped: fitglme/categorical not available in Octave.\n');
else
    fprintf('\n--- LONGITUDINAL MIXED-EFFECTS MODEL (GLME) ---\n');
    % Construct the long-format data table required by fitglme.  Each row
    % is one patient-timepoint observation.  The wide-format arrays
    % (patients x timepoints) are "melted" into columns so that the mixed-
    % effects model can estimate within-patient correlation via the random
    % effects specification.
    patient_indices = find(valid_pts);
    max_obs = length(patient_indices) * nTp;

    % Pre-allocate long-format columns (will be trimmed after population)
    long_PatientID = nan(max_obs, 1);
    long_Timepoint = nan(max_obs, 1);
    long_ADC = nan(max_obs, 1);
    long_D = nan(max_obs, 1);
    long_f = nan(max_obs, 1);
    long_Dstar = nan(max_obs, 1);
    long_LF = nan(max_obs, 1);

    obs_idx = 0;
    n_competing_excluded = 0;
    for i = 1:length(patient_indices)
        p_idx = patient_indices(i);
        % Exclude competing-risk patients (lf==2) from the GLME.  Unlike
        % the CSH Cox model where competing events are correctly treated as
        % censored, the GLME models an observed outcome state.  Recoding
        % non-cancer deaths as "local control" would conflate patients who
        % achieved LC with those who died before failure could be observed.
        if lf_group(i) == 2
            n_competing_excluded = n_competing_excluded + 1;
            continue;
        end
        for t = 1:nTp
            if ~isnan(ADC_abs(p_idx, t)) || ~isnan(D_abs(p_idx, t)) || ~isnan(f_abs(p_idx, t)) || ~isnan(Dstar_abs(p_idx, t))
                obs_idx = obs_idx + 1;
                long_PatientID(obs_idx) = i;
                long_Timepoint(obs_idx) = t;
                long_ADC(obs_idx) = ADC_abs(p_idx, t);
                long_D(obs_idx) = D_abs(p_idx, t);
                long_f(obs_idx) = f_abs(p_idx, t);
                long_Dstar(obs_idx) = Dstar_abs(p_idx, t);
                long_LF(obs_idx) = lf_group(i);
            end
        end
    end
    if n_competing_excluded > 0
        n_total_glme = length(patient_indices);
        pct_excluded = 100 * n_competing_excluded / n_total_glme;
        fprintf('  ⚠️  GLME: Excluded %d/%d (%.1f%%) competing-risk patients (lf==2).\n', ...
            n_competing_excluded, n_total_glme, pct_excluded);
        fprintf('      These patients experienced non-cancer death without documented LF.\n');
        if pct_excluded > 20
            fprintf('  ⚠️  >20%% of cohort excluded — interpret GLME results with caution.\n');
            fprintf('      Consider Fine-Gray subdistribution hazard model for sensitivity analysis.\n');
        end
    end

    long_PatientID = long_PatientID(1:obs_idx);
    long_Timepoint = long_Timepoint(1:obs_idx);
    long_ADC = long_ADC(1:obs_idx);
    long_D = long_D(1:obs_idx);
    long_f = long_f(1:obs_idx);
    long_Dstar = long_Dstar(1:obs_idx);
    long_LF = long_LF(1:obs_idx);

    % Filter to complete cases before table construction.  GLME requires
    % all predictors to be non-NaN for each observation.  Rows with any
    % NaN biomarker are dropped.  This reduces sample size but avoids
    % bias from asymmetric missingness (e.g., IVIM D* missing more often
    % than ADC due to bi-exponential fitting failures).
    clean_idx = ~isnan(long_ADC) & ~isnan(long_D) & ~isnan(long_f) & ~isnan(long_Dstar);
    long_PatientID = long_PatientID(clean_idx);
    long_Timepoint = long_Timepoint(clean_idx);
    long_ADC = long_ADC(clean_idx);
    long_D = long_D(clean_idx);
    long_f = long_f(clean_idx);
    long_Dstar = long_Dstar(clean_idx);
    long_LF = long_LF(clean_idx);

    % Z-SCORE STANDARDISATION FOR GLME
    % Biomarkers are z-scored using BASELINE (Fx1) mean and SD to make
    % coefficient magnitudes comparable across parameters with different
    % physiological scales (ADC ~0.001 mm^2/s vs f ~0.1).  Using baseline-
    % only statistics (rather than pooling all timepoints) prevents
    % treatment-induced changes from contaminating the standardisation,
    % which would attenuate the LF*Timepoint interaction effect.
    baseline_idx = long_Timepoint == 1;

    if sum(baseline_idx) < 2
        fprintf('  ⚠️  Fewer than 2 baseline observations after NaN cleaning — skipping GLME.\n');
        if ~isempty(output_folder), diary off; end
        return;
    end

    mean_ADC_base = mean(long_ADC(baseline_idx));
    std_ADC_base = std(long_ADC(baseline_idx));
    mean_D_base = mean(long_D(baseline_idx));
    std_D_base = std(long_D(baseline_idx));
    mean_f_base = mean(long_f(baseline_idx));
    std_f_base = std(long_f(baseline_idx));
    mean_Dstar_base = mean(long_Dstar(baseline_idx));
    std_Dstar_base = std(long_Dstar(baseline_idx));
    % Guard against zero std (constant baseline) producing Inf z-scores.
    % When std is zero, z-scoring is undefined — set to NaN so the GLME
    % for that biomarker is skipped, and log which biomarkers are affected.
    zero_std_names = {};
    if std_ADC_base == 0, std_ADC_base = NaN; zero_std_names{end+1} = 'ADC'; end
    if std_D_base == 0, std_D_base = NaN; zero_std_names{end+1} = 'D'; end
    if std_f_base == 0, std_f_base = NaN; zero_std_names{end+1} = 'f'; end
    if std_Dstar_base == 0, std_Dstar_base = NaN; zero_std_names{end+1} = 'Dstar'; end
    if ~isempty(zero_std_names)
        fprintf('  ⚠️  Zero baseline std for: %s — z-scoring undefined, GLME will be skipped for these.\n', ...
            strjoin(zero_std_names, ', '));
    end

    long_ADC_z = (long_ADC - mean_ADC_base) / std_ADC_base;
    long_D_z = (long_D - mean_D_base) / std_D_base;
    long_f_z = (long_f - mean_f_base) / std_f_base;
    long_Dstar_z = (long_Dstar - mean_Dstar_base) / std_Dstar_base;

    % Build the table with all columns at once (including z-scores)
    glme_table_clean = table(categorical(long_PatientID), categorical(long_Timepoint), ...
        long_ADC, long_D, long_f, long_Dstar, categorical(long_LF), ...
        long_ADC_z, long_D_z, long_f_z, long_Dstar_z, ...
        'VariableNames', {'PatientID', 'Timepoint', 'ADC', 'D', 'f', 'Dstar', 'LF', ...
                          'ADC_z', 'D_z', 'f_z', 'Dstar_z'});

    % Fit one GLME per biomarker rather than a multivariate model because:
    %   1. Collinearity between ADC and D (both reflect cellularity) would
    %      make multivariate coefficient interpretation unreliable
    %   2. Missing data patterns differ across biomarkers (IVIM parameters
    %      have more missingness than ADC due to bi-exponential fit failures)
    %   3. Separate models allow direct comparison of which biomarker's
    %      trajectory most strongly discriminates LC vs LF
    biomarkers = {'ADC_z', 'D_z', 'f_z', 'Dstar_z'};
    n_biomarkers = length(biomarkers);
    glme_pvals = nan(n_biomarkers, 1);
    warning('off', 'all');
    for b = 1:n_biomarkers
        text_progress_bar(b, n_biomarkers, 'Fitting GLME models');
        bm = biomarkers{b};
        % GLME Formula: Biomarker ~ 1 + LF * Timepoint + (1 + Timepoint|PatientID)
        %
        % Fixed effects:
        %   1           — global intercept (mean baseline biomarker value)
        %   LF          — baseline difference between LC and LF groups
        %   Timepoint   — overall time trend (average treatment response)
        %   LF*Timepoint — THE KEY INTERACTION: does the treatment response
        %                  trajectory differ between LC and LF?  A significant
        %                  interaction means the groups diverge over time,
        %                  e.g., LC patients show rising ADC while LF patients
        %                  show stable ADC — the biological basis for early
        %                  treatment response prediction.
        %
        % Random effects:
        %   (1 + Timepoint|PatientID) — random intercept (each patient starts
        %     at a different baseline) AND random slope (each patient's
        %     trajectory over time differs).  This accounts for within-patient
        %     correlation across repeated DWI scans.
        %
        % Fall back to random intercept only if the richer model fails to
        % converge (common with small N where the covariance matrix of
        % random effects becomes singular).
        formula_rs = sprintf('%s ~ 1 + LF * Timepoint + (1 + Timepoint|PatientID)', bm);
        formula_ri = sprintf('%s ~ 1 + LF * Timepoint + (1|PatientID)', bm);
        glme = [];
        used_formula = '';
        try
            glme = fitglme(glme_table_clean, formula_rs, 'OptimizerOptions', statset('MaxIter', 10000));
            used_formula = formula_rs;
        catch
            try
                glme = fitglme(glme_table_clean, formula_ri, 'OptimizerOptions', statset('MaxIter', 10000));
                used_formula = formula_ri;
                fprintf('  (random slope did not converge for %s; using random intercept only)\n', bm);
            catch ME
                fprintf('GLME model for %s failed to converge: %s\n', bm, ME.message);
            end
        end
        if ~isempty(glme)
            fprintf('\n--- %s ---\n', used_formula);

            anova_res = anova(glme);
            % MATLAB's fitglme may encode the interaction term with level
            % suffixes (e.g. 'LF_1:Timepoint_2') depending on categorical
            % encoding.  Use contains() for robust matching.
            row_idx = find(contains(anova_res.Term, 'LF') & contains(anova_res.Term, 'Timepoint') & contains(anova_res.Term, ':'));
            if numel(row_idx) > 1, row_idx = row_idx(1); end
            if ~isempty(row_idx)
                glme_pvals(b) = anova_res.pValue(row_idx);
                fprintf('Interaction P-Value (LF * Timepoint): %.4f\n', glme_pvals(b));
            else
                disp(anova_res);
            end
        end
    end
    warning('on', 'all');

    % Apply Holm-Bonferroni correction across the 4 GLME interaction tests
    % to control family-wise error rate (FWER).  Holm-Bonferroni is used
    % here instead of BH-FDR because the GLME tests are confirmatory
    % (testing a specific biological hypothesis about trajectory divergence)
    % rather than exploratory, so stronger FWER control is appropriate.
    % The step-down procedure is uniformly more powerful than Bonferroni
    % while maintaining the same FWER guarantee.
    valid_glme = ~isnan(glme_pvals);
    n_glme_tests = sum(valid_glme);
    if n_glme_tests > 0
        fprintf('\n--- GLME Interaction: Holm-Bonferroni Correction (%d tests) ---\n', n_glme_tests);
        [p_sort_glme, sort_order] = sort(glme_pvals(valid_glme));
        bm_valid = biomarkers(valid_glme);
        holm_rejected = false(n_glme_tests, 1);
        for gi = 1:n_glme_tests
            adj_alpha = 0.05 / (n_glme_tests - gi + 1);
            if gi > 1 && ~holm_rejected(gi-1)
                % Holm step-down: once a test is not rejected, all
                % subsequent (larger p-value) tests are also not rejected.
                is_sig = false;
            else
                is_sig = p_sort_glme(gi) < adj_alpha;
            end
            holm_rejected(gi) = is_sig;
            if is_sig, sig_str = 'SIGNIFICANT'; else, sig_str = 'not significant'; end
            fprintf('  %s: p=%.4f, adj_alpha=%.4f => %s\n', ...
                bm_valid{sort_order(gi)}, p_sort_glme(gi), adj_alpha, sig_str);
        end
    end
end

diary off;
end
