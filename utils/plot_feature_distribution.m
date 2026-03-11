function plot_feature_distribution(vals, lf_group, metric_name, metric_unit, plot_type)
% PLOT_FEATURE_DISTRIBUTION Plots feature distribution as histogram or boxplot.
%
% Syntax:
%   plot_feature_distribution(vals, lf_group, metric_name, metric_unit, plot_type)
%
% Inputs:
%   vals        - [N x 1] numeric vector of feature values (e.g., median ADC,
%                 perfusion fraction f, or percent change from baseline).
%   lf_group    - [N x 1] numeric vector indicating clinical outcome:
%                   0 = Local Control (LC): tumor did not recur locally
%                   1 = Local Failure (LF): tumor recurred at the primary site
%                   2 = Competing risk (death without local recurrence); excluded
%   metric_name - Character vector or string for the metric name (used in title).
%   metric_unit - Character vector or string for the metric unit (used in ylabel or xlabel).
%   plot_type   - 'histogram' or 'boxplot' indicating the type of plot.
%
% Outputs:
%   None. Generates a histogram or boxplot visualization.
%
% --- Analytical Rationale ---
% Visualizing the distribution of DWI-derived biomarkers (ADC, D, f, D*)
% stratified by treatment outcome is a foundational step in identifying
% predictive imaging signatures. If a diffusion parameter shows clear
% separation between LC and LF groups, it may serve as an early biomarker
% for adaptive radiotherapy -- allowing clinicians to escalate dose or
% modify treatment for patients showing unfavorable diffusion patterns.
%
% Histograms reveal the full distributional shape (skewness, multi-modality)
% which is common in IVIM parameters due to biological heterogeneity.
% Boxplots provide a compact summary with statistical annotation, making
% them suitable for multi-panel comparisons across many features.
%
% The Wilcoxon rank-sum test (non-parametric) is used instead of a t-test
% because IVIM parameters -- especially the perfusion fraction f and
% pseudo-diffusion coefficient D* -- are frequently non-normally distributed
% with heavy right tails, violating the normality assumption of parametric
% tests.
%

    % --- Shape Normalization ---
    % Ensure column vectors to prevent implicit expansion issues. Row vectors
    % combined with logical indexing can produce unexpected matrix results
    % in MATLAB when the index vector orientation does not match.
    vals = vals(:);
    lf_group = lf_group(:);

    % Remove NaN entries and competing risk patients (lf==2) for accurate
    % processing.  Competing risk patients are excluded because boxplot
    % labels assume exactly 2 groups (LC/LF), and ranksum requires exactly
    % 2 groups.  Callers like metrics_stats_comparisons already filter, but
    % this guard prevents crashes when called with unfiltered data.
    has_data = ~isnan(vals) & (lf_group <= 1);
    vals_clean = vals(has_data);
    lf_clean = lf_group(has_data);

    if strcmpi(plot_type, 'histogram')
        % --- Histogram: Overlaid Outcome-Stratified Distributions ---
        % Splitting by outcome allows visual assessment of whether the
        % biomarker's distribution differs between LC and LF. Overlapping
        % distributions suggest weak discriminative power; well-separated
        % distributions suggest the feature may be clinically useful.
        vals_lc = vals_clean(lf_clean == 0);   % Local Control patients
        vals_lf = vals_clean(lf_clean == 1);   % Local Failure patients

        % --- Bin Edge Computation ---
        % Use 15 equally-spaced bins (16 edges) spanning the combined range
        % of both groups. Shared bin edges are essential so that the LC and
        % LF histograms are directly comparable -- different binning would
        % make visual overlap assessment misleading.
        if exist('OCTAVE_VERSION', 'builtin')
            min_val = min(vals_clean(~isnan(vals_clean)));
            max_val = max(vals_clean(~isnan(vals_clean)));
        else
            min_val = min(vals_clean, [], 'omitnan');
            max_val = max(vals_clean, [], 'omitnan');
        end

        if isempty(vals_clean)
            % No data: use a dummy [0,1] range to avoid linspace errors.
            edges = linspace(0, 1, 16);
        elseif min_val == max_val
            % All values identical (zero variance): create a 1-unit-wide
            % range centered on the value so the histogram is still renderable.
            edges = linspace(min_val - 0.5, max_val + 0.5, 16);
        else
            edges = linspace(min_val, max_val, 16);
        end

        % --- Render Overlaid Histograms ---
        % Semi-transparent (FaceAlpha=0.6) overlaid histograms allow visual
        % inspection of distributional overlap. The blue/orange color scheme
        % follows MATLAB's default colororder and is colorblind-accessible.
        if exist('OCTAVE_VERSION', 'builtin')
            % Octave's hist() treats 2nd arg as bin centers, not edges.
            % Use histc() with edges for consistent binning with MATLAB.
            if length(unique(edges)) == length(edges) && length(edges) > 1
                counts_lc = histc(vals_lc, edges);
                counts_lf = histc(vals_lf, edges);
                % histc puts exact-edge-max counts in the last bin;
                % merge last two bins to match MATLAB histogram behavior
                % (n_edges produces n_edges-1 bins).
                counts_lc(end-1) = counts_lc(end-1) + counts_lc(end);
                counts_lc(end) = [];
                counts_lf(end-1) = counts_lf(end-1) + counts_lf(end);
                counts_lf(end) = [];
                centers = (edges(1:end-1) + edges(2:end)) / 2;
                if length(unique(centers)) == length(centers)
                    bar(centers, [counts_lc(:), counts_lf(:)], 'stacked');
                end
            end
            colormap([0 0.4470 0.7410; 0.8500 0.3250 0.0980]);
            legend(sprintf('LC (n=%d)', numel(vals_lc)), sprintf('LF (n=%d)', numel(vals_lf)), ...
                'Location', 'best');
        else
            histogram(vals_lc, edges, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.6, ...
                'EdgeColor', 'none', 'DisplayName', sprintf('LC (n=%d)', numel(vals_lc))); hold on;
            histogram(vals_lf, edges, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.6, ...
                'EdgeColor', 'none', 'DisplayName', sprintf('LF (n=%d)', numel(vals_lf)));
            hold off;
            legend('Location', 'best', 'FontSize', 9);
        end

        % Show metric name + unit together so the axis is self-describing.
        if isempty(metric_unit)
            xlabel(metric_name);
        else
            xlabel([metric_name ' (' metric_unit ')']);
        end
        ylabel('Count');
        title(metric_name, 'FontSize', 11, 'FontWeight', 'bold');
        grid on;

    elseif strcmpi(plot_type, 'boxplot')
        % --- Boxplot: Compact Group Comparison with Statistical Annotation ---
        % Boxplots show median, IQR, and outliers for each outcome group,
        % making them ideal for multi-panel feature screening where dozens
        % of DWI-derived metrics are compared simultaneously.
        if sum(has_data) > 1
            if exist('OCTAVE_VERSION', 'builtin')
                try
                    boxplot(vals_clean, lf_clean);
                catch
                    % Octave's boxplot can fail with very small groups;
                    % fall back to scatter-style display.
                    plot(lf_clean + 1, vals_clean, 'ko', 'MarkerSize', 8, ...
                        'MarkerFaceColor', [0 0.4470 0.7410]);
                    xlim([0.5, 2.5]);
                    set(gca, 'XTick', [1, 2]);
                    set(gca, 'XTickLabel', {'LC (0)', 'LF (1)'});
                end
            else
                n_lc_bp = sum(lf_clean == 0);
                n_lf_bp = sum(lf_clean == 1);
                n_groups = numel(unique(lf_clean));
                if n_groups == 2
                    boxplot(vals_clean, lf_clean, 'Labels', ...
                        {sprintf('LC (n=%d)', n_lc_bp), sprintf('LF (n=%d)', n_lf_bp)});
                else
                    % Single group present: omit Labels to avoid MATLAB error
                    boxplot(vals_clean, lf_clean);
                end
            end
        else
            % Just plot a single point if only 1 patient
            plot(lf_clean + 1, vals_clean, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.4470 0.7410]);
            xlim([0.5, 2.5]);
            set(gca, 'XTick', [1, 2]);
            set(gca, 'XTickLabel', {'LC', 'LF'});
        end

        if isempty(metric_unit)
            ylabel(metric_name);
        else
            ylabel([metric_name ' (' metric_unit ')']);
        end
        title(metric_name, 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        
        % --- Statistical Annotation ---
        % Annotate with Wilcoxon rank-sum p-value (non-parametric) for
        % consistency with the formal analysis in metrics_stats_comparisons.
        % ANOVA/t-tests assume normality, which is routinely violated by
        % IVIM parameter distributions (especially D* and f, which exhibit
        % strong right skew due to vascular heterogeneity in pancreatic tumors).
        % Require >2 data points and both groups present to avoid degenerate
        % rank-sum tests.
        if sum(has_data) > 2 && numel(unique(lf_clean)) > 1
            p = perform_statistical_test(vals_clean, lf_clean, 'ranksum');
            yl = ylim;
            text(1.5, yl(1) + 0.95*(yl(2) - yl(1)), format_p_value(p), ...
                'HorizontalAlignment', 'center', 'FontSize', 10);
        end
    else
        error('Invalid plot_type. Must be ''histogram'' or ''boxplot''.');
    end
end
