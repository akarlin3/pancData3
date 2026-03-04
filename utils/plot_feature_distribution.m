function plot_feature_distribution(vals, lf_group, metric_name, metric_unit, plot_type)
% PLOT_FEATURE_DISTRIBUTION Plots feature distribution as histogram or boxplot.
%
% Syntax:
%   plot_feature_distribution(vals, lf_group, metric_name, metric_unit, plot_type)
%
% Inputs:
%   vals        - [N x 1] numeric vector of feature values.
%   lf_group    - [N x 1] numeric vector indicating clinical outcome (0 for LC, 1 for LF).
%   metric_name - Character vector or string for the metric name (used in title).
%   metric_unit - Character vector or string for the metric unit (used in ylabel or xlabel).
%   plot_type   - 'histogram' or 'boxplot' indicating the type of plot.
%
% Outputs:
%   None. Generates a histogram or boxplot visualization.
%

    % Remove NaN entries for accurate processing
    has_data = ~isnan(vals);
    vals_clean = vals(has_data);
    lf_clean = lf_group(has_data);

    if strcmpi(plot_type, 'histogram')
        % Split values by clinical outcome
        vals_lc = vals_clean(lf_clean == 0);   % Local Control patients
        vals_lf = vals_clean(lf_clean == 1);   % Local Failure patients

        % Create 15 equally-spaced bins spanning the combined value range
        if exist('OCTAVE_VERSION', 'builtin')
            min_val = min(vals_clean(~isnan(vals_clean)));
            max_val = max(vals_clean(~isnan(vals_clean)));
        else
            min_val = min(vals_clean, [], 'omitnan');
            max_val = max(vals_clean, [], 'omitnan');
        end

        if isempty(vals_clean)
            edges = linspace(0, 1, 16);
        elseif min_val == max_val
            edges = linspace(min_val - 0.5, max_val + 0.5, 16);
        else
            edges = linspace(min_val, max_val, 16);
        end

        % Overlay semi-transparent histograms for each outcome group
        if exist('OCTAVE_VERSION', 'builtin')
            if length(unique(edges)) == length(edges) && length(edges) > 1
                [counts_lc, centers] = hist(vals_lc, edges);
                [counts_lf, ~] = hist(vals_lf, edges);
                % handle edge case where centers is not unique
                if length(unique(centers)) == length(centers)
                    bar(centers, [counts_lc(:), counts_lf(:)], 'stacked');
                end
            end
            colormap([0 0.4470 0.7410; 0.8500 0.3250 0.0980]);
            legend('Local Control', 'Local Failure', 'Location', 'best');
        else
            histogram(vals_lc, edges, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.6, ...
                'EdgeColor', 'none', 'DisplayName', 'Local Control'); hold on;
            histogram(vals_lf, edges, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.6, ...
                'EdgeColor', 'none', 'DisplayName', 'Local Failure');
            hold off;
            legend('Location', 'best', 'FontSize', 8);
        end

        xlabel(metric_unit); ylabel('Count');
        title(metric_name, 'FontSize', 11, 'FontWeight', 'bold');
        grid on;

    elseif strcmpi(plot_type, 'boxplot')
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
                boxplot(vals_clean, lf_clean, 'Labels', {'LC (0)', 'LF (1)'});
            end
        else
            % Just plot a single point if only 1 patient
            plot(lf_clean + 1, vals_clean, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.4470 0.7410]);
            xlim([0.5, 2.5]);
            set(gca, 'XTick', [1, 2]);
            set(gca, 'XTickLabel', {'LC (0)', 'LF (1)'});
        end
        
        ylabel(metric_unit);
        title(metric_name, 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        
        % Annotate with Wilcoxon rank-sum p-value (non-parametric) for
        % consistency with metrics_stats_comparisons formal analysis.
        % ANOVA assumes normality which is violated by skewed IVIM distributions.
        if sum(has_data) > 2 && numel(unique(lf_clean)) > 1
            p = perform_statistical_test(vals_clean, lf_clean, 'ranksum');
            yl = ylim;
            text(1.5, yl(2)*0.95, format_p_value(p), ...
                'HorizontalAlignment', 'center', 'FontSize', 10);
        end
    else
        error('Invalid plot_type. Must be ''histogram'' or ''boxplot''.');
    end
end
