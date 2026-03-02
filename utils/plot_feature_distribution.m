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
        edges = linspace(min(vals_clean), max(vals_clean), 16);

        % Overlay semi-transparent histograms for each outcome group
        if exist('OCTAVE_VERSION', 'builtin')
            if length(unique(edges)) == length(edges) && length(edges) > 1
                [counts_lc, centers] = hist(vals_lc, edges);
                [counts_lf, ~] = hist(vals_lf, edges);
                % handle edge case where centers is not unique
                if length(unique(centers)) == length(centers)
                    bar(centers, [counts_lc(:), counts_lf(:)], 'stacked', 'EdgeColor', 'none');
                end
            end
            colormap([0.2 0.4 0.8; 0.8 0.2 0.2]);
            legend('Local Control', 'Local Failure', 'Location', 'best');
        else
            histogram(vals_lc, edges, 'FaceColor', [0.2 0.4 0.8], 'FaceAlpha', 0.6, ...
                'EdgeColor', 'none', 'DisplayName', 'Local Control'); hold on;
            histogram(vals_lf, edges, 'FaceColor', [0.8 0.2 0.2], 'FaceAlpha', 0.6, ...
                'EdgeColor', 'none', 'DisplayName', 'Local Failure');
            hold off;
            legend('Location', 'best', 'FontSize', 8);
        end

        xlabel(metric_unit); ylabel('Count');
        title(metric_name, 'FontSize', 11, 'FontWeight', 'bold');
        grid on;

    elseif strcmpi(plot_type, 'boxplot')
        boxplot(vals_clean, lf_clean, 'Labels', {'LC (0)', 'LF (1)'});

        ylabel(metric_unit);
        title(metric_name, 'FontSize', 11, 'FontWeight', 'bold');
        grid on;

        % Annotate with a one-way ANOVA p-value comparing LC vs LF
        if sum(has_data) > 2 && numel(unique(lf_clean)) > 1
            p = anova1(vals_clean, lf_clean, 'off');
            yl = ylim;
            text(1.5, yl(2)*0.95, sprintf('p = %.3f', p), ...
                'HorizontalAlignment', 'center', 'FontSize', 10);
        end
    else
        error('Invalid plot_type. Must be ''histogram'' or ''boxplot''.');
    end
end
