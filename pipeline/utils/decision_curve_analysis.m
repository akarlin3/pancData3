function results = decision_curve_analysis(y_true, y_pred_prob, thresholds, output_folder, dtype_label, fx_label)
% DECISION_CURVE_ANALYSIS  Compute net benefit across threshold probabilities.
%
%   Performs decision curve analysis for binary outcome prediction.
%   Computes net benefit for the model, treat-all, and treat-none
%   strategies across a range of threshold probabilities.
%
% Inputs:
%   y_true       - [n x 1] binary outcome (0/1)
%   y_pred_prob  - [n x 1] predicted probabilities
%   thresholds   - [1 x m] threshold probabilities (default: 0:0.01:1)
%   output_folder - Where to save the DCA plot (optional)
%   dtype_label   - DWI type label for figure naming (optional)
%   fx_label      - Fraction label for figure naming (optional)
%
% Outputs:
%   results      - Struct with fields:
%                    thresholds           - threshold grid
%                    net_benefit_model    - net benefit for the model
%                    net_benefit_treat_all - net benefit for treat-all
%                    net_benefit_treat_none - net benefit for treat-none (zeros)

    if nargin < 3 || isempty(thresholds)
        thresholds = 0:0.01:1;
    end
    if nargin < 4, output_folder = ''; end
    if nargin < 5, dtype_label = ''; end
    if nargin < 6, fx_label = ''; end

    % Remove NaN observations
    valid = ~isnan(y_true) & ~isnan(y_pred_prob);
    y_true = y_true(valid);
    y_pred_prob = y_pred_prob(valid);
    n = numel(y_true);

    n_thresh = numel(thresholds);
    net_benefit_model = nan(1, n_thresh);
    net_benefit_treat_all = nan(1, n_thresh);
    net_benefit_treat_none = zeros(1, n_thresh);

    prevalence = mean(y_true);

    for i = 1:n_thresh
        pt = thresholds(i);

        % Treat-all strategy
        if pt < 1
            net_benefit_treat_all(i) = prevalence - (1 - prevalence) * pt / (1 - pt);
        else
            net_benefit_treat_all(i) = 0;
        end

        % Model strategy
        pred_pos = y_pred_prob >= pt;
        tp = sum(pred_pos & y_true == 1);
        fp = sum(pred_pos & y_true == 0);
        if pt < 1
            net_benefit_model(i) = tp / n - fp / n * pt / (1 - pt);
        else
            net_benefit_model(i) = 0;
        end
    end

    % Find range where model > both defaults
    model_better = (net_benefit_model > net_benefit_treat_all) & ...
                   (net_benefit_model > net_benefit_treat_none);
    useful_range = thresholds(model_better);

    results = struct();
    results.thresholds = thresholds;
    results.net_benefit_model = net_benefit_model;
    results.net_benefit_treat_all = net_benefit_treat_all;
    results.net_benefit_treat_none = net_benefit_treat_none;

    % Generate DCA plot
    if ~isempty(output_folder)
        try
            fig = figure('Visible', 'off', 'Position', [100 100 700 500]);
            plot(thresholds, net_benefit_model, 'b-', 'LineWidth', 2); hold on;
            plot(thresholds, net_benefit_treat_all, 'r--', 'LineWidth', 1.5);
            plot(thresholds, net_benefit_treat_none, 'k-', 'LineWidth', 1);
            xlabel('Threshold Probability');
            ylabel('Net Benefit');
            title_str = 'Decision Curve Analysis';
            if ~isempty(dtype_label)
                title_str = sprintf('%s (%s', title_str, dtype_label);
                if ~isempty(fx_label)
                    title_str = sprintf('%s, %s)', title_str, fx_label);
                else
                    title_str = [title_str ')'];
                end
            end
            title(title_str);
            legend({'Model', 'Treat All', 'Treat None'}, 'Location', 'best');
            grid on;
            xlim([0 1]);

            % Annotate useful range
            if ~isempty(useful_range)
                text(0.5, max(net_benefit_model) * 0.9, ...
                    sprintf('Model > defaults: [%.2f, %.2f]', ...
                    min(useful_range), max(useful_range)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 9);
            end

            fig_name = 'dca';
            if ~isempty(dtype_label)
                fig_name = sprintf('%s_%s', fig_name, dtype_label);
            end
            if ~isempty(fx_label)
                fig_name = sprintf('%s_%s', fig_name, fx_label);
            end
            fig_path = fullfile(output_folder, [fig_name '.png']);
            saveas(fig, fig_path);
            close(fig);
            fprintf('  📁 DCA plot saved: %s\n', fig_path);
        catch ME_fig
            fprintf('  ⚠️  DCA plot failed: %s\n', ME_fig.message);
        end
    end
end
