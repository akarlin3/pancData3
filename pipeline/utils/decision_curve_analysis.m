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

function tests = test_decision_curve_analysis()
    tests = functiontests(localfunctions);
end

function test_perfect_classifier(testCase)
    % Test with perfect classifier - all positives have high probability, all negatives low
    n_pos = 30;
    n_neg = 70;
    y_true = [ones(n_pos, 1); zeros(n_neg, 1)];
    y_pred_prob = [ones(n_pos, 1) * 0.9; ones(n_neg, 1) * 0.1];
    
    results = decision_curve_analysis(y_true, y_pred_prob, 0:0.1:1);
    
    % Verify output structure
    verifyTrue(testCase, isstruct(results));
    verifyTrue(testCase, all(isfield(results, {'thresholds', 'net_benefit_model', ...
                                               'net_benefit_treat_all', 'net_benefit_treat_none'})));
    
    % Verify dimensions
    verifyEqual(testCase, length(results.thresholds), 11);
    verifyEqual(testCase, length(results.net_benefit_model), 11);
    
    % Perfect classifier should have very high net benefit at optimal threshold
    [max_nb, max_idx] = max(results.net_benefit_model);
    verifyGreaterThan(testCase, max_nb, 0.25);
    
    % Treat-none is always zero
    verifyTrue(testCase, all(results.net_benefit_treat_none == 0));
end

function test_random_classifier(testCase)
    % Test with random classifier
    rng(42);  % For reproducibility
    n = 100;
    y_true = rand(n, 1) > 0.3;  % 30% prevalence
    y_pred_prob = rand(n, 1);   % Random predictions
    
    results = decision_curve_analysis(y_true, y_pred_prob, 0:0.05:1);
    
    % Random classifier should not consistently beat treat-all strategy
    prevalence = mean(y_true);
    model_advantage_count = sum(results.net_benefit_model > results.net_benefit_treat_all);
    
    % Should not dominate across most thresholds
    verifyLessThan(testCase, model_advantage_count, length(results.thresholds) * 0.7);
    
    % Net benefits should be reasonable (not extremely negative or positive)
    verifyTrue(testCase, all(results.net_benefit_model >= -1));
    verifyTrue(testCase, all(results.net_benefit_model <= 1));
end

function test_extreme_thresholds(testCase)
    % Test behavior at extreme threshold values
    n = 50;
    y_true = [ones(20, 1); zeros(30, 1)];
    y_pred_prob = linspace(0.1, 0.9, n)';
    
    % Include extreme thresholds
    thresholds = [0, 0.001, 0.5, 0.999, 1.0];
    results = decision_curve_analysis(y_true, y_pred_prob, thresholds);
    
    % At threshold = 0, everyone is classified as positive
    verifyGreaterThan(testCase, results.net_benefit_model(1), 0);
    
    % At threshold = 1, no one is classified as positive
    verifyEqual(testCase, results.net_benefit_model(end), 0);
    verifyEqual(testCase, results.net_benefit_treat_all(end), 0);
    
    % All values should be finite
    verifyTrue(testCase, all(isfinite(results.net_benefit_model)));
    verifyTrue(testCase, all(isfinite(results.net_benefit_treat_all)));
end

function test_nan_handling(testCase)
    % Test handling of NaN values
    n = 60;
    y_true = [ones(25, 1); zeros(25, 1); NaN(10, 1)];
    y_pred_prob = [linspace(0.6, 0.9, 25)'; linspace(0.1, 0.4, 25)'; NaN(10, 1)];
    
    results = decision_curve_analysis(y_true, y_pred_prob, 0:0.2:1);
    
    % Should handle NaN values without error
    verifyTrue(testCase, all(isfinite(results.net_benefit_model)));
    verifyTrue(testCase, all(isfinite(results.net_benefit_treat_all)));
    verifyTrue(testCase, all(isfinite(results.net_benefit_treat_none)));
end

function test_single_class_data(testCase)
    % Test with all positive or all negative outcomes
    n = 40;
    
    % All positive cases
    y_true_pos = ones(n, 1);
    y_pred_prob_pos = rand(n, 1);
    results_pos = decision_curve_analysis(y_true_pos, y_pred_prob_pos, 0:0.25:1);
    
    % Should not crash and produce reasonable results
    verifyTrue(testCase, all(isfinite(results_pos.net_benefit_model)));
    
    % All negative cases  
    y_true_neg = zeros(n, 1);
    y_pred_prob_neg = rand(n, 1);
    results_neg = decision_curve_analysis(y_true_neg, y_pred_prob_neg, 0:0.25:1);
    
    % Should not crash and produce reasonable results
    verifyTrue(testCase, all(isfinite(results_neg.net_benefit_model)));
    
    % Treat-all with all negatives should be very negative at high thresholds
    verifyLessThan(testCase, results_neg.net_benefit_treat_all(end-1), 0);
end

function test_net_benefit_calculations(testCase)
    % Test specific net benefit calculations with known values
    y_true = [1; 1; 1; 0; 0; 0; 0; 0];  % 3 pos, 5 neg
    y_pred_prob = [0.8; 0.7; 0.4; 0.6; 0.3; 0.2; 0.1; 0.05];
    threshold = 0.5;
    
    results = decision_curve_analysis(y_true, y_pred_prob, threshold);
    
    % Manual calculation for verification
    % At threshold 0.5: predictions >= 0.5 are [0.8, 0.7, 0.6] -> indices [1,2,4]
    % True positives: indices 1,2 (both have y_true=1) -> TP = 2
    % False positives: index 4 (has y_true=0) -> FP = 1
    % Expected net benefit = TP/n - FP/n * pt/(1-pt) = 2/8 - 1/8 * 0.5/0.5 = 0.25 - 0.125 = 0.125
    
    expected_model_nb = 2/8 - 1/8 * (0.5/(1-0.5));
    verifyEqual(testCase, results.net_benefit_model, expected_model_nb, 'AbsTol', 1e-10);
    
    % Treat-all net benefit = prevalence - (1-prevalence) * pt/(1-pt)
    prevalence = 3/8;
    expected_treat_all_nb = prevalence - (1-prevalence) * (0.5/(1-0.5));
    verifyEqual(testCase, results.net_benefit_treat_all, expected_treat_all_nb, 'AbsTol', 1e-10);
end

function test_empty_input(testCase)
    % Test with empty inputs
    verifyError(testCase, @() decision_curve_analysis([], []), 'MATLAB:badsubscript');
end

function test_default_thresholds(testCase)
    % Test that default thresholds are used when not provided
    n = 30;
    y_true = rand(n, 1) > 0.5;
    y_pred_prob = rand(n, 1);
    
    results = decision_curve_analysis(y_true, y_pred_prob);
    
    % Default should be 0:0.01:1 (101 points)
    expected_thresholds = 0:0.01:1;
    verifyEqual(testCase, results.thresholds, expected_thresholds);
    verifyEqual(testCase, length(results.net_benefit_model), 101);
end