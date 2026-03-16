function cal = compute_calibration_metrics(risk_scores, outcomes, n_bins, output_folder, dtype_label, fx_label)
% COMPUTE_CALIBRATION_METRICS  Calibration assessment for predictive models.
%
%   Computes calibration plot data, Hosmer-Lemeshow test, Brier score with
%   decomposition, and calibration-in-the-large statistic.
%
% Inputs:
%   risk_scores   - Predicted risk scores (continuous, n x 1)
%   outcomes      - Binary outcomes (0 or 1, n x 1)
%   n_bins        - Number of quantile bins for calibration (default: 5)
%   output_folder - Directory for saving calibration plot (optional)
%   dtype_label   - DWI type label for file naming (optional)
%   fx_label      - Fraction label for file naming (optional)
%
% Outputs:
%   cal           - Struct with calibration metrics

    if nargin < 3 || isempty(n_bins), n_bins = 5; end
    if nargin < 4, output_folder = ''; end
    if nargin < 5, dtype_label = ''; end
    if nargin < 6, fx_label = ''; end

    cal = struct();
    cal.brier_score = NaN;
    cal.hl_chi2 = NaN;
    cal.hl_p = NaN;
    cal.hl_df = NaN;
    cal.calibration_slope = NaN;
    cal.calibration_intercept = NaN;
    cal.mean_predicted = NaN;
    cal.mean_observed = NaN;
    cal.reliability = NaN;
    cal.resolution = NaN;
    cal.uncertainty = NaN;

    % Remove NaN
    valid = isfinite(risk_scores) & isfinite(outcomes);
    scores = risk_scores(valid);
    y = outcomes(valid);
    n = length(scores);

    if n < 10 || length(unique(y)) < 2
        fprintf('  Calibration: insufficient data (n=%d) or single class.\n', n);
        return;
    end

    % Convert scores to probabilities via logistic transform if needed
    if min(scores) < 0 || max(scores) > 1
        probs = 1 ./ (1 + exp(-scores));
    else
        probs = scores;
    end
    probs = max(0.001, min(0.999, probs));

    % --- Brier Score ---
    cal.brier_score = mean((probs - y).^2);

    % Brier decomposition: reliability - resolution + uncertainty
    cal.mean_predicted = mean(probs);
    cal.mean_observed = mean(y);
    cal.uncertainty = mean(y) * (1 - mean(y));

    % --- Hosmer-Lemeshow Test ---
    % Divide into n_bins quantile groups
    [~, sort_idx] = sort(probs);
    bin_size = floor(n / n_bins);
    hl_chi2 = 0;
    bin_means_pred = nan(n_bins, 1);
    bin_means_obs = nan(n_bins, 1);
    bin_sizes = nan(n_bins, 1);

    reliability_sum = 0;
    resolution_sum = 0;

    for g = 1:n_bins
        if g < n_bins
            idx_g = sort_idx((g-1)*bin_size + 1 : g*bin_size);
        else
            idx_g = sort_idx((g-1)*bin_size + 1 : end);
        end
        ng = length(idx_g);
        bin_sizes(g) = ng;
        obs_g = sum(y(idx_g));
        exp_g = sum(probs(idx_g));
        bin_means_pred(g) = mean(probs(idx_g));
        bin_means_obs(g) = mean(y(idx_g));

        % HL chi-squared contribution
        if exp_g > 0 && ng - exp_g > 0
            hl_chi2 = hl_chi2 + (obs_g - exp_g)^2 / (exp_g * (1 - exp_g/ng));
        end

        % Brier decomposition
        reliability_sum = reliability_sum + ng * (bin_means_pred(g) - bin_means_obs(g))^2;
        resolution_sum = resolution_sum + ng * (bin_means_obs(g) - mean(y))^2;
    end

    cal.hl_chi2 = hl_chi2;
    cal.hl_df = max(1, n_bins - 2);
    cal.hl_p = 1 - chi2cdf(hl_chi2, cal.hl_df);
    cal.reliability = reliability_sum / n;
    cal.resolution = resolution_sum / n;

    % --- Calibration slope and intercept ---
    % Logistic recalibration: fit logit(p) = a + b * logit(predicted)
    logit_pred = log(probs ./ (1 - probs));
    if std(logit_pred) > 0
        X_cal = [ones(n, 1), logit_pred];
        try
            b_cal = glmfit(logit_pred, y, 'binomial', 'link', 'logit');
            cal.calibration_intercept = b_cal(1);
            cal.calibration_slope = b_cal(2);
        catch
            % Fallback to simple linear regression
            b_cal = X_cal \ y;
            cal.calibration_intercept = b_cal(1);
            cal.calibration_slope = b_cal(2);
        end
    end

    % --- Print calibration metrics ---
    fprintf('\n  --- Calibration Assessment (%s %s) ---\n', dtype_label, fx_label);
    fprintf('  Brier Score:             %.4f', cal.brier_score);
    if cal.brier_score < 0.25
        fprintf(' (useful discrimination)\n');
    else
        fprintf(' (poor discrimination)\n');
    end
    fprintf('    Reliability:           %.4f\n', cal.reliability);
    fprintf('    Resolution:            %.4f\n', cal.resolution);
    fprintf('    Uncertainty:           %.4f\n', cal.uncertainty);
    fprintf('  Hosmer-Lemeshow:         chi2(%d) = %.2f, p = %.4f', ...
        cal.hl_df, cal.hl_chi2, cal.hl_p);
    if cal.hl_p > 0.05
        fprintf(' (adequate calibration)\n');
    else
        fprintf(' (poor calibration)\n');
    end
    fprintf('  Calibration slope:       %.3f\n', cal.calibration_slope);
    fprintf('  Calibration intercept:   %.3f\n', cal.calibration_intercept);
    fprintf('  Mean predicted:          %.3f\n', cal.mean_predicted);
    fprintf('  Mean observed:           %.3f\n', cal.mean_observed);

    % --- Generate calibration plot ---
    if ~isempty(output_folder)
        try
            fig = figure('Visible', 'off', 'Position', [100 100 600 500]);
            scatter(bin_means_pred, bin_means_obs, 80, 'filled');
            hold on;
            plot([0 1], [0 1], 'k--', 'LineWidth', 1);  % perfect calibration line
            % Add error bars (Wilson CI for proportions)
            for g = 1:n_bins
                if isfinite(bin_means_obs(g)) && isfinite(bin_means_pred(g))
                    p_hat = bin_means_obs(g);
                    n_g = bin_sizes(g);
                    se = sqrt(p_hat * (1 - p_hat) / max(1, n_g));
                    errorbar(bin_means_pred(g), bin_means_obs(g), se, 'b', 'LineWidth', 1);
                end
            end
            xlabel('Mean Predicted Probability');
            ylabel('Observed Frequency');
            title(sprintf('Calibration Plot (%s %s)\nBrier=%.3f, H-L p=%.3f', ...
                dtype_label, fx_label, cal.brier_score, cal.hl_p));
            xlim([0 1]); ylim([0 1]);
            grid on;

            fig_name = sprintf('calibration_plot_%s_%s.png', dtype_label, strrep(fx_label, ' ', '_'));
            saveas(fig, fullfile(output_folder, fig_name));
            close(fig);
            fprintf('  📁 Calibration plot saved: %s\n', fig_name);
        catch ME_cal
            fprintf('  ⚠️  Calibration plot failed: %s\n', ME_cal.message);
        end
    end
end
