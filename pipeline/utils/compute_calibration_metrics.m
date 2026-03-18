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
    cal.hl_warning = '';

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

    % --- Hosmer-Lemeshow Test with frequency validation ---
    % Divide into n_bins quantile groups
    [~, sort_idx] = sort(probs);
    
    % Initial bin allocation
    bin_size = floor(n / n_bins);
    hl_chi2 = 0;
    bin_means_pred = [];
    bin_means_obs = [];
    bin_sizes = [];
    
    % Collect initial bins with frequency validation
    temp_bins = struct('indices', {}, 'obs_freq', {}, 'exp_freq', {});
    
    for g = 1:n_bins
        if g < n_bins
            idx_g = sort_idx((g-1)*bin_size + 1 : g*bin_size);
        else
            idx_g = sort_idx((g-1)*bin_size + 1 : end);
        end
        ng = length(idx_g);
        obs_g = sum(y(idx_g));
        exp_g = sum(probs(idx_g));
        
        temp_bins(g).indices = idx_g;
        temp_bins(g).obs_freq = obs_g;
        temp_bins(g).exp_freq = exp_g;
    end
    
    % Combine bins with low expected frequencies
    min_expected_freq = 5;
    final_bins = struct('indices', {}, 'obs_freq', {}, 'exp_freq', {});
    current_bin = temp_bins(1);
    
    for g = 2:length(temp_bins)
        % Check if current bin has sufficient expected frequency
        if current_bin.exp_freq >= min_expected_freq && ...
           (length(current_bin.indices) - current_bin.exp_freq) >= min_expected_freq
            % Current bin is valid, save it and start new bin
            final_bins(end+1) = current_bin;
            current_bin = temp_bins(g);
        else
            % Combine with next bin
            current_bin.indices = [current_bin.indices; temp_bins(g).indices];
            current_bin.obs_freq = current_bin.obs_freq + temp_bins(g).obs_freq;
            current_bin.exp_freq = current_bin.exp_freq + temp_bins(g).exp_freq;
        end
    end
    
    % Add the last bin
    final_bins(end+1) = current_bin;
    
    % Check if we have enough valid bins
    n_final_bins = length(final_bins);
    if n_final_bins < 2
        fprintf('  ⚠️  H-L test: insufficient bins after frequency validation\n');
        cal.hl_warning = 'Insufficient bins after frequency validation';
        cal.hl_chi2 = NaN;
        cal.hl_p = NaN;
        cal.hl_df = NaN;
    else
        % Check for low expected frequencies in final bins
        low_freq_bins = 0;
        for g = 1:n_final_bins
            exp_g = final_bins(g).exp_freq;
            ng = length(final_bins(g).indices);
            if exp_g < min_expected_freq || (ng - exp_g) < min_expected_freq
                low_freq_bins = low_freq_bins + 1;
            end
        end
        
        if low_freq_bins > 0
            cal.hl_warning = sprintf('Warning: %d bins still have expected frequencies < 5', low_freq_bins);
            fprintf('  ⚠️  H-L test: %s\n', cal.hl_warning);
        end
        
        % Compute H-L test statistic with final bins
        reliability_sum = 0;
        resolution_sum = 0;
        
        bin_means_pred = nan(n_final_bins, 1);
        bin_means_obs = nan(n_final_bins, 1);
        bin_sizes = nan(n_final_bins, 1);
        
        for g = 1:n_final_bins
            idx_g = final_bins(g).indices;
            ng = length(idx_g);
            obs_g = final_bins(g).obs_freq;
            exp_g = final_bins(g).exp_freq;
            
            bin_sizes(g) = ng;
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
        cal.hl_df = max(1, n_final_bins - 2);
        cal.hl_p = 1 - chi2cdf(hl_chi2, cal.hl_df);
        cal.reliability = reliability_sum / n;
        cal.resolution = resolution_sum / n;
        
        if n_final_bins ~= n_bins
            fprintf('  📊 H-L test: Combined bins from %d to %d due to low expected frequencies\n', n_bins, n_final_bins);
        end
    end

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
    
    if isfinite(cal.hl_p)
        fprintf('  Hosmer-Lemeshow:         chi2(%d) = %.2f, p = %.4f', ...
            cal.hl_df, cal.hl_chi2, cal.hl_p);
        if cal.hl_p > 0.05
            fprintf(' (adequate calibration)\n');
        else
            fprintf(' (poor calibration)\n');
        end
        if ~isempty(cal.hl_warning)
            fprintf('    %s\n', cal.hl_warning);
        end
    else
        fprintf('  Hosmer-Lemeshow:         Not computed - %s\n', cal.hl_warning);
    end
    
    fprintf('  Calibration slope:       %.3f\n', cal.calibration_slope);
    fprintf('  Calibration intercept:   %.3f\n', cal.calibration_intercept);
    fprintf('  Mean predicted:          %.3f\n', cal.mean_predicted);
    fprintf('  Mean observed:           %.3f\n', cal.mean_observed);

    % --- Generate calibration plot ---
    if ~isempty(output_folder) && ~isempty(bin_means_pred)
        try
            fig = figure('Visible', 'off', 'Position', [100 100 600 500]);
            scatter(bin_means_pred, bin_means_obs, 80, 'filled');
            hold on;
            plot([0 1], [0 1], 'k--', 'LineWidth', 1);  % perfect calibration line
            % Add error bars (Wilson CI for proportions)
            for g = 1:length(bin_means_obs)
                if isfinite(bin_means_obs(g)) && isfinite(bin_means_pred(g))
                    p_hat = bin_means_obs(g);
                    n_g = bin_sizes(g);
                    se = sqrt(p_hat * (1 - p_hat) / max(1, n_g));
                    errorbar(bin_means_pred(g), bin_means_obs(g), se, 'b', 'LineWidth', 1);
                end
            end
            xlabel('Mean Predicted Probability');
            ylabel('Observed Frequency');
            hl_p_str = isfinite(cal.hl_p) ? sprintf('%.3f', cal.hl_p) : 'N/A';
            title(sprintf('Calibration Plot (%s %s)\nBrier=%.3f, H-L p=%s', ...
                dtype_label, fx_label, cal.brier_score, hl_p_str));
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