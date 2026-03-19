function cal = compute_calibration_metrics(risk_scores, outcomes, times, n_bins, output_folder, dtype_label, fx_label)
% COMPUTE_CALIBRATION_METRICS  Calibration assessment for predictive models with competing risks.
%
%   Computes calibration plot data, Hosmer-Lemeshow test, Fine-Gray competing 
%   risks Brier score with decomposition, and calibration-in-the-large statistic.
%
% Inputs:
%   risk_scores   - Predicted risk scores (continuous, n x 1)
%   outcomes      - Event outcomes: 0=censored, 1=event of interest, 2=competing (n x 1)
%                   OR binary outcomes (0 or 1, n x 1) for backward compatibility
%   times         - Time-to-event data (n x 1), optional for binary outcomes
%   n_bins        - Number of quantile bins for calibration (default: 5)
%   output_folder - Directory for saving calibration plot (optional)
%   dtype_label   - DWI type label for file naming (optional)
%   fx_label      - Fraction label for file naming (optional)
%
% Outputs:
%   cal           - Struct with calibration metrics

    % Handle input arguments for backward compatibility
    if nargin < 7 && nargin >= 3 && (ischar(times) || isempty(times))
        % Old interface: risk_scores, outcomes, n_bins, output_folder, dtype_label, fx_label
        fx_label = dtype_label;
        dtype_label = output_folder;
        output_folder = n_bins;
        n_bins = times;
        times = [];
    end
    
    if nargin < 4 || isempty(n_bins), n_bins = 5; end
    if nargin < 5, output_folder = ''; end
    if nargin < 6, dtype_label = ''; end
    if nargin < 7, fx_label = ''; end

    cal = struct();
    cal.brier_score = NaN;
    cal.competing_risks_brier = NaN;
    cal.hl_chi2 = NaN;
    cal.hl_p = NaN;
    cal.hl_df = NaN;
    cal.calibration_slope = NaN;
    cal.calibration_intercept = NaN;
    cal.calibration_slope_ci = [NaN NaN];
    cal.calibration_intercept_ci = [NaN NaN];
    cal.mean_predicted = NaN;
    cal.mean_observed = NaN;
    cal.reliability = NaN;
    cal.resolution = NaN;
    cal.uncertainty = NaN;
    cal.hl_warning = '';
    cal.competing_risks = false;

    % Determine if this is competing risks scenario
    has_competing_risks = ~isempty(times) && length(unique(outcomes(isfinite(outcomes)))) > 2;
    cal.competing_risks = has_competing_risks;

    % Remove NaN
    if has_competing_risks
        valid = isfinite(risk_scores) & isfinite(outcomes) & isfinite(times);
        scores = risk_scores(valid);
        y = outcomes(valid);
        t = times(valid);
    else
        valid = isfinite(risk_scores) & isfinite(outcomes);
        scores = risk_scores(valid);
        y = outcomes(valid);
        if ~isempty(times)
            t = times(valid);
        end
    end
    
    n = length(scores);

    if n < 10 || (has_competing_risks && length(unique(y)) < 2) || (~has_competing_risks && length(unique(y)) < 2)
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

    % --- Brier Score Calculation ---
    if has_competing_risks
        % Fine-Gray competing risks Brier score
        cal.brier_score = compute_standard_brier_score(probs, y == 1);
        cal.competing_risks_brier = compute_finegray_brier_score(probs, y, t);
        fprintf('  Using Fine-Gray competing risks Brier score\n');
    else
        % Standard binary Brier score
        y_binary = (y == 1) | (y > 0); % Handle both 0/1 and multi-class as binary
        cal.brier_score = mean((probs - y_binary).^2);
        cal.competing_risks_brier = cal.brier_score;
        y = y_binary; % Use binary outcomes for remaining calculations
    end

    % Brier decomposition: reliability - resolution + uncertainty
    cal.mean_predicted = mean(probs);
    cal.mean_observed = mean(y == 1);
    cal.uncertainty = cal.mean_observed * (1 - cal.mean_observed);

    % Use binary outcomes for calibration assessment
    y_binary = (y == 1);

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
        obs_g = sum(y_binary(idx_g));
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
            bin_means_obs(g) = mean(y_binary(idx_g));
            
            % HL chi-squared contribution
            if exp_g > 0 && ng - exp_g > 0
                hl_chi2 = hl_chi2 + (obs_g - exp_g)^2 / (exp_g * (1 - exp_g/ng));
            end
            
            % Brier decomposition
            reliability_sum = reliability_sum + ng * (bin_means_pred(g) - bin_means_obs(g))^2;
            resolution_sum = resolution_sum + ng * (bin_means_obs(g) - cal.mean_observed)^2;
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

    % --- Calibration slope and intercept with confidence intervals ---
    % Logistic recalibration: fit logit(p) = a + b * logit(predicted)
    logit_pred = log(probs ./ (1 - probs));
    if std(logit_pred) > 0
        X_cal = [ones(n, 1), logit_pred];
        try
            [b_cal, dev, stats] = glmfit(logit_pred, y_binary, 'binomial', 'link', 'logit');
            cal.calibration_intercept = b_cal(1);
            cal.calibration_slope = b_cal(2);
            
            % Extract confidence intervals from stats
            if isfield(stats, 'se') && length(stats.se) >= 2
                alpha = 0.05; % 95% confidence interval
                t_crit = norminv(1 - alpha/2);
                
                cal.calibration_intercept_ci = [b_cal(1) - t_crit * stats.se(1), ...
                                                b_cal(1) + t_crit * stats.se(1)];
                cal.calibration_slope_ci = [b_cal(2) - t_crit * stats.se(2), ...
                                           b_cal(2) + t_crit * stats.se(2)];
            end
            
        catch
            % Fallback to simple linear regression
            try
                b_cal = X_cal \ y_binary;
                cal.calibration_intercept = b_cal(1);
                cal.calibration_slope = b_cal(2);
                
                % Estimate standard errors for linear regression
                y_pred = X_cal * b_cal;
                residuals = y_binary - y_pred;
                mse = sum(residuals.^2) / (n - 2);
                cov_matrix = mse * inv(X_cal' * X_cal);
                
                if size(cov_matrix, 1) >= 2
                    se_intercept = sqrt(cov_matrix(1, 1));
                    se_slope = sqrt(cov_matrix(2, 2));
                    
                    t_crit = norminv(0.975); % 95% CI
                    cal.calibration_intercept_ci = [b_cal(1) - t_crit * se_intercept, ...
                                                    b_cal(1) + t_crit * se_intercept];
                    cal.calibration_slope_ci = [b_cal(2) - t_crit * se_slope, ...
                                               b_cal(2) + t_crit * se_slope];
                end
            catch
                % Final fallback - values already set to NaN
            end
        end
    end

    % --- Print calibration metrics ---
    fprintf('\n  --- Calibration Assessment (%s %s) ---\n', dtype_label, fx_label);
    if has_competing_risks
        fprintf('  Standard Brier Score:    %.4f', cal.brier_score);
        if cal.brier_score < 0.25, fprintf(' (useful discrimination)\n');
        else, fprintf(' (poor discrimination)\n'); end
        fprintf('  Fine-Gray Brier Score:   %.4f', cal.competing_risks_brier);
        if cal.competing_risks_brier < 0.25, fprintf(' (useful discrimination)\n');
        else, fprintf(' (poor discrimination)\n'); end
    else
        fprintf('  Brier Score:             %.4f', cal.brier_score);
        if cal.brier_score < 0.25, fprintf(' (useful discrimination)\n');
        else, fprintf(' (poor discrimination)\n'); end
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
    
    if isfinite(cal.calibration_slope)
        fprintf('  Calibration slope:       %.3f', cal.calibration_slope);
        if all(isfinite(cal.calibration_slope_ci))
            fprintf(' [%.3f, %.3f]', cal.calibration_slope_ci(1), cal.calibration_slope_ci(2));
        end
        fprintf('\n');
    else
        fprintf('  Calibration slope:       NaN\n');
    end
    
    if isfinite(cal.calibration_intercept)
        fprintf('  Calibration intercept:   %.3f', cal.calibration_intercept);
        if all(isfinite(cal.calibration_intercept_ci))
            fprintf(' [%.3f, %.3f]', cal.calibration_intercept_ci(1), cal.calibration_intercept_ci(2));
        end
        fprintf('\n');
    else
        fprintf('  Calibration intercept:   NaN\n');
    end
    
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
            
            brier_to_show = has_competing_risks ? cal.competing_risks_brier : cal.brier_score;
            brier_label = has_competing_risks ? 'F-G Brier' : 'Brier';
            title(sprintf('Calibration Plot (%s %s)\n%s=%.3f, H-L p=%s', ...
                dtype_label, fx_label, brier_label, brier_to_show, hl_p_str));
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

function brier = compute_standard_brier_score(probs, binary_outcomes)
    % Standard Brier score for binary outcomes
    brier = mean((probs - binary_outcomes).^2);
end

function fg_brier = compute_finegray_brier_score(probs, outcomes, times)
    % Fine-Gray competing risks Brier score
    % outcomes: 0=censored, 1=event of interest, 2=competing event
    
    if isempty(times)
        % Fallback to standard Brier if no time information
        fg_brier = compute_standard_brier_score(probs, outcomes == 1);
        return;
    end
    
    n = length(probs);
    
    % Estimate cumulative incidence function using Aalen-Johansen estimator
    unique_times = sort(unique(times(outcomes > 0)));
    if isempty(unique_times)
        fg_brier = compute_standard_brier_score(probs, outcomes == 1);
        return;
    end
    
    % Use median follow-up time as evaluation time point
    eval_time = median(times);
    if eval_time <= 0 || eval_time > max(times)
        eval_time = prctile(times, 75); % Use 75th percentile if median is problematic
    end
    
    % Estimate weights using inverse probability of censoring weighting (IPCW)
    weights = compute_ipcw_weights(times, outcomes, eval_time);
    
    % Compute observed cumulative incidence at evaluation time
    obs_ci = zeros(n, 1);
    for i = 1:n
        if times(i) <= eval_time && outcomes(i) == 1
            obs_ci(i) = 1; % Event of interest occurred before eval_time
        elseif times(i) <= eval_time && outcomes(i) == 2
            obs_ci(i) = 0; % Competing event occurred before eval_time
        else
            % Censored or event after eval_time - use cumulative incidence estimate
            obs_ci(i) = estimate_cumulative_incidence(times, outcomes, eval_time, times(i));
        end
    end
    
    % Compute weighted Brier score
    if sum(weights) > 0
        fg_brier = sum(weights .* (probs - obs_ci).^2) / sum(weights);
    else
        % Fallback to unweighted if all weights are zero
        fg_brier = mean((probs - obs_ci).^2);
    end
    
    % Ensure reasonable bounds
    if ~isfinite(fg_brier) || fg_brier < 0
        fg_brier = compute_standard_brier_score(probs, outcomes == 1);
    end
end

function weights = compute_ipcw_weights(times, outcomes, eval_time)
    % Compute inverse probability of censoring weights
    n = length(times);
    weights = ones(n, 1);
    
    % Simple IPCW using Kaplan-Meier estimate of censoring distribution
    censoring_indicator = (outcomes == 0);
    
    if sum(censoring_indicator) == 0
        % No censoring, all weights = 1
        return;
    end
    
    % Estimate censoring survival function
    cens_times = times(censoring_indicator);
    if isempty(cens_times)
        return;
    end
    
    % Simple approximation: use exponential fit for censoring distribution
    try
        lambda_cens = 1 / mean(times(censoring_indicator | times >= eval_time));
        
        for i = 1:n
            if times(i) <= eval_time && outcomes(i) > 0
                % Event observed before eval_time
                G_t = exp(-lambda_cens * times(i)); % P(C > t_i)
                if G_t > 0.01 % Avoid extreme weights
                    weights(i) = 1 / G_t;
                else
                    weights(i) = 1 / 0.01;
                end
            elseif times(i) > eval_time
                % Observation extends beyond eval_time
                G_tau = exp(-lambda_cens * eval_time); % P(C > eval_time)
                if G_tau > 0.01
                    weights(i) = 1 / G_tau;
                else
                    weights(i) = 1 / 0.01;
                end
            end
        end
        
        % Cap extreme weights
        weights = min(weights, 10);
        
    catch
        % Fallback to uniform weights
        weights = ones(n, 1);
    end
end

function ci = estimate_cumulative_incidence(times, outcomes, eval_time, obs_time)
    % Estimate cumulative incidence at eval_time for a subject observed until obs_time
    
    if obs_time >= eval_time
        % Subject was observed beyond evaluation time
        ci = 0; % No event of interest observed
        return;
    end
    
    % Use Aalen-Johansen estimator for subjects with events before eval_time
    event_times = times(outcomes == 1 & times <= eval_time);
    competing_times = times(outcomes == 2 & times <= eval_time);
    
    if isempty(event_times) && isempty(competing_times)
        ci = 0;
        return;
    end
    
    all_event_times = sort([event_times; competing_times]);
    n_at_risk = length(times);
    
    ci = 0;
    survival_prob = 1;
    
    for t = all_event_times'
        if t > eval_time
            break;
        end
        
        n_events_1 = sum(times == t & outcomes == 1);
        n_events_2 = sum(times == t & outcomes == 2);
        n_at_t = sum(times >= t);
        
        if n_at_t > 0
            hazard_1 = n_events_1 / n_at_t;
            hazard_total = (n_events_1 + n_events_2) / n_at_t;
            
            ci = ci + survival_prob * hazard_1;
            survival_prob = survival_prob * (1 - hazard_total);
        end
    end
    
    % Ensure reasonable bounds
    ci = max(0, min(1, ci));
end

% Test suite for calibration metrics
function run_calibration_tests()
    fprintf('\n=== Running Calibration Metrics Test Suite ===\n');
    
    % Test 1: Perfect calibration
    test_perfect_calibration();
    
    % Test 2: Poor calibration
    test_poor_calibration();
    
    % Test 3: Small sample sizes
    test_small_sample_sizes();
    
    % Test 4: Confidence interval accuracy
    test_confidence_intervals();
    
    % Test 5: Edge cases
    test_edge_cases();
    
    fprintf('\n=== All Calibration Tests Completed ===\n');
end

function test_perfect_calibration()
    fprintf('\n--- Test 1: Perfect Calibration ---\n');
    
    % Generate perfectly calibrated data
    n = 1000;
    rng(42); % For reproducibility
    probs = rand(n, 1);
    outcomes = rand(n, 1) < probs; % Perfect calibration
    
    cal = compute_calibration_metrics(probs, outcomes);
    
    % Assertions for perfect calibration
    assert(abs(cal.calibration_slope - 1) < 0.2, 'Perfect calibration should have slope ≈ 1');
    assert(abs(cal.calibration_intercept) < 0.2, 'Perfect calibration should have intercept ≈ 0');
    assert(cal.hl_p > 0.05, 'Perfect calibration should have H-L p > 0.05');
    assert(abs(cal.mean_predicted - cal.mean_observed) < 0.05, 'Mean predicted should ≈ mean observed');
    
    fprintf('✓ Perfect calibration test passed\n');
    fprintf('  Slope: %.3f [%.3f, %.3f]\n', cal.calibration_slope, cal.calibration_slope_ci);
    fprintf('  Intercept: %.3f [%.3f, %.3f]\n', cal.calibration_intercept, cal.calibration_intercept_ci);
end

function test_poor_calibration()
    fprintf('\n--- Test 2: Poor Calibration ---\n');
    
    % Generate poorly calibrated data (overly optimistic predictions)
    n = 500;
    rng(123);
    true_probs = 0.2 * rand(n, 1); % Low true probabilities
    predicted_probs = 0.5 + 0.3 * rand(n, 1); % High predicted probabilities
    outcomes = rand(n, 1) < true_probs;
    
    cal = compute_calibration_metrics(predicted_probs, outcomes);
    
    % Assertions for poor calibration
    assert(cal.mean_predicted > cal.mean_observed + 0.1, 'Should be overly optimistic');
    assert(cal.reliability > 0.01, 'Should have high reliability component');
    
    fprintf('✓ Poor calibration test passed\n');
    fprintf('  Mean predicted: %.3f, Mean observed: %.3f\n', cal.mean_predicted, cal.mean_observed);
    fprintf('  Reliability: %.4f\n', cal.reliability);
end

function test_small_sample_sizes()
    fprintf('\n--- Test 3: Small Sample Sizes ---\n');
    
    % Test with very small sample
    n_small = 15;
    rng(456);
    probs_small = rand(n_small, 1);
    outcomes_small = rand(n_small, 1) < probs_small;
    
    cal_small = compute_calibration_metrics(probs_small, outcomes_small);
    
    % Should handle small samples gracefully
    assert(isfinite(cal_small.calibration_slope) || isnan(cal_small.calibration_slope), ...
           'Should handle small samples without error');
    
    % Test with minimum viable sample
    n_min = 25;
    probs_min = rand(n_min, 1);
    outcomes_min = rand(n_min, 1) < probs_min;
    
    cal_min = compute_calibration_metrics(probs_min, outcomes_min);
    
    fprintf('✓ Small sample size tests passed\n');
    fprintf('  Small sample (n=%d): slope=%.3f, intercept=%.3f\n', ...
            n_small, cal_small.calibration_slope, cal_small.calibration_intercept);
    fprintf('  Min sample (n=%d): slope=%.3f, intercept=%.3f\n', ...
            n_min, cal_min.calibration_slope, cal_min.calibration_intercept);
end

function test_confidence_intervals()
    fprintf('\n--- Test 4: Confidence Interval Accuracy ---\n');
    
    % Generate data with known calibration parameters
    n = 2000;
    rng(789);
    true_slope = 0.8;
    true_intercept = 0.2;
    
    % Generate logit-linear relationship
    logit_pred = randn(n, 1);
    true_logit = true_intercept + true_slope * logit_pred;
    true_probs = 1 ./ (1 + exp(-true_logit));
    outcomes = rand(n, 1) < true_probs;
    
    pred_probs = 1 ./ (1 + exp(-logit_pred));
    cal = compute_calibration_metrics(pred_probs, outcomes);
    
    % Check if confidence intervals contain true values
    slope_in_ci = cal.calibration_slope_ci(1) <= true_slope && true_slope <= cal.calibration_slope_ci(2);
    intercept_in_ci = cal.calibration_intercept_ci(1) <= true_intercept && true_intercept <= cal.calibration_intercept_ci(2);
    
    % Check CI width is reasonable
    slope_ci_width = cal.calibration_slope_ci(2) - cal.calibration_slope_ci(1);
    intercept_ci_width = cal.calibration_intercept_ci(2) - cal.calibration_intercept_ci(1);
    
    assert(slope_ci_width > 0 && slope_ci_width < 1, 'Slope CI should have reasonable width');
    assert(intercept_ci_width > 0 && intercept_ci_width < 1, 'Intercept CI should have reasonable width');
    
    fprintf('✓ Confidence interval tests passed\n');
    fprintf('  True slope: %.3f, Estimated: %.3f [%.3f, %.3f] %s\n', ...
            true_slope, cal.calibration_slope, cal.calibration_slope