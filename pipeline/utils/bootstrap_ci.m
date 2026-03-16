function [ci_lo, ci_hi, boot_dist] = bootstrap_ci(data, metric_fn, n_boot, alpha)
% BOOTSTRAP_CI  Bias-corrected and accelerated (BCa) bootstrap confidence intervals.
%
%   Computes bootstrap percentile confidence intervals for an arbitrary
%   scalar metric function using the BCa method, which corrects for both
%   bias and skewness in the bootstrap distribution.
%
% Inputs:
%   data       - Input data (vector or matrix) passed to metric_fn
%   metric_fn  - Function handle: scalar = metric_fn(data)
%   n_boot     - Number of bootstrap resamples (default: 2000)
%   alpha      - Significance level (default: 0.05 for 95% CI)
%
% Outputs:
%   ci_lo      - Lower confidence bound
%   ci_hi      - Upper confidence bound
%   boot_dist  - Vector of n_boot bootstrap metric values

    if nargin < 3 || isempty(n_boot), n_boot = 2000; end
    if nargin < 4 || isempty(alpha), alpha = 0.05; end

    n = size(data, 1);
    if n == 0
        ci_lo = NaN; ci_hi = NaN; boot_dist = [];
        return;
    end

    % Original statistic
    theta_hat = metric_fn(data);

    % Bootstrap resampling
    boot_dist = nan(n_boot, 1);
    for b = 1:n_boot
        idx = randi(n, n, 1);
        try
            boot_dist(b) = metric_fn(data(idx, :));
        catch
            boot_dist(b) = NaN;
        end
    end

    % Remove NaN resamples
    boot_valid = boot_dist(~isnan(boot_dist));
    if isempty(boot_valid) || all(boot_valid == boot_valid(1))
        % Degenerate case: all identical or all NaN
        ci_lo = theta_hat;
        ci_hi = theta_hat;
        return;
    end

    % --- BCa correction ---

    % Bias correction factor z0
    z0 = norminv(mean(boot_valid < theta_hat));
    if ~isfinite(z0), z0 = 0; end

    % Acceleration factor via jackknife
    jack_vals = nan(n, 1);
    for i = 1:n
        idx_jack = [1:i-1, i+1:n]';
        try
            jack_vals(i) = metric_fn(data(idx_jack, :));
        catch
            jack_vals(i) = NaN;
        end
    end
    jack_valid = jack_vals(~isnan(jack_vals));
    if length(jack_valid) < 3
        a_hat = 0;
    else
        jack_mean = mean(jack_valid);
        num = sum((jack_mean - jack_valid).^3);
        den = 6 * (sum((jack_mean - jack_valid).^2))^1.5;
        if den == 0
            a_hat = 0;
        else
            a_hat = num / den;
        end
    end

    % BCa adjusted percentiles
    z_alpha_lo = norminv(alpha / 2);
    z_alpha_hi = norminv(1 - alpha / 2);

    % Lower percentile
    num_lo = z0 + z_alpha_lo;
    adj_lo = normcdf(z0 + num_lo / (1 - a_hat * num_lo));
    % Upper percentile
    num_hi = z0 + z_alpha_hi;
    adj_hi = normcdf(z0 + num_hi / (1 - a_hat * num_hi));

    % Clamp to [0.005, 0.995] to avoid extreme quantiles
    adj_lo = max(0.005, min(0.995, adj_lo));
    adj_hi = max(0.005, min(0.995, adj_hi));

    boot_sorted = sort(boot_valid);
    n_valid = length(boot_sorted);
    ci_lo = boot_sorted(max(1, round(adj_lo * n_valid)));
    ci_hi = boot_sorted(min(n_valid, round(adj_hi * n_valid)));
end
