function q_values = benjamini_hochberg_fdr(p_values)
% BENJAMINI_HOCHBERG_FDR — Compute FDR-corrected q-values via Benjamini-Hochberg.
%
%   The step-up procedure controls the expected proportion of false
%   discoveries (FDR) at alpha=0.05.  This is less conservative than
%   Bonferroni but appropriate for exploratory biomarker screening where
%   some false positives are acceptable as long as the overall discovery
%   rate is controlled.
%
%   Parameters
%   ----------
%   p_values : numeric vector
%       Raw p-values from multiple hypothesis tests.
%
%   Returns
%   -------
%   q_values : numeric vector (same size as p_values)
%       BH-adjusted q-values.  A test with q < alpha is significant
%       at FDR level alpha.
%
%   References
%   ----------
%   Benjamini Y, Hochberg Y. Controlling the false discovery rate: a
%   practical and powerful approach to multiple testing. J R Stat Soc
%   Series B. 1995;57(1):289-300.
%
%   See also: metrics_stats_comparisons

    n = length(p_values);
    if n == 0
        q_values = [];
        return;
    end

    % Preserve NaN positions — NaN p-values propagate as NaN q-values.
    nan_mask = isnan(p_values);
    p_clean = p_values;
    p_clean(nan_mask) = 1;  % placeholder so sort/rank is valid

    [p_sort, sort_id] = sort(p_clean);
    q_sorted = zeros(n, 1);
    q_sorted(n) = p_sort(n);  % largest p-value: q = p
    for ii = n-1:-1:1
        % Step-up: q(i) = min(q(i+1), p(i) * n/i)
        % The n/i scaling adjusts for rank position — earlier ranks
        % (smaller p) get less aggressive correction.
        q_sorted(ii) = min(q_sorted(ii+1), p_sort(ii) * (n / ii));
    end
    q_sorted = min(q_sorted, 1);  % cap at 1.0

    % Map q-values back to original order
    q_values = zeros(n, 1);
    q_values(sort_id) = q_sorted;

    % Restore NaN for inputs that were NaN
    q_values(nan_mask) = NaN;
end
