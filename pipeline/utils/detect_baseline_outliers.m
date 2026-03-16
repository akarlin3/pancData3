function [is_outlier, n_total_outliers] = detect_baseline_outliers(baseline_metrics, metric_names, lf)
% DETECT_BASELINE_OUTLIERS — Outcome-blinded IQR outlier detection.
%
% Applies 3x IQR fencing on each baseline metric independently. The detection
% is outcome-BLINDED: thresholds are derived solely from the metric distribution
% without reference to the lf outcome variable. However, the outcome distribution
% of flagged outliers is logged for informative censoring assessment.
%
% Inputs:
%   baseline_metrics - cell array of column vectors (one per metric)
%   metric_names     - cell array of metric name strings (same length)
%   lf               - local failure status vector (0=LC, 1=LF, 2=CR)
%
% Outputs:
%   is_outlier        - logical vector (true = outlier in any metric)
%   n_total_outliers  - scalar count of total outlier patients

is_outlier = false(size(lf));
n_baseline_metrics = numel(baseline_metrics);
for metric_idx = 1:n_baseline_metrics
    text_progress_bar(metric_idx, n_baseline_metrics, 'Detecting outliers');
    col = baseline_metrics{metric_idx};
    col_clean = col(~isnan(col));
    if numel(col_clean) < 3, continue; end
    med_val = median(col_clean);
    iqr_val = iqr(col_clean);
    if iqr_val == 0, continue; end
    lower_fence = med_val - 3 * iqr_val;
    upper_fence = med_val + 3 * iqr_val;
    outlier_flags = (col < lower_fence | col > upper_fence) & ~isnan(col);
    if any(outlier_flags)
        n_out = sum(outlier_flags);
        n_out_lf = sum(outlier_flags & (lf == 1));
        n_out_lc = sum(outlier_flags & (lf == 0));
        n_out_cr = sum(outlier_flags & (lf == 2));
        fprintf('  Outlier flag (%s): %d flagged (LF=%d, LC=%d, CR=%d)\n', ...
            metric_names{metric_idx}, n_out, n_out_lf, n_out_lc, n_out_cr);
    end
    is_outlier = is_outlier | outlier_flags;
end
n_total_outliers = sum(is_outlier);
if n_total_outliers > 0
    fprintf('  Total outliers removed: %d / %d (%.1f%%)\n', ...
        n_total_outliers, numel(lf), 100*n_total_outliers/numel(lf));
end
end
