function [best_thresh, best_idx, best_p, pvalues, n_lc_at_best, n_lf_at_best] = ...
    select_significance_threshold(thresholds, per_patient_vol_frac, per_patient_lf, min_per_group)
% SELECT_SIGNIFICANCE_THRESHOLD  Tactic-3 min-p threshold selection (pure).
%
%   [best_thresh, best_idx, best_p, pvalues, n_lc, n_lf] = ...
%       select_significance_threshold(thresholds, per_patient_vol_frac, ...
%                                     per_patient_lf, min_per_group)
%
%   Per-threshold Wilcoxon rank-sum p-value of per-patient sub-volume
%   fraction between LC (LF==0) and LF (LF==1) groups, then SELECT the
%   threshold with the smallest valid p-value.  This is the exact selection
%   logic of Tactic 3 in optimize_adc_threshold, factored into a standalone
%   pure function so that the bias-measurement instrument
%   (optimize_adc_threshold_permutation_test) and the nested-CV correction
%   provably reuse the *identical* selection rather than a divergent copy.
%
%   ** This function embodies the selection-on-outcome step. **  The same
%   data both chooses the threshold (argmin over `thresholds`) and is then
%   tested at it, so `best_p` is a minimum over up-to-13 correlated tests
%   reported as if it were a single pre-specified test.  It is therefore
%   optimistically biased; quantify with optimize_adc_threshold_permutation_test
%   and correct with optimize_adc_threshold_nested_cv.
%
%   Returns NaNs when no threshold has both LC>=min_per_group and
%   LF>=min_per_group finite values (Wilcoxon undefined / underpowered
%   below this floor).  Calls the Statistics-Toolbox `ranksum` directly
%   rather than routing through `perform_statistical_test`, because the
%   latter enforces the project's primary-inference sample floor of n>=5
%   per group -- too strict for an exploratory 13-threshold sweep where the
%   goal is ranking, not confirmatory inference.
%
%   Inputs:
%       thresholds           - [1 x T] candidate threshold values
%       per_patient_vol_frac - [nPatients x T] sub-volume fraction per
%                              patient per threshold (NaN where undefined)
%       per_patient_lf       - [nPatients x 1] outcome (0=LC, 1=LF, NaN=unknown)
%       min_per_group        - minimum per-group count to test a threshold
%                              (optional; default 3)
%
%   Outputs:
%       best_thresh   - threshold with smallest valid p (NaN if none valid)
%       best_idx      - its index into `thresholds` (NaN if none valid)
%       best_p        - the minimum valid p-value (the naive, biased statistic)
%       pvalues       - [1 x T] per-threshold p-values (NaN where inestimable)
%       n_lc_at_best  - LC group size at the selected threshold
%       n_lf_at_best  - LF group size at the selected threshold

    if nargin < 4 || isempty(min_per_group)
        min_per_group = 3;
    end

    n_thresh = numel(thresholds);
    pvalues  = nan(1, n_thresh);

    best_thresh   = NaN;
    best_idx      = NaN;
    best_p        = NaN;
    n_lc_at_best  = 0;
    n_lf_at_best  = 0;

    n_lc_arr = zeros(1, n_thresh);
    n_lf_arr = zeros(1, n_thresh);

    for ti = 1:n_thresh
        col = per_patient_vol_frac(:, ti);
        finite_mask = ~isnan(col) & ~isnan(per_patient_lf);
        lc_mask = finite_mask & per_patient_lf == 0;
        lf_mask = finite_mask & per_patient_lf == 1;
        n_lc = sum(lc_mask);
        n_lf = sum(lf_mask);
        n_lc_arr(ti) = n_lc;
        n_lf_arr(ti) = n_lf;
        if n_lc < min_per_group || n_lf < min_per_group
            continue;
        end
        lc_data = col(lc_mask);
        lf_data = col(lf_mask);
        try
            pvalues(ti) = ranksum(lc_data, lf_data);
        catch
            % Test failed (e.g. all-equal data, ranksum package missing
            % in a stripped Octave setup).  Leave NaN at this threshold.
        end
    end

    valid = ~isnan(pvalues);
    if ~any(valid)
        return;
    end
    [best_p, best_idx] = min(pvalues);
    best_thresh   = thresholds(best_idx);
    n_lc_at_best  = n_lc_arr(best_idx);
    n_lf_at_best  = n_lf_arr(best_idx);
end
