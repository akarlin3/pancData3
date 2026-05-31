function perm_results = optimize_adc_threshold_permutation_test( ...
    thresholds, per_patient_vol_frac, per_patient_lf, n_perm, rng_seed, min_per_group)
% OPTIMIZE_ADC_THRESHOLD_PERMUTATION_TEST  Measure Tactic-3 selection bias.
%
%   perm_results = optimize_adc_threshold_permutation_test( ...
%       thresholds, per_patient_vol_frac, per_patient_lf, ...
%       n_perm, rng_seed, min_per_group)
%
%   Quantifies the selection-on-outcome (multiple-comparisons) bias of
%   Tactic 3 in optimize_adc_threshold.  Tactic 3 sweeps T candidate ADC
%   thresholds, computes a Wilcoxon rank-sum p-value of sub-volume fraction
%   (LC vs LF) at each, and reports the MINIMUM p as `significance_pvalue`.
%   Because the same data both selects the threshold and is tested at it,
%   that minimum is optimistically biased: it is the smallest of up to T
%   correlated tests reported as a single pre-specified test.
%
%   This instrument wraps the ENTIRE Tactic-3 selection in a label-
%   permutation test.  Under each random permutation of the LF labels it
%   re-runs the full "sweep T thresholds, take the minimum p" procedure
%   (via select_significance_threshold -- the identical selection used in
%   production) and builds the null distribution of that *minimum*
%   statistic.  The selection-adjusted p-value is the fraction of
%   permutations whose minimum p is <= the observed minimum p.  The gap
%   between the naive minimum p and the selection-adjusted p IS the bias,
%   quantified.
%
%   The 13 thresholds are nested/monotone (each sub-volume is roughly a
%   superset of the next), so the per-threshold p-values are highly
%   correlated and a naive Bonferroni (xT) would be too conservative.  A
%   permutation test captures the cohort's ACTUAL correlation structure and
%   so neither over- nor under-corrects.
%
%   A per-threshold single-step (Westfall-Young minP / "maxT") adjustment is
%   also returned: for each threshold ti, the family-wise-adjusted p is the
%   fraction of permutations whose minimum p <= the observed p at ti.  This
%   yields the count of thresholds that remain notable after adjusting for
%   the 13-way selection, accounting for the correlation structure.
%
%   Inputs:
%       thresholds           - [1 x T] candidate threshold values
%       per_patient_vol_frac - [nPatients x T] sub-volume fraction per
%                              patient per threshold (NaN where undefined)
%       per_patient_lf       - [nPatients x 1] outcome (0=LC, 1=LF, NaN=unknown)
%       n_perm               - number of label permutations (default 5000)
%       rng_seed             - RNG seed for reproducibility (default 13)
%       min_per_group        - per-group floor passed through (default 3)
%
%   Output struct fields:
%       observed_min_p           - naive minimum p (the biased statistic)
%       observed_thresh          - selected threshold
%       observed_idx             - selected threshold index
%       observed_pvalues         - [1 x T] per-threshold p-values (naive)
%       perm_adjusted_min_p      - selection-adjusted p for the selection
%                                  (fraction of perms with null min-p <=
%                                  observed min-p, +1 plug-in estimator)
%       per_threshold_adjusted_p - [1 x T] Westfall-Young minP-adjusted p
%       n_notable_naive          - #thresholds with naive p < 0.05
%       n_notable_adjusted       - #thresholds with adjusted p < 0.05
%       null_min_p_dist          - [1 x n_perm_valid] null minimum-p draws
%       n_perm_requested         - n_perm
%       n_perm_valid             - permutations that yielded a valid min-p
%       n_labeled_patients       - patients with a defined LF label
%       n_lc_total, n_lf_total   - cohort LC / LF counts
%       rng_seed                 - seed used
%       method                   - description string
%
%   When no threshold is estimable on the observed data (insufficient
%   per-group counts everywhere), all numeric fields are returned NaN and
%   the instrument reports n_perm_valid = 0 rather than fabricating a value.

    if nargin < 4 || isempty(n_perm),       n_perm = 5000;       end
    if nargin < 5 || isempty(rng_seed),      rng_seed = 13;       end
    if nargin < 6 || isempty(min_per_group), min_per_group = 3;   end

    n_thresh = numel(thresholds);

    perm_results = struct();
    perm_results.observed_min_p           = NaN;
    perm_results.observed_thresh          = NaN;
    perm_results.observed_idx             = NaN;
    perm_results.observed_pvalues         = nan(1, n_thresh);
    perm_results.perm_adjusted_min_p      = NaN;
    perm_results.per_threshold_adjusted_p = nan(1, n_thresh);
    perm_results.n_notable_naive          = 0;
    perm_results.n_notable_adjusted       = 0;
    perm_results.null_min_p_dist          = [];
    perm_results.n_perm_requested         = n_perm;
    perm_results.n_perm_valid             = 0;
    perm_results.n_labeled_patients       = 0;
    perm_results.n_lc_total               = 0;
    perm_results.n_lf_total               = 0;
    perm_results.rng_seed                 = rng_seed;
    perm_results.method = ['label-permutation min-p test (Westfall-Young ' ...
        'single-step minP) over Tactic-3 selection'];

    % --- Observed (naive) selection on the real labels ---
    [obs_thresh, obs_idx, obs_min_p, obs_pvalues] = ...
        select_significance_threshold(thresholds, per_patient_vol_frac, ...
                                      per_patient_lf, min_per_group);

    perm_results.observed_min_p   = obs_min_p;
    perm_results.observed_thresh  = obs_thresh;
    perm_results.observed_idx     = obs_idx;
    perm_results.observed_pvalues = obs_pvalues;
    perm_results.n_notable_naive  = sum(obs_pvalues < 0.05);

    % Only patients with a defined LF label participate in the permutation.
    labeled_mask = ~isnan(per_patient_lf);
    labels = per_patient_lf(labeled_mask);
    vf_labeled = per_patient_vol_frac(labeled_mask, :);
    perm_results.n_labeled_patients = numel(labels);
    perm_results.n_lc_total = sum(labels == 0);
    perm_results.n_lf_total = sum(labels == 1);

    % Nothing to permute / nothing estimable -> report NaN honestly.
    if isnan(obs_min_p) || numel(labels) < 2 || ...
            perm_results.n_lc_total < min_per_group || ...
            perm_results.n_lf_total < min_per_group
        return;
    end

    % --- Permutation null of the MINIMUM-p statistic ---
    % Reseed for reproducibility without disturbing the caller's RNG stream.
    try
        prev_rng = rng(rng_seed);
        restore_rng = onCleanup(@() rng(prev_rng));
    catch
        % Octave: rng may be shimmed/limited; fall back to rand('seed', ...).
        try, rand('seed', rng_seed); randn('seed', rng_seed); end %#ok<TRYNC>
    end

    n_lab = numel(labels);
    null_min_p = nan(1, n_perm);
    % For the Westfall-Young per-threshold adjustment we also need, per
    % permutation, the same single null minimum statistic compared against
    % each observed per-threshold p -- so null_min_p is sufficient.
    for b = 1:n_perm
        perm_labels = labels(randperm(n_lab));
        [~, ~, mp] = select_significance_threshold( ...
            thresholds, vf_labeled, perm_labels, min_per_group);
        null_min_p(b) = mp;
    end

    valid = ~isnan(null_min_p);
    null_valid = null_min_p(valid);
    n_valid = numel(null_valid);
    perm_results.null_min_p_dist = null_valid;
    perm_results.n_perm_valid    = n_valid;

    if n_valid == 0
        return;  % degenerate: no permutation was estimable
    end

    % Selection-adjusted p for the selected threshold (+1 plug-in estimator,
    % keeps the p-value valid and strictly > 0).
    perm_results.perm_adjusted_min_p = ...
        (1 + sum(null_valid <= obs_min_p + 1e-15)) / (1 + n_valid);

    % Westfall-Young single-step minP adjustment, per threshold.
    adj = nan(1, n_thresh);
    for ti = 1:n_thresh
        if isnan(obs_pvalues(ti))
            continue;
        end
        adj(ti) = (1 + sum(null_valid <= obs_pvalues(ti) + 1e-15)) / (1 + n_valid);
    end
    perm_results.per_threshold_adjusted_p = adj;
    perm_results.n_notable_adjusted = sum(adj < 0.05);
end
