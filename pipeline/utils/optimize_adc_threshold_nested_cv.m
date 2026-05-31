function cv_results = optimize_adc_threshold_nested_cv( ...
    thresholds, per_patient_vol_frac, per_patient_lf, min_per_group, ...
    k_folds, n_repeats, rng_seed, id_list)
% OPTIMIZE_ADC_THRESHOLD_NESTED_CV  Honest cross-validated Tactic 3.
%
%   cv_results = optimize_adc_threshold_nested_cv( ...
%       thresholds, per_patient_vol_frac, per_patient_lf, min_per_group, ...
%       k_folds, n_repeats, rng_seed, id_list)
%
%   Corrects the selection-on-outcome bias of Tactic 3 (see CP0 and
%   docs/THRESHOLD_VALIDATION.md) by selecting the threshold INSIDE each
%   training split and evaluating the LC-vs-LF association on OUT-OF-FOLD
%   patients whose threshold was never chosen using their own data.
%
%   Two cross-validation schemes are supported through one entry point:
%     * Leave-one-patient-out (LOO): k_folds >= N (or empty).  Deterministic.
%       NOTE: at N~=42 LOO is a WEAK corrector -- removing one patient rarely
%       changes the argmin threshold, so the out-of-fold p only partially
%       relaxes toward the in-sample min-p.  Reported honestly; the
%       permutation-adjusted p (optimize_adc_threshold_permutation_test) is
%       the primary bias-corrected number.
%     * Repeated stratified grouped K-fold: k_folds < N with n_repeats > 1.
%       Leaving out ~N/k patients perturbs the selection much more, giving a
%       genuinely out-of-fold estimate.  Folds are built with
%       make_grouped_folds (the leakage-safe patient-grouped splitter).
%
%   In both schemes the inner selection reuses select_significance_threshold
%   so it is IDENTICAL to production Tactic 3, and the held-out patients never
%   participate in their own threshold selection.
%
%   Threshold STABILITY -- the distribution of inner-selected thresholds
%   across folds (x repeats) -- is reported as a first-class finding: a cut
%   that swings fold-to-fold is not robustly identifiable at this N.
%
%   Inputs:
%       thresholds           - [1 x T] candidate threshold values
%       per_patient_vol_frac - [N x T] sub-volume fraction per patient/threshold
%       per_patient_lf       - [N x 1] outcome (0=LC, 1=LF, NaN=unknown)
%       min_per_group        - per-group floor (optional; default 3)
%       k_folds              - folds; empty/Inf/>=N -> LOO (optional; default LOO)
%       n_repeats            - repeats for K-fold (optional; default 1; forced
%                              to 1 for LOO since LOO is deterministic)
%       rng_seed             - RNG seed for K-fold partitioning (default 13)
%       id_list              - [N x 1] cell of patient IDs for grouped folds
%                              (optional; defaults to one group per row)
%
%   Output struct fields:
%       thresholds, oof_vol_frac [N x1], oof_label [N x1]
%       n_folds, n_repeats, cv_scheme
%       n_inestimable_folds      - fold(-repeat)s whose inner selection failed
%       oof_pvalue               - HONEST Wilcoxon p of out-of-fold vol_frac
%       oof_auc                  - out-of-fold AUC P[LF vf > LC vf]
%       oof_n_lc, oof_n_lf
%       selected_threshold_counts- [1 x T] histogram of inner selections
%       recommended_thresh/_idx  - modal inner selection (bias-aware single cut)
%       modal_fraction, n_unique_selected
%       threshold_stability_std, threshold_stability_range
%       min_per_group, method

    if nargin < 4 || isempty(min_per_group), min_per_group = 3; end
    if nargin < 5, k_folds = []; end
    if nargin < 6 || isempty(n_repeats), n_repeats = 1; end
    if nargin < 7 || isempty(rng_seed), rng_seed = 13; end
    if nargin < 8, id_list = {}; end

    n_thresh = numel(thresholds);

    % Restrict to patients with a defined LF label.
    labeled_mask = ~isnan(per_patient_lf);
    vf = per_patient_vol_frac(labeled_mask, :);
    lf = per_patient_lf(labeled_mask);
    N = numel(lf);
    if ~isempty(id_list)
        ids = id_list(labeled_mask);
        ids = ids(:);
    else
        ids = arrayfun(@(x) sprintf('row%04d', x), (1:N)', 'UniformOutput', false);
    end

    is_loo = isempty(k_folds) || isinf(k_folds) || k_folds >= N;
    if is_loo
        n_repeats = 1;
        scheme = 'leave-one-patient-out';
    else
        scheme = sprintf('repeated grouped %d-fold x%d', k_folds, n_repeats);
    end

    cv_results = struct();
    cv_results.thresholds                = thresholds;
    cv_results.oof_vol_frac              = nan(N, 1);
    cv_results.oof_idx                   = nan(N, 1);
    cv_results.oof_label                 = lf;
    cv_results.n_folds                   = N;
    cv_results.n_repeats                 = n_repeats;
    cv_results.cv_scheme                 = scheme;
    cv_results.n_inestimable_folds       = 0;
    cv_results.oof_pvalue                = NaN;
    cv_results.oof_auc                   = NaN;
    cv_results.oof_n_lc                  = 0;
    cv_results.oof_n_lf                  = 0;
    cv_results.selected_threshold_counts = zeros(1, n_thresh);
    cv_results.recommended_thresh        = NaN;
    cv_results.recommended_idx           = NaN;
    cv_results.modal_fraction            = NaN;
    cv_results.n_unique_selected         = 0;
    cv_results.threshold_stability_std   = NaN;
    cv_results.threshold_stability_range = NaN;
    cv_results.min_per_group             = min_per_group;
    cv_results.method = ['nested CV; inner Tactic-3 min-p selection on ' ...
        'training split only (' scheme ')'];

    if N < 2
        return;
    end

    % Accumulate per-patient out-of-fold vol_frac across repeats, the
    % per-patient threshold indices applied to them, and the full list of
    % inner-selected threshold indices (folds x repeats).
    oof_accum = cell(N, 1);
    oof_idx_accum = cell(N, 1);
    selected_idx_all = [];
    n_inestimable = 0;

    for r = 1:n_repeats
        if is_loo
            % Each patient is its own fold.
            fold_id = (1:N)';
            n_fold = N;
        else
            % Reseed so each repeat is a different (reproducible) partition.
            try
                prev_rng = rng(rng_seed + r);
                restore_rng = onCleanup(@() rng(prev_rng)); %#ok<NASGU>
            catch
                try, rand('seed', rng_seed + r); end %#ok<TRYNC>
            end
            fold_id = make_grouped_folds(ids, lf, k_folds);
            n_fold = max(fold_id);
            clear restore_rng;  % restore RNG before the heavy inner loop
        end

        for f = 1:n_fold
            te_mask = (fold_id == f);
            tr_mask = ~te_mask;
            if ~any(te_mask) || ~any(tr_mask)
                continue;
            end
            [sel_thresh, sel_idx] = select_significance_threshold( ...
                thresholds, vf(tr_mask, :), lf(tr_mask), min_per_group);
            if isnan(sel_idx)
                n_inestimable = n_inestimable + 1;
                continue;
            end
            selected_idx_all(end+1) = sel_idx; %#ok<AGROW>
            % Apply the training-selected threshold to each held-out patient.
            te_rows = find(te_mask);
            for tt = 1:numel(te_rows)
                row = te_rows(tt);
                oof_accum{row}(end+1) = vf(row, sel_idx); %#ok<AGROW>
                oof_idx_accum{row}(end+1) = sel_idx; %#ok<AGROW>
            end
        end
    end

    cv_results.n_inestimable_folds = n_inestimable;

    % Average each patient's out-of-fold vol_frac across repeats; record the
    % modal threshold index applied to each patient (for LOO this is exactly
    % that patient's single fold selection -- used by the leakage test).
    oof_vol_frac = nan(N, 1);
    oof_idx = nan(N, 1);
    for i = 1:N
        a = oof_accum{i};
        a = a(~isnan(a));
        if ~isempty(a)
            oof_vol_frac(i) = mean(a);
        end
        ix = oof_idx_accum{i};
        if ~isempty(ix)
            oof_idx(i) = mode(ix);
        end
    end
    cv_results.oof_vol_frac = oof_vol_frac;
    cv_results.oof_idx = oof_idx;

    % --- Threshold stability across folds (x repeats) ---
    if ~isempty(selected_idx_all)
        counts = zeros(1, n_thresh);
        for k = 1:n_thresh
            counts(k) = sum(selected_idx_all == k);
        end
        cv_results.selected_threshold_counts = counts;
        [modal_count, modal_idx] = max(counts);
        cv_results.recommended_idx    = modal_idx;
        cv_results.recommended_thresh = thresholds(modal_idx);
        cv_results.modal_fraction     = modal_count / numel(selected_idx_all);
        cv_results.n_unique_selected  = sum(counts > 0);
        sel_thresh_all = thresholds(selected_idx_all);
        cv_results.threshold_stability_std   = std(sel_thresh_all);
        cv_results.threshold_stability_range = max(sel_thresh_all) - min(sel_thresh_all);
    end

    % --- Honest out-of-fold LC-vs-LF association ---
    valid = ~isnan(oof_vol_frac);
    lc_vals = oof_vol_frac(valid & lf == 0);
    lf_vals = oof_vol_frac(valid & lf == 1);
    cv_results.oof_n_lc = numel(lc_vals);
    cv_results.oof_n_lf = numel(lf_vals);
    if numel(lc_vals) >= min_per_group && numel(lf_vals) >= min_per_group
        try
            cv_results.oof_pvalue = ranksum(lc_vals, lf_vals);
        catch
            cv_results.oof_pvalue = NaN;
        end
        cv_results.oof_auc = local_auc(lf_vals, lc_vals);
    end
end


function auc = local_auc(pos_vals, neg_vals)
% LOCAL_AUC  AUC = P(random positive > random negative), ties = 0.5.
% Rank-based (Mann-Whitney U); no Statistics-Toolbox dependency so it runs
% identically under Octave.  NaN if either group is empty.
    n1 = numel(pos_vals);
    n0 = numel(neg_vals);
    if n1 == 0 || n0 == 0
        auc = NaN;
        return;
    end
    all_vals = [pos_vals(:); neg_vals(:)];
    r = local_tiedrank(all_vals);
    u = sum(r(1:n1)) - n1 * (n1 + 1) / 2;
    auc = u / (n1 * n0);
end


function r = local_tiedrank(x)
% LOCAL_TIEDRANK  Average ranks with tie handling (tiedrank replacement).
    x = x(:);
    n = numel(x);
    [~, order] = sort(x);
    r = zeros(n, 1);
    i = 1;
    while i <= n
        j = i;
        while j < n && x(order(j + 1)) == x(order(i))
            j = j + 1;
        end
        avg = (i + j) / 2;
        r(order(i:j)) = avg;
        i = j + 1;
    end
end
