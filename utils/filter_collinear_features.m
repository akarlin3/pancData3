function [keep_idx] = filter_collinear_features(X, y, frac_vec)
% FILTER_COLLINEAR_FEATURES Filters collinear features using Spearman |rho| > 0.8.
%
%   [keep_idx] = filter_collinear_features(X, y)
%   [keep_idx] = filter_collinear_features(X, y, frac_vec)
%
%   Inputs:
%       X        - Feature matrix [n_observations x n_features]
%       y        - Target vector [n_observations x 1] (used to break ties based on predictive power)
%       frac_vec - (Optional) Time-stratified fraction vector [n_observations x 1] to limit
%                  collinearity calculation to the baseline (fraction 1) rows only.
%
%   Outputs:
%       keep_idx - Array of column indices for features that should be retained (not pruned)
%
%   Analytical Rationale — Why Collinearity Filtering is Necessary:
%   ---------------------------------------------------------------
%   IVIM model parameters are inherently correlated: ADC and D measure
%   overlapping aspects of tissue diffusivity (ADC includes a perfusion
%   contribution that D does not), and summary statistics derived from the
%   same parameter map (e.g., mean_ADC, median_ADC, p25_ADC) are highly
%   correlated by construction.  Including redundant features in a Cox
%   regression or logistic model inflates coefficient variance (variance
%   inflation factor), destabilizes parameter estimates, and impairs
%   interpretability — a coefficient for mean_ADC becomes uninterpretable
%   when median_ADC is also in the model.
%
%   The 0.8 Spearman threshold is a standard choice that removes severe
%   redundancy while retaining features with complementary information.
%   Spearman (rank) correlation is preferred over Pearson because IVIM
%   parameters often have skewed, heavy-tailed distributions (especially
%   D* and f in heterogeneous tumors) where rank-based measures are more
%   robust.
%
%   When two features are highly correlated (|r| > 0.8), the one with the higher
%   univariate Wilcoxon rank-sum p-value (less significant) is dropped.
%
%   Strictly Non-Leaking: This function should be called on training data (X_train, y_train)
%   to obtain keep_idx. The same keep_idx MUST then be applied to the test data
%   without further correlation assessment on the test set.
%
%   Time-Stratified Mode (frac_vec provided):
%   When frac_vec is supplied, the correlation matrix and Wilcoxon significance
%   tests are computed exclusively on the Fraction 1 (baseline) rows
%   (frac_vec == 1), preventing late-stage radiation-response collinearity from
%   incorrectly pruning baseline features.  Radiation can induce artificial
%   correlations between parameters that are independent at baseline (e.g.,
%   ADC and f may both increase post-RT due to cell death and vascular
%   disruption, creating a treatment-induced correlation that is absent in
%   untreated tissue).  The resulting Boolean pruning mask is then applied
%   uniformly across ALL longitudinal rows by the caller.
%   If no Fraction 1 rows are present, the function falls back to using the
%   full matrix for correlation computation.

    if isempty(X)
        keep_idx = [];
        return;
    end

    % --- Time-stratified subset: use only Fraction 1 (baseline) rows ---
    if nargin >= 3 && ~isempty(frac_vec)
        baseline_mask = (frac_vec == 1);
        if any(baseline_mask)
            X_corr = X(baseline_mask, :);
            y_corr = y(baseline_mask);
        else
            % No Fraction 1 rows present: fall back to full matrix
            warning('filter_collinear_features:noBaseline', ...
                'No Fraction 1 (baseline) rows found. Using full matrix for collinearity detection (correlation structure may differ from baseline).');
            X_corr = X;
            y_corr = y;
        end
    else
        X_corr = X;
        y_corr = y;
    end

    % Spearman (rank) correlation is more robust than Pearson for the
    % skewed, heavy-tailed distributions typical of IVIM parameters.
    % 'Rows','pairwise' handles NaN values by computing each pairwise
    % correlation using only rows where both features are non-NaN,
    % maximizing data usage in the presence of missing scans.
    R = corr(X_corr, 'Type', 'Spearman', 'Rows', 'pairwise');
    n_feats = size(X, 2);
    drop_idx = false(1, n_feats);

    % Initialize C-indices (AUC values) with a sentinel value (-1) to
    % indicate "not yet computed".  Lazy evaluation avoids computing AUC
    % for features that are never involved in a high-correlation pair,
    % which can be many in a large feature set (e.g., 50+ histogram bins
    % and summary statistics per IVIM parameter).
    c_indices = -1 * ones(1, n_feats);
    
    % Find highly collinear pairs and sort by descending |rho| so that the
    % strongest collinearities are resolved first.  This makes the greedy
    % pruning deterministic regardless of column ordering.  Using only the
    % lower triangle of the correlation matrix avoids processing each pair
    % twice (R is symmetric) and excludes the diagonal (self-correlation = 1).
    R_lower = abs(tril(R, -1));
    [row_idx, col_idx] = find(R_lower > 0.8);
    pair_rho = zeros(length(row_idx), 1);
    for pi = 1:length(row_idx)
        pair_rho(pi) = R_lower(row_idx(pi), col_idx(pi));
    end
    [~, sort_order] = sort(pair_rho, 'descend');
    row_idx = row_idx(sort_order);
    col_idx = col_idx(sort_order);

    % Process pairs sequentially in order of decreasing |rho|, dropping
    % the weaker feature in each pair.  This greedy strategy resolves the
    % most severe collinearities first; once a feature is dropped, any
    % remaining pairs involving it are skipped (the collinearity is already
    % resolved).  The AUC-based tie-breaking ensures we retain the feature
    % with stronger univariate discriminative power for the clinical
    % outcome, preserving the most informative signal for downstream models.
    for idx = 1:length(row_idx)
        fi = row_idx(idx);
        fj = col_idx(idx);

        if drop_idx(fi) || drop_idx(fj)
            continue; % One of them is already dropped — pair resolved
        end

        % Calculate AUC lazily if not already done.  AUC (area under the
        % ROC curve via Mann-Whitney U) measures each feature's univariate
        % ability to discriminate local failures from censored patients.
        % It is preferred over p-values for tie-breaking because AUC is
        % scale-invariant and directly interpretable as a probability of
        % correct ranking.
        if c_indices(fi) == -1
            c_indices(fi) = compute_auc(X_corr(:, fi), y_corr);
        end
        if c_indices(fj) == -1
            c_indices(fj) = compute_auc(X_corr(:, fj), y_corr);
        end

        % Drop the feature that has weaker predictive power (lower AUC).
        % When AUCs are equal, feature j is dropped (arbitrary but
        % deterministic tie-breaking).
        if c_indices(fj) <= c_indices(fi)
            drop_idx(fj) = true;
        else
            drop_idx(fi) = true;
        end
    end
    keep_idx = find(~drop_idx);
end

function auc_val = compute_auc(feat_col, y_col)
    % Helper function to compute AUC (event y==1 vs censored y==0).
    %
    % AUC is equivalent to the probability that a randomly selected event
    % patient has a higher (or lower) feature value than a randomly
    % selected non-event patient.  AUC = 0.5 means no discriminative
    % power (random guessing); AUC = 1.0 means perfect separation.
    %
    % Exclude competing-risk rows (y==2) so they do not bias the
    % rank-sum denominator.  In the Cause-Specific Hazard framework,
    % competing events are censored for the cause of interest, so they
    % should not influence feature selection for local failure prediction.
    binary_mask = (y_col == 0 | y_col == 1);
    valid_idx = binary_mask & ~isnan(feat_col) & ~isnan(y_col);
    if sum(valid_idx) >= 10
        y_valid = y_col(valid_idx);
        feat_valid = feat_col(valid_idx);

        n1 = sum(y_valid == 1);  % number of events (local failures)
        n0 = sum(y_valid == 0);  % number of non-events (censored)

        if n1 == 0 || n0 == 0
            % Cannot compute AUC without both classes; return 0.5 (no
            % discrimination) so this feature is treated as uninformative
            % for tie-breaking purposes.
            auc_val = 0.5;
            return;
        end

        % Fast AUC computation via Mann-Whitney U statistic.  The U
        % statistic equals the sum of ranks of the event class minus the
        % minimum possible sum, normalized by (n1 * n0).  This avoids
        % the O(n1*n0) pairwise comparison and uses O(n*log(n)) sorting.
        ranks = tiedrank(feat_valid);
        auc = (sum(ranks(y_valid == 1)) - n1 * (n1 + 1) / 2) / (n1 * n0);

        % Take abs(auc - 0.5) + 0.5 to handle features where LOWER values
        % predict events (AUC < 0.5 raw → reflected to > 0.5).  This
        % makes the metric direction-agnostic: we care about discriminative
        % power regardless of which direction predicts the outcome.
        auc_val = abs(auc - 0.5) + 0.5;
    else
        % Insufficient data (< 10 observations): return 0.5 to avoid
        % unreliable AUC estimates dominating the tie-breaking.
        auc_val = 0.5;
    end
end
