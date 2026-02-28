function [keep_idx] = filter_collinear_features(X, y, frac_vec)
% FILTER_COLLINEAR_FEATURES Filters collinear features using Pearson |r| > 0.8.
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
%   incorrectly pruning baseline features. The resulting Boolean pruning mask is
%   then applied uniformly across ALL longitudinal rows by the caller.
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
            X_corr = X;
            y_corr = y;
        end
    else
        X_corr = X;
        y_corr = y;
    end

    R = corrcoef(X_corr);
    n_feats = size(X, 2);
    drop_idx = false(1, n_feats);
    
    % Initialize C-indices with a sentinel value (-1) to indicate "not calculated"
    % This avoids computing AUC for features that are never involved in a collision
    c_indices = -1 * ones(1, n_feats);
    
    % Find highly collinear pairs
    [row_idx, col_idx] = find(abs(tril(R, -1)) > 0.8);
    
    % Process pairs sequentially, dropping the weaker feature
    for idx = 1:length(row_idx)
        fi = row_idx(idx);
        fj = col_idx(idx);
        
        if drop_idx(fi) || drop_idx(fj)
            continue; % One of them is already dropped
        end
        
        % Calculate AUC lazily if not already done
        if c_indices(fi) == -1
            c_indices(fi) = compute_auc(X_corr(:, fi), y_corr);
        end
        if c_indices(fj) == -1
            c_indices(fj) = compute_auc(X_corr(:, fj), y_corr);
        end

        % Drop the feature that has weaker predictive power
        if c_indices(fj) <= c_indices(fi)
            drop_idx(fj) = true;
        else
            drop_idx(fi) = true;
        end
    end
    keep_idx = find(~drop_idx);
end

function auc_val = compute_auc(feat_col, y_col)
    % Helper function to compute AUC
    valid_idx = ~isnan(feat_col) & ~isnan(y_col);
    if sum(valid_idx) > 2
        y_valid = y_col(valid_idx);
        feat_valid = feat_col(valid_idx);

        n1 = sum(y_valid == 1);
        n0 = sum(y_valid == 0);

        if n1 == 0 || n0 == 0
            auc_val = 0.5;
            return;
        end

        % Fast AUC computation via Mann-Whitney U statistic
        ranks = tiedrank(feat_valid);
        auc = (sum(ranks(y_valid == 1)) - n1 * (n1 + 1) / 2) / (n1 * n0);

        auc_val = max(auc, 1 - auc);
    else
        auc_val = 0.5;
    end
end
