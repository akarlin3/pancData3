function [keep_idx] = filter_collinear_features(X, y, frac_vec)
% FILTER_COLLINEAR_FEATURES Filters collinear features using Pearson |r| > 0.8.
%
%   [keep_idx] = filter_collinear_features(X, y)
%   [keep_idx] = filter_collinear_features(X, y, frac_vec)
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
    
    for fi = 1:n_feats
        if drop_idx(fi), continue; end
        for fj = fi+1:n_feats
            if drop_idx(fj), continue; end
            if abs(R(fi, fj)) > 0.8
                % Compute univariate C-index against binary target for both features.
                % We use ROC AUC for binary targets, taking max(AUC, 1-AUC) 
                % to measure predictive magnitude irrespective of direction.
                
                valid_fi = ~isnan(X_corr(:, fi)) & ~isnan(y_corr);
                if sum(valid_fi) > 2
                    [~, ~, ~, auc_fi] = perfcurve(y_corr(valid_fi), X_corr(valid_fi, fi), 1);
                    c_idx_fi = max(auc_fi, 1 - auc_fi);
                else
                    c_idx_fi = 0.5;
                end
                
                valid_fj = ~isnan(X_corr(:, fj)) & ~isnan(y_corr);
                if sum(valid_fj) > 2
                    [~, ~, ~, auc_fj] = perfcurve(y_corr(valid_fj), X_corr(valid_fj, fj), 1);
                    c_idx_fj = max(auc_fj, 1 - auc_fj);
                else
                    c_idx_fj = 0.5;
                end
                
                % Drop the feature that has weaker predictive power (lower C-index)
                % Ties drop the latter feature by default
                if c_idx_fj <= c_idx_fi
                    drop_idx(fj) = true;
                    % Keep looking for other features to drop for fi
                else
                    drop_idx(fi) = true;
                    % fi is dropped, so we stop looking for its correlates
                    break;
                end
            end
        end
    end
    keep_idx = find(~drop_idx);
end
