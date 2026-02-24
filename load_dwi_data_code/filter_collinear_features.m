function [keep_idx] = filter_collinear_features(X, y)
% FILTER_COLLINEAR_FEATURES Filters collinear features using Pearson |r| > 0.8.
%
%   [keep_idx] = filter_collinear_features(X, y)
%
%   When two features are highly correlated (|r| > 0.8), the one with the higher 
%   univariate Wilcoxon rank-sum p-value (less significant) is dropped.
%
%   Strictly Non-Leaking: This function should be called on training data (X_train, y_train)
%   to obtain keep_idx. The same keep_idx MUST then be applied to the test data
%   without further correlation assessment on the test set.

    if isempty(X)
        keep_idx = [];
        return;
    end

    R = corrcoef(X);
    n_feats = size(X, 2);
    drop_idx = false(1, n_feats);
    
    for fi = 1:n_feats
        if drop_idx(fi), continue; end
        for fj = fi+1:n_feats
            if drop_idx(fj), continue; end
            if abs(R(fi, fj)) > 0.8
                % Compute univariate significance for both features
                p_fi = ranksum(X(y==0, fi), X(y==1, fi));
                p_fj = ranksum(X(y==0, fj), X(y==1, fj));
                
                % Drop the feature that is less significant (higher p-value)
                if p_fj >= p_fi
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
