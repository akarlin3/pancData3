function [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, k, pat_id_tr, pat_id_te)
% knn_impute_train_test Imputes missing data safely for cross-validation
%
% Integrates robust K-Nearest Neighbors (KNN) imputation ensuring that the
% imputation model is fitted strictly on the training fold (X_tr) and then
% applied to the test fold (X_te) to prevent data leakage.
%
% Also prevents temporal data leakage by excluding rows from the same patient 
% during distance calculation if pat_id_tr/pat_id_te are provided.
    if nargin < 3, k = 5; end
    if nargin < 2, X_te = []; end
    if nargin < 4, pat_id_tr = []; end
    if nargin < 5, pat_id_te = []; end
    
    [n_tr, p] = size(X_tr);
    X_tr_imp = X_tr;
    
    % Z-score standardize features based on training fold to compute distances
    mu_tr = mean(X_tr, 1, 'omitnan');
    sd_tr = std(X_tr, 0, 1, 'omitnan');
    sd_tr(sd_tr == 0 | isnan(sd_tr)) = 1; % Prevent division by zero
    Z_tr = (X_tr - mu_tr) ./ sd_tr;
    
    % --- 1. Impute Training Set (X_tr) ---
    for i = 1:n_tr
        missing_idx = isnan(X_tr(i, :));
        if any(missing_idx)
            dist = inf(n_tr, 1);
            for j = 1:n_tr
                % Condition for neighbor selection:
                % 1. Must not be the same row (i ~= j)
                % 2. If patient IDs are provided, must not be the same patient
                is_same_row = (i == j);
                is_same_patient = ~isempty(pat_id_tr) && (pat_id_tr(i) == pat_id_tr(j));
                
                if ~is_same_row && ~is_same_patient
                    valid_feat = ~missing_idx & ~isnan(X_tr(j, :));
                    if sum(valid_feat) > 0
                        % Euclidean distance normalized by mutually observed features
                        dist(j) = sqrt(sum((Z_tr(i, valid_feat) - Z_tr(j, valid_feat)).^2) / sum(valid_feat));
                    end
                end
            end
            [~, sorted_idx] = sort(dist);
            neighbors = sorted_idx(1:min(k, sum(dist < inf))); % Only use valid neighbors
            
            for m = find(missing_idx)
                vals = X_tr(neighbors, m);
                vals = vals(~isnan(vals));
                if ~isempty(vals)
                    X_tr_imp(i, m) = mean(vals);
                end
            end
        end
    end
    
    % --- 2. Impute Test/Validation Set (X_te) using Training Data ---
    if ~isempty(X_te)
        X_te_imp = X_te;
        Z_te = (X_te - mu_tr) ./ sd_tr;
        n_te = size(X_te, 1);
        for i = 1:n_te
            missing_idx = isnan(X_te(i, :));
            if any(missing_idx)
                dist = inf(n_tr, 1);
                for j = 1:n_tr
                    % Condition for neighbor selection:
                    % If patient IDs are provided, must not be the same patient
                    is_same_patient = ~isempty(pat_id_tr) && ~isempty(pat_id_te) && ...
                                      (pat_id_te(i) == pat_id_tr(j));
                    
                    if ~is_same_patient
                        valid_feat = ~missing_idx & ~isnan(X_tr(j, :));
                        if sum(valid_feat) > 0
                            dist(j) = sqrt(sum((Z_te(i, valid_feat) - Z_tr(j, valid_feat)).^2) / sum(valid_feat));
                        end
                    end
                end
                [~, sorted_idx] = sort(dist);
                neighbors = sorted_idx(1:min(k, sum(dist < inf)));
                
                for m = find(missing_idx)
                    vals = X_tr(neighbors, m);
                    vals = vals(~isnan(vals));
                    if ~isempty(vals)
                        X_te_imp(i, m) = mean(vals);
                    end
                end
            end
        end
    else
        X_te_imp = [];
    end
    
    % --- 3. Final Fallback ---
    tr_mean = mean(X_tr_imp, 1, 'omitnan');
    tr_mean(isnan(tr_mean)) = 0; % Ultimate fallback if entirely NaN
    
    for m = 1:p
        nan_tr = isnan(X_tr_imp(:, m));
        if any(nan_tr)
            X_tr_imp(nan_tr, m) = tr_mean(m);
        end
        
        if ~isempty(X_te)
            nan_te = isnan(X_te_imp(:, m));
            if any(nan_te)
                X_te_imp(nan_te, m) = tr_mean(m);
            end
        end
    end
end
