function [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, k, pat_id_tr, pat_id_te, target_cols)
% knn_impute_train_test Imputes missing data safely for cross-validation
%
% Integrates robust K-Nearest Neighbors (KNN) imputation ensuring that the
% imputation model is fitted strictly on the training fold (X_tr) and then
% applied to the test fold (X_te) to prevent data leakage.
%
% Also prevents temporal data leakage by excluding rows from the same patient 
% during distance calculation if pat_id_tr/pat_id_te are provided.
%
% target_cols (optional): Indices or logical mask of columns to exclude from
%                         distance computation (e.g., survival outcomes).
    if nargin < 6, target_cols = []; end
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
    % Create a boolean mask of valid search space coordinates
    search_space = Z_tr;
    if ~isempty(target_cols)
        search_space(:, target_cols) = nan; % mask out target columns from distance metric
    end
    
    for i = 1:n_tr
        missing_idx = isnan(X_tr(i, :));
        if any(missing_idx)
            % Extract the precise query point, masking valid features
            query_pt = search_space(i, :);
            valid_feat = ~isnan(query_pt);
            
            % We cannot search if the query has no valid features left
            if sum(valid_feat) == 0
                continue;
            end
            
            % Build temporary KD-tree strictly on mutually observed features for this query
            % Note: Incomplete data requires querying against subsets where this query's valid
            % features are non-NaN in the reference set.
            ref_idx = find(~any(isnan(search_space(:, valid_feat)), 2));
            
            % Remove self and same-patient from reference index pool
            if ~isempty(pat_id_tr)
                is_same_patient = (pat_id_tr(ref_idx) == pat_id_tr(i));
                ref_idx(is_same_patient) = [];
            else
                ref_idx(ref_idx == i) = []; % just remove self
            end
            
            if isempty(ref_idx)
                continue;
            end
            
            % Extract pure coordinates
            Y_coords = search_space(ref_idx, valid_feat);
            X_coord = query_pt(valid_feat);
            
            % KDTree distance calculation (O(N log N))
            [neighbor_indices_local, ~] = knnsearch(Y_coords, X_coord, 'K', min(k, length(ref_idx)), 'NSMethod', 'kdtree');
            
            % Map back to absolute indices
            neighbors = ref_idx(neighbor_indices_local);
            
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
        
        search_space_te = Z_te;
        if ~isempty(target_cols)
            search_space_te(:, target_cols) = nan;
        end
        
        for i = 1:n_te
            missing_idx = isnan(X_te(i, :));
            if any(missing_idx)
                query_pt = search_space_te(i, :);
                valid_feat = ~isnan(query_pt);
                
                if sum(valid_feat) == 0
                    continue;
                end
                
                ref_idx = find(~any(isnan(search_space(:, valid_feat)), 2));
                
                if ~isempty(pat_id_tr) && ~isempty(pat_id_te)
                    is_same_patient = (pat_id_tr(ref_idx) == pat_id_te(i));
                    ref_idx(is_same_patient) = [];
                end
                
                if isempty(ref_idx)
                    continue;
                end
                
                Y_coords = search_space(ref_idx, valid_feat);
                X_coord = query_pt(valid_feat);
                
                [neighbor_indices_local, ~] = knnsearch(Y_coords, X_coord, 'K', min(k, length(ref_idx)), 'NSMethod', 'kdtree');
                neighbors = ref_idx(neighbor_indices_local);
                
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
