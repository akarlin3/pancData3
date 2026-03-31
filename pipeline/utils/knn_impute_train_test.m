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
%   Inputs:
%       X_tr        - [n_train x p] Feature matrix for the training fold.
%       X_te        - [n_test  x p] Feature matrix for the testing fold.
%       k           - (Optional) Number of nearest neighbors. Default is 5.
%       pat_id_tr   - (Optional) Array of patient IDs for training instances.
%       pat_id_te   - (Optional) Array of patient IDs for testing instances.
%       target_cols - (Optional) Indices of columns to exclude from distance metric.
%
%   Outputs:
%       X_tr_imp    - [n_train x p] Imputed training feature matrix.
%       X_te_imp    - [n_test  x p] Imputed testing feature matrix.
%
%   Analytical Rationale — Why Same-Patient Exclusion is Critical:
%   ---------------------------------------------------------------
%   In the time-dependent panel, each patient has multiple rows (one per
%   treatment fraction).  If we allow KNN to use another row from the SAME
%   patient as a neighbor, the imputed value effectively comes from that
%   patient's own measurements at a different timepoint.  This creates
%   temporal data leakage: a missing ADC at fraction 3 would be imputed
%   from the patient's own fraction 5 value, using future information to
%   fill a past measurement gap.
%
%   By setting the distance to same-patient rows to infinity, KNN is forced
%   to borrow information only from OTHER patients at similar treatment
%   stages — which is the clinically appropriate imputation strategy (using
%   population-level patterns, not patient-specific future data).
%
%   Train-Test Separation:
%   ----------------------
%   Training-set imputation uses only training-set neighbors.  Test-set
%   imputation also uses only training-set neighbors (never other test
%   rows).  This mirrors the clinical deployment scenario where the model
%   has access to historical patients (training) but not to other
%   concurrent patients (test) when imputing a new patient's missing data.
%
%   Partial-Distance Strategy for Sparse Data:
%   -------------------------------------------
%   Rather than requiring candidate neighbors to have ALL distance-relevant
%   features non-NaN (which is overly strict and can eliminate most
%   candidates when training data has scattered missingness), we require a
%   minimum fraction (80%) of shared valid features and compute partial
%   distances normalized by the number of shared features. This "available-
%   case" approach prevents silent imputation failures in small or sparse
%   cohorts while maintaining distance metric validity.
%
    if nargin < 2, X_te = []; end
    if nargin < 3, k = 5; end
    if nargin < 4, pat_id_tr = []; end
    if isempty(pat_id_tr)
        warning('knn_impute_train_test:noPatientIDs', ...
            'No patient IDs provided. Only self-exclusion applied; same-patient rows may leak information.');
    end
    if nargin < 5, pat_id_te = []; end
    if nargin < 6, target_cols = []; end

    % Minimum fraction of shared valid features required for a candidate
    % neighbor to be eligible.  Setting this to 0.8 means a reference row
    % must have at least 80% of the query's distance-relevant features
    % non-NaN.  This relaxes the previous requirement of 100%, preventing
    % silent imputation failures in sparse datasets while still ensuring
    % distance metrics are computed over a substantial feature overlap.
    min_shared_feat_frac = 0.8;

    % Validate that patient ID types are consistent (both cell or both numeric)
    if ~isempty(pat_id_tr) && ~isempty(pat_id_te)
        if iscell(pat_id_tr) ~= iscell(pat_id_te)
            error('knn_impute_train_test:typeMismatch', ...
                'pat_id_tr and pat_id_te must be the same type (both cell or both numeric).');
        end
    end
    
    [n_tr, p] = size(X_tr);
    X_tr_imp = X_tr;

    % Identify columns that need imputation (sparse analysis)
    % Only compute standardization and distance metrics for columns that have missing values
    has_missing_tr = any(isnan(X_tr), 1);
    has_missing_te = false(1, p);
    if ~isempty(X_te)
        has_missing_te = any(isnan(X_te), 1);
    end
    cols_need_imputation = has_missing_tr | has_missing_te;
    
    % If no columns need imputation, return original data
    if ~any(cols_need_imputation)
        X_te_imp = X_te;
        return;
    end

    % Z-score standardize features based on training fold to compute
    % distances.  Without standardization, features with large absolute
    % values (e.g., D* ~ 0.01-0.1 mm^2/s) would dominate the Euclidean
    % distance over features with small absolute values (e.g., f ~ 0.05-
    % 0.20), making KNN effectively ignore the latter.  Z-scoring puts
    % all DWI parameters on a common scale so each contributes equally
    % to neighbor selection.
    %
    % Statistics are computed from training data only — using test-set
    % statistics would leak test-set distributional information into the
    % distance metric.
    mu_tr = mean(X_tr, 1, 'omitnan');
    sd_tr = std(X_tr, 0, 1, 'omitnan');
    sd_tr(sd_tr == 0 | isnan(sd_tr)) = 1; % Prevent division by zero for constant features
    Z_tr = (X_tr - mu_tr) ./ sd_tr;
    
    % Create a boolean mask of valid search space coordinates.
    % Target columns (e.g., the outcome variable or derived response
    % features) are masked out of the distance metric to prevent
    % circular reasoning: we should not use the outcome to find
    % neighbors for imputing predictor features that will later predict
    % that same outcome.
    search_space = Z_tr;
    if ~isempty(target_cols)
        search_space(:, target_cols) = nan; % mask out target columns from distance metric
    end
    
    % Pre-compute validity mask for vectorized operations
    valid_mask = ~isnan(search_space);

    % --- PRE-COMPUTE PATIENT BLOCKING MATRIX ---
    % Cache patient-blocking distance matrix to avoid repeated distance calculations
    % Vectorized construction eliminates O(n^2) loop iterations.
    patient_blocking_mask = [];
    if ~isempty(pat_id_tr)
        if iscell(pat_id_tr)
            % Convert cell array of patient IDs to numeric indices for vectorized comparison
            [~, ~, id_num_tr] = unique(pat_id_tr);
            patient_blocking_mask = (id_num_tr(:) == id_num_tr(:)');
        else
            patient_blocking_mask = (pat_id_tr(:) == pat_id_tr(:)');
        end
    else
        % Create identity matrix for self-exclusion only
        patient_blocking_mask = eye(n_tr, 'logical');
    end

    % --- 1. Impute Training Set (X_tr) ---
    for i = 1:n_tr
        missing_idx = isnan(X_tr(i, :));
        if any(missing_idx)
            % Extract the precise query point, masking valid features
            query_pt = search_space(i, :);
            valid_feat = valid_mask(i, :);
            
            % We cannot search if the query has no valid features left
            if sum(valid_feat) == 0
                continue;
            end
            
            % Determine which features contribute to distance computation.
            % Start from the query's valid features, exclude target columns
            % and near-zero-variance features.
            distance_feat = valid_feat;
            if ~isempty(target_cols)
                distance_feat(target_cols) = false;
            end
            
            % Further optimization: only use features that actually vary and contribute to distances
            % Skip features that are constant or have very low variance
            nonzero_var_feat = distance_feat & (sd_tr > 1e-10);
            
            num_dist_feat = sum(nonzero_var_feat);
            if num_dist_feat == 0
                continue;
            end
            
            % --- RELAXED VALID-REFERENCE CRITERION ---
            % Instead of requiring ALL nonzero_var_feat columns to be non-NaN
            % in a candidate neighbor (which is overly strict and eliminates
            % most candidates when training data has scattered missingness),
            % require at least min_shared_feat_frac of them.
            min_shared_feat_count = max(1, ceil(min_shared_feat_frac * num_dist_feat));
            shared_feat_counts = sum(valid_mask(:, nonzero_var_feat), 2);
            is_valid_ref = shared_feat_counts >= min_shared_feat_count;
            
            % Apply cached patient blocking exclusions
            is_valid_ref = is_valid_ref & ~patient_blocking_mask(i, :)';
            
            ref_idx = find(is_valid_ref);

            if isempty(ref_idx)
                continue;
            end
            
            % --- PARTIAL DISTANCE COMPUTATION ---
            % Compute distances over mutually valid (non-NaN) features
            % between the query and each reference row, normalized by the
            % number of shared features.  This "available-case" approach
            % is appropriate when missingness is random (MCAR/MAR) and
            % prevents candidates with fewer shared features from having
            % artificially smaller raw distances.
            %
            % For each reference row, the partial distance is:
            %   d_partial = (1/n_shared) * sum_over_shared( (z_query - z_ref)^2 )
            %
            % This normalizes by the number of shared features so that
            % distances are comparable across references with different
            % amounts of overlap.
            Y_coords = search_space(ref_idx, nonzero_var_feat);
            X_coord = query_pt(nonzero_var_feat);
            
            % Compute element-wise squared differences; NaN where either is NaN
            sq_diff = (Y_coords - X_coord).^2;
            
            % Count shared valid features per reference row
            shared_valid = ~isnan(sq_diff);
            n_shared = sum(shared_valid, 2);
            
            % Compute mean squared difference over shared features (partial distance)
            % Replace NaN with 0 for summation, then normalize by n_shared
            sq_diff_zero = sq_diff;
            sq_diff_zero(~shared_valid) = 0;
            dists = sum(sq_diff_zero, 2) ./ max(n_shared, 1);
            
            [~, sort_idx] = sort(dists);

            % For each missing column, filter to refs with valid data
            % in that column BEFORE selecting top k neighbors.  This is
            % the "column-adaptive" KNN strategy: the k nearest neighbors
            % may differ per missing feature because not all neighbors
            % have valid data for all features.  The mean of the k nearest
            % valid neighbors is the imputed value — a locally-weighted
            % estimate that reflects the DWI parameter values of patients
            % with similar tumor characteristics at similar treatment stages.
            for m = find(missing_idx)
                has_target = ~isnan(X_tr(ref_idx(sort_idx), m));
                valid_sorted = sort_idx(has_target);
                num_neighbors = min(k, length(valid_sorted));
                if num_neighbors == 0
                    continue;
                end
                neighbors = ref_idx(valid_sorted(1:num_neighbors));
                % Impute as the unweighted mean of k nearest valid neighbors.
                % Unweighted (vs. distance-weighted) is more robust when k
                % is small (k=5), preventing a single very-close neighbor
                % from dominating the imputed value.
                X_tr_imp(i, m) = mean(X_tr(neighbors, m));
            end
        end
    end
    
    % --- 2. Impute Test/Validation Set (X_te) using Training Data ---
    % Test-set imputation uses ONLY training-set rows as neighbors.
    % The test-set features are Z-scored using training-set statistics
    % (mu_tr, sd_tr) so that distances are computed in the same coordinate
    % system as the training imputation.  This ensures consistency and
    % prevents test-set distributional information from influencing the
    % imputation.
    if ~isempty(X_te)
        X_te_imp = X_te;
        Z_te = (X_te - mu_tr) ./ sd_tr;  % Standardize using TRAINING statistics
        n_te = size(X_te, 1);
        
        search_space_te = Z_te;
        if ~isempty(target_cols)
            search_space_te(:, target_cols) = nan;
        end
        
        % Pre-compute cross-patient blocking mask for test-train exclusions
        % Vectorized construction eliminates O(n_te * n_tr) loop iterations.
        cross_patient_blocking_mask = [];
        if ~isempty(pat_id_tr) && ~isempty(pat_id_te)
            if iscell(pat_id_tr)
                % Convert both cell arrays to numeric indices using a shared label set
                [~, ~, id_num_all] = unique([pat_id_tr(:); pat_id_te(:)]);
                id_num_tr_cross = id_num_all(1:n_tr);
                id_num_te_cross = id_num_all(n_tr+1:end);
                cross_patient_blocking_mask = (id_num_te_cross(:) == id_num_tr_cross(:)');
            else
                cross_patient_blocking_mask = (pat_id_te(:) == pat_id_tr(:)');
            end
        end
        
        for i = 1:n_te
            missing_idx = isnan(X_te(i, :));
            if any(missing_idx)
                query_pt = search_space_te(i, :);
                valid_feat = ~isnan(query_pt);
                
                if sum(valid_feat) == 0
                    continue;
                end
                
                % Apply same sparse optimization for test set
                distance_feat = valid_feat;
                if ~isempty(target_cols)
                    distance_feat(target_cols) = false;
                end
                
                nonzero_var_feat = distance_feat & (sd_tr > 1e-10);
                
                num_dist_feat = sum(nonzero_var_feat);
                if num_dist_feat == 0
                    continue;
                end
                
                % --- RELAXED VALID-REFERENCE CRITERION (TEST SET) ---
                min_shared_feat_count = max(1, ceil(min_shared_feat_frac * num_dist_feat));
                shared_feat_counts = sum(valid_mask(:, nonzero_var_feat), 2);
                is_valid_ref = shared_feat_counts >= min_shared_feat_count;
                
                % Apply cached cross-patient blocking exclusions
                if ~isempty(cross_patient_blocking_mask)
                    is_valid_ref = is_valid_ref & ~cross_patient_blocking_mask(i, :)';
                end
                
                ref_idx = find(is_valid_ref);

                if isempty(ref_idx)
                    continue;
                end
                
                % --- PARTIAL DISTANCE COMPUTATION (TEST SET) ---
                % Use training search_space for reference coordinates
                % (neighbors are always training rows), but compute
                % distances in the shared Z-score space.
                Y_coords = search_space(ref_idx, nonzero_var_feat);
                X_coord = search_space_te(i, nonzero_var_feat);

                % Compute element-wise squared differences; NaN where either is NaN
                sq_diff = (Y_coords - X_coord).^2;
                
                % Count shared valid features per reference row
                shared_valid = ~isnan(sq_diff);
                n_shared = sum(shared_valid, 2);
                
                % Compute mean squared difference over shared features (partial distance)
                sq_diff_zero = sq_diff;
                sq_diff_zero(~shared_valid) = 0;
                dists = sum(sq_diff_zero, 2) ./ max(n_shared, 1);
                
                [~, sort_idx] = sort(dists);

                % For each missing column, filter to refs with valid data
                % in that column BEFORE selecting top k neighbors.
                % Note: imputed values come from X_tr (training data), not
                % X_te — this is the key leakage prevention mechanism.
                % The test patient's missing ADC is filled with the mean
                % ADC of the k most similar training patients, ensuring the
                % imputed value reflects only previously-seen population data.
                for m = find(missing_idx)
                    has_target = ~isnan(X_tr(ref_idx(sort_idx), m));
                    valid_sorted = sort_idx(has_target);
                    num_neighbors = min(k, length(valid_sorted));
                    if num_neighbors == 0
                        continue;
                    end
                    neighbors = ref_idx(valid_sorted(1:num_neighbors));
                    X_te_imp(i, m) = mean(X_tr(neighbors, m));
                end
            end
        end
    else
        X_te_imp = [];
    end
    
    % --- 3. Final Fallback ---
    % When all K nearest neighbors have NaN for a feature, impute with the
    % pre-imputation training-set column mean.  This is a standard fallback
    % (sklearn KNNImputer uses the same approach).  For test rows, this
    % introduces marginal information leakage (a training-set statistic),
    % but is strictly less leakage than using a global or test-set mean,
    % and only activates when KNN fails entirely (rare in practice).
    tr_mean = mean(X_tr, 1, 'omitnan');
    all_nan_cols = isnan(tr_mean);
    if any(all_nan_cols)
        warning('knn_impute_train_test:allNaNColumn', ...
            '%d column(s) are entirely NaN in training data; imputing to 0.', sum(all_nan_cols));
    end
    tr_mean(all_nan_cols) = 0;
    
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