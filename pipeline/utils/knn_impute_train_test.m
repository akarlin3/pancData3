function test_knn_impute()
% test_knn_impute  Unit tests for knn_impute_train_test edge cases.
%
%   Covers:
%     1. All-NaN neighbor column (should fall back to column-mean or 0)
%     2. k larger than available non-blocked neighbors (should use fewer neighbors)
%     3. Single-patient training set (all rows blocked, imputation impossible)
%     4. Mixed cell/numeric patient IDs (should error on mismatch, work on match)
%     5. Empty X_te
%     6. Basic correctness sanity check
%     7. No missing data passthrough
%     8. Partial NaN neighbor column
%     9. Target column exclusion from distance metric
%    10. Multiple missing columns in single row
%    11. Training imputation with patient blocking
%    12. Single training row (degenerate case)
%    13. Patient blocking changes imputed value (blocked vs unblocked)
%    14. Constant column (zero variance) produces no NaN/Inf
%    15. Negative feature values are correctly handled
%    16. Wide matrix (more columns than rows)
%    17. Distance-weighted vs uniform-weighted imputation
%    18. Scale invariance (standardization before distance computation)
%
%   Usage:
%     >> test_knn_impute   % runs all tests, prints PASS/FAIL summary

    n_pass = 0;
    n_fail = 0;
    test_names = {};
    test_results = {};

    % =====================================================================
    % Test 1: All-NaN neighbor column fallback
    %   When every training row has NaN in a particular column, KNN cannot
    %   find any valid neighbor for that column.  The function should fall
    %   back gracefully (impute to 0 per the all-NaN-column warning path)
    %   rather than error.
    % =====================================================================
    try
        X_tr = [1 NaN 3; 4 NaN 6; 7 NaN 9];
        X_te = [2 NaN 5];
        pat_id_tr = [1; 2; 3];
        pat_id_te = [4];
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, 5, pat_id_tr, pat_id_te);

        % Column 2 is entirely NaN in training data, so fallback should produce 0
        assert(~any(isnan(X_tr_imp(:))), 'Training output should have no NaN');
        assert(~any(isnan(X_te_imp(:))), 'Test output should have no NaN');
        assert(all(X_tr_imp(:,2) == 0), 'All-NaN column in training should be imputed to 0');
        assert(X_te_imp(1,2) == 0, 'All-NaN column in test should be imputed to 0');
        record_pass('AllNaN_neighbor_column_fallback');
    catch me
        record_fail('AllNaN_neighbor_column_fallback', me.message);
    end

    % =====================================================================
    % Test 2: k larger than available non-blocked neighbors
    %   With 3 training rows from different patients and k=10, the function
    %   should silently use fewer neighbors (up to 3) without error.
    % =====================================================================
    try
        X_tr = [1 2; 3 4; 5 6];
        X_te = [2 NaN];
        pat_id_tr = [1; 2; 3];
        pat_id_te = [4];
        k_large = 10;  % much larger than n_tr
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, k_large, pat_id_tr, pat_id_te);

        assert(~any(isnan(X_te_imp(:))), 'Test output should have no NaN with k > n_neighbors');
        % The imputed value should be the mean of column 2 from all 3 training rows
        expected_val = mean([2, 4, 6]);
        assert(abs(X_te_imp(1,2) - expected_val) < 1e-10, ...
            sprintf('Expected %.4f but got %.4f', expected_val, X_te_imp(1,2)));
        record_pass('k_larger_than_available_neighbors');
    catch me
        record_fail('k_larger_than_available_neighbors', me.message);
    end

    % =====================================================================
    % Test 3: Single-patient training set — all rows blocked
    %   When all training rows belong to the same patient as the test row,
    %   patient blocking excludes all neighbors.  KNN imputation is
    %   impossible, so the function should fall back to column-mean
    %   imputation without error.
    % =====================================================================
    try
        X_tr = [1 2; 3 4; 5 6];
        X_te = [2 NaN];
        % All rows (train and test) belong to patient 1
        pat_id_tr = [1; 1; 1];
        pat_id_te = [1];
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, 5, pat_id_tr, pat_id_te);

        assert(~any(isnan(X_te_imp(:))), 'Test output should have no NaN even when all neighbors blocked');
        % Fallback should be column mean of X_tr column 2 = mean([2,4,6]) = 4
        expected_val = mean([2, 4, 6]);
        assert(abs(X_te_imp(1,2) - expected_val) < 1e-10, ...
            sprintf('Expected column-mean fallback %.4f but got %.4f', expected_val, X_te_imp(1,2)));
        record_pass('single_patient_all_blocked');
    catch me
        record_fail('single_patient_all_blocked', me.message);
    end

    % =====================================================================
    % Test 3b: Single unique patient in training set — training self-imputation
    %   All training rows are from the same patient.  Patient blocking within
    %   training means no row can use another row as a neighbor.  KNN
    %   imputation is impossible; fallback should handle gracefully.
    % =====================================================================
    try
        X_tr = [1 NaN; 3 4; NaN 6];
        X_te = [];
        pat_id_tr = [1; 1; 1];
        pat_id_te = [];
        [X_tr_imp, ~] = knn_impute_train_test(X_tr, X_te, 5, pat_id_tr, pat_id_te);

        assert(~any(isnan(X_tr_imp(:))), 'Training output should have no NaN even with single patient');
        % Row 1 col 2 is NaN; no neighbors available, fallback = mean of col 2 = mean([4,6]) = 5
        assert(abs(X_tr_imp(1,2) - 5) < 1e-10, 'Expected column-mean fallback for training NaN');
        % Row 3 col 1 is NaN; fallback = mean of col 1 = mean([1,3]) = 2
        assert(abs(X_tr_imp(3,1) - 2) < 1e-10, 'Expected column-mean fallback for training NaN');
        record_pass('single_patient_training_self_imputation');
    catch me
        record_fail('single_patient_training_self_imputation', me.message);
    end

    % =====================================================================
    % Test 4a: Mixed cell/numeric patient IDs should error
    %   Providing cell pat_id_tr with numeric pat_id_te (or vice versa)
    %   should produce a clear error, not silent incorrect behavior.
    % =====================================================================
    try
        X_tr = [1 2; 3 4];
        X_te = [2 NaN];
        pat_id_tr = {'P1'; 'P2'};  % cell
        pat_id_te = [1];            % numeric
        errored = false;
        try
            [~, ~] = knn_impute_train_test(X_tr, X_te, 5, pat_id_tr, pat_id_te);
        catch inner_me
            if contains(inner_me.identifier, 'typeMismatch')
                errored = true;
            else
                rethrow(inner_me);
            end
        end
        assert(errored, 'Should have thrown typeMismatch error for mixed cell/numeric IDs');
        record_pass('mixed_cell_numeric_IDs_error');
    catch me
        record_fail('mixed_cell_numeric_IDs_error', me.message);
    end

    % =====================================================================
    % Test 4b: Cell patient IDs should work correctly
    %   Both pat_id_tr and pat_id_te as cell arrays of strings should
    %   produce correct patient blocking behavior.
    % =====================================================================
    try
        X_tr = [1 2; 3 4; 5 6];
        X_te = [2 NaN];
        pat_id_tr = {'P1'; 'P2'; 'P3'};
        pat_id_te = {'P4'};
        [~, X_te_imp] = knn_impute_train_test(X_tr, X_te, 5, pat_id_tr, pat_id_te);

        assert(~any(isnan(X_te_imp(:))), 'Cell patient IDs should work without error');
        expected_val = mean([2, 4, 6]);
        assert(abs(X_te_imp(1,2) - expected_val) < 1e-10, ...
            'Cell patient IDs should produce same result as numeric');
        record_pass('cell_patient_IDs_work');
    catch me
        record_fail('cell_patient_IDs_work', me.message);
    end

    % =====================================================================
    % Test 4c: Cell patient IDs with blocking
    %   Verify that cell patient IDs correctly block same-patient neighbors.
    % =====================================================================
    try
        X_tr = [1 10; 3 30; 5 50];
        X_te = [2 NaN];
        % Test patient matches training patient P2
        pat_id_tr = {'P1'; 'P2'; 'P3'};
        pat_id_te = {'P2'};
        [~, X_te_imp] = knn_impute_train_test(X_tr, X_te, 5, pat_id_tr, pat_id_te);

        assert(~any(isnan(X_te_imp(:))), 'Cell IDs with blocking should not produce NaN');
        % P2 (row 2, col2=30) should be blocked; neighbors are P1 (col2=10) and P3 (col2=50)
        expected_val = mean([10, 50]);
        assert(abs(X_te_imp(1,2) - expected_val) < 1e-10, ...
            sprintf('Expected %.4f with P2 blocked, got %.4f', expected_val, X_te_imp(1,2)));
        record_pass('cell_patient_IDs_blocking');
    catch me
        record_fail('cell_patient_IDs_blocking', me.message);
    end

    % =====================================================================
    % Test 5: Empty X_te
    %   Passing an empty test matrix should return an empty matrix without
    %   error.  Training imputation should still proceed normally.
    % =====================================================================
    try
        X_tr = [1 NaN; 3 4; 5 6];
        X_te = [];
        pat_id_tr = [1; 2; 3];
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, 5, pat_id_tr);

        assert(isempty(X_te_imp), 'Empty X_te should return empty X_te_imp');
        assert(~any(isnan(X_tr_imp(:))), 'Training imputation should still work with empty X_te');
        record_pass('empty_X_te');
    catch me
        record_fail('empty_X_te', me.message);
    end

    % =====================================================================
    % Test 5b: Empty X_te with no patient IDs
    % =====================================================================
    try
        X_tr = [1 NaN; 3 4; 5 6];
        X_te = [];
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, 5);

        assert(isempty(X_te_imp), 'Empty X_te should return empty X_te_imp (no pat IDs)');
        assert(~any(isnan(X_tr_imp(:))), 'Training imputation should work without pat IDs');
        record_pass('empty_X_te_no_pat_ids');
    catch me
        record_fail('empty_X_te_no_pat_ids', me.message);
    end

    % =====================================================================
    % Test 6: Basic correctness sanity check
    %   With clean training data and a single missing test value, verify
    %   that the imputed value is the mean of the k nearest neighbors'
    %   values in the missing column.
    % =====================================================================
    try
        % 5 training rows, no missing data, well-separated
        X_tr = [1 10; 2 20; 3 30; 4 40; 5 50];
        X_te = [2.1 NaN];  % closest to row 2 (value 2), then row 3, then row 1
        pat_id_tr = [1; 2; 3; 4; 5];
        pat_id_te = [6];
        k = 3;
        [~, X_te_imp] = knn_impute_train_test(X_tr, X_te, k, pat_id_tr, pat_id_te);

        assert(~isnan(X_te_imp(1,2)), 'Basic imputation should produce non-NaN');
        % After z-scoring col 1: mean=3, std=sqrt(2.5)
        % Z-scores: [-1.2649, -0.6325, 0, 0.6325, 1.2649]
        % Query z-score for 2.1: (2.1-3)/sqrt(2.5) = -0.5692
        % Distances to z-scored query: closest are row2, row1, row3
        % k=3 neighbors: rows 2,1,3 -> col2 values: 20,10,30 -> mean=20
        expected_val = mean([20, 10, 30]);
        assert(abs(X_te_imp(1,2) - expected_val) < 1e-10, ...
            sprintf('Expected %.4f, got %.4f', expected_val, X_te_imp(1,2)));
        record_pass('basic_correctness');
    catch me
        record_fail('basic_correctness', me.message);
    end

    % =====================================================================
    % Test 7: No missing data — should pass through unchanged
    % =====================================================================
    try
        X_tr = [1 2; 3 4; 5 6];
        X_te = [2 3];
        pat_id_tr = [1; 2; 3];
        pat_id_te = [4];
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, 5, pat_id_tr, pat_id_te);

        assert(isequal(X_tr_imp, X_tr), 'No-missing training data should be unchanged');
        assert(isequal(X_te_imp, X_te), 'No-missing test data should be unchanged');
        record_pass('no_missing_data_passthrough');
    catch me
        record_fail('no_missing_data_passthrough', me.message);
    end

    % =====================================================================
    % Test 8: Partial NaN neighbor column
    %   Some neighbors have NaN in the target column, some don't.  The
    %   function should use only neighbors with valid data in that column,
    %   potentially fewer than k.
    % =====================================================================
    try
        % Row 1 and 3 have NaN in col 2; row 2 and 4 have valid col 2
        X_tr = [1 NaN; 2 20; 3 NaN; 4 40];
        X_te = [1.5 NaN];
        pat_id_tr = [1; 2; 3; 4];
        pat_id_te = [5];
        k = 3;
        [~, X_te_imp] = knn_impute_train_test(X_tr, X_te, k, pat_id_tr, pat_id_te);

        assert(~isnan(X_te_imp(1,2)), 'Partial NaN column should still impute');
        % Among k=3 nearest by col 1 (after z-scoring), only those with valid col 2 used
        % The imputed value should come from the valid neighbors only
        record_pass('partial_NaN_neighbor_column');
    catch me
        record_fail('partial_NaN_neighbor_column', me.message);
    end

    % =====================================================================
    % Test 9: Target column exclusion from distance metric
    %   Verify that target_cols are not used in distance computation.
    % =====================================================================
    try
        % Col 3 is the target; cols 1-2 are features
        X_tr = [1 1 100; 2 2 200; 3 3 300; 4 4 400; 5 5 500];
        X_te = [2.5 2.5 NaN];
        pat_id_tr = [1; 2; 3; 4; 5];
        pat_id_te = [6];
        target_cols = 3;
        k = 2;
        [~, X_te_imp] = knn_impute_train_test(X_tr, X_te, k, pat_id_tr, pat_id_te, target_cols);

        assert(~isnan(X_te_imp(1,3)), 'Target column imputation should work');
        % Nearest by cols 1-2 are rows 2 and 3 (values 2,2 and 3,3)
        % Their col 3 values: 200 and 300 -> mean = 250
        expected_val = mean([200, 300]);
        assert(abs(X_te_imp(1,3) - expected_val) < 1e-10, ...
            sprintf('Expected %.4f with target exclusion, got %.4f', expected_val, X_te_imp(1,3)));
        record_pass('target_column_exclusion');
    catch me
        record_fail('target_column_exclusion', me.message);
    end

    % =====================================================================
    % Test 10: Multiple missing columns in single row
    %   A test row with multiple NaN values should have each imputed
    %   independently (potentially from different neighbor sets).
    % =====================================================================
    try
        X_tr = [1 10 100; 2 20 200; 3 30 300];
        X_te = [1.5 NaN NaN];
        pat_id_tr = [1; 2; 3];
        pat_id_te = [4];
        k = 2;
        [~, X_te_imp] = knn_impute_train_test(X_tr, X_te, k, pat_id_tr, pat_id_te);

        assert(~any(isnan(X_te_imp(:))), 'Multiple missing columns should all be imputed');
        record_pass('multiple_missing_columns');
    catch me
        record_fail('multiple_missing_columns', me.message);
    end

    % =====================================================================
    % Test 11: Training imputation with patient blocking
    %   Verify that training rows from the same patient are excluded when
    %   imputing within the training set.
    % =====================================================================
    try
        % Patient 1 has rows 1,2; patient 2 has row 3; patient 3 has row 4
        X_tr = [1 NaN; 1.1 20; 5 50; 6 60];
        X_te = [];
        pat_id_tr = [1; 1; 2; 3];
        % Row 1 is missing col 2.  Without blocking, row 2 (same patient,
        % col2=20) would be the nearest neighbor.  With blocking, rows 3,4
        % are used instead.
        [X_tr_imp, ~] = knn_impute_train_test(X_tr, X_te, 5, pat_id_tr);

        assert(~any(isnan(X_tr_imp(:))), 'Training imputation with blocking should produce no NaN');
        % With patient blocking, row 1 cannot use row 2 as neighbor
        % Neighbors are rows 3 (col2=50) and 4 (col2=60) -> mean = 55
        expected_val = mean([50, 60]);
        assert(abs(X_tr_imp(1,2) - expected_val) < 1e-10, ...
            sprintf('Expected %.4f with patient blocking, got %.4f', expected_val, X_tr_imp(1,2)));
        record_pass('training_patient_blocking');
    catch me
        record_fail('training_patient_blocking', me.message);
    end

    % =====================================================================
    % Test 12: Single training row (degenerate case)
    %   With only one training row that has a NaN, KNN has no neighbors.
    %   Fallback should handle this without error.
    % =====================================================================
    try
        X_tr = [1 NaN];
        X_te = [2 NaN];
        pat_id_tr = [1];
        pat_id_te = [2];
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, 5, pat_id_tr, pat_id_te);

        assert(~any(isnan(X_tr_imp(:))), 'Single-row training should not produce NaN');
        assert(~any(isnan(X_te_imp(:))), 'Single-row test should not produce NaN');
        % Col 2 entirely NaN -> fallback to 0
        assert(X_tr_imp(1,2) == 0, 'All-NaN col with single row should impute to 0');
        assert(X_te_imp(1,2) == 0, 'All-NaN col with single row test should impute to 0');
        record_pass('single_training_row');
    catch me
        record_fail('single_training_row', me.message);
    end

    % =====================================================================
    % Test 13: Patient blocking changes imputed value
    %   With 3 patients, verify that blocking patient P1 produces a
    %   different imputed value than if P1 were not blocked.  This
    %   confirms the blocking mechanism actually affects the result.
    %
    %   Setup: 3 training rows from patients P1, P2, P3.
    %   Test row belongs to P1.  P1's training row (row 1) is the nearest
    %   neighbor by feature distance, so blocking it should change the
    %   imputed value.
    % =====================================================================
    try
        % Col 1 = feature, Col 2 = value to impute
        % Row 1 (P1): feature=1, col2=100  (nearest to test query feature=1.1)
        % Row 2 (P2): feature=5, col2=200
        % Row 3 (P3): feature=6, col2=300
        X_tr = [1 100; 5 200; 6 300];
        X_te = [1.1 NaN];
        pat_id_tr = [1; 2; 3];
        k = 2;

        % Run WITH blocking (test patient = P1, so row 1 is blocked)
        pat_id_te_blocked = [1];
        [~, X_te_blocked] = knn_impute_train_test(X_tr, X_te, k, pat_id_tr, pat_id_te_blocked);

        % Run WITHOUT blocking (test patient = P4, no rows blocked)
        pat_id_te_unblocked = [4];
        [~, X_te_unblocked] = knn_impute_train_test(X_tr, X_te, k, pat_id_tr, pat_id_te_unblocked);

        assert(~isnan(X_te_blocked(1,2)), 'Blocked imputation should not be NaN');
        assert(~isnan(X_te_unblocked(1,2)), 'Unblocked imputation should not be NaN');

        % The blocked and unblocked values should differ because:
        % Unblocked k=2 neighbors: row 1 (P1, col2=100) and row 2 (P2, col2=200)
        % Blocked k=2 neighbors: row 2 (P2, col2=200) and row 3 (P3, col2=300)
        assert(abs(X_te_blocked(1,2) - X_te_unblocked(1,2)) > 1e-6, ...
            sprintf('Blocking should change imputed value: blocked=%.4f, unblocked=%.4f', ...
                X_te_blocked(1,2), X_te_unblocked(1,2)));

        % Verify the blocked result uses only P2 and P3
        % With k=2 and P1 blocked, neighbors are P2 (col2=200) and P3 (col2=300)
        expected_blocked = mean([200, 300]);
        assert(abs(X_te_blocked(1,2) - expected_blocked) < 1e-10, ...
            sprintf('Blocked imputation expected %.4f, got %.4f', expected_blocked, X_te_blocked(1,2)));

        record_pass('patient_blocking_changes_value');
    catch me
        record_fail('patient_blocking_changes_value', me.message);
    end

    % =====================================================================
    % Test 14: Constant column (zero variance)
    %   A column with zero variance would produce NaN when standardized
    %   (division by zero std).  The function should handle this gracefully,
    %   producing no NaN or Inf in distances or imputed values.
    % =====================================================================
    try
        % Col 1 varies, Col 2 is constant (zero variance), Col 3 has missing
        X_tr = [1 5 10; 2 5 20; 3 5 30; 4 5 40];
        X_te = [2.5 5 NaN];
        pat_id_tr = [1; 2; 3; 4];
        pat_id_te = [5];
        k = 2;
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, k, pat_id_tr, pat_id_te);

        assert(~any(isnan(X_te_imp(:))), 'Constant column should not cause NaN in output');
        assert(~any(isinf(X_te_imp(:))), 'Constant column should not cause Inf in output');
        assert(~any(isnan(X_tr_imp(:))), 'Constant column should not cause NaN in training output');
        assert(~any(isinf(X_tr_imp(:))), 'Constant column should not cause Inf in training output');

        % The imputed value in col 3 should be based on nearest neighbors
        % by col 1 (since col 2 is constant and contributes nothing to distance)
        % Nearest to 2.5 in col 1: rows 2 (val=2) and 3 (val=3) -> col3: 20, 30
        expected_val = mean([20, 30]);
        assert(abs(X_te_imp(1,3) - expected_val) < 1e-10, ...
            sprintf('Expected %.4f with constant column, got %.4f', expected_val, X_te_imp(1,3)));

        record_pass('constant_column_zero_variance');
    catch me
        record_fail('constant_column_zero_variance', me.message);
    end

    % =====================================================================
    % Test 15: Negative feature values
    %   Percent-change features (e.g., D_pct, f_delta) can be negative.
    %   Verify that negative values in the feature matrix are handled
    %   correctly during distance computation and imputation.
    % =====================================================================
    try
        % Features include negative values (like percent changes)
        X_tr = [-10 -5 100; -2 3 200; 5 8 300; 12 15 400];
        X_te = [-8 -3 NaN];  % closest to row 1 by feature distance
        pat_id_tr = [1; 2; 3; 4];
        pat_id_te = [5];
        k = 2;
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, k, pat_id_tr, pat_id_te);

        assert(~any(isnan(X_te_imp(:))), 'Negative feature values should not cause NaN');
        assert(~any(isinf(X_te_imp(:))), 'Negative feature values should not cause Inf');
        assert(~any(isnan(X_tr_imp(:))), 'Negative feature training values should not cause NaN');

        % Test query [-8, -3] should be nearest to row 1 [-10, -5] then row 2 [-2, 3]
        % After z-scoring cols 1-2 using training stats:
        %   col1: mean=1.25, std=9.2157 -> z-scores: [-1.22, -0.35, 0.41, 1.17], query: -1.00
        %   col2: mean=5.25, std=8.54   -> z-scores: [-1.20, -0.26, 0.32, 1.14], query: -0.97
        % Row 1 is nearest, then row 2.  k=2 neighbors -> col3 values: 100, 200
        expected_val = mean([100, 200]);
        assert(abs(X_te_imp(1,3) - expected_val) < 1e-10, ...
            sprintf('Negative features: expected %.4f, got %.4f', expected_val, X_te_imp(1,3)));

        record_pass('negative_feature_values');
    catch me
        record_fail('negative_feature_values', me.message);
    end

    % =====================================================================
    % Test 16: Wide matrix (more columns than rows)
    %   With 3 rows and 10 columns (common when assembling 22 features for
    %   small cohorts of 5-10 patients), Euclidean distances can become
    %   unreliable due to the curse of dimensionality.  The function should
    %   still produce valid (no NaN, no Inf) imputed values.
    % =====================================================================
    try
        rng(42);  % reproducibility
        % 3 training rows x 10 columns (wide matrix)
        X_tr = randn(3, 10);
        % Introduce a NaN in the last column of the test row
        X_te = randn(1, 10);
        X_te(1, 10) = NaN;
        % Also introduce a NaN in training to test