classdef test_perf_knn < matlab.unittest.TestCase
    % TEST_PERF_KNN Performance benchmark for the KNN imputation utility.
    %
    % Validates knn_impute_train_test.m (in utils/) with a moderately sized
    % synthetic dataset (500 rows x 50 features, ~10% missing). This serves
    % as both a correctness check (no NaNs remain after imputation) and a
    % runtime benchmark for regression detection.
    %
    % The test uses training-only imputation (empty test set) with k=5
    % neighbors and 100 unique patient IDs. The patient ID vector enforces
    % the temporal leakage safeguard: rows from the same patient receive
    % infinite distance so they cannot serve as neighbors for each other.

    methods (Test)
        function testKnnImputePerformance(testCase)
            % Benchmark: 500 training rows, 50 features, k=5 nearest neighbors.
            % Approximately 10% of entries are set to NaN (missing at random).
            n_tr = 500; p = 50; k = 5;
            X_tr = randn(n_tr, p);

            % Inject ~10% missing values uniformly at random
            missing_idx = rand(n_tr, p) < 0.1;
            X_tr(missing_idx) = NaN;

            % Assign each row to one of 100 patients (multiple rows per patient
            % simulates longitudinal data). The imputer must respect these IDs
            % to prevent intra-patient leakage.
            pat_id_tr = randi(100, n_tr, 1);

            % Time the imputation (empty test set: training-only mode)
            tic;
            X_tr_imp = knn_impute_train_test(X_tr, [], k, pat_id_tr, []);
            t = toc;

            fprintf('KNN impute time: %f s\n', t);

            % Verify output shape is preserved (no rows/columns dropped)
            testCase.verifyEqual(size(X_tr_imp), size(X_tr));
            % Verify all NaN values have been successfully imputed
            testCase.verifyEqual(sum(isnan(X_tr_imp(:))), 0, ...
                'All NaN values should be imputed');
        end
    end
end
