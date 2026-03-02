classdef test_perf_knn < matlab.unittest.TestCase
    % Performance benchmark for knn_impute_train_test

    methods (Test)
        function testKnnImputePerformance(testCase)
            n_tr = 500; p = 50; k = 5;
            X_tr = randn(n_tr, p);
            missing_idx = rand(n_tr, p) < 0.1;
            X_tr(missing_idx) = NaN;
            pat_id_tr = randi(100, n_tr, 1);

            tic;
            X_tr_imp = knn_impute_train_test(X_tr, [], k, pat_id_tr, []);
            t = toc;

            fprintf('KNN impute time: %f s\n', t);
            testCase.verifyEqual(size(X_tr_imp), size(X_tr));
            testCase.verifyEqual(sum(isnan(X_tr_imp(:))), 0, ...
                'All NaN values should be imputed');
        end
    end
end
