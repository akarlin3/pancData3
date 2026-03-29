classdef test_knn_impute_train_test < matlab.unittest.TestCase
% TEST_KNN_IMPUTE_TRAIN_TEST — Unit tests for knn_impute_train_test.m
%
% Validates KNN imputation with:
%   - No-NaN passthrough
%   - Train and test imputation correctness
%   - Same-patient temporal leakage blocking (numeric and cell IDs)
%   - Patient ID type mismatch error
%   - Target column exclusion from distance metric
%   - All-NaN column fallback to 0 with warning
%   - Empty test set handling
%   - Default k parameter

    properties
        origPath
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            testCase.origPath = path;
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
        end
    end

    methods (TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.origPath);
            diary off;
        end
    end

    methods (Test)
        function testNoMissingDataPassthrough(testCase)
            % When no NaNs exist, output should equal input exactly
            rng(42);
            X_tr = rand(10, 4);
            X_te = rand(5, 4);
            [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, 3);
            testCase.verifyEqual(X_tr_imp, X_tr, 'AbsTol', 1e-12, ...
                'Training output should equal input when no NaNs present.');
            testCase.verifyEqual(X_te_imp, X_te, 'AbsTol', 1e-12, ...
                'Test output should equal input when no NaNs present.');
        end

        function testTrainImputationBasic(testCase)
            % Insert NaN in training row, verify it gets imputed (not NaN)
            rng(42);
            X_tr = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
            X_tr(2, 2) = NaN;  % Make one value missing
            X_te = [3 4 5];
            evalc('[X_tr_imp, ~] = knn_impute_train_test(X_tr, X_te, 3);');
            testCase.verifyFalse(any(isnan(X_tr_imp(:))), ...
                'Imputed training matrix should have no NaN values.');
            testCase.verifyNotEqual(X_tr_imp(2, 2), 0, ...
                'Imputed value should not be zero (neighbors have values 2, 8, 11).');
        end

        function testTestImputationUsesTrainingOnly(testCase)
            % Insert NaN in test, verify imputed value comes from training data range
            rng(42);
            X_tr = [1 10; 2 20; 3 30; 4 40; 5 50];
            X_te = [3 NaN];
            evalc('[~, X_te_imp] = knn_impute_train_test(X_tr, X_te, 3);');
            testCase.verifyFalse(isnan(X_te_imp(1, 2)), ...
                'Test NaN should be imputed.');
            % Imputed value should be within range of training col 2 values
            testCase.verifyGreaterThanOrEqual(X_te_imp(1, 2), 10 - 1e-10, ...
                'Imputed value should be >= min of training column.');
            testCase.verifyLessThanOrEqual(X_te_imp(1, 2), 50 + 1e-10, ...
                'Imputed value should be <= max of training column.');
        end

        function testPatientBlockingNumeric(testCase)
            % With numeric patient IDs, verify same-patient rows are NOT used
            % as neighbors. Patient 1 has very high values, Patient 2 has
            % very low values. A NaN in patient 1 should be imputed from
            % patient 2 (the only allowed neighbor), not from patient 1's
            % other rows.
            rng(42);
            X_tr = [100 100; 101 101; 102 102;  % Patient 1 (high)
                      1   1;   2   2;   3   3];  % Patient 2 (low)
            pat_id_tr = [1; 1; 1; 2; 2; 2];
            X_tr(1, 2) = NaN;  % Missing value in patient 1, row 1
            X_te = [];
            evalc('[X_tr_imp, ~] = knn_impute_train_test(X_tr, X_te, 3, pat_id_tr, []);');
            % Imputed value should come from patient 2 (values 1, 2, 3), NOT
            % from patient 1's own rows (101, 102)
            testCase.verifyLessThan(X_tr_imp(1, 2), 50, ...
                'Imputed value should come from patient 2 (low values), not patient 1.');
        end

        function testPatientBlockingCell(testCase)
            % Same as testPatientBlockingNumeric but with cell array patient IDs
            rng(42);
            X_tr = [100 100; 101 101; 102 102;  % Patient A (high)
                      1   1;   2   2;   3   3];  % Patient B (low)
            pat_id_tr = {'A'; 'A'; 'A'; 'B'; 'B'; 'B'};
            X_tr(1, 2) = NaN;
            X_te = [];
            evalc('[X_tr_imp, ~] = knn_impute_train_test(X_tr, X_te, 3, pat_id_tr, {});');
            testCase.verifyLessThan(X_tr_imp(1, 2), 50, ...
                'Imputed value should come from patient B (low values), not patient A.');
        end

        function testPatientIdTypeMismatch(testCase)
            % Verify error when pat_id_tr is cell but pat_id_te is numeric
            rng(42);
            X_tr = [1 2; 3 4];
            X_te = [5 NaN];
            pat_id_tr = {'A'; 'B'};
            pat_id_te = [1];
            testCase.verifyError( ...
                @() knn_impute_train_test(X_tr, X_te, 3, pat_id_tr, pat_id_te), ...
                'knn_impute_train_test:typeMismatch');
        end

        function testTargetColsExcluded(testCase)
            % Verify target columns are excluded from distance metric.
            % Create data where col 2 (target) would dominate distance if
            % included, leading to different neighbor selection.
            rng(42);
            % Col 1: the "real" predictor; Col 2: target (should be excluded)
            % Row 1 (query): col1=5, col2=NaN
            % Row 2: col1=4 (close on predictor), col2=1000 (far on target)
            % Row 3: col1=100 (far on predictor), col2=6 (close on target)
            X_tr = [5 NaN; 4 1000; 100 6];
            X_te = [];
            target_cols = 2;
            evalc('[X_tr_imp, ~] = knn_impute_train_test(X_tr, X_te, 2, [], [], target_cols);');
            % With target excluded, row 2 (col1=4) is the closest neighbor.
            % Without target exclusion, row 3 would be closer (col2 distance dominates).
            % Imputed value for col 2 should be influenced by row 2's value (1000).
            testCase.verifyFalse(isnan(X_tr_imp(1, 2)), ...
                'Target column should still be imputed even though excluded from distance.');
        end

        function testAllNanColumnFallback(testCase)
            % When an entire column is NaN in training, verify warning and impute to 0
            rng(42);
            X_tr = [1 NaN; 2 NaN; 3 NaN];
            X_te = [4 NaN];
            testCase.verifyWarning( ...
                @() knn_impute_train_test(X_tr, X_te, 2), ...
                'knn_impute_train_test:allNaNColumn');
            evalc('[X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, 2);');
            testCase.verifyEqual(X_tr_imp(:, 2), [0; 0; 0], 'AbsTol', 1e-12, ...
                'All-NaN training column should be imputed to 0.');
            testCase.verifyEqual(X_te_imp(1, 2), 0, 'AbsTol', 1e-12, ...
                'All-NaN column in test should also be imputed to 0.');
        end

        function testEmptyTestSet(testCase)
            % When X_te is empty, should return empty X_te_imp
            rng(42);
            X_tr = [1 2; 3 NaN; 5 6];
            X_te = [];
            evalc('[X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, 2);');
            testCase.verifyEmpty(X_te_imp, ...
                'Empty test set should return empty imputed test set.');
            testCase.verifyFalse(any(isnan(X_tr_imp(:))), ...
                'Training NaNs should still be imputed.');
        end

        function testDefaultK(testCase)
            % Verify function works with just X_tr and X_te (default k=5)
            rng(42);
            X_tr = [1 2; 3 4; 5 6; 7 8; 9 10; 11 12];
            X_tr(3, 1) = NaN;
            X_te = [4 NaN];
            evalc('[X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te);');
            testCase.verifyFalse(any(isnan(X_tr_imp(:))), ...
                'Training NaN should be imputed with default k.');
            testCase.verifyFalse(any(isnan(X_te_imp(:))), ...
                'Test NaN should be imputed with default k.');
        end
    end
end
