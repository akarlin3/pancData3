classdef test_remove_constant_columns < matlab.unittest.TestCase
    % TEST_REMOVE_CONSTANT_COLUMNS Unit tests for remove_constant_columns.
    %
    % Validates that constant, all-NaN, and zero-variance columns are
    % correctly removed while non-constant columns are preserved.

    methods(Test)
        function test_removes_constant_column(testCase)
            % A column with identical values has zero range and should be removed.
            X = [1, 5; 2, 5; 3, 5];
            [X_clean, keep] = remove_constant_columns(X);
            testCase.verifyEqual(X_clean, [1; 2; 3]);
            testCase.verifyEqual(keep, [true, false]);
        end

        function test_keeps_varying_columns(testCase)
            % All columns have variance — nothing should be removed.
            X = [1, 4; 2, 5; 3, 6];
            [X_clean, keep] = remove_constant_columns(X);
            testCase.verifyEqual(X_clean, X);
            testCase.verifyEqual(keep, [true, true]);
        end

        function test_removes_all_nan_column(testCase)
            % A column of all NaN has no finite values and should be removed.
            X = [1, NaN; 2, NaN; 3, NaN];
            [X_clean, keep] = remove_constant_columns(X);
            testCase.verifyEqual(X_clean, [1; 2; 3]);
            testCase.verifyEqual(keep, [true, false]);
        end

        function test_constant_with_nans(testCase)
            % A column with one finite value and NaNs is constant (range=0).
            X = [1, NaN; 2, 5; 3, NaN];
            [~, keep] = remove_constant_columns(X);
            testCase.verifyEqual(keep, [true, false]);
        end

        function test_mixed_columns(testCase)
            % Mix of constant, varying, and all-NaN columns.
            X = [1, 5, NaN, 10;
                 2, 5, NaN, 20;
                 3, 5, NaN, 30];
            [X_clean, keep] = remove_constant_columns(X);
            testCase.verifyEqual(keep, [true, false, false, true]);
            testCase.verifyEqual(X_clean, [1, 10; 2, 20; 3, 30]);
        end

        function test_single_row(testCase)
            % A single-row matrix: all columns are trivially constant.
            X = [1, 2, 3];
            [X_clean, keep] = remove_constant_columns(X);
            testCase.verifyTrue(isempty(X_clean) || size(X_clean, 2) == 0);
            testCase.verifyEqual(keep, [false, false, false]);
        end

        function test_empty_matrix(testCase)
            % Empty matrix should return empty.
            X = zeros(0, 3);
            [X_clean, keep] = remove_constant_columns(X);
            testCase.verifyTrue(isempty(X_clean));
            testCase.verifyEqual(numel(keep), 3);
        end

        function test_varying_with_nans(testCase)
            % A column that varies among its finite values should be kept.
            X = [1, NaN; 2, 5; 3, 10];
            [X_clean, keep] = remove_constant_columns(X);
            testCase.verifyEqual(X_clean, X);
            testCase.verifyEqual(keep, [true, true]);
        end
    end
end
