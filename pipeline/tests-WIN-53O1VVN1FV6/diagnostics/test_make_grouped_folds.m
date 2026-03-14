classdef test_make_grouped_folds < matlab.unittest.TestCase
% TEST_MAKE_GROUPED_FOLDS  Minimal diagnostic sanity check for
%   make_grouped_folds.
%
%   Verifies basic fold generation: all rows get a valid fold assignment
%   (between 1 and n_folds), and the output length matches the input.
%   This is a lightweight complement to the more thorough tests in
%   tests/test_grouped_folds.m.

    methods (Test)
        function testBasicFoldGeneration(testCase)
        %TESTBASICFOLDGENERATION Verify make_grouped_folds returns valid
        %   fold IDs for a small 4-patient, 2-fold scenario.

            % 4 unique patients (A, B, C, D) with some patients having
            % multiple rows (A has 2 rows, C has 2 rows)
            id_list_cell = {'A'; 'A'; 'B'; 'C'; 'C'; 'D'};
            y = [1; 0; 0; 1; 1; 0];  % Binary outcome per row
            n_folds = 2;

            fold_id = make_grouped_folds(id_list_cell, y, n_folds);

            % Every row must receive a fold assignment
            testCase.verifyEqual(numel(fold_id), numel(y));
            % All fold IDs must be in the valid range [1, n_folds]
            testCase.verifyTrue(all(fold_id >= 1 & fold_id <= n_folds));
        end
    end
end
