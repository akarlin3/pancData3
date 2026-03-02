classdef test_make_grouped_folds < matlab.unittest.TestCase
    % Minimal sanity check for make_grouped_folds

    methods (Test)
        function testBasicFoldGeneration(testCase)
            id_list_cell = {'A'; 'A'; 'B'; 'C'; 'C'; 'D'};
            y = [1; 0; 0; 1; 1; 0];
            n_folds = 2;
            fold_id = make_grouped_folds(id_list_cell, y, n_folds);
            testCase.verifyEqual(numel(fold_id), numel(y));
            testCase.verifyTrue(all(fold_id >= 1 & fold_id <= n_folds));
        end
    end
end
