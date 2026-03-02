classdef test_opt2 < matlab.unittest.TestCase
    % Benchmark comparing fold-assignment optimization approaches

    methods (Test)
        function testOptimizedMatchesBaseline(testCase)
            n_rows = 100000;
            k = 5;
            id_list_cell = cell(n_rows, 1);
            for i = 1:n_rows
                id_list_cell{i} = sprintf('PT%04d', randi(1000));
            end
            [unique_ids, ~, ic] = unique(id_list_cell);
            n_unique = numel(unique_ids);
            pt_fold_assignments = randi(k, n_unique, 1);
            test_mock = @(f) pt_fold_assignments == f;

            % Baseline: ismember approach
            tic;
            fold_id1 = zeros(n_rows, 1);
            for f = 1:k
                pt_idx = find(test_mock(f));
                fold_id1(ismember(id_list_cell, unique_ids(pt_idx))) = f;
            end
            t1 = toc;

            % Optimized: index-based approach
            tic;
            pt_fold = zeros(n_unique, 1);
            for f = 1:k
                pt_idx = find(test_mock(f));
                pt_fold(pt_idx) = f;
            end
            fold_id2 = pt_fold(ic);
            t2 = toc;

            fprintf('Baseline: %f s, Optimized: %f s\n', t1, t2);
            testCase.verifyEqual(fold_id1, fold_id2);
        end
    end
end
