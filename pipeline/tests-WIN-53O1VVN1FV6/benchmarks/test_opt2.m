classdef test_opt2 < matlab.unittest.TestCase
% TEST_OPT2  Benchmark comparing fold-assignment broadcast strategies
%   used in make_grouped_folds.
%
%   After assigning each unique patient to a CV fold, the fold label must be
%   propagated to every row belonging to that patient. This benchmark compares:
%     1. Baseline: ismember() string matching per fold (O(k * n_rows))
%     2. Optimized: build a patient-level fold vector, then index via ic
%        from unique() (O(n_rows), no string matching)

    methods (Test)
        function testOptimizedMatchesBaseline(testCase)
        %TESTOPTIMIZEDMATCHESBASELINE Verify both fold-broadcast approaches
        %   produce identical row-level fold assignments, and report timing.

            % Generate 100k rows with ~1000 unique patients and 5 CV folds
            n_rows = 100000;
            k = 5;
            id_list_cell = cell(n_rows, 1);
            for i = 1:n_rows
                id_list_cell{i} = sprintf('PT%04d', randi(1000));
            end
            [unique_ids, ~, ic] = unique(id_list_cell);
            n_unique = numel(unique_ids);
            % Random fold assignment per unique patient
            pt_fold_assignments = randi(k, n_unique, 1);
            test_mock = @(f) pt_fold_assignments == f;

            % Approach 1 (Baseline): For each fold, find assigned patients
            % then use ismember() to locate all their rows. String matching
            % on every iteration makes this O(k * n_rows * n_unique).
            tic;
            fold_id1 = zeros(n_rows, 1);
            for f = 1:k
                pt_idx = find(test_mock(f));
                fold_id1(ismember(id_list_cell, unique_ids(pt_idx))) = f;
            end
            t1 = toc;

            % Approach 2 (Optimized): Build a patient-level fold vector,
            % then use the index-code vector ic from unique() to broadcast
            % fold labels to all rows in O(n_rows) with no string ops.
            tic;
            pt_fold = zeros(n_unique, 1);
            for f = 1:k
                pt_idx = find(test_mock(f));
                pt_fold(pt_idx) = f;
            end
            fold_id2 = pt_fold(ic);
            t2 = toc;

            fprintf('Baseline: %f s, Optimized: %f s\n', t1, t2);
            % Both must produce identical row-level fold IDs
            testCase.verifyEqual(fold_id1, fold_id2);
        end
    end
end
