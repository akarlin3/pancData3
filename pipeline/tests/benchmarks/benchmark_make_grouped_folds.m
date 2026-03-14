classdef benchmark_make_grouped_folds < matlab.unittest.TestCase
% BENCHMARK_MAKE_GROUPED_FOLDS  Compares loop-based vs accumarray-based
%   patient-level outcome aggregation.
%
%   In make_grouped_folds, each patient's outcome must be reduced to a single
%   value (e.g., "did any row for this patient have y > 0?"). This benchmark
%   compares a naive strcmp-loop approach against a vectorized accumarray
%   approach on a 10,000-row dataset with 1,000 patients, verifying that both
%   produce identical results.

    methods (Test)
        function testAggregationApproaches(testCase)
        %TESTAGGREGATIONAPPROACHES Time loop-based vs accumarray aggregation
        %   and verify they produce identical patient-level outcomes.

            % Generate synthetic patient data: 10k rows, ~1000 unique patients
            n_rows = 10000;
            n_patients = 1000;
            ids = cell(n_rows, 1);
            for i = 1:n_rows
                ids{i} = sprintf('Patient_%04d', randi(n_patients));
            end
            y = rand(n_rows, 1);

            % Approach 1: Loop-based — for each unique patient, use strcmp to
            % find matching rows and check if any outcome is positive.
            % This is O(n_unique * n_rows) and slow for large datasets.
            tic;
            unique_ids = unique(ids);
            n_unique = numel(unique_ids);
            pt_y_old = zeros(n_unique, 1);
            for i = 1:n_unique
                pt_y_old(i) = double(any(y(strcmp(ids, unique_ids{i})) > 0));
            end
            t_old = toc;

            % Approach 2: accumarray — uses unique()'s index vector to group
            % rows by patient in O(n_rows) time, then applies @max per group.
            tic;
            [unique_ids_new, ~, ic] = unique(ids);
            n_unique_new = numel(unique_ids_new);
            pt_y_new = double(accumarray(ic, double(y > 0), [n_unique_new, 1], @max));
            t_new = toc;

            fprintf('Loop: %f s, accumarray: %f s\n', t_old, t_new);
            % Both approaches must produce the same patient-level outcomes
            testCase.verifyEqual(norm(pt_y_old - pt_y_new), 0, 'AbsTol', 1e-12);
        end
    end
end
