classdef test_opt < matlab.unittest.TestCase
% TEST_OPT  Benchmark comparing three patient-level outcome aggregation
%   strategies at scale (100k rows, ~1000 patients).
%
%   Approaches compared:
%     1. Baseline loop   — strcmp per unique patient (O(n_unique * n_rows))
%     2. Index-code loop — unique() returns index codes; loop over codes
%     3. accumarray       — fully vectorized grouping with custom function
%
%   All three must produce identical results; timing shows the speedup.

    methods (Test)
        function testAccumarrayMatchesLoop(testCase)
        %TESTACCUMARRAYMATCHESLOOP Verify that all three aggregation
        %   approaches produce identical patient-level outcomes, and report
        %   their wall-clock times.

            % Generate 100k rows with ~1000 unique patient IDs and sparse outcomes
            n_rows = 100000;
            id_list_cell = cell(n_rows, 1);
            y = rand(n_rows, 1) > 0.9;  % ~10% positive rate
            for i = 1:n_rows
                id_list_cell{i} = sprintf('PT%04d', randi(1000));
            end

            % Approach 1: Loop-based with strcmp matching (slowest)
            tic;
            unique_ids = unique(id_list_cell);
            n_unique = numel(unique_ids);
            pt_y = zeros(n_unique, 1);
            for i = 1:n_unique
                pt_y(i) = double(any(y(strcmp(id_list_cell, unique_ids{i})) > 0));
            end
            t1 = toc;

            % Approach 2: Loop using unique()'s index-code vector (ic)
            % Avoids strcmp by using numeric group indices
            tic;
            [unique_ids2, ~, ic] = unique(id_list_cell);
            n_unique2 = numel(unique_ids2);
            pt_y2 = zeros(n_unique2, 1);
            for i = 1:n_unique2
                pt_y2(i) = double(any(y(ic == i) > 0));
            end
            t2 = toc;

            % Approach 3: accumarray — fully vectorized, no explicit loop
            tic;
            [~, ~, ic3] = unique(id_list_cell);
            n_unique3 = numel(unique(id_list_cell));
            pt_y3 = accumarray(ic3, y, [n_unique3, 1], @(x) double(any(x > 0)));
            t3 = toc;

            fprintf('Baseline: %f s, ic loop: %f s, accumarray: %f s\n', t1, t2, t3);
            % All approaches must yield identical patient-level outcomes
            testCase.verifyEqual(pt_y, pt_y2, 'AbsTol', 1e-12);
            testCase.verifyEqual(pt_y, pt_y3, 'AbsTol', 1e-12);
        end
    end
end
