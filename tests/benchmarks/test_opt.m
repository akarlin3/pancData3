classdef test_opt < matlab.unittest.TestCase
    % Benchmark comparing patient-level outcome aggregation (loop vs accumarray)

    methods (Test)
        function testAccumarrayMatchesLoop(testCase)
            n_rows = 100000;
            id_list_cell = cell(n_rows, 1);
            y = rand(n_rows, 1) > 0.9;
            for i = 1:n_rows
                id_list_cell{i} = sprintf('PT%04d', randi(1000));
            end

            % Loop-based
            tic;
            unique_ids = unique(id_list_cell);
            n_unique = numel(unique_ids);
            pt_y = zeros(n_unique, 1);
            for i = 1:n_unique
                pt_y(i) = double(any(y(strcmp(id_list_cell, unique_ids{i})) > 0));
            end
            t1 = toc;

            % ic loop
            tic;
            [unique_ids2, ~, ic] = unique(id_list_cell);
            n_unique2 = numel(unique_ids2);
            pt_y2 = zeros(n_unique2, 1);
            for i = 1:n_unique2
                pt_y2(i) = double(any(y(ic == i) > 0));
            end
            t2 = toc;

            % accumarray
            tic;
            [~, ~, ic3] = unique(id_list_cell);
            n_unique3 = numel(unique(id_list_cell));
            pt_y3 = accumarray(ic3, y, [n_unique3, 1], @(x) double(any(x > 0)));
            t3 = toc;

            fprintf('Baseline: %f s, ic loop: %f s, accumarray: %f s\n', t1, t2, t3);
            testCase.verifyEqual(pt_y, pt_y2, 'AbsTol', 1e-12);
            testCase.verifyEqual(pt_y, pt_y3, 'AbsTol', 1e-12);
        end
    end
end
