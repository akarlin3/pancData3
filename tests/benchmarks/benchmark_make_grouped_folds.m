classdef benchmark_make_grouped_folds < matlab.unittest.TestCase
    % Benchmark comparing patient-level outcome aggregation strategies

    methods (Test)
        function testAggregationApproaches(testCase)
            n_rows = 10000;
            n_patients = 1000;
            ids = cell(n_rows, 1);
            for i = 1:n_rows
                ids{i} = sprintf('Patient_%04d', randi(n_patients));
            end
            y = rand(n_rows, 1);

            % Loop-based approach
            tic;
            unique_ids = unique(ids);
            n_unique = numel(unique_ids);
            pt_y_old = zeros(n_unique, 1);
            for i = 1:n_unique
                pt_y_old(i) = double(any(y(strcmp(ids, unique_ids{i})) > 0));
            end
            t_old = toc;

            % accumarray approach
            tic;
            [unique_ids_new, ~, ic] = unique(ids);
            n_unique_new = numel(unique_ids_new);
            pt_y_new = double(accumarray(ic, double(y > 0), [n_unique_new, 1], @max));
            t_new = toc;

            fprintf('Loop: %f s, accumarray: %f s\n', t_old, t_new);
            testCase.verifyEqual(norm(pt_y_old - pt_y_new), 0, 'AbsTol', 1e-12);
        end
    end
end
