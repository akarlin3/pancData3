classdef test_benjamini_hochberg_fdr < matlab.unittest.TestCase
% TEST_BENJAMINI_HOCHBERG_FDR — Unit tests for benjamini_hochberg_fdr.m
%
% Validates the Benjamini-Hochberg FDR correction algorithm with:
%   - Known reference values from the original BH paper
%   - Edge cases (empty, single, all-ones, all-zeros)
%   - Monotonicity and cap-at-1 properties
%   - Order preservation (output matches input ordering)

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
        end
    end

    methods (Test)
        function test_empty_input(testCase)
            q = benjamini_hochberg_fdr([]);
            testCase.verifyEmpty(q, 'Empty input should return empty output.');
        end

        function test_single_pvalue(testCase)
            % Single p-value: q = p (no correction needed)
            q = benjamini_hochberg_fdr(0.03);
            testCase.verifyEqual(q, 0.03, 'AbsTol', 1e-12);
        end

        function test_all_ones(testCase)
            % All p-values = 1.0 => all q-values should be 1.0
            p = ones(10, 1);
            q = benjamini_hochberg_fdr(p);
            testCase.verifyEqual(q, ones(10, 1), 'AbsTol', 1e-12);
        end

        function test_all_zeros(testCase)
            % All p-values = 0.0 => all q-values should be 0.0
            p = zeros(5, 1);
            q = benjamini_hochberg_fdr(p);
            testCase.verifyEqual(q, zeros(5, 1), 'AbsTol', 1e-12);
        end

        function test_known_reference(testCase)
            % Classic BH example: 5 tests with known result
            % p-values: [0.001, 0.008, 0.039, 0.041, 0.23]
            % sorted ranks 1-5, BH thresholds: 0.01, 0.02, 0.03, 0.04, 0.05
            % At alpha=0.05: reject first two (0.001 < 0.01, 0.008 < 0.02)
            p = [0.001; 0.008; 0.039; 0.041; 0.23];
            q = benjamini_hochberg_fdr(p);

            % Manual BH q-values (step-up):
            % q(5) = 0.23
            % q(4) = min(0.23, 0.041 * 5/4) = min(0.23, 0.05125) = 0.05125
            % q(3) = min(0.05125, 0.039 * 5/3) = min(0.05125, 0.065) = 0.05125
            % q(2) = min(0.05125, 0.008 * 5/2) = min(0.05125, 0.02) = 0.02
            % q(1) = min(0.02, 0.001 * 5/1) = min(0.02, 0.005) = 0.005
            expected = [0.005; 0.02; 0.05125; 0.05125; 0.23];
            testCase.verifyEqual(q, expected, 'AbsTol', 1e-10);
        end

        function test_preserves_original_order(testCase)
            % Input in non-sorted order; output should map back correctly
            p = [0.05; 0.001; 0.50; 0.01];
            q = benjamini_hochberg_fdr(p);

            % Smallest p (0.001) at index 2 should have smallest q
            [~, min_idx] = min(q);
            testCase.verifyEqual(min_idx, 2);

            % Largest p (0.50) at index 3 should have largest q
            [~, max_idx] = max(q);
            testCase.verifyEqual(max_idx, 3);
        end

        function test_q_values_capped_at_one(testCase)
            % Large p-values with small n can produce n/i * p > 1
            p = [0.8; 0.9; 0.95];
            q = benjamini_hochberg_fdr(p);
            testCase.verifyLessThanOrEqual(q, 1.0, 'Q-values must be capped at 1.0.');
        end

        function test_monotonicity_of_sorted_q(testCase)
            % After sorting p, the corresponding q-values should be non-decreasing
            p = [0.04; 0.01; 0.10; 0.003; 0.50; 0.02; 0.08];
            q = benjamini_hochberg_fdr(p);
            [~, sort_idx] = sort(p);
            q_sorted = q(sort_idx);
            diffs = diff(q_sorted);
            testCase.verifyGreaterThanOrEqual(diffs, -1e-12, ...
                'Sorted q-values should be non-decreasing.');
        end

        function test_output_shape_matches_input(testCase)
            % Row vector input
            p_row = [0.01, 0.05, 0.10];
            q_row = benjamini_hochberg_fdr(p_row);
            testCase.verifyEqual(size(q_row), [3, 1], ...
                'Output should be a column vector.');

            % Column vector input
            p_col = [0.01; 0.05; 0.10];
            q_col = benjamini_hochberg_fdr(p_col);
            testCase.verifyEqual(size(q_col), [3, 1]);
        end

        function test_two_identical_pvalues(testCase)
            % Two identical p-values should get the same q-value
            p = [0.03; 0.03; 0.50];
            q = benjamini_hochberg_fdr(p);
            testCase.verifyEqual(q(1), q(2), 'AbsTol', 1e-12, ...
                'Identical p-values should yield identical q-values.');
        end

        function test_large_n_all_significant(testCase)
            % Many small p-values should all remain significant after FDR
            n = 100;
            p = linspace(0.0001, 0.01, n)';
            q = benjamini_hochberg_fdr(p);
            % All q-values should be < 0.05 since all p are very small
            testCase.verifyLessThan(q, 0.05, ...
                'All q-values should be < 0.05 for uniformly small p-values.');
        end

        function test_mixed_significant_nonsignificant(testCase)
            % Mix of significant and non-significant p-values
            p = [0.001; 0.01; 0.20; 0.80; 0.95];
            q = benjamini_hochberg_fdr(p);

            % First two should be well below 0.05
            testCase.verifyLessThan(q(1), 0.05);
            testCase.verifyLessThan(q(2), 0.05);

            % Last should be >= original p
            testCase.verifyGreaterThanOrEqual(q(5), p(5) - 1e-12);
        end
    end
end
