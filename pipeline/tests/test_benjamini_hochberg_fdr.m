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
            % Single p-value: q = p * (n/rank) = p * (1/1) = p
            q = benjamini_hochberg_fdr(0.03);
            testCase.verifyEqual(q, 0.03, 'AbsTol', 1e-12, ...
                'Single p-value should be returned unchanged.');
        end

        function test_single_pvalue_one(testCase)
            % Single p-value of 1.0
            q = benjamini_hochberg_fdr(1.0);
            testCase.verifyEqual(q, 1.0, 'AbsTol', 1e-12);
        end

        function test_single_pvalue_zero(testCase)
            % Single p-value of 0.0
            q = benjamini_hochberg_fdr(0.0);
            testCase.verifyEqual(q, 0.0, 'AbsTol', 1e-12);
        end

        function test_all_ones(testCase)
            % All p-values = 1.0 => all q-values should be 1.0
            p = ones(10, 1);
            q = benjamini_hochberg_fdr(p);
            testCase.verifyEqual(q, ones(10, 1), 'AbsTol', 1e-12, ...
                'All p=1 input should return all q=1.');
        end

        function test_all_zeros(testCase)
            % All p-values = 0.0 => all q-values should be 0.0
            p = zeros(5, 1);
            q = benjamini_hochberg_fdr(p);
            testCase.verifyEqual(q, zeros(5, 1), 'AbsTol', 1e-12, ...
                'All p=0 input should return all q=0.');
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
            testCase.verifyEqual(q, expected, 'AbsTol', 1e-10, ...
                'Q-values must match hand-calculated BH step-up values.');
        end

        function test_known_reference_unsorted_input(testCase)
            % Same classic BH example but with p-values in non-sorted order
            % to verify both numeric correctness AND order preservation together
            p = [0.039; 0.001; 0.23; 0.041; 0.008];
            q = benjamini_hochberg_fdr(p);

            % Expected q-values mapped back to original positions:
            % p=0.039 -> q=0.05125, p=0.001 -> q=0.005, p=0.23 -> q=0.23,
            % p=0.041 -> q=0.05125, p=0.008 -> q=0.02
            expected = [0.05125; 0.005; 0.23; 0.05125; 0.02];
            testCase.verifyEqual(q, expected, 'AbsTol', 1e-10, ...
                'Q-values for unsorted input must match hand-calculated values at original positions.');
        end

        function test_known_reference_three_tests(testCase)
            % Another hand-calculated example with 3 tests
            % p = [0.01, 0.04, 0.30]
            % Sorted: rank1=0.01, rank2=0.04, rank3=0.30
            % adjusted(3) = 0.30 * 3/3 = 0.30
            % adjusted(2) = min(0.30, 0.04 * 3/2) = min(0.30, 0.06) = 0.06
            % adjusted(1) = min(0.06, 0.01 * 3/1) = min(0.06, 0.03) = 0.03
            p = [0.01; 0.04; 0.30];
            q = benjamini_hochberg_fdr(p);
            expected = [0.03; 0.06; 0.30];
            testCase.verifyEqual(q, expected, 'AbsTol', 1e-10, ...
                'Three-test hand-calculated q-values must match.');
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

        function test_q_values_capped_at_one_extreme(testCase)
            % Two p-values where n/i * p would exceed 1 without clamping
            % p = [0.60, 0.90], n=2
            % adjusted(2) = min(1, 0.90 * 2/2) = 0.90
            % adjusted(1) = min(0.90, 0.60 * 2/1) = min(0.90, 1.20) = 0.90
            p = [0.60; 0.90];
            q = benjamini_hochberg_fdr(p);
            testCase.verifyLessThanOrEqual(q, 1.0, ...
                'All q-values must be clamped to [0, 1].');
            testCase.verifyGreaterThanOrEqual(q, 0.0, ...
                'All q-values must be non-negative.');
        end

        function test_q_values_nonnegative(testCase)
            % Verify output is always >= 0
            p = [0.0; 0.001; 0.5; 1.0];
            q = benjamini_hochberg_fdr(p);
            testCase.verifyGreaterThanOrEqual(q, 0.0, ...
                'Q-values must be non-negative.');
            testCase.verifyLessThanOrEqual(q, 1.0, ...
                'Q-values must be at most 1.');
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

        function test_monotonicity_random_input(testCase)
            % Monotonicity with a wider range of values
            rng(42);  % reproducible
            p = rand(50, 1);
            q = benjamini_hochberg_fdr(p);
            [~, sort_idx] = sort(p);
            q_sorted = q(sort_idx);
            diffs = diff(q_sorted);
            testCase.verifyGreaterThanOrEqual(diffs, -1e-12, ...
                'Sorted q-values should be non-decreasing for random input.');
        end

        function test_monotonicity_already_sorted(testCase)
            % p-values already sorted — q should be non-decreasing as-is
            p = [0.001; 0.005; 0.01; 0.05; 0.10; 0.50; 0.90];
            q = benjamini_hochberg_fdr(p);
            diffs = diff(q);
            testCase.verifyGreaterThanOrEqual(diffs, -1e-12, ...
                'Q-values for sorted p should be non-decreasing.');
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

        function test_large_n_clamped_and_monotone(testCase)
            % Large array: verify clamping and monotonicity together
            n = 500;
            rng(123);
            p = rand(n, 1);
            q = benjamini_hochberg_fdr(p);

            % All in [0, 1]
            testCase.verifyGreaterThanOrEqual(q, 0.0);
            testCase.verifyLessThanOrEqual(q, 1.0);

            % Monotonicity in sorted-p order
            [~, sort_idx] = sort(p);
            q_sorted = q(sort_idx);
            diffs = diff(q_sorted);
            testCase.verifyGreaterThanOrEqual(diffs, -1e-12, ...
                'Sorted q-values must be non-decreasing for large input.');
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

        function test_q_geq_p(testCase)
            % BH correction can only increase (or keep equal) p-values
            p = [0.001; 0.01; 0.05; 0.10; 0.50; 0.99];
            q = benjamini_hochberg_fdr(p);
            testCase.verifyGreaterThanOrEqual(q, p - 1e-12, ...
                'Q-values should always be >= corresponding p-values.');
        end

        function test_p_at_fdr_threshold(testCase)
            % p-values exactly at BH threshold boundaries
            % For n=4, alpha=0.05: thresholds are i/n * alpha = [0.0125, 0.025, 0.0375, 0.05]
            % p-values sitting exactly on these thresholds
            p = [0.0125; 0.025; 0.0375; 0.05];
            q = benjamini_hochberg_fdr(p);

            % All should be clamped in [0, 1]
            testCase.verifyGreaterThanOrEqual(q, 0.0);
            testCase.verifyLessThanOrEqual(q, 1.0);

            % Verify monotonicity (already sorted)
            diffs = diff(q);
            testCase.verifyGreaterThanOrEqual(diffs, -1e-12, ...
                'Q-values at FDR thresholds should be non-decreasing.');

            % For these specific values, adjusted = p_i * n/i = 0.05 for all
            expected_all = 0.05 * ones(4, 1);
            testCase.verifyEqual(q, expected_all, 'AbsTol', 1e-10, ...
                'P-values at exact BH thresholds should all yield q = alpha.');
        end
    end
end