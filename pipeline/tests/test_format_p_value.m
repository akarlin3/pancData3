classdef test_format_p_value < matlab.unittest.TestCase
    % TEST_FORMAT_P_VALUE Unit tests for the format_p_value utility.
    %
    % Validates formatting of p-values for biomedical reporting conventions:
    %   - NaN inputs produce 'p = NaN'
    %   - Values below 0.001 produce 'p < 0.001' (standard threshold)
    %   - Values >= 0.001 produce 'p = X.XXX' (3 decimal places)
    %   - Boundary cases at exactly 0.001, 0, and 1.0

    methods(Test)
        function test_nan_returns_nan_string(testCase)
            % NaN p-values (e.g., from insufficient data) should display
            % as 'p = NaN' rather than causing a formatting error.
            result = format_p_value(NaN);
            testCase.verifyEqual(result, 'p = NaN');
        end

        function test_very_small_p_returns_inequality(testCase)
            % P-values below 0.001 should use inequality notation per
            % biomedical convention (e.g., 'p < 0.001' instead of exact).
            result = format_p_value(0.0001);
            testCase.verifyEqual(result, 'p < 0.001');
        end

        function test_zero_p_returns_inequality(testCase)
            % Edge case: p = 0 (exact) should also produce 'p < 0.001'
            result = format_p_value(0);
            testCase.verifyEqual(result, 'p < 0.001');
        end

        function test_boundary_below_threshold(testCase)
            % Just below 0.001 should still use inequality notation
            result = format_p_value(0.0009);
            testCase.verifyEqual(result, 'p < 0.001');
        end

        function test_boundary_at_threshold(testCase)
            % Exactly 0.001 should use equality notation (not inequality)
            result = format_p_value(0.001);
            testCase.verifyEqual(result, 'p = 0.001');
        end

        function test_typical_significant_p(testCase)
            % A typical significant p-value (< 0.05) formatted to 3 decimals
            result = format_p_value(0.032);
            testCase.verifyEqual(result, 'p = 0.032');
        end

        function test_typical_nonsignificant_p(testCase)
            % A typical non-significant p-value formatted to 3 decimals
            result = format_p_value(0.456);
            testCase.verifyEqual(result, 'p = 0.456');
        end

        function test_p_value_of_one(testCase)
            % Upper bound: p = 1.0 should display as 'p = 1.000'
            result = format_p_value(1.0);
            testCase.verifyEqual(result, 'p = 1.000');
        end

        function test_three_decimal_precision(testCase)
            % Verify correct rounding to 3 decimal places:
            % 0.04567 rounds to 0.046 (not 0.045)
            result = format_p_value(0.04567);
            testCase.verifyEqual(result, 'p = 0.046');
        end
    end
end
