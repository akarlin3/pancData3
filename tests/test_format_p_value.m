classdef test_format_p_value < matlab.unittest.TestCase
    % TEST_FORMAT_P_VALUE Unit tests for the format_p_value utility.
    % Validates formatting of p-values for biomedical reporting conventions.

    methods(Test)
        function test_nan_returns_nan_string(testCase)
            result = format_p_value(NaN);
            testCase.verifyEqual(result, 'p = NaN');
        end

        function test_very_small_p_returns_inequality(testCase)
            result = format_p_value(0.0001);
            testCase.verifyEqual(result, 'p < 0.001');
        end

        function test_zero_p_returns_inequality(testCase)
            result = format_p_value(0);
            testCase.verifyEqual(result, 'p < 0.001');
        end

        function test_boundary_below_threshold(testCase)
            result = format_p_value(0.0009);
            testCase.verifyEqual(result, 'p < 0.001');
        end

        function test_boundary_at_threshold(testCase)
            result = format_p_value(0.001);
            testCase.verifyEqual(result, 'p = 0.001');
        end

        function test_typical_significant_p(testCase)
            result = format_p_value(0.032);
            testCase.verifyEqual(result, 'p = 0.032');
        end

        function test_typical_nonsignificant_p(testCase)
            result = format_p_value(0.456);
            testCase.verifyEqual(result, 'p = 0.456');
        end

        function test_p_value_of_one(testCase)
            result = format_p_value(1.0);
            testCase.verifyEqual(result, 'p = 1.000');
        end

        function test_three_decimal_precision(testCase)
            % Verify rounding to 3 decimal places
            result = format_p_value(0.04567);
            testCase.verifyEqual(result, 'p = 0.046');
        end
    end
end
