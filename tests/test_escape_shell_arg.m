classdef test_escape_shell_arg < matlab.unittest.TestCase
    % TEST_ESCAPE_SHELL_ARG Unit tests for shell argument escaping function
    methods (Test)

        function test_unix_simple(testCase)
            % Test simple Unix escaping
            arg = 'simple';
            expected = "'simple'";
            actual = escape_shell_arg(arg, 'unix');
            testCase.verifyEqual(actual, expected);
        end

        function test_unix_single_quote(testCase)
            % Test escaping single quotes in Unix
            arg = "O'Connor";
            expected = "'O'\''Connor'";
            actual = escape_shell_arg(arg, 'unix');
            testCase.verifyEqual(actual, expected);
        end

        function test_unix_spaces(testCase)
            % Test escaping spaces in Unix
            arg = 'path with spaces';
            expected = "'path with spaces'";
            actual = escape_shell_arg(arg, 'unix');
            testCase.verifyEqual(actual, expected);
        end

        function test_pc_simple(testCase)
            % Test simple PC escaping
            arg = 'simple';
            expected = '"simple"';
            actual = escape_shell_arg(arg, 'pc');
            testCase.verifyEqual(actual, expected);
        end

        function test_pc_double_quote(testCase)
            % Test escaping double quotes in PC
            arg = 'quote"mark';
            expected = '"quote\"mark"';
            actual = escape_shell_arg(arg, 'pc');
            testCase.verifyEqual(actual, expected);
        end

        function test_pc_spaces(testCase)
            % Test escaping spaces in PC
            arg = 'path with spaces';
            expected = '"path with spaces"';
            actual = escape_shell_arg(arg, 'pc');
            testCase.verifyEqual(actual, expected);
        end

        function test_pc_trailing_backslash(testCase)
            % Test escaping trailing backslash in PC
            arg = 'path\to\dir\';
            % On Windows, trailing backslash before closing quote must be doubled
            expected = '"path\to\dir\\"';
            actual = escape_shell_arg(arg, 'pc');
            testCase.verifyEqual(actual, expected);
        end

        function test_default_os(testCase)
            % Verify that it runs without style argument (uses ispc)
            arg = 'test';
            actual = escape_shell_arg(arg);
            % Cannot assert value without knowing OS, but verify no crash
            testCase.verifyNotEmpty(actual);
        end

    end

end
