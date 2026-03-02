classdef test_escape_shell_arg < matlab.unittest.TestCase
    % TEST_ESCAPE_SHELL_ARG Unit tests for shell argument escaping function
    methods (Test)

        function test_unix_simple(testCase)
            % Test simple Unix escaping
            arg = 'simple';
            expected = '''simple''';
            actual = escape_shell_arg(arg, 'unix');
            testCase.verifyEqual(actual, expected);
        end

        function test_unix_single_quote(testCase)
            % Test escaping single quotes in Unix
            arg = 'O''Connor';
            expected = '''O''\''''Connor''';
            actual = escape_shell_arg(arg, 'unix');
            testCase.verifyEqual(actual, expected);
        end

        function test_unix_spaces(testCase)
            % Test escaping spaces in Unix
            arg = 'path with spaces';
            expected = '''path with spaces''';
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

        function test_empty_string_unix(testCase)
            % Empty string should be wrapped in quotes
            actual = escape_shell_arg('', 'unix');
            testCase.verifyEqual(actual, '''''');
        end

        function test_empty_string_pc(testCase)
            % Empty string should be wrapped in double quotes
            actual = escape_shell_arg('', 'pc');
            testCase.verifyEqual(actual, '""');
        end

        function test_unicode_unix(testCase)
            % Non-ASCII characters should pass through correctly
            arg = 'path/to/Pat_Mueller';
            actual = escape_shell_arg(arg, 'unix');
            testCase.verifyEqual(actual, ['''path/to/Pat_Mueller''']);
        end

        function test_string_with_newline_unix(testCase)
            % Newline characters inside a single-quoted string are preserved by the shell
            arg = sprintf('line1\nline2');
            actual = escape_shell_arg(arg, 'unix');
            testCase.verifyTrue(contains(actual, sprintf('\n')), ...
                'Newline should be preserved within the escaped string.');
        end

        function test_string_with_newline_pc(testCase)
            % Newline characters should be preserved in PC escaping
            arg = sprintf('line1\nline2');
            actual = escape_shell_arg(arg, 'pc');
            testCase.verifyTrue(contains(actual, sprintf('\n')), ...
                'Newline should be preserved within the escaped string.');
        end

        function test_special_shell_chars_unix(testCase)
            % Characters like $, `, !, etc. should be safely escaped
            arg = 'file$name`test';
            actual = escape_shell_arg(arg, 'unix');
            expected = '''file$name`test''';
            testCase.verifyEqual(actual, expected);
        end

        function test_numeric_input_throws(testCase)
            % Non-string input should throw an error
            testCase.verifyError(@() escape_shell_arg(123, 'unix'), ...
                'escape_shell_arg:invalidInput');
        end

        function test_pc_multiple_double_quotes(testCase)
            % Multiple double quotes should all be escaped
            arg = 'say "hello" and "bye"';
            actual = escape_shell_arg(arg, 'pc');
            expected = '"say \"hello\" and \"bye\""';
            testCase.verifyEqual(actual, expected);
        end

    end

end
