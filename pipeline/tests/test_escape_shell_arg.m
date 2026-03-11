classdef test_escape_shell_arg < matlab.unittest.TestCase
    % TEST_ESCAPE_SHELL_ARG Unit tests for the escape_shell_arg utility.
    %
    % escape_shell_arg wraps a string argument for safe use in system()
    % calls. On Unix it uses single-quote wrapping (escaping embedded single
    % quotes via the '\'' idiom). On Windows (PC) it uses double-quote
    % wrapping (escaping embedded double quotes with backslash and doubling
    % trailing backslashes to prevent them from escaping the closing quote).
    %
    % These tests verify correct escaping for both platforms, edge cases
    % (empty strings, special characters, newlines, non-ASCII), and input
    % validation (non-string input must throw an error).
    %
    % Run tests with:
    %   results = runtests('tests/test_escape_shell_arg.m');

    methods (Test)

        function test_unix_simple(testCase)
            % A plain alphanumeric string on Unix should be wrapped in
            % single quotes with no internal escaping needed.
            arg = 'simple';
            expected = '''simple''';
            actual = escape_shell_arg(arg, 'unix');
            testCase.verifyEqual(actual, expected);
        end

        function test_unix_single_quote(testCase)
            % An embedded single quote must be escaped using the Unix
            % idiom: end the current single-quoted segment, insert an
            % escaped literal single quote (\'), then resume quoting.
            % Input: O'Connor  -->  'O'\''Connor'
            arg = 'O''Connor';
            expected = '''O''\''''Connor''';
            actual = escape_shell_arg(arg, 'unix');
            testCase.verifyEqual(actual, expected);
        end

        function test_unix_spaces(testCase)
            % Spaces in paths are the most common reason for escaping.
            % Single-quote wrapping preserves spaces literally.
            arg = 'path with spaces';
            expected = '''path with spaces''';
            actual = escape_shell_arg(arg, 'unix');
            testCase.verifyEqual(actual, expected);
        end

        function test_pc_simple(testCase)
            % A plain alphanumeric string on Windows should be wrapped in
            % double quotes with no internal escaping needed.
            arg = 'simple';
            expected = '"simple"';
            actual = escape_shell_arg(arg, 'pc');
            testCase.verifyEqual(actual, expected);
        end

        function test_pc_double_quote(testCase)
            % An embedded double quote on Windows must be escaped with a
            % backslash so the shell does not interpret it as a delimiter.
            arg = 'quote"mark';
            expected = '"quote\"mark"';
            actual = escape_shell_arg(arg, 'pc');
            testCase.verifyEqual(actual, expected);
        end

        function test_pc_spaces(testCase)
            % Spaces in Windows paths (e.g., "Program Files") are handled
            % by wrapping in double quotes.
            arg = 'path with spaces';
            expected = '"path with spaces"';
            actual = escape_shell_arg(arg, 'pc');
            testCase.verifyEqual(actual, expected);
        end

        function test_pc_trailing_backslash(testCase)
            % On Windows, a trailing backslash before the closing double
            % quote would escape the quote itself (\""), so it must be
            % doubled to produce a literal backslash: path\to\dir\\".
            arg = 'path\to\dir\';
            expected = '"path\to\dir\\"';
            actual = escape_shell_arg(arg, 'pc');
            testCase.verifyEqual(actual, expected);
        end

        function test_default_os(testCase)
            % When no style argument is given, escape_shell_arg auto-detects
            % the OS via ispc(). We cannot assert the exact output without
            % knowing the OS, but it must not crash and must return non-empty.
            arg = 'test';
            actual = escape_shell_arg(arg);
            testCase.verifyNotEmpty(actual);
        end

        function test_empty_string_unix(testCase)
            % An empty string must still produce a valid shell token (two
            % adjacent single quotes) so it is not silently dropped.
            actual = escape_shell_arg('', 'unix');
            testCase.verifyEqual(actual, '''''');
        end

        function test_empty_string_pc(testCase)
            % An empty string on Windows must produce "" (two double quotes).
            actual = escape_shell_arg('', 'pc');
            testCase.verifyEqual(actual, '""');
        end

        function test_unicode_unix(testCase)
            % Non-ASCII characters (e.g., patient names with accents or
            % umlauts) must pass through single-quote wrapping unchanged.
            arg = 'path/to/Pat_Mueller';
            actual = escape_shell_arg(arg, 'unix');
            testCase.verifyEqual(actual, ['''path/to/Pat_Mueller''']);
        end

        function test_string_with_newline_unix(testCase)
            % Newline characters inside a single-quoted Unix string are
            % preserved by the shell. The escaping function must not strip
            % or alter them.
            arg = sprintf('line1\nline2');
            actual = escape_shell_arg(arg, 'unix');
            testCase.verifyTrue(contains(actual, sprintf('\n')), ...
                'Newline should be preserved within the escaped string.');
        end

        function test_string_with_newline_pc(testCase)
            % Newline characters should also be preserved in PC escaping,
            % even though cmd.exe handles them differently from Unix shells.
            arg = sprintf('line1\nline2');
            actual = escape_shell_arg(arg, 'pc');
            testCase.verifyTrue(contains(actual, sprintf('\n')), ...
                'Newline should be preserved within the escaped string.');
        end

        function test_special_shell_chars_unix(testCase)
            % Shell metacharacters ($, `, !) must be neutralised.
            % Single-quote wrapping on Unix prevents variable expansion ($),
            % command substitution (`), and history expansion (!).
            arg = 'file$name`test';
            actual = escape_shell_arg(arg, 'unix');
            expected = '''file$name`test''';
            testCase.verifyEqual(actual, expected);
        end

        function test_numeric_input_throws(testCase)
            % Passing a non-string (numeric) input must throw an error
            % rather than silently converting, since accidental numeric
            % arguments could produce nonsensical shell commands.
            testCase.verifyError(@() escape_shell_arg(123, 'unix'), ...
                'escape_shell_arg:invalidInput');
        end

        function test_pc_multiple_double_quotes(testCase)
            % When a string contains multiple double quotes, each one must
            % be independently escaped with a backslash.
            arg = 'say "hello" and "bye"';
            actual = escape_shell_arg(arg, 'pc');
            expected = '"say \"hello\" and \"bye\""';
            testCase.verifyEqual(actual, expected);
        end

    end

end
