classdef test_has_table_var < matlab.unittest.TestCase
    % TEST_HAS_TABLE_VAR  Table-variable existence test.
    %
    % Locks in behaviour for the input variants the helper has to handle:
    %   * MATLAB tables (uses T.Properties.VariableNames)
    %   * Octave-fallback structs (uses isfield)
    %   * Empty / non-tabular inputs (returns false)
    %
    % Regression: load_dwi_data previously used isfield(T, 'col') which
    % returned false for tables on some MATLAB configurations even when
    % the column existed, silently dropping LF/Immuno/Pat reads.

    methods(TestMethodSetup)
        function addUtilsToPath(testCase) %#ok<INUSD>
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(Test)
        function testTableHasNamedVariable(testCase)
            T = table([1; 2; 3], [0; 1; 0], 'VariableNames', {'Pat', 'LF'});
            testCase.verifyTrue(has_table_var(T, 'Pat'));
            testCase.verifyTrue(has_table_var(T, 'LF'));
        end

        function testTableLacksNamedVariable(testCase)
            T = table([1; 2; 3], [0; 1; 0], 'VariableNames', {'Pat', 'LF'});
            testCase.verifyFalse(has_table_var(T, 'Immuno'));
            testCase.verifyFalse(has_table_var(T, 'Pat ')); % trailing space
        end

        function testStructFallback(testCase)
            S = struct('Pat', {{'P01'}}, 'LF', 0);
            testCase.verifyTrue(has_table_var(S, 'Pat'));
            testCase.verifyTrue(has_table_var(S, 'LF'));
            testCase.verifyFalse(has_table_var(S, 'Immuno'));
        end

        function testEmptyInputReturnsFalse(testCase)
            testCase.verifyFalse(has_table_var([], 'Pat'));
            testCase.verifyFalse(has_table_var({}, 'Pat'));
        end

        function testNonTabularInputReturnsFalse(testCase)
            testCase.verifyFalse(has_table_var(42, 'Pat'));
            testCase.verifyFalse(has_table_var('not a table', 'Pat'));
        end

        function testCaseSensitivityMatchesNative(testCase)
            % MATLAB column names are case-sensitive; helper must mirror.
            T = table([1; 2], 'VariableNames', {'Pat'});
            testCase.verifyTrue(has_table_var(T, 'Pat'));
            testCase.verifyFalse(has_table_var(T, 'pat'));
            testCase.verifyFalse(has_table_var(T, 'PAT'));
        end
    end
end
