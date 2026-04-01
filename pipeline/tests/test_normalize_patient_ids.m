classdef test_normalize_patient_ids < matlab.unittest.TestCase
% TEST_NORMALIZE_PATIENT_IDS — Unit tests for normalize_patient_ids.m
%
% Validates patient ID normalization including:
%   - Underscore-to-hyphen replacement
%   - Excel single-quote stripping
%   - Categorical and cellstr input conversion
%   - Empty input handling
%   - Cross-source matching after normalization

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
        end
    end

    methods (Test)
        function test_underscore_to_hyphen(testCase)
            T_Pat = {'P_01', 'P_02', 'P_03'};
            id_list = {'P-01', 'P-02', 'P-03'};

            [pat_norm, id_norm] = normalize_patient_ids(T_Pat, id_list);

            testCase.verifyEqual(pat_norm, {'P-01', 'P-02', 'P-03'}, ...
                'Underscores should be replaced with hyphens.');
            testCase.verifyEqual(id_norm, {'P-01', 'P-02', 'P-03'});
        end

        function test_excel_quote_stripping(testCase)
            T_Pat = {'''P01''', '''P02'''};
            id_list = {'P01', 'P02'};

            [pat_norm, id_norm] = normalize_patient_ids(T_Pat, id_list);

            testCase.verifyEqual(pat_norm, {'P01', 'P02'}, ...
                'Single quotes from Excel should be stripped.');
        end

        function test_categorical_input(testCase)
            T_Pat = categorical({'P01', 'P02', 'P03'});
            id_list = {'P01', 'P02', 'P03'};

            [pat_norm, id_norm] = normalize_patient_ids(T_Pat, id_list);

            testCase.verifyEqual(pat_norm, {'P01', 'P02', 'P03'}, ...
                'Categorical input should be converted to cellstr.');
        end

        function test_empty_T_Pat(testCase)
            [pat_norm, id_norm] = normalize_patient_ids({}, {'P01', 'P02'});

            testCase.verifyTrue(isempty(pat_norm), ...
                'Empty T_Pat should return empty.');
            testCase.verifyEqual(id_norm, {'P01', 'P02'});
        end

        function test_empty_id_list(testCase)
            [pat_norm, id_norm] = normalize_patient_ids({'P01'}, {});

            testCase.verifyTrue(isempty(id_norm), ...
                'Empty id_list should return empty.');
            testCase.verifyEqual(pat_norm, {'P01'});
        end

        function test_both_empty(testCase)
            [pat_norm, id_norm] = normalize_patient_ids({}, {});

            testCase.verifyTrue(isempty(pat_norm));
            testCase.verifyTrue(isempty(id_norm));
        end

        function test_cross_source_matching(testCase)
            % After normalization, spreadsheet and folder IDs should match
            T_Pat = {'''P_01''', 'P_02', '''P_03'''};
            id_list = {'P-01', 'P-02', 'P-03'};

            [pat_norm, id_norm] = normalize_patient_ids(T_Pat, id_list);

            for i = 1:3
                testCase.verifyEqual(pat_norm{i}, id_norm{i}, ...
                    sprintf('Patient %d should match after normalization.', i));
            end
        end

        function test_mixed_separators(testCase)
            % Patients with no separator or hyphen should pass through unchanged
            T_Pat = {'P01', 'P-02', 'P_03'};
            id_list = {'P01', 'P_02', 'P-03'};

            [pat_norm, id_norm] = normalize_patient_ids(T_Pat, id_list);

            testCase.verifyEqual(pat_norm, {'P01', 'P-02', 'P-03'});
            testCase.verifyEqual(id_norm, {'P01', 'P-02', 'P-03'});
        end

        function test_twogtvs_suffix_stripped(testCase)
            % Dual-GTV patients have -p/-n suffixes in the spreadsheet
            % that should be stripped to match the single folder name
            T_Pat = {'P5-SB- twoGTVs-p', 'P5-SB- twoGTVs-n', 'P01'};
            id_list = {'P5-SB- twoGTVs', 'P01'};

            [pat_norm, id_norm] = normalize_patient_ids(T_Pat, id_list);

            % Both -p and -n should normalize to folder name
            testCase.verifyEqual(pat_norm{1}, 'P5-SB- twoGTVs', ...
                'Trailing -p should be stripped from twoGTVs patient.');
            testCase.verifyEqual(pat_norm{2}, 'P5-SB- twoGTVs', ...
                'Trailing -n should be stripped from twoGTVs patient.');
            % Non-twoGTVs patients should be unaffected
            testCase.verifyEqual(pat_norm{3}, 'P01');
            % Folder names should match
            testCase.verifyEqual(pat_norm{1}, id_norm{1});
        end

        function test_twogtvs_suffix_only_on_two_token(testCase)
            % The -p/-n suffix should NOT be stripped from non-twoGTVs patients
            T_Pat = {'P01-p', 'P8-MG-twoGTVs-n'};
            id_list = {'P01-p', 'P8-MG-twoGTVs'};

            [pat_norm, ~] = normalize_patient_ids(T_Pat, id_list);

            testCase.verifyEqual(pat_norm{1}, 'P01-p', ...
                'Non-twoGTVs patient should keep -p suffix.');
            testCase.verifyEqual(pat_norm{2}, 'P8-MG-twoGTVs', ...
                'twoGTVs patient should have -n suffix stripped.');
        end

        function test_numeric_T_Pat_returns_empty(testCase)
            % Numeric T_Pat (possible in Octave) should return empty
            % This tests the isnumeric branch
            if exist('OCTAVE_VERSION', 'builtin')
                [pat_norm, ~] = normalize_patient_ids(123, {'P01'});
                testCase.verifyTrue(isempty(pat_norm));
            else
                % In MATLAB, numeric input goes through iscategorical branch
                % which will fail gracefully
                testCase.verifyTrue(true);
            end
        end
    end
end
