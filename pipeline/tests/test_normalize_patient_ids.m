classdef test_normalize_patient_ids < matlab.unittest.TestCase
% TEST_NORMALIZE_PATIENT_IDS — Unit tests for normalize_patient_ids.m
%
% Validates patient ID normalization including:
%   - Categorical to cellstr conversion
%   - Single-quote stripping (Excel artefact)
%   - Underscore to hyphen replacement
%   - Empty input handling
%   - Numeric T_Pat returns empty
%   - Char array (multi-row) conversion

    methods (TestMethodSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            addpath(fullfile(testDir, '..', 'utils'));
        end
    end

    methods (Test)
        function test_underscore_to_hyphen(testCase)
            T_Pat = categorical({'P_01', 'P_02', 'P_03'});
            id_list = {'P_01', 'P_02', 'P_03'};

            [pat_norm, id_norm] = normalize_patient_ids(T_Pat, id_list);

            testCase.verifyEqual(pat_norm, {'P-01', 'P-02', 'P-03'});
            testCase.verifyEqual(id_norm, {'P-01', 'P-02', 'P-03'});
        end

        function test_single_quote_stripping(testCase)
            T_Pat = categorical({'''P01''', '''P02'''});
            id_list = {'P01', 'P02'};

            [pat_norm, ~] = normalize_patient_ids(T_Pat, id_list);

            testCase.verifyEqual(pat_norm, {'P01', 'P02'});
        end

        function test_cellstr_input_passthrough(testCase)
            T_Pat = {'P-01', 'P-02'};
            id_list = {'P-01', 'P-02'};

            [pat_norm, id_norm] = normalize_patient_ids(T_Pat, id_list);

            testCase.verifyEqual(pat_norm, {'P-01', 'P-02'});
            testCase.verifyEqual(id_norm, {'P-01', 'P-02'});
        end

        function test_empty_input(testCase)
            [pat_norm, id_norm] = normalize_patient_ids({}, {});

            testCase.verifyEmpty(pat_norm);
            testCase.verifyEmpty(id_norm);
        end

        function test_numeric_T_Pat_returns_empty(testCase)
            % Numeric T_Pat (rare edge case) should return empty
            if exist('OCTAVE_VERSION', 'builtin')
                [pat_norm, ~] = normalize_patient_ids(42, {'P01'});
                testCase.verifyEmpty(pat_norm);
            end
        end

        function test_mixed_formats_match(testCase)
            % After normalization, spreadsheet and folder IDs should match
            T_Pat = categorical({'''P_01''', 'P_02'});
            id_list = {'P-01', 'P_02'};

            [pat_norm, id_norm] = normalize_patient_ids(T_Pat, id_list);

            testCase.verifyEqual(pat_norm{1}, id_norm{1}, ...
                'Quoted+underscore should match hyphen after normalization.');
            testCase.verifyEqual(pat_norm{2}, id_norm{2});
        end

        function test_no_modification_needed(testCase)
            T_Pat = categorical({'P01', 'P02'});
            id_list = {'P01', 'P02'};

            [pat_norm, id_norm] = normalize_patient_ids(T_Pat, id_list);

            testCase.verifyEqual(pat_norm, {'P01', 'P02'});
            testCase.verifyEqual(id_norm, {'P01', 'P02'});
        end

        function test_id_list_only_empty(testCase)
            T_Pat = categorical({'P01'});
            [~, id_norm] = normalize_patient_ids(T_Pat, {});
            testCase.verifyEmpty(id_norm);
        end
    end
end
