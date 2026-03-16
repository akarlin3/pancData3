classdef test_load_auxiliary_biomarkers < matlab.unittest.TestCase
% TEST_LOAD_AUXILIARY_BIOMARKERS  Tests for auxiliary biomarker loading.
%
%   Verifies: (1) well-formed CSV parses correctly, (2) duplicate rows
%   warn and take last value, (3) empty CSV returns empty struct,
%   (4) patient ID mismatch produces a warning.

    properties
        OriginalPath
        TempDir
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.OriginalPath);
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testWellFormedCSV(testCase)
            % A properly formatted CSV should parse into a struct with
            % one field per biomarker.
            csv_path = fullfile(testCase.TempDir, 'biomarkers.csv');
            fid = fopen(csv_path, 'w');
            fprintf(fid, 'patient_id,biomarker_name,timepoint,value\n');
            fprintf(fid, 'P01,CA199,1,100.5\n');
            fprintf(fid, 'P01,CA199,2,80.3\n');
            fprintf(fid, 'P02,CA199,1,200.0\n');
            fprintf(fid, 'P01,CEA,1,5.2\n');
            fclose(fid);

            id_list = {'P01', 'P02', 'P03'};
            aux = load_auxiliary_biomarkers(csv_path, id_list);

            testCase.verifyTrue(isfield(aux, 'CA199'), 'Should have CA199 field.');
            testCase.verifyTrue(isfield(aux, 'CEA'), 'Should have CEA field.');

            % P01 CA199 tp1 = 100.5
            testCase.verifyEqual(aux.CA199(1, 1), 100.5, 'AbsTol', 1e-6, ...
                'P01 CA199 tp1 should be 100.5.');
            % P01 CA199 tp2 = 80.3
            testCase.verifyEqual(aux.CA199(1, 2), 80.3, 'AbsTol', 1e-6, ...
                'P01 CA199 tp2 should be 80.3.');
            % P02 CA199 tp1 = 200
            testCase.verifyEqual(aux.CA199(2, 1), 200.0, 'AbsTol', 1e-6, ...
                'P02 CA199 tp1 should be 200.');
            % P03 should be NaN (no data)
            testCase.verifyTrue(isnan(aux.CA199(3, 1)), ...
                'P03 CA199 tp1 should be NaN (no data).');
        end

        function testDuplicateEntryWarns(testCase)
            % Duplicate patient/timepoint rows should warn and take last.
            csv_path = fullfile(testCase.TempDir, 'dups.csv');
            fid = fopen(csv_path, 'w');
            fprintf(fid, 'patient_id,biomarker_name,timepoint,value\n');
            fprintf(fid, 'P01,CA199,1,100\n');
            fprintf(fid, 'P01,CA199,1,150\n');  % duplicate
            fclose(fid);

            id_list = {'P01'};
            testCase.verifyWarning( ...
                @() load_auxiliary_biomarkers(csv_path, id_list), ...
                'load_auxiliary_biomarkers:duplicateEntry', ...
                'Duplicate entry should produce a warning.');

            aux = load_auxiliary_biomarkers(csv_path, id_list);
            testCase.verifyEqual(aux.CA199(1, 1), 150, ...
                'Duplicate should take last value (150).');
        end

        function testEmptyCSV(testCase)
            % An empty CSV (header only) should return empty struct.
            csv_path = fullfile(testCase.TempDir, 'empty.csv');
            fid = fopen(csv_path, 'w');
            fprintf(fid, 'patient_id,biomarker_name,timepoint,value\n');
            fclose(fid);

            id_list = {'P01'};
            aux = load_auxiliary_biomarkers(csv_path, id_list);
            testCase.verifyTrue(isempty(fieldnames(aux)), ...
                'Empty CSV should return struct with no fields.');
        end

        function testPatientMismatchWarns(testCase)
            % Patient IDs in CSV not matching pipeline cohort should warn.
            csv_path = fullfile(testCase.TempDir, 'mismatch.csv');
            fid = fopen(csv_path, 'w');
            fprintf(fid, 'patient_id,biomarker_name,timepoint,value\n');
            fprintf(fid, 'P99,CA199,1,100\n');
            fclose(fid);

            id_list = {'P01', 'P02'};
            testCase.verifyWarning( ...
                @() load_auxiliary_biomarkers(csv_path, id_list), ...
                'load_auxiliary_biomarkers:unmatchedPatients', ...
                'Unmatched patient IDs should produce a warning.');
        end
    end
end
