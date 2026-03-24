classdef test_data_integrity_check < matlab.unittest.TestCase
    % TEST_DATA_INTEGRITY_CHECK Integration tests for patient_data_check.m
    %
    % Creates synthetic patient directory structures and verifies that the
    % pre-pipeline data integrity scanner correctly identifies missing data,
    % categorizes issue severity, and produces valid report structs.
    %
    % Run tests with:
    %   results = runtests('tests/test_data_integrity_check.m');

    properties
        TempDir
        ConfigFile
        DataDir
        OrigPath
    end

    methods (TestMethodSetup)
        function setupTemp(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            testCase.DataDir = fullfile(testCase.TempDir, 'patient_data');
            mkdir(testCase.DataDir);
            testCase.ConfigFile = fullfile(testCase.TempDir, 'config.json');

            testCase.OrigPath = path;
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(baseDir);
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'dependencies'));
        end
    end

    methods (TestMethodTeardown)
        function cleanupTemp(testCase)
            diary off;
            if isfield(testCase, 'OrigPath') && ~isempty(testCase.OrigPath)
                entries = strsplit(testCase.OrigPath, pathsep);
                entries = entries(cellfun(@(p) isempty(p) || isfolder(p), entries));
                path(strjoin(entries, pathsep));
            end
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods (Access = private)
        function writeConfig(testCase, cfg)
            % Write a config struct to the config file as JSON.
            fid = fopen(testCase.ConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);
        end

        function createPatientDir(testCase, dataDir, patName, fractions, varargin) %#ok<INUSL>
            % Create a synthetic patient directory with optional data.
            %   fractions: cell array of fraction names, e.g. {'Fx1_date', 'Fx2_date'}
            %   Name-value pairs:
            %     'dwi'  - logical array matching fractions (create DWI DICOM files)
            %     'gtv'  - logical array matching fractions (create GTV mask files)
            %     'dose' - logical array matching fractions (create RT dose files)
            options = struct('dwi', true(size(fractions)), ...
                             'gtv', true(size(fractions)), ...
                             'dose', true(size(fractions)));
            for vi = 1:2:numel(varargin)
                options.(varargin{vi}) = varargin{vi+1};
            end

            patDir = fullfile(dataDir, patName);
            mkdir(patDir);

            for i = 1:numel(fractions)
                fxDir = fullfile(patDir, fractions{i});
                mkdir(fxDir);

                % DWI DICOM folder
                if options.dwi(i)
                    dwiDir = fullfile(fxDir, 'DWI_series');
                    mkdir(dwiDir);
                    % Create a dummy .dcm file
                    fid = fopen(fullfile(dwiDir, 'test.dcm'), 'w');
                    fwrite(fid, uint8(zeros(1, 10)));
                    fclose(fid);
                end

                % GTV mask
                if options.gtv(i)
                    fid = fopen(fullfile(fxDir, 'GTVp_mask.mat'), 'w');
                    fwrite(fid, uint8(zeros(1, 10)));
                    fclose(fid);
                end

                % RT dose folder (only for Fx fractions, not post)
                if options.dose(i) && ~contains(fractions{i}, 'post')
                    doseDir = fullfile(fxDir, 'rtdose_folder');
                    mkdir(doseDir);
                    fid = fopen(fullfile(doseDir, 'dose.dcm'), 'w');
                    fwrite(fid, uint8(zeros(1, 10)));
                    fclose(fid);
                end
            end
        end
    end

    methods (Test)

        function testCompletePatientNoErrors(testCase)
            % A patient with all 6 fractions, DWI, GTV, and dose should
            % produce zero errors and zero warnings.
            fractions = {'Fx1_20230101', 'Fx2_20230102', 'Fx3_20230103', ...
                         'Fx4_20230104', 'Fx5_20230105', 'post_20230401'};
            testCase.createPatientDir(testCase.DataDir, '001-TestPat', fractions);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyEqual(report.n_patients, 1, ...
                'Should detect 1 patient.');
            testCase.verifyEqual(report.n_errors, 0, ...
                'Complete patient should have no errors.');
        end

        function testMissingBaselineFx1IsError(testCase)
            % Missing Fx1 folder should be classified as ERROR.
            fractions = {'Fx2_20230102', 'Fx3_20230103'};
            testCase.createPatientDir(testCase.DataDir, '002-NoBL', fractions);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyGreaterThan(report.n_errors, 0, ...
                'Missing Fx1 should produce at least one error.');
            % Check for the specific error message
            errorMsgs = report.issues(strcmp(report.issues(:,1), 'ERROR'), 3);
            hasBaselineError = any(cellfun(@(m) contains(m, 'Fx1'), errorMsgs));
            testCase.verifyTrue(hasBaselineError, ...
                'Should report missing Fx1 folder as error.');
        end

        function testMissingDwiAtBaselineIsError(testCase)
            % Fx1 folder exists but has no DWI DICOM files -> ERROR.
            fractions = {'Fx1_20230101'};
            testCase.createPatientDir(testCase.DataDir, '003-NoDWI', fractions, ...
                'dwi', false);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyGreaterThan(report.n_errors, 0, ...
                'Missing DWI at Fx1 should produce an error.');
        end

        function testMissingGtvAtBaselineIsError(testCase)
            % Fx1 folder exists with DWI but no GTV mask -> ERROR.
            fractions = {'Fx1_20230101'};
            testCase.createPatientDir(testCase.DataDir, '004-NoGTV', fractions, ...
                'gtv', false);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyGreaterThan(report.n_errors, 0, ...
                'Missing GTV at Fx1 should produce an error.');
        end

        function testMissingDwiAtLaterFractionIsWarning(testCase)
            % Missing DWI at Fx2 (not baseline) should be WARNING, not ERROR.
            fractions = {'Fx1_20230101', 'Fx2_20230102'};
            testCase.createPatientDir(testCase.DataDir, '005-LateDWI', fractions, ...
                'dwi', [true, false]);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyEqual(report.n_errors, 0, ...
                'Missing DWI at Fx2 should not be an error.');
            testCase.verifyGreaterThan(report.n_warnings, 0, ...
                'Missing DWI at Fx2 should produce a warning.');
        end

        function testNonPatientFoldersSkipped(testCase)
            % Folders without numeric ID prefix should be skipped with INFO.
            mkdir(fullfile(testCase.DataDir, 'config_backup'));
            mkdir(fullfile(testCase.DataDir, 'scripts'));
            fractions = {'Fx1_20230101'};
            testCase.createPatientDir(testCase.DataDir, '006-RealPat', fractions);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyEqual(report.n_patients, 1, ...
                'Should detect only 1 real patient, skipping non-patient folders.');
            % Check for INFO about skipped folders
            infoMsgs = report.issues(strcmp(report.issues(:,1), 'INFO'), 3);
            hasSkipInfo = any(cellfun(@(m) contains(m, 'non-patient'), infoMsgs));
            testCase.verifyTrue(hasSkipInfo, ...
                'Should report skipped non-patient folders as INFO.');
        end

        function testMultiplePatientsWithMixedCompleteness(testCase)
            % Multiple patients: one complete, one missing data.
            fractions_complete = {'Fx1_20230101', 'Fx2_20230102', 'Fx3_20230103', ...
                                  'Fx4_20230104', 'Fx5_20230105', 'post_20230401'};
            testCase.createPatientDir(testCase.DataDir, '007-Complete', fractions_complete);

            fractions_incomplete = {'Fx2_20230102'}; % Missing Fx1
            testCase.createPatientDir(testCase.DataDir, '008-Incomplete', fractions_incomplete);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyEqual(report.n_patients, 2, ...
                'Should detect 2 patients.');
            testCase.verifyGreaterThan(report.n_errors, 0, ...
                'Incomplete patient should produce errors.');
        end

        function testPatientIdFiltering(testCase)
            % Verify that patient_ids config filters the scan.
            fractions = {'Fx1_20230101'};
            testCase.createPatientDir(testCase.DataDir, '010-IncludedPat', fractions);
            testCase.createPatientDir(testCase.DataDir, '011-ExcludedPat', fractions);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '', ...
                         'patient_ids', [10]);
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyEqual(report.n_patients, 1, ...
                'Should only scan the filtered patient.');
        end

        function testMissingDataDirectory(testCase)
            % Missing dataloc should produce an error and return early.
            cfg = struct('dataloc', fullfile(testCase.TempDir, 'nonexistent'), ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyGreaterThan(report.n_errors, 0, ...
                'Missing data directory should produce an error.');
            testCase.verifyEqual(report.n_patients, 0, ...
                'Should report 0 patients when data dir is missing.');
        end

        function testEmptyDatalocConfig(testCase)
            % Empty dataloc in config should produce an error.
            cfg = struct('dataloc', '', ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyGreaterThan(report.n_errors, 0, ...
                'Empty dataloc should produce an error.');
        end

        function testReportStructFields(testCase)
            % Verify the report struct has all expected fields.
            fractions = {'Fx1_20230101'};
            testCase.createPatientDir(testCase.DataDir, '012-Fields', fractions);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyTrue(isfield(report, 'n_patients'), ...
                'Report should have n_patients field.');
            testCase.verifyTrue(isfield(report, 'issues'), ...
                'Report should have issues field.');
            testCase.verifyTrue(isfield(report, 'n_errors'), ...
                'Report should have n_errors field.');
            testCase.verifyTrue(isfield(report, 'n_warnings'), ...
                'Report should have n_warnings field.');
            testCase.verifyTrue(isfield(report, 'patients'), ...
                'Report should have patients table.');
        end

        function testPatientsTableStructure(testCase)
            % Verify the patients table has correct columns.
            fractions = {'Fx1_20230101', 'Fx2_20230102'};
            testCase.createPatientDir(testCase.DataDir, '013-Table', fractions);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyTrue(isfield(report, 'patients'), ...
                'Report should contain patients table.');

            if exist('OCTAVE_VERSION', 'builtin')
                % Octave table shim may not support VariableNames
                return;
            end

            expected_cols = {'Patient', 'Fractions', 'DWI', 'GTV', 'Dose', 'Issues'};
            actual_cols = report.patients.Properties.VariableNames;
            for i = 1:numel(expected_cols)
                testCase.verifyTrue(ismember(expected_cols{i}, actual_cols), ...
                    sprintf('Patients table should have column: %s', expected_cols{i}));
            end
        end

        function testEmptyPatientDirectory(testCase)
            % Data directory exists but contains no patient folders.
            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyEqual(report.n_patients, 0, ...
                'Empty data directory should report 0 patients.');
            testCase.verifyEqual(report.n_errors, 0, ...
                'Empty data directory should have no errors.');
        end

        function testMissingDoseIsInfo(testCase)
            % Missing RT dose at a fraction should be INFO, not error/warning.
            fractions = {'Fx1_20230101'};
            testCase.createPatientDir(testCase.DataDir, '014-NoDose', fractions, ...
                'dose', false);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            % Missing dose should be INFO
            infoMsgs = report.issues(strcmp(report.issues(:,1), 'INFO'), 3);
            hasDoseInfo = any(cellfun(@(m) contains(m, 'RT dose'), infoMsgs));
            testCase.verifyTrue(hasDoseInfo, ...
                'Missing RT dose should be reported as INFO.');
        end

        function testMissingClinicalSpreadsheetIsWarning(testCase)
            % Missing clinical spreadsheet should produce a warning.
            fractions = {'Fx1_20230101'};
            testCase.createPatientDir(testCase.DataDir, '015-NoSheet', fractions);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', 'nonexistent_sheet.xlsx');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            testCase.verifyGreaterThan(report.n_warnings, 0, ...
                'Missing clinical spreadsheet should produce a warning.');
        end

        function testIssuesSeverityCategories(testCase)
            % Verify issues cell array uses only valid severity levels.
            fractions = {'Fx2_20230102'}; % Missing Fx1 -> ERROR
            testCase.createPatientDir(testCase.DataDir, '016-Severity', fractions, ...
                'dose', false); % Missing dose -> INFO

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            valid_severities = {'ERROR', 'WARNING', 'INFO'};
            for i = 1:size(report.issues, 1)
                testCase.verifyTrue(ismember(report.issues{i,1}, valid_severities), ...
                    sprintf('Issue severity "%s" should be ERROR, WARNING, or INFO.', ...
                    report.issues{i,1}));
            end
        end

        function testIssuesCellArrayFormat(testCase)
            % Verify each issue row has exactly 3 columns: severity, patient, message.
            fractions = {'Fx1_20230101'};
            testCase.createPatientDir(testCase.DataDir, '017-Format', fractions);

            cfg = struct('dataloc', testCase.DataDir, ...
                         'dcm2nii_call', '', ...
                         'clinical_data_sheet', '');
            testCase.writeConfig(cfg);

            report = patient_data_check(testCase.ConfigFile);

            if ~isempty(report.issues)
                testCase.verifyEqual(size(report.issues, 2), 3, ...
                    'Issues cell array should have exactly 3 columns.');
            end
        end

    end
end
