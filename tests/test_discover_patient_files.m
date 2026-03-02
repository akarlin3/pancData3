classdef test_discover_patient_files < matlab.unittest.TestCase
    % TEST_DISCOVER_PATIENT_FILES Formal unit test suite for discover_patient_files.m

    properties
        MockDataDir
        MockPathDir
        OrigPath
    end

    methods(TestMethodSetup)
        function createMockData(testCase)
            % 1. Create a temporary directory for mock data
            testCase.MockDataDir = fullfile(tempdir, 'mock_patient_data');
            if exist(testCase.MockDataDir, 'dir')
                rmdir(testCase.MockDataDir, 's');
            end
            mkdir(testCase.MockDataDir);

            % 2. Create a temporary directory for mock path (to mock dicominfo)
            testCase.MockPathDir = fullfile(tempdir, 'mock_path_dir');
            if exist(testCase.MockPathDir, 'dir')
                rmdir(testCase.MockPathDir, 's');
            end
            mkdir(testCase.MockPathDir);

            % Write a mock dicominfo function
            mock_dicominfo_file = fullfile(testCase.MockPathDir, 'dicominfo.m');
            fid = fopen(mock_dicominfo_file, 'w');
            fprintf(fid, 'function info = dicominfo(filename)\n');
            fprintf(fid, '    info.PatientID = ''MOCK-MRN'';\n');
            fprintf(fid, '    info.StudyDate = ''20240101'';\n');
            fprintf(fid, 'end\n');
            fclose(fid);

            % Add to path, saving the original path
            testCase.OrigPath = path;
            addpath(testCase.MockPathDir);

            % Add necessary paths for tests
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'dependencies'));

            % --- Create Mock Patient Directory Structure ---

            % Patient 1: Standard patient (P01-ABC)
            % Has Fx1 and post.
            % Fx1 has DWI (with 5+ .dcm files), rtdose (with .dcm), and GTV
            p1_dir = fullfile(testCase.MockDataDir, 'P01-ABC');
            mkdir(p1_dir);
            mkdir(fullfile(p1_dir, 'Fx1'));
            mkdir(fullfile(p1_dir, 'post'));

            p1_fx1_dwi = fullfile(p1_dir, 'Fx1', '1_DWI_dir');
            mkdir(p1_fx1_dwi);
            for i=1:5
                fclose(fopen(fullfile(p1_fx1_dwi, sprintf('im%d.dcm', i)), 'w'));
            end

            p1_fx1_rtdose = fullfile(p1_dir, 'Fx1', '2_rtdose');
            mkdir(p1_fx1_rtdose);
            fclose(fopen(fullfile(p1_fx1_rtdose, 'dose.dcm'), 'w'));

            % GTV file for P01-ABC in Fx1
            fclose(fopen(fullfile(p1_dir, 'Fx1', 'GTV_1.mat'), 'w'));

            % Patient 2: Missing Data (P02_DEF)
            % Has Fx1, but missing DWI and rtdose folders
            p2_dir = fullfile(testCase.MockDataDir, 'P02_DEF');
            mkdir(p2_dir);
            mkdir(fullfile(p2_dir, 'Fx1'));

            % Patient 3: Two GTVs (P03_twoGTV)
            % Has Fx1 with DWI (DICOM subfolder fallback) and both GTVp and GTVn
            p3_dir = fullfile(testCase.MockDataDir, 'P03_twoGTV');
            mkdir(p3_dir);
            mkdir(fullfile(p3_dir, 'Fx1'));

            p3_fx1_dwi = fullfile(p3_dir, 'Fx1', 'DICOM_dir');
            mkdir(p3_fx1_dwi);
            % Fallback: DICOM nested folder
            mkdir(fullfile(p3_fx1_dwi, 'DICOM'));
            for i=1:5
                fclose(fopen(fullfile(p3_fx1_dwi, 'DICOM', sprintf('im%d.dcm', i)), 'w'));
            end

            % GTV files for P03_twoGTV
            fclose(fopen(fullfile(p3_dir, 'Fx1', 'GTVp1.mat'), 'w'));
            fclose(fopen(fullfile(p3_dir, 'Fx1', 'GTV_LN1.mat'), 'w'));

            % Non-patient folder to ignore
            mkdir(fullfile(testCase.MockDataDir, 'template'));
        end
    end

    methods(TestMethodTeardown)
        function removeMockData(testCase)
            % Restore original path
            path(testCase.OrigPath);

            if exist(testCase.MockDataDir, 'dir')
                rmdir(testCase.MockDataDir, 's');
            end
            if exist(testCase.MockPathDir, 'dir')
                rmdir(testCase.MockPathDir, 's');
            end
        end
    end

    methods(Test)
        function testBasicDiscovery(testCase)
            % Run the function
            [id_list, mrn_list, fx_dates, dwi_locations, rtdose_locations, gtv_locations, gtvn_locations] = discover_patient_files(testCase.MockDataDir);

            % Should find 3 patients, template should be ignored
            testCase.verifyEqual(length(id_list), 3, 'Should discover exactly 3 valid patients.');

            % Verify IDs are sorted ascending (1, 2, 3)
            testCase.verifyEqual(id_list{1}, 'P01-ABC');
            testCase.verifyEqual(id_list{2}, 'P02_DEF');
            testCase.verifyEqual(id_list{3}, 'P03_twoGTV');
        end

        function testStandardPatient(testCase)
            [~, mrn_list, fx_dates, dwi_locations, rtdose_locations, gtv_locations, gtvn_locations] = discover_patient_files(testCase.MockDataDir);

            % P01-ABC is index 1
            % MRN and Date from mock dicominfo
            testCase.verifyEqual(mrn_list{1}, 'MOCK-MRN');
            testCase.verifyEqual(fx_dates{1, 1}, '20240101');

            % DWI and rtdose found in Fx1
            testCase.verifyTrue(~isempty(dwi_locations{1, 1, 1}), 'DWI location should be populated for P01.');
            testCase.verifyTrue(contains(dwi_locations{1, 1, 1}, '1_DWI_dir'));

            testCase.verifyTrue(~isempty(rtdose_locations{1, 1}), 'RTdose location should be populated for P01.');
            testCase.verifyTrue(contains(rtdose_locations{1, 1}, '2_rtdose'));

            % GTV found, no GTVn
            testCase.verifyTrue(~isempty(gtv_locations{1, 1, 1}), 'GTV location should be found.');
            testCase.verifyTrue(isempty(gtvn_locations{1, 1, 1}), 'GTVn should not be found for P01.');
        end

        function testMissingDataPatient(testCase)
            [~, mrn_list, ~, dwi_locations, rtdose_locations, gtv_locations, ~] = discover_patient_files(testCase.MockDataDir);

            % P02_DEF is index 2
            testCase.verifyTrue(isempty(mrn_list{2}), 'MRN should be empty for missing data patient.');

            % No DWI or rtdose
            testCase.verifyTrue(isempty(dwi_locations{2, 1, 1}), 'DWI location should be empty for P02.');
            testCase.verifyTrue(isempty(rtdose_locations{2, 1}), 'RTdose location should be empty for P02.');
            testCase.verifyTrue(isempty(gtv_locations{2, 1, 1}), 'GTV location should be empty for P02.');
        end

        function testTwoGTVPatient(testCase)
            [~, ~, ~, dwi_locations, ~, gtv_locations, gtvn_locations] = discover_patient_files(testCase.MockDataDir);

            % P03_twoGTV is index 3
            % It uses nested DICOM fallback, which changes the DWI location name
            testCase.verifyTrue(~isempty(dwi_locations{3, 1, 1}), 'DWI location should be found via DICOM fallback.');
            testCase.verifyTrue(contains(dwi_locations{3, 1, 1}, 'DICOM'), 'Should use nested DICOM fallback path.');

            % Both GTVp and GTVn should be found
            testCase.verifyTrue(~isempty(gtv_locations{3, 1, 1}), 'Primary GTV should be found for P03.');
            testCase.verifyTrue(contains(gtv_locations{3, 1, 1}, 'GTVp1.mat'));

            testCase.verifyTrue(~isempty(gtvn_locations{3, 1, 1}), 'Nodal GTV should be found for P03.');
            testCase.verifyTrue(contains(gtvn_locations{3, 1, 1}, 'GTV_LN1.mat'));
        end
    end
end
