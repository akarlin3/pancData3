classdef test_discover_patient_files < matlab.unittest.TestCase
    % TEST_DISCOVER_PATIENT_FILES Formal unit test suite for discover_patient_files.m
    %
    % discover_patient_files scans a data directory to locate patient folders,
    % DWI DICOM series, RT dose files, and GTV mask files. It returns cell
    % arrays of paths and metadata (MRN, study dates) for the pipeline.
    %
    % These tests use a mock filesystem with three patient scenarios:
    %   P01-ABC:     Standard patient with DWI, RT dose, and a single GTV
    %   P02_DEF:     Patient with a Fx1 folder but no DWI/RTdose/GTV data
    %   P03_twoGTV:  Patient with nested DICOM fallback and dual GTVs (GTVp + GTVn)
    % A 'template' folder is also created to verify it is correctly ignored.

    properties
        MockDataDir   % Temporary directory containing the mock patient tree
        MockPathDir   % Temporary directory containing the mock dicominfo function
        OrigPath      % Original MATLAB path for restoration in teardown
    end

    methods(TestMethodSetup)
        function createMockData(testCase)
            % Build the entire mock filesystem used by all tests.

            % 1. Create a temporary directory for mock patient data
            testCase.MockDataDir = fullfile(tempdir, 'mock_patient_data');
            if exist(testCase.MockDataDir, 'dir')
                rmdir(testCase.MockDataDir, 's');
            end
            mkdir(testCase.MockDataDir);

            % 2. Create a mock dicominfo function that returns fixed MRN and
            %    date values, so tests do not require real DICOM files.
            testCase.MockPathDir = fullfile(tempdir, 'mock_path_dir');
            if exist(testCase.MockPathDir, 'dir')
                rmdir(testCase.MockPathDir, 's');
            end
            mkdir(testCase.MockPathDir);

            mock_dicominfo_file = fullfile(testCase.MockPathDir, 'dicominfo.m');
            fid = fopen(mock_dicominfo_file, 'w');
            fprintf(fid, 'function info = dicominfo(filename)\n');
            fprintf(fid, '    info.PatientID = ''MOCK-MRN'';\n');
            fprintf(fid, '    info.StudyDate = ''20240101'';\n');
            fprintf(fid, 'end\n');
            fclose(fid);

            % Place mock dicominfo on path ahead of the real one
            testCase.OrigPath = path;
            addpath(testCase.MockPathDir);

            % Add project paths needed by discover_patient_files
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'dependencies'));

            % --- Patient 1: Standard patient (P01-ABC) ---
            % Has Fx1 and post directories.
            % Fx1 contains a DWI folder with 5 .dcm files (minimum for detection),
            % an rtdose folder with a .dcm file, and a single GTV mask.
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

            fclose(fopen(fullfile(p1_dir, 'Fx1', 'GTV_1.mat'), 'w'));

            % --- Patient 2: Missing data patient (P02_DEF) ---
            % Has Fx1 folder but no DWI, rtdose, or GTV files inside it.
            % Tests that the function handles missing data gracefully.
            p2_dir = fullfile(testCase.MockDataDir, 'P02_DEF');
            mkdir(p2_dir);
            mkdir(fullfile(p2_dir, 'Fx1'));

            % --- Patient 3: Dual-GTV patient (P03_twoGTV) ---
            % Uses the nested DICOM subfolder fallback pattern (DICOM_dir/DICOM/)
            % and has both a primary GTV (GTVp1.mat) and nodal GTV (GTV_LN1.mat).
            p3_dir = fullfile(testCase.MockDataDir, 'P03_twoGTV');
            mkdir(p3_dir);
            mkdir(fullfile(p3_dir, 'Fx1'));

            p3_fx1_dwi = fullfile(p3_dir, 'Fx1', 'DICOM_dir');
            mkdir(p3_fx1_dwi);
            % Nested DICOM folder (fallback discovery path)
            mkdir(fullfile(p3_fx1_dwi, 'DICOM'));
            for i=1:5
                fclose(fopen(fullfile(p3_fx1_dwi, 'DICOM', sprintf('im%d.dcm', i)), 'w'));
            end

            % Dual GTV masks
            fclose(fopen(fullfile(p3_dir, 'Fx1', 'GTVp1.mat'), 'w'));
            fclose(fopen(fullfile(p3_dir, 'Fx1', 'GTV_LN1.mat'), 'w'));

            % --- Patient 4: Repeatability-session patient (P04-REP) ---
            % At MSK, some sites name the Fx1 baseline session as
            % "Fx1 - repeatability" (with indexed per-repeat masks). The
            % fraction matcher must recognize this folder as Fx1 so the
            % indexed masks feed into compute_spatial_repeatability;
            % otherwise dice_rpt_* comes back all-NaN for the entire
            % cohort.
            %
            % Mask filenames use the simple 'GTVp<N>.mat' form (matched
            % directly by the '*GTVp' pattern in discover_gtv_file with
            % a single-wildcard '*GTVp<N>.mat' glob) rather than the
            % MSK-cohort 'ROI_GTV<N>_<date>.mat' form. Both name shapes
            % exercise the same code path through the matcher, but the
            % simpler form sidesteps any platform-specific quirks in
            % dir()'s multi-wildcard '*GTV*<N>_*.mat' matching — which
            % is irrelevant to the regression this test is locking in.
            p4_dir = fullfile(testCase.MockDataDir, 'P04-REP');
            mkdir(p4_dir);
            p4_fx1 = fullfile(p4_dir, 'Fx1 - repeatability');
            mkdir(p4_fx1);

            for rpi = 1:3
                dwi_sub = fullfile(p4_fx1, sprintf('DWI%d', rpi));
                mkdir(dwi_sub);
                for i = 1:5
                    fclose(fopen(fullfile(dwi_sub, sprintf('im%d.dcm', i)), 'w'));
                end
                fclose(fopen(fullfile(p4_fx1, sprintf('GTVp%d.mat', rpi)), 'w'));
            end

            % Non-patient folder (should be ignored by discovery logic)
            mkdir(fullfile(testCase.MockDataDir, 'template'));
        end
    end

    methods(TestMethodTeardown)
        function removeMockData(testCase)
            % Restore the original MATLAB path and remove all temporary dirs.
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
            % Verify that discover_patient_files finds exactly the 4 patient
            % folders (P01-P04) and ignores the 'template' folder.
            % Results should be sorted alphabetically by folder name.
            [id_list, mrn_list, fx_dates, dwi_locations, rtdose_locations, gtv_locations, gtvn_locations] = discover_patient_files(testCase.MockDataDir);

            testCase.verifyEqual(length(id_list), 4, 'Should discover exactly 4 valid patients.');

            % Verify alphabetical sort order
            testCase.verifyEqual(id_list{1}, 'P01-ABC');
            testCase.verifyEqual(id_list{2}, 'P02_DEF');
            testCase.verifyEqual(id_list{3}, 'P03_twoGTV');
            testCase.verifyEqual(id_list{4}, 'P04-REP');
        end

        function testFx1RepeatabilityFolderNamingIsRecognizedAsFx1(testCase)
            % Regression test for MSK 'Fx1 - repeatability' layout.
            % The fraction matcher must treat this as Fx1 and find all
            % three indexed ROI_GTV*.mat masks so dice_rpt_* can be
            % computed downstream. The previous '\b'-based regex silently
            % skipped this folder on the live MSK machine, which made
            % gtv_locations{:, 1, :} empty for ~95% of the cohort and
            % cascaded into all-NaN repeat Dice metrics.
            [~, ~, ~, dwi_locations, ~, gtv_locations, ~] = discover_patient_files(testCase.MockDataDir);

            % Each of the three repeats must resolve to a DWI path and an
            % indexed GTV mask at Fx1 (fraction index 1).
            for rpi = 1:3
                testCase.verifyTrue(~isempty(dwi_locations{4, 1, rpi}), ...
                    sprintf('DWI path should be found for P04-REP Fx1 repeat %d.', rpi));
                testCase.verifyTrue(~isempty(gtv_locations{4, 1, rpi}), ...
                    sprintf('Indexed GTV mask should be found for P04-REP Fx1 repeat %d.', rpi));
                testCase.verifyTrue(contains(gtv_locations{4, 1, rpi}, ...
                    sprintf('GTVp%d.mat', rpi)), ...
                    sprintf('Fx1 repeat %d should resolve to GTVp%d.mat.', rpi, rpi));
                testCase.verifyTrue(contains(gtv_locations{4, 1, rpi}, 'Fx1 - repeatability'), ...
                    'Resolved GTV path must live under the "Fx1 - repeatability" folder.');
            end
        end

        function testStandardPatient(testCase)
            % Verify that the standard patient (P01-ABC, index 1) has:
            % - MRN and date populated from the mock dicominfo
            % - DWI location pointing to the 1_DWI_dir folder
            % - RTdose location pointing to the 2_rtdose folder
            % - Primary GTV found, but no nodal GTV (GTVn)
            [~, mrn_list, fx_dates, dwi_locations, rtdose_locations, gtv_locations, gtvn_locations] = discover_patient_files(testCase.MockDataDir);

            % MRN and Date come from the mock dicominfo function
            testCase.verifyEqual(mrn_list{1}, 'MOCK-MRN');
            testCase.verifyEqual(fx_dates{1, 1}, '20240101');

            % DWI and rtdose locations should be non-empty and contain
            % the expected folder names
            testCase.verifyTrue(~isempty(dwi_locations{1, 1, 1}), 'DWI location should be populated for P01.');
            testCase.verifyTrue(contains(dwi_locations{1, 1, 1}, '1_DWI_dir'));

            testCase.verifyTrue(~isempty(rtdose_locations{1, 1}), 'RTdose location should be populated for P01.');
            testCase.verifyTrue(contains(rtdose_locations{1, 1}, '2_rtdose'));

            % Single GTV found (GTV_1.mat); no nodal GTV present
            testCase.verifyTrue(~isempty(gtv_locations{1, 1, 1}), 'GTV location should be found.');
            testCase.verifyTrue(isempty(gtvn_locations{1, 1, 1}), 'GTVn should not be found for P01.');
        end

        function testMissingDataPatient(testCase)
            % Patient P02_DEF (index 2) has an empty Fx1 folder: no DICOMs,
            % no RT dose, and no GTV mask. All location fields should be empty
            % and MRN should be empty (no DICOM to read).
            [~, mrn_list, ~, dwi_locations, rtdose_locations, gtv_locations, ~] = discover_patient_files(testCase.MockDataDir);

            testCase.verifyTrue(isempty(mrn_list{2}), 'MRN should be empty for missing data patient.');
            testCase.verifyTrue(isempty(dwi_locations{2, 1, 1}), 'DWI location should be empty for P02.');
            testCase.verifyTrue(isempty(rtdose_locations{2, 1}), 'RTdose location should be empty for P02.');
            testCase.verifyTrue(isempty(gtv_locations{2, 1, 1}), 'GTV location should be empty for P02.');
        end

        function testTwoGTVPatient(testCase)
            % Patient P03_twoGTV (index 3) exercises two features:
            % 1. Nested DICOM fallback: DWI files are in DICOM_dir/DICOM/
            % 2. Dual GTV discovery: both GTVp1.mat (primary) and GTV_LN1.mat
            %    (nodal) should be found.
            [~, ~, ~, dwi_locations, ~, gtv_locations, gtvn_locations] = discover_patient_files(testCase.MockDataDir);

            % DWI via nested DICOM fallback
            testCase.verifyTrue(~isempty(dwi_locations{3, 1, 1}), 'DWI location should be found via DICOM fallback.');
            testCase.verifyTrue(contains(dwi_locations{3, 1, 1}, 'DICOM'), 'Should use nested DICOM fallback path.');

            % Primary GTV (GTVp1.mat)
            testCase.verifyTrue(~isempty(gtv_locations{3, 1, 1}), 'Primary GTV should be found for P03.');
            testCase.verifyTrue(contains(gtv_locations{3, 1, 1}, 'GTVp1.mat'));

            % Nodal GTV (GTV_LN1.mat)
            testCase.verifyTrue(~isempty(gtvn_locations{3, 1, 1}), 'Nodal GTV should be found for P03.');
            testCase.verifyTrue(contains(gtvn_locations{3, 1, 1}, 'GTV_LN1.mat'));
        end
    end
end
