classdef test_convert_dicom < matlab.unittest.TestCase
    % TEST_CONVERT_DICOM Unit tests for the convert_dicom function.
    %
    % convert_dicom wraps the external dcm2niix tool to convert DICOM
    % directories into NIfTI + bval/bvec files.  These tests use mock
    % shell scripts in place of the real dcm2niix binary to control what
    % output files are created, allowing verification of:
    %   - Success when exactly the expected 3 output files are generated
    %   - Failure detection when the wrong number of files appear
    %   - Short-circuit skip when the output NIfTI already exists

    properties
        TestDir      % Root temp directory for this test run
        DicomLoc     % Simulated DICOM input directory
        OutLoc       % NIfTI output directory
        MockScript   % Path to the mock dcm2niix shell script
    end

    methods (TestMethodSetup)
        function createTempDirs(testCase)
            % Set up isolated temp directories and a mock dcm2niix script.
            % The mock script is a minimal shell/batch script that can be
            % overwritten by individual tests to produce specific file counts.
            testCase.TestDir = tempname;
            mkdir(testCase.TestDir);

            testCase.DicomLoc = fullfile(testCase.TestDir, 'dicom');
            mkdir(testCase.DicomLoc);

            testCase.OutLoc = fullfile(testCase.TestDir, 'out');
            mkdir(testCase.OutLoc);

            % Create a platform-appropriate mock script that exits cleanly
            if ispc
                testCase.MockScript = fullfile(testCase.TestDir, 'mock_dcm2niix.bat');
                % A batch script that we can overwrite dynamically
                fid = fopen(testCase.MockScript, 'w');
                fprintf(fid, '@echo off\nexit /b 0\n');
                fclose(fid);
            else
                testCase.MockScript = fullfile(testCase.TestDir, 'mock_dcm2niix.sh');
                fid = fopen(testCase.MockScript, 'w');
                fprintf(fid, '#!/bin/bash\nexit 0\n');
                fclose(fid);
                system(['chmod +x ' testCase.MockScript]);
            end
        end
    end

    methods (TestMethodTeardown)
        function removeTempDirs(testCase)
            % Remove the entire temp tree to avoid leftover mock scripts/files.
            if exist(testCase.TestDir, 'dir')
                rmdir(testCase.TestDir, 's');
            end
        end
    end

    methods (Test)

        function test_success_three_files(testCase)
            % Verifies the success path: the mock script creates the 3
            % expected output files (.nii.gz, .bval, .bvec).  convert_dicom
            % counts the directory listing delta (3 new files + 1 dir entry
            % for '.' = 4) and should return bad_dwi_found = 0.
            scanID = 'scan_001';
            fx_id = 'fx_test';

            % Update mock script to generate the 3 expected output files (.nii.gz, .bval, .bvec)
            fid = fopen(testCase.MockScript, 'w');
            if ispc
                fprintf(fid, '@echo off\n');
                fprintf(fid, 'echo dummy > "%s\\%s.nii.gz"\n', testCase.OutLoc, scanID);
                fprintf(fid, 'echo dummy > "%s\\%s.bval"\n', testCase.OutLoc, scanID);
                fprintf(fid, 'echo dummy > "%s\\%s.bvec"\n', testCase.OutLoc, scanID);
            else
                fprintf(fid, '#!/bin/bash\n');
                fprintf(fid, 'touch "%s/%s.nii.gz"\n', testCase.OutLoc, scanID);
                fprintf(fid, 'touch "%s/%s.bval"\n', testCase.OutLoc, scanID);
                fprintf(fid, 'touch "%s/%s.bvec"\n', testCase.OutLoc, scanID);
            end
            fclose(fid);
            if ~ispc
                system(['chmod +x ' testCase.MockScript]);
            end

            bad_dwi_found = convert_dicom(testCase.DicomLoc, testCase.OutLoc, scanID, testCase.MockScript, fx_id);

            testCase.verifyEqual(bad_dwi_found, 0, 'Should return 0 when dir length increases by exactly 4.');
        end

        function test_failure_wrong_files(testCase)
            % Verifies failure detection: the mock script creates 3 generic
            % .txt files instead of the expected .nii.gz/.bval/.bvec triple.
            % The directory delta is 3 (not the required 4), so
            % convert_dicom should return bad_dwi_found = 1.
            scanID = 'scan_002';
            fx_id = 'fx_test';

            % Update mock script to generate only 3 files delta in OutLoc
            fid = fopen(testCase.MockScript, 'w');
            if ispc
                fprintf(fid, '@echo off\n');
                fprintf(fid, 'echo dummy > "%s\\file1.txt"\n', testCase.OutLoc);
                fprintf(fid, 'echo dummy > "%s\\file2.txt"\n', testCase.OutLoc);
                fprintf(fid, 'echo dummy > "%s\\file3.txt"\n', testCase.OutLoc);
            else
                fprintf(fid, '#!/bin/bash\n');
                fprintf(fid, 'touch "%s/file1.txt"\n', testCase.OutLoc);
                fprintf(fid, 'touch "%s/file2.txt"\n', testCase.OutLoc);
                fprintf(fid, 'touch "%s/file3.txt"\n', testCase.OutLoc);
            end
            fclose(fid);
            if ~ispc
                system(['chmod +x ' testCase.MockScript]);
            end

            warnState = warning('off', 'convert_dicom:missingFiles');
            restoreWarn = onCleanup(@() warning(warnState));
            bad_dwi_found = convert_dicom(testCase.DicomLoc, testCase.OutLoc, scanID, testCase.MockScript, fx_id);

            testCase.verifyEqual(bad_dwi_found, 1, 'Should return 1 when expected files are not created.');
        end

        function test_skip_existing_file(testCase)
            % Test when the primary .nii.gz file already exists
            scanID = 'scan_003';
            fx_id = 'fx_test';

            % Create the output file beforehand
            fid = fopen(fullfile(testCase.OutLoc, [scanID '.nii.gz']), 'w');
            fprintf(fid, 'dummy content');
            fclose(fid);

            % Update mock script to generate 10 files (which would fail if executed)
            fid = fopen(testCase.MockScript, 'w');
            if ispc
                fprintf(fid, '@echo off\n');
                for i = 1:10
                    fprintf(fid, 'echo dummy > "%s\\file%d.txt"\n', testCase.OutLoc, i);
                end
            else
                fprintf(fid, '#!/bin/bash\n');
                for i = 1:10
                    fprintf(fid, 'touch "%s/file%d.txt"\n', testCase.OutLoc, i);
                end
            end
            fclose(fid);
            if ~ispc
                system(['chmod +x ' testCase.MockScript]);
            end

            bad_dwi_found = convert_dicom(testCase.DicomLoc, testCase.OutLoc, scanID, testCase.MockScript, fx_id);

            % Should return 0 and skip execution because file exists
            testCase.verifyEqual(bad_dwi_found, 0, 'Should return 0 and skip generation if file already exists.');
        end

    end

end
