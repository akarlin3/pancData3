classdef test_find_gtv_files < matlab.unittest.TestCase
    % TEST_FIND_GTV_FILES Unit tests for find_gtv_files.
    %
    % find_gtv_files delegates to discover_gtv_file with different pattern
    % sets depending on whether the patient name contains "two".
    %
    % Tests:
    %   - Single-GTV patient (no "two"): gtvp_path found, gtvn_path empty
    %   - Dual-GTV patient ("two" in name): both GTVp and GTVn paths found
    %   - Dual-GTV patient: only GTVp file present, GTVn is empty string
    %   - Single-GTV patient: no matching file → gtvp_path is empty string
    %   - dwii index is forwarded to discover_gtv_file correctly

    properties
        TempDir        % Isolated temp directory for dummy GTV mask files
        OriginalPath   % Saved MATLAB path, restored in teardown
    end

    methods(TestMethodSetup)
        function setup(testCase)
            % Create isolated temp dir and add utils/ to path.
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            % Remove temp files and restore the original MATLAB path to
            % prevent cross-test side effects.
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testSingleGtvPatientFindsGtvp(testCase)
            % Patient name without "two" → single-GTV branch.
            % Pattern used: '*GTV*' + dwii.
            % File 'GTV_panc1.mat' matches '*GTV*1*.mat' (dwii = 1).
            fid = fopen(fullfile(testCase.TempDir, 'GTV_panc1.mat'), 'w');
            fclose(fid);

            [gtvp_path, gtvn_path] = find_gtv_files(testCase.TempDir, 1, 'Patient_A');

            testCase.verifyNotEmpty(gtvp_path, ...
                'gtvp_path should be found for a single-GTV patient.');
            testCase.verifyEmpty(gtvn_path, ...
                'gtvn_path should be empty for a single-GTV patient.');
        end

        function testSingleGtvPatientNoFileReturnsEmpty(testCase)
            % No GTV file in the folder → gtvp_path should be ''.
            [gtvp_path, gtvn_path] = find_gtv_files(testCase.TempDir, 1, 'PatientB');

            testCase.verifyEqual(gtvp_path, '', ...
                'gtvp_path should be empty string when no file matches.');
            testCase.verifyEqual(gtvn_path, '', ...
                'gtvn_path should be empty string for single-GTV patient.');
        end

        function testDualGtvPatientFindsBothPaths(testCase)
            % Patient name contains "two" → dual-GTV branch.
            % GTVp patterns: '*GTV_MR', '*GTVp', '*GTV_panc*'
            % GTVn patterns: '*GTV*LN', '*GTVn', '*GTV_node*'
            %
            % 'GTV_panc1.mat' matches '*GTV_panc*1*.mat' (GTVp, dwii=1)
            % 'GTVn1.mat'     matches '*GTVn1*.mat'       (GTVn, dwii=1)
            fid = fopen(fullfile(testCase.TempDir, 'GTV_panc1.mat'), 'w'); fclose(fid);
            fid = fopen(fullfile(testCase.TempDir, 'GTVn1.mat'), 'w');     fclose(fid);

            [gtvp_path, gtvn_path] = find_gtv_files(testCase.TempDir, 1, 'Patient_two_lesions');

            testCase.verifyNotEmpty(gtvp_path, ...
                'gtvp_path should be found for a dual-GTV patient.');
            testCase.verifyNotEmpty(gtvn_path, ...
                'gtvn_path should be found for a dual-GTV patient when GTVn file exists.');
        end

        function testDualGtvPatientMissingGtvnReturnsEmptyGtvn(testCase)
            % Only GTVp file exists; GTVn file absent → gtvn_path = ''.
            fid = fopen(fullfile(testCase.TempDir, 'GTV_panc1.mat'), 'w'); fclose(fid);

            [gtvp_path, gtvn_path] = find_gtv_files(testCase.TempDir, 1, 'two_region_patient');

            testCase.verifyNotEmpty(gtvp_path, 'gtvp_path should be found.');
            testCase.verifyEqual(gtvn_path, '', ...
                'gtvn_path should be empty string when no GTVn file exists.');
        end

        function testDwiIndexForwardedCorrectly(testCase)
            % For a single-GTV patient with two GTV files at different indices,
            % only the file matching the requested dwii is returned.
            fid = fopen(fullfile(testCase.TempDir, 'GTV_panc1.mat'), 'w'); fclose(fid);
            fid = fopen(fullfile(testCase.TempDir, 'GTV_panc2.mat'), 'w'); fclose(fid);

            [gtvp1, ~] = find_gtv_files(testCase.TempDir, 1, 'PatientC');
            [gtvp2, ~] = find_gtv_files(testCase.TempDir, 2, 'PatientC');

            testCase.verifyTrue(contains(gtvp1, '1'), ...
                'dwii=1 should select the file with index 1.');
            testCase.verifyTrue(contains(gtvp2, '2'), ...
                'dwii=2 should select the file with index 2.');
        end

        function testPatientNameSubstringTwoIsRespected(testCase)
            % "two" as a substring (not just the exact word) triggers the
            % dual-GTV branch.  Verify the GTVn patterns are searched.
            fid = fopen(fullfile(testCase.TempDir, 'GTV_LN1.mat'), 'w'); fclose(fid);

            [~, gtvn_path] = find_gtv_files(testCase.TempDir, 1, 'pat_two_sites');

            testCase.verifyNotEmpty(gtvn_path, ...
                'GTVn should be found when "two" appears anywhere in pat_name.');
        end

        function testGtvInNonIndexedSubfolder(testCase)
            % Some sites nest GTV masks in a "GTV/" subfolder shared across
            % repeats.  find_gtv_files should locate masks there when the Fx
            % folder itself is empty.
            subfolder = fullfile(testCase.TempDir, 'GTV');
            mkdir(subfolder);
            fid = fopen(fullfile(subfolder, 'ROI_GTV1_20210212.mat'), 'w'); fclose(fid);

            [gtvp_path, ~] = find_gtv_files(testCase.TempDir, 1, 'Patient_X');

            testCase.verifyNotEmpty(gtvp_path, ...
                'gtvp_path should be found inside a non-indexed GTV/ subfolder.');
            testCase.verifyTrue(contains(gtvp_path, 'GTV'), ...
                'Returned path should be inside the GTV subfolder.');
        end

        function testGtvInIndexedTimepointSubfolder(testCase)
            % Some sites nest masks in timepoint-specific subfolders like
            % GTVtimepoint1/, GTVtp2/, GTVfx3/.  find_gtv_files should only
            % pick up the folder whose trailing index matches dwii.
            mkdir(fullfile(testCase.TempDir, 'GTVtimepoint1'));
            mkdir(fullfile(testCase.TempDir, 'GTVtimepoint2'));
            fid = fopen(fullfile(testCase.TempDir, 'GTVtimepoint1', ...
                'ROI_GTV1_20210212.mat'), 'w'); fclose(fid);
            fid = fopen(fullfile(testCase.TempDir, 'GTVtimepoint2', ...
                'ROI_GTV2_20210212.mat'), 'w'); fclose(fid);

            [gtvp1, ~] = find_gtv_files(testCase.TempDir, 1, 'Patient_Y');
            [gtvp2, ~] = find_gtv_files(testCase.TempDir, 2, 'Patient_Y');

            testCase.verifyTrue(contains(gtvp1, 'timepoint1'), ...
                'dwii=1 should resolve to GTVtimepoint1/.');
            testCase.verifyTrue(contains(gtvp2, 'timepoint2'), ...
                'dwii=2 should resolve to GTVtimepoint2/.');
        end

        function testGtvOriPatternIsMatchedByBroadFallback(testCase)
            % ROI_GTVori1_<date>.mat is the most common naming variant in the
            % MSK cohort.  It must be caught by the '*GTV*' broad fallback
            % via the '*GTV*1_*.mat' glob, not missed.
            fid = fopen(fullfile(testCase.TempDir, 'ROI_GTVori1_20210212.mat'), 'w');
            fclose(fid);

            [gtvp_path, ~] = find_gtv_files(testCase.TempDir, 1, 'Patient_Z');

            testCase.verifyNotEmpty(gtvp_path, ...
                'Broad *GTV* fallback should match ROI_GTVori1_<date>.mat.');
        end

        function testFxFolderWinsOverSubfolder(testCase)
            % When matching masks exist both in fxfolder and in a GTV/
            % subfolder, fxfolder must take priority (candidate list order).
            fid = fopen(fullfile(testCase.TempDir, 'GTV_panc1.mat'), 'w'); fclose(fid);
            subfolder = fullfile(testCase.TempDir, 'GTV');
            mkdir(subfolder);
            fid = fopen(fullfile(subfolder, 'GTV_panc1.mat'), 'w'); fclose(fid);

            [gtvp_path, ~] = find_gtv_files(testCase.TempDir, 1, 'Patient_Q');

            testCase.verifyTrue(strcmp(fileparts(gtvp_path), testCase.TempDir), ...
                'fxfolder match should win over subfolder match.');
        end

    end
end
