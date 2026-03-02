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
        TempDir
        OriginalPath
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
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

    end
end
