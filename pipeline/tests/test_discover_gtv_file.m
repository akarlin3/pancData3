classdef test_discover_gtv_file < matlab.unittest.TestCase
    % TEST_DISCOVER_GTV_FILE Unit tests for discover_gtv_file.
    %
    % discover_gtv_file locates a GTV (Gross Tumor Volume) mask .mat file
    % in a patient directory using pattern matching and a DWI scan index.
    % It supports indexed filenames (e.g., GTVp1_MR.mat, GTVp2_MR.mat),
    % fallback to a single unindexed file, multiple search patterns tried
    % in order, and char vs cell pattern input.
    %
    % Tests:
    %   - Exact match with numeric index suffix
    %   - Fallback when only one file matches (regardless of requested index)
    %   - Empty directory returns empty string
    %   - Multiple patterns tried in order until a match is found
    %   - Char input (not cell) for patterns works correctly

    properties
        TempDir   % Isolated temp directory populated with dummy .mat files
    end

    methods(TestMethodSetup)
        function setup(testCase)
            % Create an empty temp directory and add utils/ to the path.
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            current_dir = fileparts(mfilename('fullpath'));
            utils_dir = fullfile(current_dir, '..', 'utils');
            addpath(utils_dir);
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            % Clean up the temp directory and all dummy files.
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testExactMatchWithIndex(testCase)
            % Two indexed GTV files exist (GTVp1_MR.mat, GTVp2_MR.mat).
            % Requesting index 1 should return the first; index 2 the second.
            fid = fopen(fullfile(testCase.TempDir, 'GTVp1_MR.mat'), 'w'); fclose(fid);
            fid = fopen(fullfile(testCase.TempDir, 'GTVp2_MR.mat'), 'w'); fclose(fid);

            % Test exact match with index 1
            filepath = discover_gtv_file(testCase.TempDir, {'GTVp'}, 1);
            expectedPath = fullfile(testCase.TempDir, 'GTVp1_MR.mat');
            testCase.verifyEqual(filepath, expectedPath);

            % Test exact match with index 2
            filepath = discover_gtv_file(testCase.TempDir, {'GTVp'}, 2);
            expectedPath = fullfile(testCase.TempDir, 'GTVp2_MR.mat');
            testCase.verifyEqual(filepath, expectedPath);
        end

        function testFallbackSingleMatch(testCase)
            % When only one file matches the pattern (no numeric index in
            % the filename), that file is returned regardless of which
            % index is requested.  This handles legacy naming conventions.
            fid = fopen(fullfile(testCase.TempDir, 'GTVp_MR.mat'), 'w'); fclose(fid);
            fid = fopen(fullfile(testCase.TempDir, 'other_file.mat'), 'w'); fclose(fid);

            % Test fallback single match when asking for index 1
            filepath = discover_gtv_file(testCase.TempDir, {'GTVp_MR'}, 1);
            expectedPath = fullfile(testCase.TempDir, 'GTVp_MR.mat');
            testCase.verifyEqual(filepath, expectedPath);

            % Test fallback single match when asking for index 2
            filepath = discover_gtv_file(testCase.TempDir, {'GTVp_MR'}, 2);
            testCase.verifyEqual(filepath, expectedPath);
        end

        function testEmptyMatchResults(testCase)
            % When the directory contains no files matching any pattern,
            % the function should return an empty string (not error).
            filepath = discover_gtv_file(testCase.TempDir, {'GTVp'}, 1);
            testCase.verifyEqual(filepath, '');
        end

        function testMultiplePatterns(testCase)
            % When multiple patterns are provided, the function tries them
            % in order.  Here only a 'GTV_node' file exists (not 'GTVp'),
            % so the second pattern should match.
            fid = fopen(fullfile(testCase.TempDir, 'GTV_node1.mat'), 'w'); fclose(fid);

            % Test multiple patterns, should match the second one
            filepath = discover_gtv_file(testCase.TempDir, {'GTVp', 'GTV_node'}, 1);
            expectedPath = fullfile(testCase.TempDir, 'GTV_node1.mat');
            testCase.verifyEqual(filepath, expectedPath);
        end

        function testCharInput(testCase)
            % Patterns can be passed as a plain char vector instead of a
            % cell array.  Verify that the function wraps it internally.
            fid = fopen(fullfile(testCase.TempDir, 'GTVp1_MR.mat'), 'w'); fclose(fid);

            filepath = discover_gtv_file(testCase.TempDir, 'GTVp', 1);
            expectedPath = fullfile(testCase.TempDir, 'GTVp1_MR.mat');
            testCase.verifyEqual(filepath, expectedPath);
        end
    end
end
