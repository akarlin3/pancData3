classdef test_discover_gtv_file < matlab.unittest.TestCase
    % TEST_DISCOVER_GTV_FILE Unit tests for discover_gtv_file logic

    properties
        TempDir
    end

    methods(TestMethodSetup)
        function setup(testCase)
            % Create a temporary directory for dummy files
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            % Add utils directory to path
            current_dir = fileparts(mfilename('fullpath'));
            utils_dir = fullfile(current_dir, '..', 'utils');
            addpath(utils_dir);
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            % Remove the temporary directory and all its contents
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testExactMatchWithIndex(testCase)
            % Create dummy files
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
            % Create dummy files without index
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
            % Directory is empty, should return empty string
            filepath = discover_gtv_file(testCase.TempDir, {'GTVp'}, 1);
            testCase.verifyEqual(filepath, '');
        end

        function testMultiplePatterns(testCase)
            % Create dummy files matching the second pattern
            fid = fopen(fullfile(testCase.TempDir, 'GTV_node1.mat'), 'w'); fclose(fid);

            % Test multiple patterns, should match the second one
            filepath = discover_gtv_file(testCase.TempDir, {'GTVp', 'GTV_node'}, 1);
            expectedPath = fullfile(testCase.TempDir, 'GTV_node1.mat');
            testCase.verifyEqual(filepath, expectedPath);
        end

        function testCharInput(testCase)
            % Test that char input for patterns works
            fid = fopen(fullfile(testCase.TempDir, 'GTVp1_MR.mat'), 'w'); fclose(fid);

            filepath = discover_gtv_file(testCase.TempDir, 'GTVp', 1);
            expectedPath = fullfile(testCase.TempDir, 'GTVp1_MR.mat');
            testCase.verifyEqual(filepath, expectedPath);
        end
    end
end
