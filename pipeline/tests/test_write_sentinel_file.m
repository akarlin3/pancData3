classdef test_write_sentinel_file < matlab.unittest.TestCase
% TEST_WRITE_SENTINEL_FILE — Unit tests for write_sentinel_file.m
%
% Validates sentinel file creation:
%   - File is created with correct name pattern
%   - File contains the expected message
%   - DWI type name is included in filename
%   - Warning on invalid output folder (no crash)

    properties
        TempDir
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            addpath(fullfile(testDir, '..', 'utils'));
        end

        function setupTempDir(testCase)
            testCase.TempDir = fullfile(tempdir, ['test_sentinel_' char(datetime('now','Format','yyyyMMddHHmmss'))]);
            mkdir(testCase.TempDir);
        end
    end

    methods (TestMethodTeardown)
        function cleanupTempDir(testCase)
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods (Test)
        function test_file_created(testCase)
            write_sentinel_file(testCase.TempDir, 'metrics_done', 'All good.', 'Standard');

            expected = fullfile(testCase.TempDir, 'metrics_done_Standard.txt');
            testCase.verifyTrue(isfile(expected), ...
                'Sentinel file should be created.');
        end

        function test_file_content(testCase)
            msg = 'Metrics computed successfully.';
            write_sentinel_file(testCase.TempDir, 'step_complete', msg, 'dnCNN');

            filepath = fullfile(testCase.TempDir, 'step_complete_dnCNN.txt');
            content = fileread(filepath);
            testCase.verifyTrue(contains(content, msg), ...
                'Sentinel file should contain the message.');
        end

        function test_different_dwi_types(testCase)
            write_sentinel_file(testCase.TempDir, 'done', 'OK', 'Standard');
            write_sentinel_file(testCase.TempDir, 'done', 'OK', 'dnCNN');
            write_sentinel_file(testCase.TempDir, 'done', 'OK', 'IVIMnet');

            testCase.verifyTrue(isfile(fullfile(testCase.TempDir, 'done_Standard.txt')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempDir, 'done_dnCNN.txt')));
            testCase.verifyTrue(isfile(fullfile(testCase.TempDir, 'done_IVIMnet.txt')));
        end

        function test_invalid_folder_warns(testCase)
            % Writing to a non-existent folder should warn, not crash
            testCase.verifyWarning(@() ...
                write_sentinel_file('/nonexistent/path/xyz', 'test', 'msg', 'Standard'), ...
                'write_sentinel_file:fileWriteFailed');
        end

        function test_overwrite_existing(testCase)
            % Writing twice should overwrite the existing file
            write_sentinel_file(testCase.TempDir, 'redo', 'First', 'Standard');
            write_sentinel_file(testCase.TempDir, 'redo', 'Second', 'Standard');

            filepath = fullfile(testCase.TempDir, 'redo_Standard.txt');
            content = fileread(filepath);
            testCase.verifyTrue(contains(content, 'Second'));
            testCase.verifyFalse(contains(content, 'First'));
        end
    end
end
