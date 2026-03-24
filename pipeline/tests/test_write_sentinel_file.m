classdef test_write_sentinel_file < matlab.unittest.TestCase
% TEST_WRITE_SENTINEL_FILE — Unit tests for write_sentinel_file.m
%
% Validates sentinel file creation including:
%   - Correct file path construction
%   - Message content written to file
%   - Different prefix and dwi_type combinations
%   - Graceful warning on write failure (invalid path)

    properties
        TempDir
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
        end

        function createTempDir(testCase)
            testCase.TempDir = fullfile(tempdir, ['test_sentinel_' strrep(tempname, tempdir, '')]);
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
        function test_creates_file_with_correct_name(testCase)
            write_sentinel_file(testCase.TempDir, 'metrics_baseline', ...
                'Baseline done.', 'Standard');

            expected = fullfile(testCase.TempDir, 'metrics_baseline_Standard.txt');
            testCase.verifyTrue(exist(expected, 'file') > 0, ...
                'Sentinel file should be created with prefix_dwitype.txt name.');
        end

        function test_file_contains_message(testCase)
            msg = 'Longitudinal metrics generated successfully.';
            write_sentinel_file(testCase.TempDir, 'metrics_longitudinal', ...
                msg, 'dnCNN');

            filepath = fullfile(testCase.TempDir, 'metrics_longitudinal_dnCNN.txt');
            content = fileread(filepath);
            testCase.verifySubstring(content, msg, ...
                'Sentinel file should contain the specified message.');
        end

        function test_different_dwi_types(testCase)
            write_sentinel_file(testCase.TempDir, 'step', 'msg1', 'Standard');
            write_sentinel_file(testCase.TempDir, 'step', 'msg2', 'dnCNN');
            write_sentinel_file(testCase.TempDir, 'step', 'msg3', 'IVIMnet');

            testCase.verifyTrue(exist(fullfile(testCase.TempDir, 'step_Standard.txt'), 'file') > 0);
            testCase.verifyTrue(exist(fullfile(testCase.TempDir, 'step_dnCNN.txt'), 'file') > 0);
            testCase.verifyTrue(exist(fullfile(testCase.TempDir, 'step_IVIMnet.txt'), 'file') > 0);
        end

        function test_overwrites_existing_sentinel(testCase)
            write_sentinel_file(testCase.TempDir, 'test_prefix', 'first', 'Standard');
            write_sentinel_file(testCase.TempDir, 'test_prefix', 'second', 'Standard');

            filepath = fullfile(testCase.TempDir, 'test_prefix_Standard.txt');
            content = fileread(filepath);
            testCase.verifySubstring(content, 'second');
            testCase.verifyTrue(~contains(content, 'first'), ...
                'Overwritten sentinel should only contain the new message.');
        end

        function test_warning_on_invalid_path(testCase)
            % Writing to a non-existent directory should warn, not error
            bad_path = fullfile(testCase.TempDir, 'nonexistent', 'subdir');

            testCase.verifyWarning( ...
                @() write_sentinel_file(bad_path, 'fail', 'msg', 'Standard'), ...
                'write_sentinel_file:fileWriteFailed');
        end

        function test_empty_message(testCase)
            write_sentinel_file(testCase.TempDir, 'empty_msg', '', 'Standard');

            filepath = fullfile(testCase.TempDir, 'empty_msg_Standard.txt');
            testCase.verifyTrue(exist(filepath, 'file') > 0, ...
                'Sentinel should be created even with empty message.');
        end
    end
end
