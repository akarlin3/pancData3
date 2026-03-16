classdef test_setup_output_folders < matlab.unittest.TestCase
    % TEST_SETUP_OUTPUT_FOLDERS Unit tests for setup_output_folders.
    %
    % Validates master output folder creation: explicit reuse, timestamped
    % auto-creation with sentinel, and non-existent parent handling.

    properties
        TmpDir
    end

    methods(TestMethodSetup)
        function setupTempDir(testCase)
            testCase.TmpDir = tempname;
            mkdir(testCase.TmpDir);
        end
    end

    methods(TestMethodTeardown)
        function cleanupTempDir(testCase)
            if isfolder(testCase.TmpDir)
                rmdir(testCase.TmpDir, 's');
            end
        end
    end

    methods(Test)
        function test_explicit_folder_reuse(testCase)
            % When an explicit folder is provided, it should be returned as-is.
            explicit = fullfile(testCase.TmpDir, 'my_output');
            mkdir(explicit);
            master = setup_output_folders(testCase.TmpDir, explicit);
            testCase.verifyEqual(master, explicit);
        end

        function test_explicit_folder_created_if_missing(testCase)
            % If the explicit folder doesn't exist, it should be created.
            explicit = fullfile(testCase.TmpDir, 'new_folder');
            master = setup_output_folders(testCase.TmpDir, explicit);
            testCase.verifyTrue(isfolder(master));
            testCase.verifyEqual(master, explicit);
        end

        function test_auto_creates_timestamped_folder(testCase)
            % Empty master_output_folder triggers auto-creation with timestamp.
            master = setup_output_folders(testCase.TmpDir, '');
            testCase.verifyTrue(isfolder(master));
            % Folder name should contain 'saved_files_'
            [~, name] = fileparts(master);
            testCase.verifySubstring(name, 'saved_files_');
        end

        function test_auto_creates_sentinel(testCase)
            % Auto-created folder should have a .pipeline_created sentinel.
            master = setup_output_folders(testCase.TmpDir, '');
            sentinel = fullfile(master, '.pipeline_created');
            testCase.verifyTrue(exist(sentinel, 'file') == 2);
        end

        function test_explicit_folder_no_sentinel(testCase)
            % Explicit folder should NOT get a sentinel (not pipeline-created).
            explicit = fullfile(testCase.TmpDir, 'explicit');
            setup_output_folders(testCase.TmpDir, explicit);
            sentinel = fullfile(explicit, '.pipeline_created');
            testCase.verifyFalse(exist(sentinel, 'file') == 2);
        end
    end
end
