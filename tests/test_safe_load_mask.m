classdef test_safe_load_mask < matlab.unittest.TestCase
    % TEST_SAFE_LOAD_MASK
    %
    % Verify the security and functionality of the safe_load_mask utility.
    %
    % This test suite creates temporary .mat files with various content types
    % (valid numeric, valid logical, invalid strings/structs/objects) and
    % asserts that safe_load_mask correctly handles each case without executing
    % unsafe code or loading dangerous payloads.
    %
    % Run tests with:
    %   results = runtests('tests/test_safe_load_mask.m');
    % or
    %   mcp_matlab-mcp_run_matlab_file tests/test_safe_load_mask.m

    properties
        TempDir
        OriginalPath
    end

    methods(TestMethodSetup)
        function createTempEnv(testCase)
            % Create temporary directory for test files
            testCase.TempDir = fullfile(tempdir, ['test_safe_load_' datestr(now, 'yyyymmddHHMMSSFFF')]);
            if ~exist(testCase.TempDir, 'dir')
                mkdir(testCase.TempDir);
            end

            % Add necessary paths
            testCase.OriginalPath = path;
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function cleanTempEnv(testCase)
            % Clean up temporary directory
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end

            % Restore path
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testValidNumericMask(testCase)
            % Should successfully load a double array
            filename = fullfile(testCase.TempDir, 'valid_numeric.mat');
            Stvol3d = rand(5,5,5);
            save(filename, 'Stvol3d');

            loaded_mask = safe_load_mask(filename, 'Stvol3d');

            testCase.verifyEqual(loaded_mask, Stvol3d, ...
                'Failed to load valid numeric mask.');
        end

        function testValidLogicalMask(testCase)
            % Should successfully load a logical array
            filename = fullfile(testCase.TempDir, 'valid_logical.mat');
            Stvol3d = rand(5,5,5) > 0.5;
            save(filename, 'Stvol3d');

            loaded_mask = safe_load_mask(filename, 'Stvol3d');

            testCase.verifyEqual(loaded_mask, Stvol3d, ...
                'Failed to load valid logical mask.');
        end

        function testInvalidTypeString(testCase)
            % Should reject a string variable and issue a warning
            filename = fullfile(testCase.TempDir, 'invalid_string.mat');
            Stvol3d = "malicious_payload";
            save(filename, 'Stvol3d');

            % Verify warning is issued
            testCase.verifyWarning(@() safe_load_mask(filename, 'Stvol3d'), ...
                'safe_load_mask:SecurityRisk', ...
                'Should trigger warning for unsafe string type.');

            % Verify output is empty
            loaded_mask = safe_load_mask(filename, 'Stvol3d');
            testCase.verifyEmpty(loaded_mask, ...
                'Should return empty for unsafe string type.');
        end

        function testInvalidTypeStruct(testCase)
            % Should reject a struct variable and issue a warning
            filename = fullfile(testCase.TempDir, 'invalid_struct.mat');
            Stvol3d = struct('field', 'value');
            save(filename, 'Stvol3d');

            % Verify warning is issued
            testCase.verifyWarning(@() safe_load_mask(filename, 'Stvol3d'), ...
                'safe_load_mask:SecurityRisk', ...
                'Should trigger warning for unsafe struct type.');

            % Verify output is empty
            loaded_mask = safe_load_mask(filename, 'Stvol3d');
            testCase.verifyEmpty(loaded_mask, ...
                'Should return empty for unsafe struct type.');
        end

        function testVariableNotFound(testCase)
            % Should return empty if variable name is incorrect
            filename = fullfile(testCase.TempDir, 'missing_var.mat');
            DifferentVar = rand(3,3);
            save(filename, 'DifferentVar');

            loaded_mask = safe_load_mask(filename, 'Stvol3d');
            testCase.verifyEmpty(loaded_mask, ...
                'Should return empty if target variable is missing.');
        end

        function testFileNotFound(testCase)
             % Should return empty if file does not exist
             filename = fullfile(testCase.TempDir, 'non_existent.mat');
             loaded_mask = safe_load_mask(filename, 'Stvol3d');
             testCase.verifyEmpty(loaded_mask, ...
                 'Should return empty if file does not exist.');
        end

    end
end
