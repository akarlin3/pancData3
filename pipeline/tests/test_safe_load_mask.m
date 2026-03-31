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
            % Verifies that a standard double-precision 3D array (the typical
            % GTV mask format) loads successfully without any warnings.
            filename = fullfile(testCase.TempDir, 'valid_numeric.mat');
            Stvol3d = rand(8,8,8);
            save(filename, 'Stvol3d');

            loaded_mask = safe_load_mask(filename, 'Stvol3d');

            testCase.verifyEqual(loaded_mask, Stvol3d, ...
                'Failed to load valid numeric mask.');
        end

        function testValidLogicalMask(testCase)
            % Verifies that logical arrays (binary masks) are accepted.
            % Many GTV masks are stored as logical to save space.
            filename = fullfile(testCase.TempDir, 'valid_logical.mat');
            Stvol3d = rand(8,8,8) > 0.5;
            save(filename, 'Stvol3d');

            loaded_mask = safe_load_mask(filename, 'Stvol3d');

            testCase.verifyEqual(loaded_mask, Stvol3d, ...
                'Failed to load valid logical mask.');
        end

        function testInvalidTypeString(testCase)
            % Security test: a string variable stored as 'Stvol3d' could be
            % crafted to exploit eval-based loaders. safe_load_mask must
            % reject non-numeric types with a SecurityRisk warning and
            % return empty instead of the string content.
            filename = fullfile(testCase.TempDir, 'invalid_string.mat');
            Stvol3d = repmat("malicious_payload", 8, 8, 8);
            save(filename, 'Stvol3d');

            % Verify warning is issued
            testCase.verifyWarning(@() safe_load_mask(filename, 'Stvol3d'), ...
                'safe_load_mask:SecurityRisk', ...
                'Should trigger warning for unsafe string type.');

            % Verify output is empty
            w = warning('off', 'safe_load_mask:SecurityRisk');
            loaded_mask = safe_load_mask(filename, 'Stvol3d');
            warning(w);
            testCase.verifyEmpty(loaded_mask, ...
                'Should return empty for unsafe string type.');
        end

        function testInvalidTypeStruct(testCase)
            % Security test: a struct stored under the expected variable name
            % could contain function handles or other executable objects.
            % safe_load_mask must reject structs with SecurityRisk warning.
            filename = fullfile(testCase.TempDir, 'invalid_struct.mat');
            Stvol3d = repmat(struct('field', 'value'), 8, 8, 8);
            save(filename, 'Stvol3d');

            % Verify warning is issued
            testCase.verifyWarning(@() safe_load_mask(filename, 'Stvol3d'), ...
                'safe_load_mask:SecurityRisk', ...
                'Should trigger warning for unsafe struct type.');

            % Verify output is empty
            w = warning('off', 'safe_load_mask:SecurityRisk');
            loaded_mask = safe_load_mask(filename, 'Stvol3d');
            warning(w);
            testCase.verifyEmpty(loaded_mask, ...
                'Should return empty for unsafe struct type.');
        end

        function testVariableNotFound(testCase)
            % Tests the case where the .mat file exists but does not contain
            % the requested variable name ('Stvol3d'). This can happen if the
            % mask was saved under a non-standard variable name.
            filename = fullfile(testCase.TempDir, 'missing_var.mat');
            DifferentVar = rand(3,3);
            save(filename, 'DifferentVar');

            w = warning('off', 'safe_load_mask:VariableNotFound');
            loaded_mask = safe_load_mask(filename, 'Stvol3d');
            warning(w);
            testCase.verifyEmpty(loaded_mask, ...
                'Should return empty if target variable is missing.');
        end

        function testFileNotFound(testCase)
             % Tests graceful handling of a non-existent file path.
             % In production, missing mask files occur when patients lack
             % GTV contours; the pipeline must continue without crashing.
             filename = fullfile(testCase.TempDir, 'non_existent.mat');
             w = warning('off', 'safe_load_mask:FileNotFound');
             loaded_mask = safe_load_mask(filename, 'Stvol3d');
             warning(w);
             testCase.verifyEmpty(loaded_mask, ...
                 'Should return empty if file does not exist.');
        end

        function testCorruptedMatFile(testCase)
            % A truncated/corrupted .mat file should not crash;
            % safe_load_mask should return empty gracefully.
            filename = fullfile(testCase.TempDir, 'corrupted.mat');
            fid = fopen(filename, 'w');
            fwrite(fid, 'NOT_A_VALID_MAT_FILE_HEADER');
            fclose(fid);

            w = warning('off', 'safe_load_mask:FileReadError');
            loaded_mask = safe_load_mask(filename, 'Stvol3d');
            warning(w);
            testCase.verifyEmpty(loaded_mask, ...
                'Should return empty for corrupted .mat file.');
        end

        function testMultiVariableMatFile(testCase)
            % In practice, .mat files may contain multiple variables (e.g.,
            % both a mask and metadata). safe_load_mask must extract only
            % the named variable and ignore others, even if some of those
            % others are unsafe types (like 'AnotherVar' which is a string).
            filename = fullfile(testCase.TempDir, 'multi_var.mat');
            OtherVar = rand(8, 8);
            Stvol3d = rand(8, 8, 8);
            AnotherVar = 'hello';
            save(filename, 'OtherVar', 'Stvol3d', 'AnotherVar');

            loaded_mask = safe_load_mask(filename, 'Stvol3d');
            testCase.verifyEqual(loaded_mask, Stvol3d, ...
                'Should load the correct variable from multi-variable .mat file.');
        end

        function testValidIntegerMask(testCase)
            % Integer types (uint8, int16, etc.) are commonly used for masks
            % from third-party segmentation tools. They must be accepted as
            % valid numeric types alongside double and logical.
            filename = fullfile(testCase.TempDir, 'valid_uint8.mat');
            Stvol3d = uint8(randi([0 1], 8, 8, 8));
            save(filename, 'Stvol3d');

            loaded_mask = safe_load_mask(filename, 'Stvol3d');
            testCase.verifyEqual(loaded_mask, Stvol3d, ...
                'Should successfully load uint8 mask.');
        end

    end
end
