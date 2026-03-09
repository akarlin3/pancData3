classdef test_parse_config < matlab.unittest.TestCase
    % TEST_PARSE_CONFIG Unit tests for configuration parsing logic.
    %
    % parse_config reads a JSON config file, validates its contents, and
    % fills in default values for any missing fields using the isfield +
    % fallback pattern. These tests verify:
    %   - Missing fields receive correct default values
    %   - User-specified values are preserved (not overwritten by defaults)
    %   - Non-existent file, malformed JSON, and empty file produce errors
    %   - dwi_type string maps to the correct numeric dwi_types_to_run index
    %   - Invalid dwi_type strings are rejected
    %   - Missing dwi_type defaults to running all three types [1:3]
    %   - Boolean fields (clear_cache, use_checkpoints) default correctly

    properties
        TempConfigFile  % Path to a temporary JSON file created per test
    end

    methods(TestMethodSetup)
        function setup(testCase)
            % Create a unique temporary file path for the test config and
            % add the utils directory (which contains parse_config) to path.
            testCase.TempConfigFile = [tempname '.json'];

            current_dir = fileparts(mfilename('fullpath'));
            utils_dir = fullfile(current_dir, '..', 'utils');
            addpath(utils_dir);
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            % Clean up the temporary config file after each test.
            if exist(testCase.TempConfigFile, 'file')
                delete(testCase.TempConfigFile);
            end
        end
    end

    methods(Test)
        function testMissingFieldsAssignedDefaults(testCase)
            % A minimal config with only 'dataloc' should still produce a
            % fully populated struct with all default values filled in by
            % parse_config's isfield + fallback logic.
            minimal_config = struct();
            minimal_config.dataloc = '/tmp';

            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(minimal_config));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);

            % Verify key numeric defaults match expected values
            testCase.verifyEqual(config.adc_thresh, 0.001);
            testCase.verifyEqual(config.high_adc_thresh, 0.00115);
            testCase.verifyEqual(config.d_thresh, 0.001);
            testCase.verifyEqual(config.f_thresh, 0.1);
            testCase.verifyEqual(config.dstar_thresh, 0.01);
            testCase.verifyEqual(config.ivim_bthr, 100);
            testCase.verifyEqual(config.min_vox_hist, 100);
            testCase.verifyEqual(config.adc_max, 0.003);
            % Boolean default: checkpointing off by default
            testCase.verifyFalse(config.use_checkpoints);
        end

        function testExistingFieldsPreserved(testCase)
            % User-specified values should not be overwritten by defaults.
            % Here adc_thresh and d_thresh are set to non-default values;
            % f_thresh is omitted and should get its default.
            custom_config = struct();
            custom_config.dataloc = '/tmp';
            custom_config.adc_thresh = 0.002;   % non-default
            custom_config.d_thresh = 0.005;      % non-default

            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(custom_config));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);

            % Custom values preserved
            testCase.verifyEqual(config.adc_thresh, 0.002);
            testCase.verifyEqual(config.d_thresh, 0.005);

            % Omitted field gets its default
            testCase.verifyEqual(config.f_thresh, 0.1);
        end

        function testNonExistentFileThrowsError(testCase)
            % Requesting a config file that does not exist should throw
            % 'parse_config:fileNotFound'.
            testCase.verifyError(@() parse_config('/nonexistent/path/config.json'), ...
                'parse_config:fileNotFound');
        end

        function testMalformedJsonThrowsError(testCase)
            % Invalid JSON syntax (trailing text after valid JSON) should
            % throw 'parse_config:invalidJSON'.
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '{"dataloc": "/tmp", BROKEN}');
            fclose(fid);

            testCase.verifyError(@() parse_config(testCase.TempConfigFile), ...
                'parse_config:invalidJSON');
        end

        function testEmptyFileThrowsError(testCase)
            % A zero-byte file contains no valid JSON and should trigger
            % 'parse_config:invalidJSON'.
            fid = fopen(testCase.TempConfigFile, 'w');
            fclose(fid);

            testCase.verifyError(@() parse_config(testCase.TempConfigFile), ...
                'parse_config:invalidJSON');
        end

        function testDwiTypeStandard(testCase)
            % dwi_type "Standard" should map to numeric index 1.
            cfg = struct('dataloc', '/tmp', 'dwi_type', 'Standard');
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);
            testCase.verifyEqual(config.dwi_types_to_run, 1);
        end

        function testDwiTypeDnCNN(testCase)
            % dwi_type "dnCNN" should map to numeric index 2.
            cfg = struct('dataloc', '/tmp', 'dwi_type', 'dnCNN');
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);
            testCase.verifyEqual(config.dwi_types_to_run, 2);
        end

        function testDwiTypeIVIMnet(testCase)
            % dwi_type "IVIMnet" should map to numeric index 3.
            cfg = struct('dataloc', '/tmp', 'dwi_type', 'IVIMnet');
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);
            testCase.verifyEqual(config.dwi_types_to_run, 3);
        end

        function testInvalidDwiTypeThrowsError(testCase)
            % An unrecognized dwi_type string (not "Standard", "dnCNN", or
            % "IVIMnet") should throw an error. The outer try-catch in
            % parse_config wraps this as 'parse_config:invalidJSON'.
            cfg = struct('dataloc', '/tmp', 'dwi_type', 'InvalidType');
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            testCase.verifyError(@() parse_config(testCase.TempConfigFile), ...
                'parse_config:invalidJSON');
        end

        function testMissingDwiTypeDefaultsToAll(testCase)
            % When dwi_type field is absent, the pipeline should run all
            % three DWI types: dwi_types_to_run = [1, 2, 3].
            cfg = struct('dataloc', '/tmp');
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);
            testCase.verifyEqual(config.dwi_types_to_run, 1:3);
        end

        function testClearCacheDefaultsFalse(testCase)
            % clear_cache should default to false when not specified,
            % preserving cached checkpoint files by default.
            cfg = struct('dataloc', '/tmp');
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);
            testCase.verifyFalse(config.clear_cache);
        end

        function testClearCachePreservedWhenSet(testCase)
            % When clear_cache is explicitly set to true in the config file,
            % that value should be preserved (not overwritten by the default).
            cfg = struct('dataloc', '/tmp', 'clear_cache', true);
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);
            testCase.verifyTrue(config.clear_cache);
        end
    end
end
