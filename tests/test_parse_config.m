classdef test_parse_config < matlab.unittest.TestCase
    % TEST_PARSE_CONFIG Unit tests for configuration parsing logic

    properties
        TempConfigFile
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.TempConfigFile = [tempname '.json'];

            % Add utils directory to path
            current_dir = fileparts(mfilename('fullpath'));
            utils_dir = fullfile(current_dir, '..', 'utils');
            addpath(utils_dir);
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            if exist(testCase.TempConfigFile, 'file')
                delete(testCase.TempConfigFile);
            end
        end
    end

    methods(Test)
        function testMissingFieldsAssignedDefaults(testCase)
            % Create a minimal config file
            minimal_config = struct();
            minimal_config.dataloc = '/tmp';

            % Write to file
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(minimal_config));
            fclose(fid);

            % Parse config
            config = parse_config(testCase.TempConfigFile);

            % Verify defaults are assigned
            testCase.verifyEqual(config.adc_thresh, 0.00115);
            testCase.verifyEqual(config.high_adc_thresh, 0.001);
            testCase.verifyEqual(config.d_thresh, 0.001);
            testCase.verifyEqual(config.f_thresh, 0.1);
            testCase.verifyEqual(config.dstar_thresh, 0.01);
            testCase.verifyEqual(config.ivim_bthr, 100);
            testCase.verifyEqual(config.min_vox_hist, 100);
            testCase.verifyEqual(config.adc_max, 0.003);
            testCase.verifyFalse(config.use_checkpoints);
        end

        function testExistingFieldsPreserved(testCase)
            % Create a config file with custom values
            custom_config = struct();
            custom_config.dataloc = '/tmp';
            custom_config.adc_thresh = 0.002;
            custom_config.d_thresh = 0.005;

            % Write to file
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(custom_config));
            fclose(fid);

            % Parse config
            config = parse_config(testCase.TempConfigFile);

            % Verify custom values are preserved
            testCase.verifyEqual(config.adc_thresh, 0.002);
            testCase.verifyEqual(config.d_thresh, 0.005);

            % Verify other defaults are still assigned
            testCase.verifyEqual(config.f_thresh, 0.1);
        end

        function testNonExistentFileThrowsError(testCase)
            % parse_config should error when the file does not exist
            testCase.verifyError(@() parse_config('/nonexistent/path/config.json'), ...
                'MATLAB:error');
        end

        function testMalformedJsonThrowsError(testCase)
            % parse_config should error on invalid JSON syntax
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '{"dataloc": "/tmp", BROKEN}');
            fclose(fid);

            testCase.verifyError(@() parse_config(testCase.TempConfigFile), ...
                'MATLAB:error');
        end

        function testEmptyFileThrowsError(testCase)
            % A zero-byte config file should trigger an error
            fid = fopen(testCase.TempConfigFile, 'w');
            fclose(fid);

            testCase.verifyError(@() parse_config(testCase.TempConfigFile), ...
                'MATLAB:error');
        end

        function testDwiTypeStandard(testCase)
            % dwi_type "Standard" should map to dwi_types_to_run = 1
            cfg = struct('dataloc', '/tmp', 'dwi_type', 'Standard');
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);
            testCase.verifyEqual(config.dwi_types_to_run, 1);
        end

        function testDwiTypeDnCNN(testCase)
            % dwi_type "dnCNN" should map to dwi_types_to_run = 2
            cfg = struct('dataloc', '/tmp', 'dwi_type', 'dnCNN');
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);
            testCase.verifyEqual(config.dwi_types_to_run, 2);
        end

        function testDwiTypeIVIMnet(testCase)
            % dwi_type "IVIMnet" should map to dwi_types_to_run = 3
            cfg = struct('dataloc', '/tmp', 'dwi_type', 'IVIMnet');
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);
            testCase.verifyEqual(config.dwi_types_to_run, 3);
        end

        function testInvalidDwiTypeDefaultsToAll(testCase)
            % An unrecognized dwi_type should default to running all types [1:3]
            cfg = struct('dataloc', '/tmp', 'dwi_type', 'InvalidType');
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);
            testCase.verifyEqual(config.dwi_types_to_run, 1:3);
        end

        function testMissingDwiTypeDefaultsToAll(testCase)
            % When dwi_type field is absent, dwi_types_to_run should be [1:3]
            cfg = struct('dataloc', '/tmp');
            fid = fopen(testCase.TempConfigFile, 'w');
            fprintf(fid, '%s', jsonencode(cfg));
            fclose(fid);

            config = parse_config(testCase.TempConfigFile);
            testCase.verifyEqual(config.dwi_types_to_run, 1:3);
        end
    end
end
