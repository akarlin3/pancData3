classdef test_load_data_from_disk < matlab.unittest.TestCase
    % TEST_LOAD_DATA_FROM_DISK Unit tests for load_data_from_disk.
    %
    % Validates data loading with type-specific files, legacy fallback for
    % Standard only, and error on missing files.

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
            rmdir(testCase.TmpDir, 's');
        end
    end

    methods(Test)
        function test_loads_type_specific_file(testCase)
            % When the type-specific file exists, it should be loaded.
            data_vectors_gtvp = struct('adc_vector', [1; 2; 3]);
            data_vectors_gtvn = struct('adc_vector', [4; 5]);
            dwi_file = fullfile(testCase.TmpDir, 'dwi_vectors_Standard.mat');
            save(dwi_file, 'data_vectors_gtvp', 'data_vectors_gtvn');

            [gtvp, gtvn, sm] = load_data_from_disk( ...
                dwi_file, '', '', 1, 'Standard');
            testCase.verifyEqual(gtvp.adc_vector, [1; 2; 3]);
            testCase.verifyEqual(gtvn.adc_vector, [4; 5]);
            testCase.verifyTrue(isempty(sm));
        end

        function test_legacy_fallback_for_standard(testCase)
            % Standard (dtype=1) should fall back to legacy dwi_vectors.mat.
            data_vectors_gtvp = struct('adc_vector', [10; 20]);
            data_vectors_gtvn = struct('adc_vector', [30]);
            fallback = fullfile(testCase.TmpDir, 'dwi_vectors.mat');
            save(fallback, 'data_vectors_gtvp', 'data_vectors_gtvn');

            nonexistent = fullfile(testCase.TmpDir, 'dwi_vectors_Standard.mat');
            [gtvp, ~, ~] = load_data_from_disk( ...
                nonexistent, fallback, '', 1, 'Standard');
            testCase.verifyEqual(gtvp.adc_vector, [10; 20]);
        end

        function test_no_legacy_fallback_for_dncnn(testCase)
            % dnCNN (dtype=2) must NOT fall back to legacy file.
            data_vectors_gtvp = struct('adc_vector', [10; 20]);
            data_vectors_gtvn = struct('adc_vector', [30]);
            fallback = fullfile(testCase.TmpDir, 'dwi_vectors.mat');
            save(fallback, 'data_vectors_gtvp', 'data_vectors_gtvn');

            nonexistent = fullfile(testCase.TmpDir, 'dwi_vectors_dnCNN.mat');
            testCase.verifyError( ...
                @() load_data_from_disk(nonexistent, fallback, '', 2, 'dnCNN'), ...
                'LoadDataFromDisk:NotFound');
        end

        function test_loads_summary_metrics(testCase)
            % Summary metrics file should be loaded when present.
            data_vectors_gtvp = struct('adc_vector', [1]);
            data_vectors_gtvn = struct('adc_vector', [2]);
            dwi_file = fullfile(testCase.TmpDir, 'dwi_vectors_Standard.mat');
            save(dwi_file, 'data_vectors_gtvp', 'data_vectors_gtvn');

            summary_metrics = struct('ADC_abs', [0.001, 0.002]);
            sm_file = fullfile(testCase.TmpDir, 'summary_metrics_Standard.mat');
            save(sm_file, 'summary_metrics');

            [~, ~, sm] = load_data_from_disk( ...
                dwi_file, '', sm_file, 1, 'Standard');
            testCase.verifyEqual(sm.ADC_abs, [0.001, 0.002]);
        end

        function test_missing_file_errors(testCase)
            % When neither file exists, should throw an error.
            testCase.verifyError( ...
                @() load_data_from_disk('/nonexistent1.mat', '/nonexistent2.mat', '', 1, 'Standard'), ...
                'LoadDataFromDisk:NotFound');
        end
    end
end
