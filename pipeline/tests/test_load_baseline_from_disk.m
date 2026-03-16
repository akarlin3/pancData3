classdef test_load_baseline_from_disk < matlab.unittest.TestCase
    % TEST_LOAD_BASELINE_FROM_DISK Unit tests for load_baseline_from_disk.
    %
    % Validates loading of persisted metrics_baseline outputs and error
    % handling when the file is missing.

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
        function test_loads_all_fields(testCase)
            % All baseline fields should be accessible in the returned struct.
            baseline_file = fullfile(testCase.TmpDir, 'metrics_baseline_results_Standard.mat');
            m_lf = [1; 0; 1];
            m_total_time = [12; 24; 6];
            m_total_follow_up_time = [24; 36; 12];
            m_gtv_vol = [10; 20; 15];
            m_adc_mean = rand(3, 5);
            m_d_mean = rand(3, 5);
            m_f_mean = rand(3, 5);
            m_dstar_mean = rand(3, 5);
            m_id_list = {'P01', 'P02', 'P03'};
            m_mrn_list = {'MRN1', 'MRN2', 'MRN3'};
            m_d95_gtvp = rand(3, 5);
            m_v50gy_gtvp = rand(3, 5);
            m_data_vectors_gtvp = struct('adc', 1);
            lf_group = [1; 0; 1];
            valid_pts = [true; true; false];
            ADC_abs = rand(3, 5);
            D_abs = rand(3, 5);
            f_abs = rand(3, 5);
            Dstar_abs = rand(3, 5);
            ADC_pct = rand(3, 5);
            D_pct = rand(3, 5);
            f_delta = rand(3, 5);
            Dstar_pct = rand(3, 5);
            nTp = 5;
            metric_sets = {{ADC_abs, D_abs}};
            set_names = {{'ADC', 'D'}};
            time_labels = {'Fx1', 'Fx2', 'Fx3', 'Fx4', 'Fx5'};
            dtype_label = 'Standard';
            dl_provenance = struct('train_ids', {{}});
            save(baseline_file, 'm_lf', 'm_total_time', 'm_total_follow_up_time', ...
                'm_gtv_vol', 'm_adc_mean', 'm_d_mean', 'm_f_mean', 'm_dstar_mean', ...
                'm_id_list', 'm_mrn_list', 'm_d95_gtvp', 'm_v50gy_gtvp', ...
                'm_data_vectors_gtvp', 'lf_group', 'valid_pts', ...
                'ADC_abs', 'D_abs', 'f_abs', 'Dstar_abs', ...
                'ADC_pct', 'D_pct', 'f_delta', 'Dstar_pct', ...
                'nTp', 'metric_sets', 'set_names', 'time_labels', ...
                'dtype_label', 'dl_provenance');

            baseline = load_baseline_from_disk(baseline_file);

            testCase.verifyEqual(baseline.m_lf, m_lf);
            testCase.verifyEqual(baseline.nTp, 5);
            testCase.verifyEqual(baseline.dtype_label, 'Standard');
            testCase.verifyEqual(baseline.valid_pts, valid_pts);
            testCase.verifyEqual(baseline.ADC_abs, ADC_abs);
        end

        function test_error_on_missing_file(testCase)
            % Should throw LoadBaseline:NotFound when file doesn't exist.
            missing = fullfile(testCase.TmpDir, 'nonexistent.mat');
            testCase.verifyError(@() load_baseline_from_disk(missing), ...
                'LoadBaseline:NotFound');
        end
    end
end
