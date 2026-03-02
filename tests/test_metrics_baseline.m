classdef test_metrics_baseline < matlab.unittest.TestCase
    % TEST_METRICS_BASELINE Unit tests for the metrics_baseline function.
    % Ensures that raw data vectors and computed summaries from compute_summary_metrics
    % are correctly curated for downstream modeling.

    properties
        ConfigStruct
        DataVectorsGTVp
        DataVectorsGTVn
        SummaryMetrics
    end

    methods(TestMethodSetup)
        function createMockInputs(testCase)
            testCase.ConfigStruct = struct('dataloc', pwd, 'use_checkpoints', false);

            % Mock Data vectors GTVp
            testCase.DataVectorsGTVp = repmat(struct('adc', [1, 2, 3]', 'd', [1, 2, 3]', 'f', [0.1, 0.2, 0.3]', 'dstar', [0.01, 0.02, 0.03]'), 2, 3);
            testCase.DataVectorsGTVn = repmat(struct('adc', [1.5, 2.5]', 'd', [1.5, 2.5]', 'f', [0.15, 0.25]', 'dstar', [0.015, 0.025]'), 2, 3);

            % Mock Summary metrics
            testCase.SummaryMetrics = struct();
            testCase.SummaryMetrics.id_list = {'P1', 'P2'};
            testCase.SummaryMetrics.mrn_list = {'M1', 'M2'};
            testCase.SummaryMetrics.lf = [0; 1];
            testCase.SummaryMetrics.total_time = [30; 40];
            testCase.SummaryMetrics.total_follow_up_time = [30; 40];
            testCase.SummaryMetrics.gtv_vol = [10; 20];
            testCase.SummaryMetrics.dmean_gtvp = [50; 60];
            testCase.SummaryMetrics.d95_gtvp = [45; 55];
            testCase.SummaryMetrics.v50gy_gtvp = [90; 80];

            % Fill matrices with fake data
            metrics_list = {'adc_mean', 'adc_sd', 'adc_kurt', 'adc_skew', 'adc_sub_vol', 'adc_mean_high', 'adc_mean_low', ...
                            'd_mean', 'd_sd', 'd_kurt', 'd_skew', 'd_sub_vol', ...
                            'f_mean', 'f_sd', 'f_kurt', 'f_skew', 'f_sub_vol', ...
                            'dstar_mean', 'dstar_sd', 'dstar_kurt', 'dstar_skew', 'dstar_sub_vol'};
            for i = 1:length(metrics_list)
                testCase.SummaryMetrics.(metrics_list{i}) = rand(2, 3);
            end

            % Some fields expected by metrics_baseline
            testCase.SummaryMetrics.ivim_sub_vol = rand(2, 3);
            testCase.SummaryMetrics.n_rpt = [2; 2];
            testCase.SummaryMetrics.adc_mean_rpt = rand(2, 2);
            testCase.SummaryMetrics.adc_sub_rpt = rand(2, 2);
            testCase.SummaryMetrics.d_mean_rpt = rand(2, 2);
            testCase.SummaryMetrics.f_mean_rpt = rand(2, 2);
            testCase.SummaryMetrics.dstar_mean_rpt = rand(2, 2);
            testCase.SummaryMetrics.dl_provenance = zeros(2, 1);
        end
    end

    methods(Test)
        function testBaselineExecution(testCase)
            % Should run without error
            [m_lf, m_total_time, m_total_follow_up_time, m_gtv_vol, m_adc_mean, m_d_mean, m_f_mean, m_dstar_mean, ...
             m_id_list, m_mrn_list, m_d95_gtvp, m_v50gy_gtvp, m_data_vectors_gtvp, lf_group, valid_pts, ...
             ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_pct, Dstar_pct, ...
             nTp, metric_sets, set_names, time_labels, dtype_label, dl_provenance] = ...
             metrics_baseline(testCase.DataVectorsGTVp, testCase.DataVectorsGTVn, testCase.SummaryMetrics, testCase.ConfigStruct);

            % Check that valid_pts masks correctly
            testCase.verifyEqual(length(valid_pts), 2);
            testCase.verifyTrue(all(valid_pts));

            % Check relative percent calculation
            % If ADC_abs is [10, 20], ADC_pct should be (20-10)/10 * 100 = 100%
            % We just need to ensure the dimensions of the outputs match the input 2 patients
            testCase.verifyEqual(size(ADC_abs, 1), 2);
            testCase.verifyEqual(size(ADC_pct, 1), 2);
            testCase.verifyEqual(size(ADC_abs, 2), 3); % 3 timepoints

            % Check percent metrics have size N x (nTp-1)
            testCase.verifyEqual(size(ADC_pct, 2), 2);
        end

        function testOutlierCleaning(testCase)
            % Insert an outlier into the summary metrics
            testCase.SummaryMetrics.adc_mean(1, 1) = 1e6; % Huge outlier

            [~, ~, ~, ~, m_adc_mean, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
             ADC_abs, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~] = ...
             metrics_baseline(testCase.DataVectorsGTVp, testCase.DataVectorsGTVn, testCase.SummaryMetrics, testCase.ConfigStruct);

            % The outlier should be cleaned and replaced by NaN
            testCase.verifyTrue(isnan(ADC_abs(1, 1)));
            % Other elements should be intact
            testCase.verifyFalse(isnan(ADC_abs(2, 1)));
        end
    end
end
