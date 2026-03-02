classdef test_metrics_baseline < matlab.unittest.TestCase
    % TEST_METRICS_BASELINE Unit tests for the metrics_baseline function.
    % Ensures that raw data vectors and computed summaries from compute_summary_metrics
    % are correctly curated for downstream modeling.

    properties
        ConfigStruct
        DataVectorsGTVp
        DataVectorsGTVn
        SummaryMetrics
        TempDir
    end

    methods(TestMethodSetup)
        function createMockInputs(testCase)
            % Create a temp dir for the mock clinical Excel file
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            testCase.ConfigStruct = struct( ...
                'dataloc',            testCase.TempDir, ...
                'use_checkpoints',    false, ...
                'dwi_types_to_run',   1, ...
                'clinical_data_sheet','mock_clinical.xlsx');

            % Write a minimal mock clinical Excel file with 4 matching patient IDs
            id_list_xls = {'P1'; 'P2'; 'P3'; 'P4'};
            lf_vals     = [0; 1; 0; 1];
            dt_event    = datetime({'2023-06-01'; '2023-06-01'; '2023-06-01'; '2023-06-01'});
            dt_censor   = datetime({'2023-06-01'; '2023-06-01'; '2023-06-01'; '2023-06-01'});
            dt_reg      = datetime({'2023-06-01'; '2023-06-01'; '2023-06-01'; '2023-06-01'});
            dt_rtstart  = datetime({'2022-01-01'; '2022-01-01'; '2022-01-01'; '2022-01-01'});
            dt_rtstop   = datetime({'2022-03-01'; '2022-03-01'; '2022-03-01'; '2022-03-01'});
            T_clin = table(id_list_xls, lf_vals, dt_event, dt_censor, dt_reg, ...
                dt_rtstart, dt_rtstop, ...
                'VariableNames', {'Pat', 'LocalOrRegionalFailure', ...
                'LocoregionalFailureDateOfLocalOrRegionalFailure', ...
                'LocalFailureDateOfLocalFailureOrCensor', ...
                'RegionalFailureDateOfRegionalFailureOrCensor', ...
                'RTStartDate', 'RTStopDate'});
            writetable(T_clin, fullfile(testCase.TempDir, 'mock_clinical.xlsx'));

            % Mock Data vectors GTVp (unused by metrics_baseline directly but required input)
            testCase.DataVectorsGTVp = repmat( ...
                struct('adc', [1, 2, 3]', 'd', [1, 2, 3]', 'f', [0.1, 0.2, 0.3]', 'dstar', [0.01, 0.02, 0.03]'), ...
                4, 3);
            testCase.DataVectorsGTVn = repmat( ...
                struct('adc', [1.5, 2.5]', 'd', [1.5, 2.5]', 'f', [0.15, 0.25]', 'dstar', [0.015, 0.025]'), ...
                4, 3);

            % Mock Summary metrics — 4 patients, 3 timepoints
            testCase.SummaryMetrics = struct();
            testCase.SummaryMetrics.id_list    = {'P1', 'P2', 'P3', 'P4'};
            testCase.SummaryMetrics.mrn_list   = {'M1', 'M2', 'M3', 'M4'};
            testCase.SummaryMetrics.lf         = [0; 1; 0; 1];
            testCase.SummaryMetrics.gtv_vol    = repmat([10; 15; 8; 12], 1, 3);
            testCase.SummaryMetrics.dmean_gtvp = repmat([50; 55; 48; 52], 1, 3);
            testCase.SummaryMetrics.d95_gtvp   = repmat([45; 50; 43; 47], 1, 3);
            testCase.SummaryMetrics.v50gy_gtvp = repmat([90; 85; 88; 92], 1, 3);

            % Deterministic DWI metric arrays to avoid non-reproducible outlier detection
            testCase.SummaryMetrics.adc_mean  = repmat([1.0; 1.2; 0.8; 1.1], 1, 3) * 1e-3;
            testCase.SummaryMetrics.d_mean    = repmat([0.8; 1.0; 0.7; 0.9], 1, 3) * 1e-3;
            testCase.SummaryMetrics.f_mean    = repmat([0.20; 0.20; 0.20; 0.20], 1, 3);
            testCase.SummaryMetrics.dstar_mean = repmat([0.01; 0.01; 0.01; 0.01], 1, 3);
            testCase.SummaryMetrics.adc_sd    = repmat([0.1; 0.1; 0.1; 0.1], 1, 3) * 1e-3;

            % Repeatability fields
            testCase.SummaryMetrics.n_rpt          = [2; 2; 2; 2];
            testCase.SummaryMetrics.adc_mean_rpt   = rand(4, 2);
            testCase.SummaryMetrics.adc_sub_rpt    = rand(4, 2);
            testCase.SummaryMetrics.d_mean_rpt     = rand(4, 2);
            testCase.SummaryMetrics.f_mean_rpt     = rand(4, 2);
            testCase.SummaryMetrics.dstar_mean_rpt = rand(4, 2);

            % Fields accessed by some metric modules
            testCase.SummaryMetrics.gtv_locations  = {};
            testCase.SummaryMetrics.dwi_locations  = {};
        end
    end

    methods(TestMethodTeardown)
        function cleanupTempDir(testCase)
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
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

            % Check that valid_pts masks correctly (4 patients, all with finite lf)
            testCase.verifyEqual(length(valid_pts), 4);
            testCase.verifyTrue(all(valid_pts));

            % Check output matrix dimensions
            testCase.verifyEqual(size(ADC_abs, 1), 4);  % 4 patients
            testCase.verifyEqual(size(ADC_pct, 1), 4);
            testCase.verifyEqual(size(ADC_abs, 2), 3);  % 3 timepoints
            testCase.verifyEqual(size(ADC_pct, 2), 3);  % same shape as ADC_abs
        end

        function testOutlierCleaning(testCase)
            % Insert an outlier into the summary metrics for patient 1 at baseline.
            % With 4 patients the IQR-based fence is computed and 1e6 is detected.
            testCase.SummaryMetrics.adc_mean(1, 1) = 1e6;

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
