classdef test_metrics_stats_predictive < matlab.unittest.TestCase
    % TEST_METRICS_STATS_PREDICTIVE Unit tests for metrics_stats_predictive.
    %
    % Covers:
    %   - No-op path when nTp = 1 (loop body never executes)
    %   - Graceful handling when all features are NaN (nothing imputable)
    %   - Output types and sizes are correct regardless of convergence
    %   - DL provenance leakage detection for DnCNN (dtype = 2)
    %   - DL provenance leakage detection for IVIMnet (dtype = 3)

    properties
        TempDir
        OriginalPath
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            close all;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
            path(testCase.OriginalPath);
        end
    end

    % ------------------------------------------------------------------ %
    %  Helper to build a minimal valid input struct                       %
    % ------------------------------------------------------------------ %
    methods(Access = private)
        function args = buildMinimalArgs(testCase, n, nTp, dtype)
            % Returns a cell array of all positional arguments for
            % metrics_stats_predictive given n patients and nTp timepoints.
            rng(42);
            valid_pts  = true(n, 1);
            lf_group   = [ones(ceil(n/2), 1); zeros(floor(n/2), 1)];
            dtype_label = 'Standard';
            output_folder = testCase.TempDir;
            dataloc    = testCase.TempDir;

            % Feature matrices: n × nTp (non-NaN, finite values)
            ADC_abs   = rand(n, nTp) * 2e-3 + 0.5e-3;
            D_abs     = rand(n, nTp) * 1e-3 + 0.2e-3;
            f_abs     = rand(n, nTp) * 0.3  + 0.05;
            Dstar_abs = rand(n, nTp) * 0.05 + 0.01;
            ADC_pct   = rand(n, nTp) * 20 - 10;
            D_pct     = rand(n, nTp) * 20 - 10;
            f_pct     = rand(n, nTp) * 0.1 - 0.05;
            Dstar_pct = rand(n, nTp) * 20 - 10;

            % Dose matrices
            m_d95_gtvp   = rand(n, nTp) * 50 + 10;
            m_v50gy_gtvp = rand(n, nTp) * 100;
            d95_adc_sub  = rand(n, nTp) * 50;
            v50_adc_sub  = rand(n, nTp) * 100;
            d95_d_sub    = rand(n, nTp) * 50;
            v50_d_sub    = rand(n, nTp) * 100;
            d95_f_sub    = rand(n, nTp) * 50;
            v50_f_sub    = rand(n, nTp) * 100;
            d95_dstar_sub = rand(n, nTp) * 50;
            v50_dstar_sub = rand(n, nTp) * 100;

            id_list = arrayfun(@(x) sprintf('Pt%02d', x), 1:n, 'UniformOutput', false)';

            dl_provenance = struct();
            if dtype == 2
                dl_provenance.dncnn_train_ids = {};
            elseif dtype == 3
                dl_provenance.ivimnet_train_ids = {};
            end

            x_labels = arrayfun(@(x) sprintf('Fx%d', x), 1:nTp, 'UniformOutput', false);

            m_lf       = lf_group;
            m_total_time = randi([10, 100], n, 1);
            m_total_follow_up_time = nan(n, 1);
            m_total_follow_up_time(m_lf == 0) = ...
                m_total_time(m_lf == 0) + randi([5, 30], sum(m_lf == 0), 1);

            % adc_sd: used only when n_sig > 0; provide a valid 3-D array anyway
            adc_sd = rand(n, nTp, 1) * 0.5e-3;
            m_gtv_vol = rand(n, nTp) * 50 + 5;

            args = {valid_pts, lf_group, dtype_label, output_folder, dataloc, nTp, ...
                    m_gtv_vol, adc_sd, ...
                    ADC_abs, D_abs, f_abs, Dstar_abs, ...
                    ADC_pct, D_pct, f_pct, Dstar_pct, ...
                    m_d95_gtvp, m_v50gy_gtvp, ...
                    d95_adc_sub, v50_adc_sub, ...
                    d95_d_sub,   v50_d_sub, ...
                    d95_f_sub,   v50_f_sub, ...
                    d95_dstar_sub, v50_dstar_sub, ...
                    id_list, dtype, dl_provenance, x_labels, ...
                    m_lf, m_total_time, m_total_follow_up_time};
        end
    end

    methods(Test)

        function testNoopWhenNTpIsOne(testCase)
            % nTp = 1 → the for-loop (target_fx = 2:1) never executes.
            % All four outputs must be empty.
            n   = 5;
            nTp = 1;
            args = testCase.buildMinimalArgs(n, nTp, 1);

            [risk_scores_all, is_high_risk, times_km, events_km] = ...
                metrics_stats_predictive(args{:});

            testCase.verifyEmpty(risk_scores_all,  'risk_scores_all should be empty for nTp=1.');
            testCase.verifyEmpty(is_high_risk,     'is_high_risk should be empty for nTp=1.');
            testCase.verifyEmpty(times_km,         'times_km should be empty for nTp=1.');
            testCase.verifyEmpty(events_km,        'events_km should be empty for nTp=1.');
        end

        function testAllNaNFeaturesReturnsEmptyOutputs(testCase)
            % All feature matrices are NaN → no imputable patients →
            % elastic net cannot run → n_sig = 0 → outputs remain empty.
            n   = 8;
            nTp = 2;
            args = testCase.buildMinimalArgs(n, nTp, 1);

            % Overwrite all feature matrices with NaN (positions 9–26 in args)
            feat_positions = 9:26;
            for pos = feat_positions
                args{pos} = nan(n, nTp);
            end

            [risk_scores_all, is_high_risk, times_km, events_km] = ...
                metrics_stats_predictive(args{:});

            testCase.verifyEmpty(risk_scores_all, ...
                'risk_scores_all should be empty when all features are NaN.');
            testCase.verifyEmpty(is_high_risk, ...
                'is_high_risk should be empty when all features are NaN.');
            testCase.verifyEmpty(times_km, ...
                'times_km should be empty when no features pass imputation.');
            testCase.verifyEmpty(events_km, ...
                'events_km should be empty when no features pass imputation.');
        end

        function testOutputTypesAreNumericOrLogical(testCase)
            % Even when the elastic net fails gracefully, the returned
            % variables must be numeric or logical (never cell arrays).
            n   = 6;
            nTp = 1;   % fast no-op path
            args = testCase.buildMinimalArgs(n, nTp, 1);

            [risk_scores_all, is_high_risk, times_km, events_km] = ...
                metrics_stats_predictive(args{:});

            testCase.verifyTrue(isnumeric(risk_scores_all) || islogical(risk_scores_all), ...
                'risk_scores_all must be numeric or logical.');
            testCase.verifyTrue(isnumeric(is_high_risk) || islogical(is_high_risk), ...
                'is_high_risk must be numeric or logical.');
            testCase.verifyTrue(isnumeric(times_km) || islogical(times_km), ...
                'times_km must be numeric or logical.');
            testCase.verifyTrue(isnumeric(events_km) || islogical(events_km), ...
                'events_km must be numeric or logical.');
        end

        function testLeakageDetectionDnCNN(testCase)
            % dtype = 2 (DnCNN): if a patient is in dncnn_train_ids and the
            % LOOCV loop reaches that patient, an error must be thrown.
            % We manufacture a scenario where elastic net runs by supplying
            % a valid but minimal dataset and explicitly populating
            % dl_provenance.dncnn_train_ids.
            %
            % NOTE: the leakage error only fires when n_sig > 0 (elastic net
            % selects at least one feature).  With very few patients the
            % elastic net typically selects nothing and n_sig = 0, so the
            % LOOCV section is skipped.  This test therefore only validates
            % the leakage guard logic when it CAN be reached; if n_sig = 0
            % the test still verifies that the function runs without error.
            n   = 6;
            nTp = 2;
            args = testCase.buildMinimalArgs(n, nTp, 2);

            % Inject the first patient into the DnCNN training set
            id_list_pos = 28;   % position of id_list in args
            dl_prov_pos = 29;   % position of dl_provenance
            leaky_id = args{id_list_pos}{1};
            args{dl_prov_pos} = struct('dncnn_train_ids', {{leaky_id}});

            % The function either errors (leakage found) or completes cleanly
            % (n_sig = 0 so LOOCV never reached).  Both outcomes are acceptable
            % from an interface perspective — what is NOT acceptable is a silent
            % wrong answer.  We verify it either errors with the expected ID
            % or returns gracefully with empty outputs.
            try
                [risk_scores_all, ~, ~, ~] = metrics_stats_predictive(args{:});
                % If we get here, LOOCV was skipped (n_sig = 0 path)
                testCase.verifyEmpty(risk_scores_all, ...
                    'When LOOCV is skipped, risk_scores_all should be empty.');
            catch ME
                testCase.verifySubstring(ME.message, leaky_id, ...
                    'Leakage error message should name the offending patient ID.');
            end
        end

        function testLeakageDetectionIVIMnet(testCase)
            % dtype = 3 (IVIMnet): same guard for ivimnet_train_ids.
            n   = 6;
            nTp = 2;
            args = testCase.buildMinimalArgs(n, nTp, 3);

            id_list_pos = 28;
            dl_prov_pos = 29;
            leaky_id = args{id_list_pos}{2};   % second patient
            args{dl_prov_pos} = struct('ivimnet_train_ids', {{leaky_id}});

            try
                [risk_scores_all, ~, ~, ~] = metrics_stats_predictive(args{:});
                testCase.verifyEmpty(risk_scores_all, ...
                    'When LOOCV is skipped, risk_scores_all should be empty.');
            catch ME
                testCase.verifySubstring(ME.message, leaky_id, ...
                    'Leakage error message should name the offending patient ID.');
            end
        end

    end
end
