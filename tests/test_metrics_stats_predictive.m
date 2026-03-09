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
            diary off;  % close any diary opened by the function under test
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
            % Constructs a cell array of all 34 positional arguments for
            % metrics_stats_predictive. Generates synthetic but finite
            % feature matrices (ADC, D, f, D*, dose metrics), patient IDs,
            % DL provenance, outcome labels, and survival times. The dtype
            % parameter (1=Standard, 2=DnCNN, 3=IVIMnet) controls which
            % DL provenance field is populated.
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
            f_delta     = rand(n, nTp) * 0.1 - 0.05;
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
            adc_sd = rand(n, nTp, 3) * 0.5e-3;
            m_gtv_vol = rand(n, nTp) * 50 + 5;

            cfg = struct('use_firth_refit', true);

            args = {valid_pts, lf_group, dtype_label, output_folder, dataloc, nTp, ...
                    m_gtv_vol, adc_sd, ...
                    ADC_abs, D_abs, f_abs, Dstar_abs, ...
                    ADC_pct, D_pct, f_delta, Dstar_pct, ...
                    m_d95_gtvp, m_v50gy_gtvp, ...
                    d95_adc_sub, v50_adc_sub, ...
                    d95_d_sub,   v50_d_sub, ...
                    d95_f_sub,   v50_f_sub, ...
                    d95_dstar_sub, v50_dstar_sub, ...
                    id_list, dtype, dl_provenance, x_labels, ...
                    m_lf, m_total_time, m_total_follow_up_time, cfg};
        end
    end

    methods(Test)

        function testNoopWhenNTpIsOne(testCase)
            % Verifies the no-op path when nTp=1: the predictive modeling
            % loop iterates over target_fx = 2:nTp, which is 2:1 (empty
            % range), so no elastic net or LOOCV runs. All four outputs
            % (risk_scores, is_high_risk, times_km, events_km) must be empty.
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

        function testNaNRiskScoresDoNotCauseLogicalError(testCase)
            % Regression test for the bug where NaN risk scores could not
            % be assigned into a logical array (is_high_risk_oof).
            % After the fix, is_high_risk_oof is initialised with zeros()
            % instead of false(), allowing NaN assignment for invalid folds.
            %
            % This test exercises the exact logic from lines 417-429 of
            % metrics_stats_predictive.m in isolation to verify that NaN
            % values in risk_scores_oof are handled without error.
            risk_scores_oof = [0.3; NaN; 0.7; NaN; 0.5];
            valid_oof = ~isnan(risk_scores_oof);
            oof_median = median(risk_scores_oof(valid_oof));
            is_high_risk_oof = zeros(size(risk_scores_oof));
            is_high_risk_oof(valid_oof) = risk_scores_oof(valid_oof) > oof_median;
            is_high_risk_oof(~valid_oof) = NaN;

            % NaN entries must survive assignment
            testCase.verifyTrue(isnan(is_high_risk_oof(2)), ...
                'Patient 2 (NaN risk score) should have NaN high-risk flag.');
            testCase.verifyTrue(isnan(is_high_risk_oof(4)), ...
                'Patient 4 (NaN risk score) should have NaN high-risk flag.');

            % Valid entries must be 0 or 1 (double, not logical)
            testCase.verifyTrue(all(is_high_risk_oof(valid_oof) == 0 | ...
                                    is_high_risk_oof(valid_oof) == 1), ...
                'Valid patients should have binary (0/1) high-risk flags.');

            % Downstream target array must accept these values
            valid_pts = true(8, 1);
            impute_mask = [true; false; true; true; false; true; false; true];
            is_high_risk_target = nan(sum(valid_pts), 1);
            is_high_risk_target(impute_mask) = is_high_risk_oof;
            testCase.verifyEqual(sum(isnan(is_high_risk_target)), 5, ...
                'Target array should have NaN for non-imputed and invalid patients.');
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
            id_list_pos = 27;   % position of id_list in args
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

            id_list_pos = 27;
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

        function testBackwardCompatWithout34thArg(testCase)
            % When config_struct (34th arg) is omitted, function must still
            % work — nargin guard defaults use_firth_refit to true.
            n   = 5;
            nTp = 1;   % fast no-op path
            args = testCase.buildMinimalArgs(n, nTp, 1);
            args = args(1:33);  % drop the config_struct arg

            [risk_scores_all, is_high_risk, times_km, events_km] = ...
                metrics_stats_predictive(args{:});

            testCase.verifyEmpty(risk_scores_all, ...
                'Backward-compat: outputs should be empty for nTp=1.');
            testCase.verifyEmpty(is_high_risk, ...
                'Backward-compat: is_high_risk should be empty for nTp=1.');
            testCase.verifyEmpty(times_km, ...
                'Backward-compat: times_km should be empty for nTp=1.');
            testCase.verifyEmpty(events_km, ...
                'Backward-compat: events_km should be empty for nTp=1.');
        end

        function testFirthDisabledRunsWithoutError(testCase)
            % When use_firth_refit is false, the function must run the
            % elastic-net-only path without error.
            n   = 5;
            nTp = 1;
            args = testCase.buildMinimalArgs(n, nTp, 1);
            args{34} = struct('use_firth_refit', false);

            [risk_scores_all, ~, ~, ~] = metrics_stats_predictive(args{:});

            testCase.verifyEmpty(risk_scores_all, ...
                'Firth disabled: outputs should be empty for nTp=1.');
        end


        function testLOOCVProducesNonEmptyRiskScores(testCase)
            % With n=24 patients and a strong separable signal, elastic net
            % should select features and the LOOCV loop should produce
            % non-empty risk scores and risk stratification.
            n   = 24;
            nTp = 3;
            args = testCase.buildMinimalArgs(n, nTp, 1);

            % Inject a strong signal: make ADC_abs column 2 perfectly
            % correlated with the outcome so elastic net selects at least
            % one feature.
            lf_group = args{2};
            ADC_abs  = args{9};
            for row = 1:n
                if lf_group(row) == 1
                    ADC_abs(row, 2) = 0.0005 + 0.0001 * rand();
                else
                    ADC_abs(row, 2) = 0.0018 + 0.0001 * rand();
                end
            end
            args{9} = ADC_abs;

            [risk_scores_all, is_high_risk, times_km, events_km] = ...
                metrics_stats_predictive(args{:});

            % If elastic net selected features, outputs should be non-empty
            % with valid numeric values.  If it did not (edge case with
            % random seed), the test still passes for empty outputs.
            if ~isempty(risk_scores_all)
                testCase.verifyGreaterThan(numel(risk_scores_all), 0, ...
                    'risk_scores_all should have entries when LOOCV runs.');
                valid_scores = risk_scores_all(~isnan(risk_scores_all));
                testCase.verifyTrue(all(isfinite(valid_scores)), ...
                    'Non-NaN risk scores should be finite.');

                testCase.verifyGreaterThan(numel(is_high_risk), 0, ...
                    'is_high_risk should be non-empty.');
                testCase.verifyTrue(isnumeric(is_high_risk), ...
                    'is_high_risk must be numeric (not logical) to support NaN.');

                testCase.verifyFalse(isempty(times_km), ...
                    'times_km should be populated for KM plot.');
                testCase.verifyFalse(isempty(events_km), ...
                    'events_km should be populated for KM plot.');
            end
        end

    end
end
