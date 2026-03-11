classdef test_metrics_baseline < matlab.unittest.TestCase
    % TEST_METRICS_BASELINE Unit tests for the metrics_baseline function.
    % Ensures that raw data vectors and computed summaries from compute_summary_metrics
    % are correctly curated for downstream modeling.

    properties
        ConfigStruct      % Minimal pipeline config for metrics_baseline
        DataVectorsGTVp   % Mock GTVp voxel-level data (unused directly by metrics_baseline but required input)
        DataVectorsGTVn   % Mock GTVn voxel-level data
        SummaryMetrics    % Mock summary metrics with deterministic values for reproducible tests
        TempDir           % Temporary directory for mock clinical Excel files and output
    end

    methods(TestMethodSetup)
        function createMockInputs(testCase)
            % Creates 8 mock patients with 3 timepoints: deterministic DWI
            % metric values (to avoid non-reproducible outlier detection),
            % a mock clinical Excel file with LF/date columns, and
            % repeatability sub-fields. The mock data is designed so that
            % all patients pass the valid_pts filter (finite LF values).
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            testCase.ConfigStruct = struct( ...
                'dataloc',            testCase.TempDir, ...
                'use_checkpoints',    false, ...
                'dwi_types_to_run',   1, ...
                'clinical_data_sheet','mock_clinical.xlsx');

            % Number of mock patients
            nPat = 8;

            % Write a minimal mock clinical Excel file with matching patient IDs
            if ~exist('OCTAVE_VERSION', 'builtin')
                id_list_xls = arrayfun(@(x) sprintf('P%d', x), (1:nPat)', 'UniformOutput', false);
                lf_vals     = repmat([0; 1], nPat/2, 1);
                dt_event    = repmat(datetime('2023-06-01'), nPat, 1);
                dt_censor   = repmat(datetime('2023-06-01'), nPat, 1);
                dt_reg      = repmat(datetime('2023-06-01'), nPat, 1);
                dt_rtstart  = repmat(datetime('2022-01-01'), nPat, 1);
                dt_rtstop   = repmat(datetime('2022-03-01'), nPat, 1);
                T_clin = table(id_list_xls, lf_vals, dt_event, dt_censor, dt_reg, ...
                    dt_rtstart, dt_rtstop, ...
                    'VariableNames', {'Pat', 'LocalOrRegionalFailure', ...
                    'LocoregionalFailureDateOfLocalOrRegionalFailure', ...
                    'LocalFailureDateOfLocalFailureOrCensor', ...
                    'RegionalFailureDateOfRegionalFailureOrCensor', ...
                    'RTStartDate', 'RTStopDate'});
                writetable(T_clin, fullfile(testCase.TempDir, 'mock_clinical.xlsx'));
            end

            % Mock Data vectors GTVp (unused by metrics_baseline directly but required input)
            testCase.DataVectorsGTVp = repmat( ...
                struct('adc', [1, 2, 3]', 'd', [1, 2, 3]', 'f', [0.1, 0.2, 0.3]', 'dstar', [0.01, 0.02, 0.03]'), ...
                nPat, 3);
            testCase.DataVectorsGTVn = repmat( ...
                struct('adc', [1.5, 2.5]', 'd', [1.5, 2.5]', 'f', [0.15, 0.25]', 'dstar', [0.015, 0.025]'), ...
                nPat, 3);

            % Mock Summary metrics — 8 patients, 3 timepoints
            testCase.SummaryMetrics = struct();
            testCase.SummaryMetrics.id_list    = arrayfun(@(x) sprintf('P%d', x), 1:nPat, 'UniformOutput', false);
            testCase.SummaryMetrics.mrn_list   = arrayfun(@(x) sprintf('M%d', x), 1:nPat, 'UniformOutput', false);
            testCase.SummaryMetrics.lf         = repmat([0; 1], nPat/2, 1);
            testCase.SummaryMetrics.gtv_vol    = repmat([10; 15; 8; 12; 11; 14; 9; 13], 1, 3);
            testCase.SummaryMetrics.dmean_gtvp = repmat([50; 55; 48; 52; 51; 54; 49; 53], 1, 3);
            testCase.SummaryMetrics.d95_gtvp   = repmat([45; 50; 43; 47; 46; 49; 44; 48], 1, 3);
            testCase.SummaryMetrics.v50gy_gtvp = repmat([90; 85; 88; 92; 89; 86; 87; 91], 1, 3);

            % Deterministic DWI metric arrays to avoid non-reproducible outlier detection
            testCase.SummaryMetrics.adc_mean  = repmat([1.0; 1.2; 0.8; 1.1; 0.9; 1.3; 0.85; 1.15], 1, 3) * 1e-3;
            testCase.SummaryMetrics.d_mean    = repmat([0.8; 1.0; 0.7; 0.9; 0.75; 1.05; 0.72; 0.95], 1, 3) * 1e-3;
            testCase.SummaryMetrics.f_mean    = repmat(0.20 * ones(nPat, 1), 1, 3);
            testCase.SummaryMetrics.dstar_mean = repmat(0.01 * ones(nPat, 1), 1, 3);
            testCase.SummaryMetrics.adc_sd    = repmat(0.1e-3 * ones(nPat, 1), 1, 3);

            % Repeatability fields
            testCase.SummaryMetrics.n_rpt          = 2 * ones(nPat, 1);
            testCase.SummaryMetrics.adc_mean_rpt   = rand(nPat, 2, 3);
            testCase.SummaryMetrics.adc_sub_rpt    = rand(nPat, 2, 3);
            testCase.SummaryMetrics.d_mean_rpt     = rand(nPat, 2, 3);
            testCase.SummaryMetrics.f_mean_rpt     = rand(nPat, 2, 3);
            testCase.SummaryMetrics.dstar_mean_rpt = rand(nPat, 2, 3);

            % Fields accessed by some metric modules
            testCase.SummaryMetrics.gtv_locations  = {};
            testCase.SummaryMetrics.dwi_locations  = {};
        end
    end

    methods(TestMethodTeardown)
        function cleanupTempDir(testCase)
            % Close any diary opened by metrics_baseline before deleting
            % the temp directory (avoids file-lock errors on Windows).
            diary off;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testBaselineExecution(testCase)
            % Smoke test: verifies metrics_baseline runs without error and
            % returns correctly dimensioned output arrays. Checks that all
            % 8 patients pass the valid_pts filter, and that absolute and
            % percent-change matrices have the expected [nPat x nTp] shape.
            [m_lf, m_total_time, m_total_follow_up_time, m_gtv_vol, m_adc_mean, m_d_mean, m_f_mean, m_dstar_mean, ...
             m_id_list, m_mrn_list, m_d95_gtvp, m_v50gy_gtvp, m_data_vectors_gtvp, lf_group, valid_pts, ...
             ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
             nTp, metric_sets, set_names, time_labels, dtype_label, dl_provenance] = ...
             metrics_baseline(testCase.DataVectorsGTVp, testCase.DataVectorsGTVn, testCase.SummaryMetrics, testCase.ConfigStruct);

            % Check that valid_pts masks correctly (8 patients, all with finite lf)
            testCase.verifyEqual(length(valid_pts), 8);
            testCase.verifyTrue(all(valid_pts));

            % Check output matrix dimensions
            testCase.verifyEqual(size(ADC_abs, 1), 8);  % 8 patients
            testCase.verifyEqual(size(ADC_pct, 1), 8);
            testCase.verifyEqual(size(ADC_abs, 2), 3);  % 3 timepoints
            testCase.verifyEqual(size(ADC_pct, 2), 3);  % same shape as ADC_abs
        end

        function testOutlierCleaning(testCase)
            % Verifies that the IQR-based outlier cleaning in metrics_baseline
            % replaces extreme values with NaN. Patient 1's baseline ADC is
            % set to 1e6 (far outside the 3*IQR fence for the 8-patient cohort),
            % and should be replaced by NaN while other patients remain intact.
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

        function testCauseOfDeathCompetingRisk(testCase)
            % Verifies competing risk classification logic:
            %   - Patient 1: no LF + lung cancer death -> lf=2 (competing risk)
            %   - Patient 2: no LF + pancreatic cancer death -> lf=0 (stays censored)
            %   - Patient 3: no LF + unknown cause -> lf=0 (stays censored)
            %   - Patient 5: LF + lung cancer death -> lf=1 (LF overrides)
            % This implements the Cause-Specific Hazards framework where
            % non-cancer-of-interest deaths are treated as competing events.
            if exist('OCTAVE_VERSION', 'builtin')
                testCase.assumeFail('Test requires MATLAB readtable/writetable.');
            end
            nPat = 8;
            id_list_xls = arrayfun(@(x) sprintf('P%d', x), (1:nPat)', 'UniformOutput', false);
            % Patient 1-4: no LF, Patient 5-8: have LF
            lf_vals     = [0; 0; 0; 0; 1; 1; 1; 1];
            dt_event    = repmat(datetime('2023-06-01'), nPat, 1);
            dt_censor   = repmat(datetime('2023-06-01'), nPat, 1);
            dt_reg      = repmat(datetime('2023-06-01'), nPat, 1);
            dt_rtstart  = repmat(datetime('2022-01-01'), nPat, 1);
            dt_rtstop   = repmat(datetime('2022-03-01'), nPat, 1);
            % Patient 1: lung cancer death (competing risk, no LF)
            % Patient 2: pancreatic cancer death (not competing)
            % Patient 3: unknown cause (not competing)
            % Patient 4: empty cause (not competing)
            % Patient 5: lung cancer death WITH LF (should stay lf==1)
            % Patient 6-8: empty cause
            cod_vals    = {'lung cancer'; 'pancreatic cancer'; 'unknown'; ''; ...
                           'lung cancer'; ''; ''; ''};
            T_clin = table(id_list_xls, lf_vals, dt_event, dt_censor, dt_reg, ...
                dt_rtstart, dt_rtstop, cod_vals, ...
                'VariableNames', {'Pat', 'LocalOrRegionalFailure', ...
                'LocoregionalFailureDateOfLocalOrRegionalFailure', ...
                'LocalFailureDateOfLocalFailureOrCensor', ...
                'RegionalFailureDateOfRegionalFailureOrCensor', ...
                'RTStartDate', 'RTStopDate', 'CauseOfDeath'});
            writetable(T_clin, fullfile(testCase.TempDir, 'mock_clinical.xlsx'));

            testCase.SummaryMetrics.lf = lf_vals;

            [m_lf, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~] = ...
             metrics_baseline(testCase.DataVectorsGTVp, testCase.DataVectorsGTVn, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);

            % Patient 1: no LF + non-pancreatic death => competing risk (2)
            testCase.verifyEqual(m_lf(1), 2, ...
                'Non-pancreatic-cancer death without LF should be coded as competing risk.');
            % Patient 2: no LF + pancreatic cancer death => stays 0
            testCase.verifyEqual(m_lf(2), 0, ...
                'Pancreatic cancer death should not be coded as competing risk.');
            % Patient 3: unknown cause => stays 0
            testCase.verifyEqual(m_lf(3), 0, ...
                'Unknown cause of death should not be coded as competing risk.');
            % Patient 5: LF + non-pancreatic death => stays 1 (LF overrides)
            testCase.verifyEqual(m_lf(5), 1, ...
                'Patient with LF should stay coded as event regardless of cause of death.');
        end

        function testCauseOfDeathNoWarningWhenPresent(testCase)
            % Verifies that when the CauseOfDeath column IS present in the
            % clinical spreadsheet, the metrics_baseline:noCauseOfDeath
            % warning is NOT emitted. This warning is only meant for cases
            % where the column is entirely absent.
            if exist('OCTAVE_VERSION', 'builtin')
                testCase.assumeFail('Test requires MATLAB readtable/writetable.');
            end
            nPat = 8;
            id_list_xls = arrayfun(@(x) sprintf('P%d', x), (1:nPat)', 'UniformOutput', false);
            lf_vals     = repmat([0; 1], nPat/2, 1);
            dt_event    = repmat(datetime('2023-06-01'), nPat, 1);
            dt_censor   = repmat(datetime('2023-06-01'), nPat, 1);
            dt_reg      = repmat(datetime('2023-06-01'), nPat, 1);
            dt_rtstart  = repmat(datetime('2022-01-01'), nPat, 1);
            dt_rtstop   = repmat(datetime('2022-03-01'), nPat, 1);
            % Use at least one non-empty value so readtable preserves the
            % column (all-empty text columns may be dropped by readtable).
            cod_vals    = [{'pancreatic cancer'}; repmat({''}, nPat-1, 1)];
            T_clin = table(id_list_xls, lf_vals, dt_event, dt_censor, dt_reg, ...
                dt_rtstart, dt_rtstop, cod_vals, ...
                'VariableNames', {'Pat', 'LocalOrRegionalFailure', ...
                'LocoregionalFailureDateOfLocalOrRegionalFailure', ...
                'LocalFailureDateOfLocalFailureOrCensor', ...
                'RegionalFailureDateOfRegionalFailureOrCensor', ...
                'RTStartDate', 'RTStopDate', 'CauseOfDeath'});
            writetable(T_clin, fullfile(testCase.TempDir, 'mock_clinical.xlsx'));

            % Capture warnings explicitly instead of verifyWarningFree to
            % avoid conflicts between the test framework's warning capture
            % and metrics_baseline's internal warning state management
            % (diary, warning('off', ...)).
            lastwarn('', '');
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~] = ...
             metrics_baseline(testCase.DataVectorsGTVp, testCase.DataVectorsGTVn, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);
            [~, warnId] = lastwarn;
            testCase.verifyNotEqual(warnId, 'metrics_baseline:noCauseOfDeath', ...
                'No noCauseOfDeath warning should be emitted when column is present.');
        end

        function testCustomCauseOfDeathColumnName(testCase)
            % Verifies that the cause_of_death_column config option allows
            % users to specify a non-default column name (e.g., 'DeathCause')
            % in their clinical spreadsheet. The function should find and
            % use this column for competing risk classification.
            if exist('OCTAVE_VERSION', 'builtin')
                testCase.assumeFail('Test requires MATLAB readtable/writetable.');
            end
            nPat = 8;
            id_list_xls = arrayfun(@(x) sprintf('P%d', x), (1:nPat)', 'UniformOutput', false);
            lf_vals     = [0; 0; 0; 0; 1; 1; 1; 1];
            dt_event    = repmat(datetime('2023-06-01'), nPat, 1);
            dt_censor   = repmat(datetime('2023-06-01'), nPat, 1);
            dt_reg      = repmat(datetime('2023-06-01'), nPat, 1);
            dt_rtstart  = repmat(datetime('2022-01-01'), nPat, 1);
            dt_rtstop   = repmat(datetime('2022-03-01'), nPat, 1);
            cod_vals    = {'lung cancer'; ''; ''; ''; ''; ''; ''; ''};
            % Use a non-standard column name
            T_clin = table(id_list_xls, lf_vals, dt_event, dt_censor, dt_reg, ...
                dt_rtstart, dt_rtstop, cod_vals, ...
                'VariableNames', {'Pat', 'LocalOrRegionalFailure', ...
                'LocoregionalFailureDateOfLocalOrRegionalFailure', ...
                'LocalFailureDateOfLocalFailureOrCensor', ...
                'RegionalFailureDateOfRegionalFailureOrCensor', ...
                'RTStartDate', 'RTStopDate', 'DeathCause'});
            writetable(T_clin, fullfile(testCase.TempDir, 'mock_clinical.xlsx'));

            testCase.SummaryMetrics.lf = lf_vals;
            cfg = testCase.ConfigStruct;
            cfg.cause_of_death_column = 'DeathCause';

            [m_lf, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~] = ...
             metrics_baseline(testCase.DataVectorsGTVp, testCase.DataVectorsGTVn, ...
                testCase.SummaryMetrics, cfg);

            % Patient 1: no LF + non-pancreatic death => competing risk (2)
            testCase.verifyEqual(m_lf(1), 2, ...
                'Custom column name should be recognized for competing risk classification.');
        end

        function testCaseInsensitiveCauseOfDeathColumn(testCase)
            % Verifies case-insensitive column matching: the config expects
            % 'CauseOfDeath' but the spreadsheet has 'causeofdeath'. The
            % function should still find the column and apply competing
            % risk classification without emitting a noCauseOfDeath warning.
            % This handles MATLAB's readtable behavior which may alter casing.
            if exist('OCTAVE_VERSION', 'builtin')
                testCase.assumeFail('Test requires MATLAB readtable/writetable.');
            end
            nPat = 8;
            id_list_xls = arrayfun(@(x) sprintf('P%d', x), (1:nPat)', 'UniformOutput', false);
            lf_vals     = [0; 0; 0; 0; 1; 1; 1; 1];
            dt_event    = repmat(datetime('2023-06-01'), nPat, 1);
            dt_censor   = repmat(datetime('2023-06-01'), nPat, 1);
            dt_reg      = repmat(datetime('2023-06-01'), nPat, 1);
            dt_rtstart  = repmat(datetime('2022-01-01'), nPat, 1);
            dt_rtstop   = repmat(datetime('2022-03-01'), nPat, 1);
            cod_vals    = {'lung cancer'; ''; ''; ''; ''; ''; ''; ''};
            % Use mismatched casing in the column name
            T_clin = table(id_list_xls, lf_vals, dt_event, dt_censor, dt_reg, ...
                dt_rtstart, dt_rtstop, cod_vals, ...
                'VariableNames', {'Pat', 'LocalOrRegionalFailure', ...
                'LocoregionalFailureDateOfLocalOrRegionalFailure', ...
                'LocalFailureDateOfLocalFailureOrCensor', ...
                'RegionalFailureDateOfRegionalFailureOrCensor', ...
                'RTStartDate', 'RTStopDate', 'causeofdeath'});
            writetable(T_clin, fullfile(testCase.TempDir, 'mock_clinical.xlsx'));

            testCase.SummaryMetrics.lf = lf_vals;

            % Config uses default 'CauseOfDeath' but spreadsheet has 'causeofdeath'
            lastwarn('', '');
            [m_lf, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~] = ...
             metrics_baseline(testCase.DataVectorsGTVp, testCase.DataVectorsGTVn, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);
            [~, warnId] = lastwarn;

            % Patient 1: no LF + non-pancreatic death => competing risk (2)
            testCase.verifyEqual(m_lf(1), 2, ...
                'Case-insensitive column match should enable competing risk classification.');
            testCase.verifyNotEqual(warnId, 'metrics_baseline:noCauseOfDeath', ...
                'No noCauseOfDeath warning should be emitted with case-insensitive match.');
        end

        function testStringScalarEmptyDisablesCauseOfDeath(testCase)
            % Verifies that setting cause_of_death_column to string("")
            % (the value jsondecode produces from JSON "") explicitly
            % disables competing risk classification without emitting a
            % warning. This is the expected behavior when the user sets
            % "cause_of_death_column": "" in config.json to opt out.
            if exist('OCTAVE_VERSION', 'builtin')
                testCase.assumeFail('Test requires MATLAB readtable/writetable.');
            end
            nPat = 8;
            id_list_xls = arrayfun(@(x) sprintf('P%d', x), (1:nPat)', 'UniformOutput', false);
            lf_vals     = [0; 0; 0; 0; 1; 1; 1; 1];
            dt_event    = repmat(datetime('2023-06-01'), nPat, 1);
            dt_censor   = repmat(datetime('2023-06-01'), nPat, 1);
            dt_reg      = repmat(datetime('2023-06-01'), nPat, 1);
            dt_rtstart  = repmat(datetime('2022-01-01'), nPat, 1);
            dt_rtstop   = repmat(datetime('2022-03-01'), nPat, 1);
            % Spreadsheet has NO CauseOfDeath column
            T_clin = table(id_list_xls, lf_vals, dt_event, dt_censor, dt_reg, ...
                dt_rtstart, dt_rtstop, ...
                'VariableNames', {'Pat', 'LocalOrRegionalFailure', ...
                'LocoregionalFailureDateOfLocalOrRegionalFailure', ...
                'LocalFailureDateOfLocalFailureOrCensor', ...
                'RegionalFailureDateOfRegionalFailureOrCensor', ...
                'RTStartDate', 'RTStopDate'});
            writetable(T_clin, fullfile(testCase.TempDir, 'mock_clinical.xlsx'));

            testCase.SummaryMetrics.lf = lf_vals;
            cfg = testCase.ConfigStruct;
            % Simulate jsondecode("") producing a string scalar
            cfg.cause_of_death_column = string("");

            lastwarn('', '');
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~] = ...
             metrics_baseline(testCase.DataVectorsGTVp, testCase.DataVectorsGTVn, ...
                testCase.SummaryMetrics, cfg);
            [~, warnId] = lastwarn;
            testCase.verifyNotEqual(warnId, 'metrics_baseline:noCauseOfDeath', ...
                'No noCauseOfDeath warning when cause_of_death_column is string("").');
        end

        function testNormalizedCauseOfDeathColumn(testCase)
            % Verifies the normalized fallback for column name matching:
            % the spreadsheet uses 'Cause_of_Death' (with underscores) but
            % the config expects 'CauseOfDeath'. After stripping underscores
            % and spaces and doing a case-insensitive compare, the column
            % should be found and competing risk classification should work.
            if exist('OCTAVE_VERSION', 'builtin')
                testCase.assumeFail('Test requires MATLAB readtable/writetable.');
            end
            nPat = 8;
            id_list_xls = arrayfun(@(x) sprintf('P%d', x), (1:nPat)', 'UniformOutput', false);
            lf_vals     = [0; 0; 0; 0; 1; 1; 1; 1];
            dt_event    = repmat(datetime('2023-06-01'), nPat, 1);
            dt_censor   = repmat(datetime('2023-06-01'), nPat, 1);
            dt_reg      = repmat(datetime('2023-06-01'), nPat, 1);
            dt_rtstart  = repmat(datetime('2022-01-01'), nPat, 1);
            dt_rtstop   = repmat(datetime('2022-03-01'), nPat, 1);
            cod_vals    = {'lung cancer'; ''; ''; ''; ''; ''; ''; ''};
            % Use underscore-separated column name
            T_clin = table(id_list_xls, lf_vals, dt_event, dt_censor, dt_reg, ...
                dt_rtstart, dt_rtstop, cod_vals, ...
                'VariableNames', {'Pat', 'LocalOrRegionalFailure', ...
                'LocoregionalFailureDateOfLocalOrRegionalFailure', ...
                'LocalFailureDateOfLocalFailureOrCensor', ...
                'RegionalFailureDateOfRegionalFailureOrCensor', ...
                'RTStartDate', 'RTStopDate', 'Cause_of_Death'});
            writetable(T_clin, fullfile(testCase.TempDir, 'mock_clinical.xlsx'));

            testCase.SummaryMetrics.lf = lf_vals;

            % Config uses default 'CauseOfDeath' but spreadsheet has 'Cause_of_Death'
            lastwarn('', '');
            [m_lf, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~, ~, ~, ...
             ~, ~, ~, ~, ~, ~] = ...
             metrics_baseline(testCase.DataVectorsGTVp, testCase.DataVectorsGTVn, ...
                testCase.SummaryMetrics, testCase.ConfigStruct);
            [~, warnId] = lastwarn;

            % Patient 1: no LF + non-pancreatic death => competing risk (2)
            testCase.verifyEqual(m_lf(1), 2, ...
                'Normalized column match should enable competing risk classification.');
            testCase.verifyNotEqual(warnId, 'metrics_baseline:noCauseOfDeath', ...
                'No noCauseOfDeath warning should be emitted with normalized match.');
        end
    end
end
