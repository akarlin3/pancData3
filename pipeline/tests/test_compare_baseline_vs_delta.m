classdef test_compare_baseline_vs_delta < matlab.unittest.TestCase
    % TEST_COMPARE_BASELINE_VS_DELTA  Tests for compare_baseline_vs_delta.m

    properties
        TempDir
        BaselineResults
        ConfigStruct
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'core'));
            addpath(fullfile(repoRoot, 'utils'));
            if exist('OCTAVE_VERSION', 'builtin')
                addpath(fullfile(repoRoot, '.octave_compat'));
            end

            set(0, 'DefaultFigureVisible', 'off');

            [testCase.BaselineResults, testCase.ConfigStruct] = buildCompareData(testCase.TempDir);
        end
    end

    methods(TestMethodTeardown)
        function cleanup(testCase)
            diary off;
            close all;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)

        function testOutputStructFields(testCase)
            % Skip if coxphfit is unavailable (no Statistics Toolbox)
            if ~exist('coxphfit', 'file')
                testCase.assumeTrue(false, 'coxphfit not available; skipping.');
            end
            comparison = compare_baseline_vs_delta(testCase.BaselineResults, testCase.ConfigStruct);

            testCase.verifyTrue(isfield(comparison, 'parameters'));
            testCase.verifyTrue(isfield(comparison, 'timepoints'));
            testCase.verifyTrue(isfield(comparison, 'results'));
            testCase.verifyEqual(comparison.parameters, {'ADC', 'D', 'f', 'Dstar'});
            testCase.verifyEqual(comparison.timepoints, [2, 3]);
        end

        function testResultRowFields(testCase)
            if ~exist('coxphfit', 'file')
                testCase.assumeTrue(false, 'coxphfit not available; skipping.');
            end
            comparison = compare_baseline_vs_delta(testCase.BaselineResults, testCase.ConfigStruct);
            testCase.verifyGreaterThan(numel(comparison.results), 0);

            row = comparison.results(1);
            expected = {'parameter', 'timepoint', 'baseline_hr', 'baseline_p', ...
                'baseline_cindex', 'delta_hr', 'delta_p', 'delta_cindex', ...
                'better_predictor', 'n_events', 'n_patients'};
            for i = 1:numel(expected)
                testCase.verifyTrue(isfield(row, expected{i}), ...
                    sprintf('Result row missing field %s.', expected{i}));
            end
        end

        function testBetterPredictorDelta(testCase)
            if ~exist('coxphfit', 'file')
                testCase.assumeTrue(false, 'coxphfit not available; skipping.');
            end
            comparison = compare_baseline_vs_delta(testCase.BaselineResults, testCase.ConfigStruct);
            better_list = {comparison.results.better_predictor};
            testCase.verifyTrue(any(strcmp(better_list, 'delta')), ...
                'Given the synthetic data (delta is the true predictor), at least one row should favor delta.');
        end

        function testFigureSaved(testCase)
            if ~exist('coxphfit', 'file')
                testCase.assumeTrue(false, 'coxphfit not available; skipping.');
            end
            compare_baseline_vs_delta(testCase.BaselineResults, testCase.ConfigStruct);
            png_path = fullfile(testCase.TempDir, ...
                sprintf('baseline_vs_delta_%s.png', testCase.ConfigStruct.dwi_type_name));
            testCase.verifyTrue(exist(png_path, 'file') == 2, ...
                'Comparison PNG should be saved.');
        end

        function testMissingDpctHandled(testCase)
            % Should not error when D_pct is missing.
            br = struct();
            br.m_lf = zeros(10, 1);
            testCase.verifyWarningFree(@() compare_baseline_vs_delta(br, testCase.ConfigStruct));
        end

    end
end


function [baseline, cfg] = buildCompareData(tempDir)
% Build 20-patient dataset: 8 events, 12 censored. Delta has better
% separation than baseline (baseline nearly identical across groups;
% delta is much more negative for event patients).

    rng(100);
    n = 20;
    n_tp = 3;
    n_dwi = 1;

    % Outcome: first 8 are events (LF), last 12 are censored (LC)
    m_lf = [ones(8, 1); zeros(12, 1)];
    % Follow-up times: events have shorter times
    td_time = [5 + 2 * rand(8, 1); 20 + 5 * rand(12, 1)];

    % Baseline values: nearly identical across LF and LC (bad separator)
    base_adc = 1.2e-3 + 1e-5 * randn(n, 1);
    base_d   = 1.1e-3 + 1e-5 * randn(n, 1);
    base_f   = 0.12   + 0.005 * randn(n, 1);
    base_dstar = 0.02 + 0.0005 * randn(n, 1);

    % Delta (% change): strongly different between groups (good separator).
    % LF group: large negative delta; LC group: small positive delta.
    delta_adc = [-30 + 3 * randn(8, 1); 10 + 3 * randn(12, 1)];
    delta_d   = [-35 + 4 * randn(8, 1); 12 + 4 * randn(12, 1)];
    delta_f   = [-0.05 + 0.01 * randn(8, 1); 0.02 + 0.01 * randn(12, 1)];
    delta_dstar = [-30 + 4 * randn(8, 1); 8 + 4 * randn(12, 1)];

    % Assemble 3D arrays
    ADC_abs   = nan(n, n_tp, n_dwi); ADC_abs(:, 1, 1)   = base_adc;
    D_abs     = nan(n, n_tp, n_dwi); D_abs(:, 1, 1)     = base_d;
    f_abs     = nan(n, n_tp, n_dwi); f_abs(:, 1, 1)     = base_f;
    Dstar_abs = nan(n, n_tp, n_dwi); Dstar_abs(:, 1, 1) = base_dstar;

    ADC_pct   = nan(n, n_tp, n_dwi);
    ADC_pct(:, 2, 1) = delta_adc;
    ADC_pct(:, 3, 1) = delta_adc + 2 * randn(n, 1);

    D_pct = nan(n, n_tp, n_dwi);
    D_pct(:, 2, 1) = delta_d;
    D_pct(:, 3, 1) = delta_d + 2 * randn(n, 1);

    f_delta = nan(n, n_tp, n_dwi);
    f_delta(:, 2, 1) = delta_f;
    f_delta(:, 3, 1) = delta_f + 0.005 * randn(n, 1);

    Dstar_pct = nan(n, n_tp, n_dwi);
    Dstar_pct(:, 2, 1) = delta_dstar;
    Dstar_pct(:, 3, 1) = delta_dstar + 2 * randn(n, 1);

    baseline = struct();
    baseline.m_lf = m_lf;
    baseline.td_time = td_time;
    baseline.m_total_time = td_time;
    baseline.ADC_abs = ADC_abs;
    baseline.D_abs = D_abs;
    baseline.f_abs = f_abs;
    baseline.Dstar_abs = Dstar_abs;
    baseline.ADC_pct = ADC_pct;
    baseline.D_pct = D_pct;
    baseline.f_delta = f_delta;
    baseline.Dstar_pct = Dstar_pct;

    cfg = struct();
    cfg.output_folder = tempDir;
    cfg.dwi_types_to_run = 1;
    cfg.dwi_type_name = 'Standard';
end
