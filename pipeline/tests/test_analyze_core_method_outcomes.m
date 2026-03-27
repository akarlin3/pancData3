classdef test_analyze_core_method_outcomes < matlab.unittest.TestCase
    % TEST_ANALYZE_CORE_METHOD_OUTCOMES  Tests for core method outcome analysis.

    properties
        TempDir
        ConfigStruct
        PerMethodDos
        BaselineResults
        ActiveMethods
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

            [testCase.ConfigStruct, testCase.PerMethodDos, ...
                testCase.BaselineResults, testCase.ActiveMethods] = buildTestData(testCase.TempDir);
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
            results = analyze_core_method_outcomes(testCase.PerMethodDos, ...
                testCase.BaselineResults, testCase.ActiveMethods, testCase.ConfigStruct);

            testCase.verifyTrue(isfield(results, 'method_results'));
            testCase.verifyTrue(isfield(results, 'ranking'));
            testCase.verifyTrue(isfield(results, 'summary_table'));
            testCase.verifyTrue(isfield(results, 'active_methods'));
        end

        function testMethodResultsCount(testCase)
            results = analyze_core_method_outcomes(testCase.PerMethodDos, ...
                testCase.BaselineResults, testCase.ActiveMethods, testCase.ConfigStruct);

            % Both adc_threshold and otsu have dosimetry data
            testCase.verifyEqual(numel(results.method_results), 2, ...
                'Should have results for 2 methods with dosimetry data.');
        end

        function testUnivariableHRPositive(testCase)
            results = analyze_core_method_outcomes(testCase.PerMethodDos, ...
                testCase.BaselineResults, testCase.ActiveMethods, testCase.ConfigStruct);

            for i = 1:numel(results.method_results)
                mr = results.method_results(i);
                if ~isempty(mr.univariable)
                    for j = 1:numel(mr.univariable)
                        testCase.verifyGreaterThan(mr.univariable(j).hr, 0, ...
                            'HR must be positive (exp of any real number).');
                    end
                end
            end
        end

        function testKMFieldsPresent(testCase)
            results = analyze_core_method_outcomes(testCase.PerMethodDos, ...
                testCase.BaselineResults, testCase.ActiveMethods, testCase.ConfigStruct);

            for i = 1:numel(results.method_results)
                mr = results.method_results(i);
                testCase.verifyTrue(isfield(mr.km, 'logrank_p'));
                testCase.verifyTrue(isfield(mr.km, 'best_metric'));
            end
        end

        function testRankingOrder(testCase)
            results = analyze_core_method_outcomes(testCase.PerMethodDos, ...
                testCase.BaselineResults, testCase.ActiveMethods, testCase.ConfigStruct);

            testCase.verifyGreaterThanOrEqual(numel(results.ranking), 1, ...
                'Should have at least one method in ranking.');
        end

        function testClearSeparationSignificant(testCase)
            % With synthetic data designed for clear dose-outcome separation,
            % the best method should have p < 0.05.
            results = analyze_core_method_outcomes(testCase.PerMethodDos, ...
                testCase.BaselineResults, testCase.ActiveMethods, testCase.ConfigStruct);

            % Find best p-value across all methods
            best_p = Inf;
            for i = 1:numel(results.method_results)
                mr = results.method_results(i);
                if ~isempty(mr.univariable)
                    p_vals = [mr.univariable.p_value];
                    best_p = min(best_p, min(p_vals));
                end
            end
            testCase.verifyLessThan(best_p, 0.05, ...
                'With clear separation, best p should be < 0.05.');
        end

        function testInsufficientEvents(testCase)
            % No events -> should return with empty results.
            br = testCase.BaselineResults;
            br.m_lf = zeros(size(br.m_lf));  % all censored

            results = analyze_core_method_outcomes(testCase.PerMethodDos, ...
                br, testCase.ActiveMethods, testCase.ConfigStruct);

            testCase.verifyEmpty(results.method_results, ...
                'Should return empty method_results with no events.');
        end

        function testEmptyPerMethodDosimetry(testCase)
            results = analyze_core_method_outcomes(struct(), ...
                testCase.BaselineResults, testCase.ActiveMethods, testCase.ConfigStruct);

            testCase.verifyEmpty(results.method_results);
        end

        function testMissingMethodInDosimetry(testCase)
            % Active methods includes 'nonexistent' which is not in dosimetry
            active = [testCase.ActiveMethods, {'nonexistent_method'}];

            results = analyze_core_method_outcomes(testCase.PerMethodDos, ...
                testCase.BaselineResults, active, testCase.ConfigStruct);

            % Should still have results for 2 real methods, skip nonexistent
            testCase.verifyEqual(numel(results.method_results), 2);
        end

    end
end


function [cfg, pmd, br, active] = buildTestData(tempDir)
% Build synthetic test data with clear dose-outcome separation.

    rng(42);
    n_patients = 20;
    nTp = 3;

    cfg = struct();
    cfg.output_folder = tempDir;
    cfg.dwi_type_name = 'Standard';
    cfg.active_core_methods = {'adc_threshold', 'otsu'};

    active = {'adc_threshold', 'otsu'};

    % --- Baseline results ---
    br = struct();
    br.valid_pts = true(n_patients, 1);
    br.nTp = nTp;
    br.m_id_list = arrayfun(@(x) sprintf('P%d', x), 1:n_patients, 'UniformOutput', false);

    % Outcomes: 10 LC (0), 8 LF (1), 2 competing (2)
    br.m_lf = zeros(n_patients, 1);
    br.m_lf(1:8) = 1;    % LF events
    br.m_lf(19:20) = 2;  % competing risk

    % Times: LF patients fail early, LC patients censored late
    br.m_total_time = zeros(n_patients, 1);
    br.m_total_time(1:8) = 60 + 60*rand(8, 1);      % LF: 60-120 days
    br.m_total_time(9:18) = 200 + 200*rand(10, 1);   % LC: 200-400 days
    br.m_total_time(19:20) = 100 + 50*rand(2, 1);    % competing: 100-150

    br.m_total_follow_up_time = br.m_total_time;
    br.m_total_follow_up_time(9:18) = br.m_total_time(9:18) + 100;  % longer follow-up for LC

    % --- Per method dosimetry ---
    % Design: LF patients (1-8) get LOW dose coverage, LC patients (9-18) get HIGH
    pmd = struct();

    % adc_threshold: clear separation
    pmd.adc_threshold.d95_adc_sub = nan(n_patients, nTp);
    pmd.adc_threshold.v50_adc_sub = nan(n_patients, nTp);
    pmd.adc_threshold.d95_d_sub = nan(n_patients, nTp);
    pmd.adc_threshold.v50_d_sub = nan(n_patients, nTp);

    % Fx1 values: LF patients get 30-40 Gy, LC patients get 45-55 Gy
    pmd.adc_threshold.d95_adc_sub(1:8, 1) = 30 + 10*rand(8, 1);
    pmd.adc_threshold.d95_adc_sub(9:20, 1) = 45 + 10*rand(12, 1);
    pmd.adc_threshold.v50_adc_sub(:, 1) = 0.5 + 0.4*rand(n_patients, 1);
    pmd.adc_threshold.d95_d_sub(:, 1) = pmd.adc_threshold.d95_adc_sub(:, 1) + 2*randn(n_patients, 1);
    pmd.adc_threshold.v50_d_sub(:, 1) = 0.6 + 0.3*rand(n_patients, 1);

    % otsu: also clear separation but slightly noisier
    pmd.otsu.d95_adc_sub = nan(n_patients, nTp);
    pmd.otsu.v50_adc_sub = nan(n_patients, nTp);
    pmd.otsu.d95_d_sub = nan(n_patients, nTp);
    pmd.otsu.v50_d_sub = nan(n_patients, nTp);

    pmd.otsu.d95_adc_sub(1:8, 1) = 32 + 10*rand(8, 1);
    pmd.otsu.d95_adc_sub(9:20, 1) = 43 + 12*rand(12, 1);
    pmd.otsu.v50_adc_sub(:, 1) = 0.5 + 0.4*rand(n_patients, 1);
    pmd.otsu.d95_d_sub(:, 1) = pmd.otsu.d95_adc_sub(:, 1) + 3*randn(n_patients, 1);
    pmd.otsu.v50_d_sub(:, 1) = 0.6 + 0.3*rand(n_patients, 1);
end
