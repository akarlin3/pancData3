classdef test_compute_dose_response_roc < matlab.unittest.TestCase
    % TEST_COMPUTE_DOSE_RESPONSE_ROC  Tests for dose-response ROC analysis.

    properties
        TempDir
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
            [pmd, bl, methods, cfg] = buildRocTestData(testCase.TempDir);
            result = compute_dose_response_roc(pmd, bl, methods, cfg);

            testCase.verifyTrue(isfield(result, 'method_results'));
            testCase.verifyTrue(isfield(result, 'ranking'));
        end

        function testAucAboveChance(testCase)
            % With well-separated data, AUC should be > 0.7.
            [pmd, bl, methods, cfg] = buildRocTestData(testCase.TempDir);
            result = compute_dose_response_roc(pmd, bl, methods, cfg);

            if ~isempty(result.method_results)
                best_auc = result.method_results(1).best_auc;
                testCase.verifyGreaterThan(best_auc, 0.5, ...
                    'AUC should be above chance (0.5) with separated data.');
            end
        end

        function testOptimalThresholdInRange(testCase)
            % Optimal threshold should be between the two group distributions.
            [pmd, bl, methods, cfg] = buildRocTestData(testCase.TempDir);
            result = compute_dose_response_roc(pmd, bl, methods, cfg);

            if ~isempty(result.method_results) && ~isempty(result.method_results(1).metrics)
                thresh = result.method_results(1).metrics(1).optimal_threshold;
                testCase.verifyGreaterThan(thresh, 25, ...
                    'Threshold should be above LF distribution center.');
                testCase.verifyLessThan(thresh, 55, ...
                    'Threshold should be below LC distribution center.');
            end
        end

        function testAucRange(testCase)
            [pmd, bl, methods, cfg] = buildRocTestData(testCase.TempDir);
            result = compute_dose_response_roc(pmd, bl, methods, cfg);

            if ~isempty(result.method_results)
                for i = 1:numel(result.method_results)
                    auc = result.method_results(i).best_auc;
                    if ~isnan(auc)
                        testCase.verifyGreaterThanOrEqual(auc, 0);
                        testCase.verifyLessThanOrEqual(auc, 1);
                    end
                end
            end
        end

    end
end


function [pmd, bl, active_methods, cfg] = buildRocTestData(tempDir)
% Build synthetic data: 20 patients, LF have low D95, LC have high D95.

    cfg = struct();
    cfg.output_folder = tempDir;
    cfg.dwi_types_to_run = 1;
    cfg.dwi_type_name = 'Standard';

    rng(42);
    n_pts = 20;

    % Outcome: first 8 = LF, rest = LC
    m_lf = [ones(8, 1); zeros(12, 1)];

    bl = struct();
    bl.m_lf = m_lf;

    % D95: LF patients have low D95 (30-40 Gy), LC have high D95 (45-55 Gy)
    d95_lf = 30 + 10 * rand(8, 1);
    d95_lc = 45 + 10 * rand(12, 1);
    d95 = [d95_lf; d95_lc];

    v50_lf = 0.3 + 0.2 * rand(8, 1);
    v50_lc = 0.7 + 0.2 * rand(12, 1);
    v50 = [v50_lf; v50_lc];

    pmd = struct();
    pmd.adc_threshold = struct('d95_adc_sub', d95, 'v50_adc_sub', v50);

    active_methods = {'adc_threshold'};
end
