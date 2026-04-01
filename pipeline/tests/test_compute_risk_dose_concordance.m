classdef test_compute_risk_dose_concordance < matlab.unittest.TestCase
    % TEST_COMPUTE_RISK_DOSE_CONCORDANCE  Tests for risk-dose concordance analysis.

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
            [pred, pmd, bl, methods, cfg] = buildConcordanceTestData(testCase.TempDir, 'concordant');
            result = compute_risk_dose_concordance(pred, pmd, bl, methods, cfg);

            testCase.verifyTrue(isfield(result, 'method_results'));
            testCase.verifyTrue(isfield(result, 'summary'));
        end

        function testPerfectConcordance(testCase)
            % When same patients flagged by both, concordance should be high.
            [pred, pmd, bl, methods, cfg] = buildConcordanceTestData(testCase.TempDir, 'concordant');
            result = compute_risk_dose_concordance(pred, pmd, bl, methods, cfg);

            if ~isempty(result.method_results) && ~isnan(result.method_results(1).concordance_pct)
                testCase.verifyGreaterThan(result.method_results(1).concordance_pct, 50, ...
                    'Concordant data should show >50% concordance.');
            end
        end

        function testComplementaryCase(testCase)
            % When no overlap, concordance should be low.
            [pred, pmd, bl, methods, cfg] = buildConcordanceTestData(testCase.TempDir, 'complementary');
            result = compute_risk_dose_concordance(pred, pmd, bl, methods, cfg);

            if ~isempty(result.method_results) && ~isnan(result.method_results(1).concordance_pct)
                testCase.verifyLessThan(result.method_results(1).concordance_pct, 80, ...
                    'Complementary data should show lower concordance.');
            end
        end

        function testKappaRange(testCase)
            [pred, pmd, bl, methods, cfg] = buildConcordanceTestData(testCase.TempDir, 'concordant');
            result = compute_risk_dose_concordance(pred, pmd, bl, methods, cfg);

            if ~isempty(result.method_results)
                kappa = result.method_results(1).cohen_kappa;
                if ~isnan(kappa)
                    testCase.verifyGreaterThanOrEqual(kappa, -1);
                    testCase.verifyLessThanOrEqual(kappa, 1);
                end
            end
        end

        function testNoPredictiveResults(testCase)
            % Missing is_high_risk should return gracefully.
            [pred, pmd, bl, methods, cfg] = buildConcordanceTestData(testCase.TempDir, 'concordant');
            pred.is_high_risk = nan(20, 1);
            result = compute_risk_dose_concordance(pred, pmd, bl, methods, cfg);

            testCase.verifyTrue(contains(result.summary, 'not available') || ...
                contains(result.summary, 'No'));
        end

    end
end


function [pred, pmd, bl, active_methods, cfg] = buildConcordanceTestData(tempDir, scenario)
% Build synthetic data for concordance tests.

    cfg = struct();
    cfg.output_folder = tempDir;
    cfg.dwi_types_to_run = 1;
    cfg.dwi_type_name = 'Standard';

    rng(42);
    n_pts = 20;

    m_lf = [ones(8, 1); zeros(12, 1)];

    bl = struct();
    bl.m_lf = m_lf;
    bl.m_total_time = 10 + 20 * rand(n_pts, 1);

    switch scenario
        case 'concordant'
            % High-risk patients = patients with low D95 (same patients)
            is_high_risk = [true(8, 1); false(12, 1)];
            risk_scores = [0.7 + 0.2*rand(8,1); 0.2 + 0.2*rand(12,1)];
            d95 = [30 + 5*rand(8,1); 50 + 5*rand(12,1)];  % low for high-risk

        case 'complementary'
            % High-risk from model doesn't match low-D95
            is_high_risk = [false(4, 1); true(4, 1); true(6, 1); false(6, 1)];
            risk_scores = [0.3*rand(4,1); 0.7+0.2*rand(4,1); 0.7+0.2*rand(6,1); 0.3*rand(6,1)];
            d95 = [30+5*rand(4,1); 50+5*rand(4,1); 30+5*rand(6,1); 50+5*rand(6,1)];
    end

    pred = struct();
    pred.is_high_risk = is_high_risk;
    pred.risk_scores_all = risk_scores;

    pmd = struct();
    pmd.adc_threshold = struct('d95_adc_sub', d95);

    active_methods = {'adc_threshold'};
end
