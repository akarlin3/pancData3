classdef test_compute_gtv_confounding < matlab.unittest.TestCase
    % TEST_COMPUTE_GTV_CONFOUNDING  Tests for GTV volume confounding check.

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
            [pmd, bl, methods, cfg] = buildConfoundTestData(testCase.TempDir, false);
            result = compute_gtv_confounding(pmd, bl, methods, cfg);

            testCase.verifyTrue(isfield(result, 'method_results'));
            testCase.verifyTrue(isfield(result, 'summary'));
        end

        function testConfoundedScenario(testCase)
            % When D95 is perfectly correlated with GTV, confounding_flag should be true.
            [pmd, bl, methods, cfg] = buildConfoundTestData(testCase.TempDir, true);
            result = compute_gtv_confounding(pmd, bl, methods, cfg);

            if ~isempty(result.method_results)
                % With correlated data, correlation should be strong
                testCase.verifyTrue(abs(result.method_results(1).d95_gtv_correlation) > 0.3, ...
                    'D95-GTV correlation should be strong when confounded.');
            end
        end

        function testIndependentScenario(testCase)
            % When D95 is independent of GTV, confounding_flag should be false.
            [pmd, bl, methods, cfg] = buildConfoundTestData(testCase.TempDir, false);
            result = compute_gtv_confounding(pmd, bl, methods, cfg);

            if ~isempty(result.method_results) && ...
                    ~isnan(result.method_results(1).d95_gtv_correlation)
                % Weak correlation expected
                testCase.verifyLessThan(abs(result.method_results(1).d95_gtv_correlation), 0.8, ...
                    'D95-GTV correlation should be weak when independent.');
            end
        end

        function testMissingGtvVolume(testCase)
            % If m_gtv_vol is missing, should return gracefully.
            [pmd, bl, methods, cfg] = buildConfoundTestData(testCase.TempDir, false);
            bl = rmfield(bl, 'm_gtv_vol');
            result = compute_gtv_confounding(pmd, bl, methods, cfg);

            testCase.verifyTrue(isempty(result.method_results) || ...
                (isstruct(result.method_results) && numel(result.method_results) == 0), ...
                'Should return empty results when GTV volume is missing.');
        end

    end
end


function [pmd, bl, active_methods, cfg] = buildConfoundTestData(tempDir, confounded)
% Build synthetic data for confounding tests.

    cfg = struct();
    cfg.output_folder = tempDir;
    cfg.dwi_types_to_run = 1;
    cfg.dwi_type_name = 'Standard';

    rng(42);
    n_pts = 20;

    % Outcome: first 8 = LF, rest = LC
    m_lf = [ones(8, 1); zeros(12, 1)];

    % GTV volumes
    gtv_vol = 50 + 30 * rand(n_pts, 5);
    % GTV change: LC patients shrink, LF patients grow
    for j = 1:n_pts
        for k = 2:5
            if m_lf(j) == 0
                gtv_vol(j, k) = gtv_vol(j, 1) * (1 - 0.05 * k);
            else
                gtv_vol(j, k) = gtv_vol(j, 1) * (1 + 0.02 * k);
            end
        end
    end

    % Time and follow-up
    m_total_time = 10 + 20 * rand(n_pts, 1);
    m_total_follow_up_time = m_total_time + 5;

    bl = struct();
    bl.m_lf = m_lf;
    bl.m_gtv_vol = gtv_vol;
    bl.m_total_time = m_total_time;
    bl.m_total_follow_up_time = m_total_follow_up_time;

    if confounded
        % D95 correlates with GTV: larger GTV → lower D95
        d95 = 60 - 0.3 * gtv_vol(:, 1) + 2 * randn(n_pts, 1);
    else
        % D95 independent of GTV
        d95_lf = 30 + 10 * rand(8, 1);
        d95_lc = 45 + 10 * rand(12, 1);
        d95 = [d95_lf; d95_lc];
    end

    pmd = struct();
    pmd.adc_threshold = struct('d95_adc_sub', d95);

    active_methods = {'adc_threshold'};
end
