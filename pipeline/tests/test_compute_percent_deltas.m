classdef test_compute_percent_deltas < matlab.unittest.TestCase
% TEST_COMPUTE_PERCENT_DELTAS — Unit tests for compute_percent_deltas.m
%
% Validates treatment-induced percent/absolute change computation:
%   - Basic percent change arithmetic
%   - Baseline column is zero change
%   - f uses absolute delta (not percent)
%   - Near-zero baselines produce NaN
%   - Winsorization at +/-500%

    properties
        origPath
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            testCase.origPath = path;
            baseDir = fullfile(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods (TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.origPath);
        end
    end

    methods (Test)
        function testBasicPercentChange(testCase)
            % ADC doubles from 0.001 to 0.002 => 100% change.
            ADC   = [0.001, 0.002];
            D     = [0.001, 0.001];
            f     = [0.10,  0.10];
            Dstar = [0.01,  0.01];
            evalc('[ADC_pct, ~, ~, ~] = compute_percent_deltas(ADC, D, f, Dstar);');
            testCase.verifyEqual(ADC_pct(1, 2), 100, 'AbsTol', 1e-6, ...
                'Doubling ADC from 0.001 to 0.002 should be 100% change.');
        end

        function testBaselineIsZero(testCase)
            % First column of all outputs should be 0 (change from self).
            ADC   = [0.001, 0.002; 0.002, 0.003];
            D     = [0.001, 0.002; 0.002, 0.003];
            f     = [0.10,  0.15;  0.12,  0.14];
            Dstar = [0.01,  0.02;  0.02,  0.03];
            evalc('[ADC_pct, D_pct, f_delta, Dstar_pct] = compute_percent_deltas(ADC, D, f, Dstar);');
            testCase.verifyEqual(ADC_pct(:, 1), [0; 0], 'AbsTol', 1e-12, ...
                'ADC baseline column should be 0.');
            testCase.verifyEqual(D_pct(:, 1), [0; 0], 'AbsTol', 1e-12, ...
                'D baseline column should be 0.');
            testCase.verifyEqual(f_delta(:, 1), [0; 0], 'AbsTol', 1e-12, ...
                'f_delta baseline column should be 0.');
            testCase.verifyEqual(Dstar_pct(:, 1), [0; 0], 'AbsTol', 1e-12, ...
                'D* baseline column should be 0.');
        end

        function testFUsesAbsoluteDelta(testCase)
            % f_delta should be absolute difference, not percent.
            ADC   = [0.001, 0.001];
            D     = [0.001, 0.001];
            f     = [0.10,  0.20];
            Dstar = [0.01,  0.01];
            evalc('[~, ~, f_delta, ~] = compute_percent_deltas(ADC, D, f, Dstar);');
            testCase.verifyEqual(f_delta(1, 2), 0.10, 'AbsTol', 1e-12, ...
                'f_delta should be absolute change (0.20 - 0.10 = 0.10).');
        end

        function testNearZeroBaselineExcluded(testCase)
            % Very small baseline (below epsilon) should produce NaN.
            ADC   = [1e-7, 0.002];   % baseline below adc_eps (1e-5)
            D     = [0.001, 0.002];
            f     = [0.10,  0.15];
            Dstar = [0.01,  0.02];
            evalc('[ADC_pct, ~, ~, ~] = compute_percent_deltas(ADC, D, f, Dstar);');
            testCase.verifyTrue(isnan(ADC_pct(1, 2)), ...
                'Near-zero ADC baseline should produce NaN percent change.');
        end

        function testWinsorization(testCase)
            % Changes exceeding 500% should be clipped to +/-500%.
            ADC   = [0.001, 0.010];   % 900% increase
            D     = [0.001, 0.001];
            f     = [0.10,  0.10];
            Dstar = [0.01,  0.01];
            evalc('[ADC_pct, ~, ~, ~] = compute_percent_deltas(ADC, D, f, Dstar);');
            testCase.verifyEqual(ADC_pct(1, 2), 500, 'AbsTol', 1e-6, ...
                'Percent change >500%% should be clipped to 500.');
        end
    end
end
