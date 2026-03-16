classdef test_bootstrap_ci < matlab.unittest.TestCase
    % TEST_BOOTSTRAP_CI  Tests for BCa bootstrap confidence intervals.
    %
    % Covers:
    %   - Known-distribution recovery (normal mean)
    %   - Degenerate input (all identical values)
    %   - metric_fn that returns NaN for some resamples
    %   - Output dimensions

    properties
        OriginalPath
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testNormalMeanCIContainsTrueMean(testCase)
            % For a large sample from N(5, 1), the 95% CI for the mean
            % should contain the true mean most of the time.
            rng(42);
            true_mean = 5;
            data = true_mean + randn(200, 1);

            [ci_lo, ci_hi, boot_dist] = bootstrap_ci(data, @mean, 1000, 0.05);

            testCase.verifyLessThan(ci_lo, true_mean, ...
                'Lower CI should be below true mean.');
            testCase.verifyGreaterThan(ci_hi, true_mean, ...
                'Upper CI should be above true mean.');
            testCase.verifyEqual(length(boot_dist), 1000, ...
                'Boot distribution should have n_boot elements.');
        end

        function testCIOrderingCorrect(testCase)
            % CI lower bound should always be <= upper bound.
            rng(7);
            data = randn(50, 1);
            [ci_lo, ci_hi] = bootstrap_ci(data, @mean, 500);

            testCase.verifyLessThanOrEqual(ci_lo, ci_hi, ...
                'Lower CI should be <= upper CI.');
        end

        function testDegenerateInput(testCase)
            % All identical values: CI should collapse to the point estimate.
            data = ones(30, 1) * 3.14;
            [ci_lo, ci_hi] = bootstrap_ci(data, @mean, 500);

            testCase.verifyEqual(ci_lo, 3.14, 'AbsTol', 1e-10);
            testCase.verifyEqual(ci_hi, 3.14, 'AbsTol', 1e-10);
        end

        function testNaNResilience(testCase)
            % metric_fn that returns NaN for some resamples should still work.
            rng(42);
            data = [1; 2; 3; 4; 5; NaN; 7; 8];
            fn = @(x) nanmean(x);

            [ci_lo, ci_hi, boot_dist] = bootstrap_ci(data, fn, 500);

            testCase.verifyTrue(isfinite(ci_lo), 'Lower CI should be finite.');
            testCase.verifyTrue(isfinite(ci_hi), 'Upper CI should be finite.');
        end

        function testEmptyInput(testCase)
            % Empty data should return NaN.
            [ci_lo, ci_hi, boot_dist] = bootstrap_ci([], @mean, 100);

            testCase.verifyTrue(isnan(ci_lo));
            testCase.verifyTrue(isnan(ci_hi));
            testCase.verifyTrue(isempty(boot_dist));
        end

        function testMedianCI(testCase)
            % Test with median as the metric function.
            rng(42);
            data = exprnd(3, 100, 1);  % skewed distribution
            true_median = 3 * log(2);  % theoretical median of Exp(3)

            [ci_lo, ci_hi] = bootstrap_ci(data, @median, 1000);

            testCase.verifyLessThan(ci_lo, true_median + 1, ...
                'Lower CI should be reasonable.');
            testCase.verifyGreaterThan(ci_hi, true_median - 1, ...
                'Upper CI should be reasonable.');
        end

        function testDefaultArguments(testCase)
            % Test with default n_boot and alpha.
            rng(42);
            data = randn(50, 1);
            [ci_lo, ci_hi, boot_dist] = bootstrap_ci(data, @mean);

            testCase.verifyEqual(length(boot_dist), 2000, ...
                'Default n_boot should be 2000.');
            testCase.verifyTrue(isfinite(ci_lo));
            testCase.verifyTrue(isfinite(ci_hi));
        end

    end
end
