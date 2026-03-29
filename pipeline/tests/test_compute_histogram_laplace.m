classdef test_compute_histogram_laplace < matlab.unittest.TestCase
% TEST_COMPUTE_HISTOGRAM_LAPLACE — Unit tests for compute_histogram_laplace.m
%
% Validates Laplace-smoothed histogram probability distribution:
%   - Uniform data produces roughly equal probabilities
%   - Empty bins get non-zero probability (Laplace smoothing)
%   - NaN exclusion
%   - Empty vector handling
%   - Probabilities sum to approximately 1

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
        function testUniformDistribution(testCase)
            % Uniform data across 5 bins should produce roughly equal probabilities.
            rng(42);
            vec = linspace(0, 1, 1000);
            bin_edges = 0:0.2:1;
            p = compute_histogram_laplace(vec, bin_edges);
            % With Laplace smoothing, each bin ~(200+1)/(1000+5) ~ 0.2
            testCase.verifyEqual(numel(p), 5, ...
                'Should have 5 bins for 6 edges.');
            for i = 1:numel(p)
                testCase.verifyGreaterThan(p(i), 0.15, ...
                    sprintf('Bin %d probability should be roughly 0.2.', i));
                testCase.verifyLessThan(p(i), 0.25, ...
                    sprintf('Bin %d probability should be roughly 0.2.', i));
            end
        end

        function testLaplaceSmoothingPreventZero(testCase)
            % Data concentrated in one bin; other bins should still be > 0.
            vec = ones(100, 1) * 0.5;
            bin_edges = [0, 0.25, 0.75, 1.0];
            p = compute_histogram_laplace(vec, bin_edges);
            testCase.verifyTrue(all(p > 0), ...
                'Laplace smoothing should prevent zero-probability bins.');
        end

        function testNaNsExcluded(testCase)
            % NaN values should be excluded from the histogram.
            vec = [1, 2, NaN, 3, NaN, 4];
            bin_edges = [0.5, 1.5, 2.5, 3.5, 4.5];
            p = compute_histogram_laplace(vec, bin_edges);
            % 4 data points + Laplace => sum should be ~1
            testCase.verifyEqual(sum(p), 1, 'AbsTol', 1e-10, ...
                'Probabilities should sum to 1 after NaN exclusion.');
        end

        function testEmptyVecReturnsZeros(testCase)
            % Empty vector should return all zeros.
            bin_edges = [0, 1, 2, 3];
            p = compute_histogram_laplace([], bin_edges);
            testCase.verifyEqual(p, zeros(1, 3), ...
                'Empty vector should return zeros.');
        end

        function testProbabilitiesSumToApproxOne(testCase)
            % Non-empty data: output should sum to approximately 1.
            rng(42);
            vec = randn(500, 1);
            bin_edges = -4:0.5:4;
            p = compute_histogram_laplace(vec, bin_edges);
            testCase.verifyEqual(sum(p), 1, 'AbsTol', 1e-10, ...
                'Laplace-smoothed probabilities should sum to 1.');
        end
    end
end
