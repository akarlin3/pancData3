classdef test_nanmean_safe < matlab.unittest.TestCase
% TEST_NANMEAN_SAFE — Unit tests for nanmean_safe.m
%
% Validates Octave-compatible NaN-ignoring mean:
%   - Normal vectors, vectors with NaNs
%   - All-NaN and empty edge cases
%   - Single element

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
        function testNormalVector(testCase)
            result = nanmean_safe([1 2 3 4]);
            testCase.verifyEqual(result, 2.5, 'AbsTol', 1e-12, ...
                'Mean of [1 2 3 4] should be 2.5.');
        end

        function testWithNaNs(testCase)
            result = nanmean_safe([1 NaN 3]);
            testCase.verifyEqual(result, 2, 'AbsTol', 1e-12, ...
                'Mean of [1 NaN 3] should be 2 (ignoring NaN).');
        end

        function testAllNaN(testCase)
            result = nanmean_safe([NaN NaN]);
            testCase.verifyTrue(isnan(result), ...
                'All-NaN input should return NaN.');
        end

        function testEmptyInput(testCase)
            result = nanmean_safe([]);
            % MATLAB nanmean([]) returns NaN; Octave may return 0.
            testCase.verifyTrue(isnan(result) || result == 0, ...
                'Empty input should return NaN or 0.');
        end

        function testSingleElement(testCase)
            result = nanmean_safe([5]);
            testCase.verifyEqual(result, 5, 'AbsTol', 1e-12, ...
                'Single element [5] should return 5.');
        end
    end
end
