classdef test_nanstd_safe < matlab.unittest.TestCase
% TEST_NANSTD_SAFE — Unit tests for nanstd_safe.m
%
% Validates Octave-compatible NaN-ignoring standard deviation:
%   - Normal vectors, vectors with NaNs
%   - All-NaN, single element, constant vector edge cases

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
            result = nanstd_safe([1 2 3 4 5]);
            expected = std([1 2 3 4 5]);
            testCase.verifyEqual(result, expected, 'AbsTol', 1e-12, ...
                'nanstd_safe should match std for vectors without NaN.');
        end

        function testWithNaNs(testCase)
            result = nanstd_safe([1 NaN 3]);
            expected = std([1 3]);
            testCase.verifyEqual(result, expected, 'AbsTol', 1e-12, ...
                'nanstd_safe([1 NaN 3]) should match std([1 3]).');
        end

        function testAllNaN(testCase)
            result = nanstd_safe([NaN NaN]);
            testCase.verifyTrue(isnan(result), ...
                'All-NaN input should return NaN.');
        end

        function testSingleElement(testCase)
            result = nanstd_safe([5]);
            testCase.verifyEqual(result, 0, 'AbsTol', 1e-12, ...
                'Single element should have std = 0.');
        end

        function testConstantVector(testCase)
            result = nanstd_safe([3 3 3]);
            testCase.verifyEqual(result, 0, 'AbsTol', 1e-12, ...
                'Constant vector should have std = 0.');
        end
    end
end
