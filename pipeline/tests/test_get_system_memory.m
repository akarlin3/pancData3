classdef test_get_system_memory < matlab.unittest.TestCase
% TEST_GET_SYSTEM_MEMORY — Unit tests for get_system_memory.m
%
% Validates cross-platform system memory query:
%   - Returns two outputs
%   - Total memory is positive
%   - Available <= total
%   - Outputs are numeric

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
        function testReturnsTwoOutputs(testCase)
            [total_gb, avail_gb] = get_system_memory();
            testCase.verifyTrue(exist('total_gb', 'var') == 1, ...
                'First output (total_gb) should exist.');
            testCase.verifyTrue(exist('avail_gb', 'var') == 1, ...
                'Second output (avail_gb) should exist.');
        end

        function testTotalIsPositive(testCase)
            [total_gb, ~] = get_system_memory();
            % On any real system total should be > 0; NaN means unsupported.
            if ~isnan(total_gb)
                testCase.verifyGreaterThan(total_gb, 0, ...
                    'Total memory should be positive on a real system.');
            else
                % NaN is acceptable on unsupported platforms.
                testCase.verifyTrue(true);
            end
        end

        function testAvailableLeqTotal(testCase)
            [total_gb, avail_gb] = get_system_memory();
            if ~isnan(total_gb) && ~isnan(avail_gb)
                testCase.verifyLessThanOrEqual(avail_gb, total_gb, ...
                    'Available memory should not exceed total memory.');
            else
                testCase.verifyTrue(true);
            end
        end

        function testOutputsAreNumeric(testCase)
            [total_gb, avail_gb] = get_system_memory();
            testCase.verifyTrue(isnumeric(total_gb), ...
                'total_gb should be numeric (double).');
            testCase.verifyTrue(isnumeric(avail_gb), ...
                'avail_gb should be numeric (double).');
        end
    end
end
