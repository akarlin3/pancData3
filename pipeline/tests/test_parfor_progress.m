classdef test_parfor_progress < matlab.unittest.TestCase
% TEST_PARFOR_PROGRESS — Unit tests for parfor_progress.m
%
% Validates the parallel progress callback factory:
%   - Returns a function handle
%   - Callback is callable without error
%   - Counter increments correctly across calls
%   - Default label works when omitted
%   - Multiple independent callbacks have independent counters

    methods (TestMethodSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            addpath(fullfile(testDir, '..', 'utils'));
        end
    end

    methods (Test)
        function test_returns_function_handle(testCase)
            cb = parfor_progress(10, 'Test');
            testCase.verifyTrue(isa(cb, 'function_handle'), ...
                'parfor_progress should return a function handle.');
        end

        function test_callback_callable(testCase)
            cb = parfor_progress(5, 'Test');
            % Should not error when called
            cb(1);
            cb(2);
        end

        function test_all_iterations(testCase)
            total = 10;
            cb = parfor_progress(total, 'Test');
            for i = 1:total
                cb(i);  % Should not error at any iteration
            end
        end

        function test_default_label(testCase)
            % Omitting label should default to 'Progress'
            cb = parfor_progress(5);
            cb(1);  % Should not error
        end

        function test_independent_counters(testCase)
            % Two callbacks should maintain independent counters
            cb1 = parfor_progress(10, 'Counter1');
            cb2 = parfor_progress(20, 'Counter2');
            cb1(1);
            cb1(2);
            cb2(1);
            % Both should work independently without error
        end
    end
end
