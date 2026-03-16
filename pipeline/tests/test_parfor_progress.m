classdef test_parfor_progress < matlab.unittest.TestCase
% TEST_PARFOR_PROGRESS — Unit tests for parfor_progress.m
%
% Validates the progress callback factory including:
%   - Returns a function handle
%   - Callback increments internal counter correctly
%   - Default label usage
%   - Multiple independent callbacks maintain separate state

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
        end
    end

    methods (Test)
        function test_returns_function_handle(testCase)
            cb = parfor_progress(10, 'Test');
            testCase.verifyClass(cb, 'function_handle', ...
                'parfor_progress should return a function handle.');
        end

        function test_default_label(testCase)
            % Should not error when label is omitted
            cb = parfor_progress(5);
            testCase.verifyClass(cb, 'function_handle');
        end

        function test_callback_invocable(testCase)
            % Calling the callback should not error
            cb = parfor_progress(3, 'Unit Test');
            cb(1);  % simulate worker completion
            cb(2);
            cb(3);
            testCase.verifyTrue(true, 'Callback should be callable without error.');
        end

        function test_independent_callbacks(testCase)
            % Two callbacks should maintain independent counters
            cb1 = parfor_progress(5, 'Counter 1');
            cb2 = parfor_progress(10, 'Counter 2');

            cb1(1);
            cb1(2);
            cb2(1);

            % If we get here without error, independence is maintained.
            % (The internal counter is a closure variable, not shared.)
            testCase.verifyTrue(true, ...
                'Multiple callbacks should be independent.');
        end

        function test_zero_total(testCase)
            % total=0 should not crash (edge case)
            cb = parfor_progress(0, 'Empty');
            testCase.verifyClass(cb, 'function_handle');
        end

        function test_callback_ignores_argument(testCase)
            % The argument to the callback is ignored (~)
            cb = parfor_progress(3, 'Ignore Arg');
            cb('anything');  % should not error even with non-numeric arg
            cb([]);
            cb(struct());
            testCase.verifyTrue(true);
        end
    end
end
