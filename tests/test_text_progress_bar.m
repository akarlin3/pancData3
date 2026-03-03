classdef test_text_progress_bar < matlab.unittest.TestCase
    % TEST_TEXT_PROGRESS_BAR Unit tests for the text progress bar utilities.

    methods (Test)

        function test_basic_call_no_error(testCase)
            % Verify basic call does not error
            text_progress_bar(1, 10, 'Test');
            fprintf('\n');  % clean up line
        end

        function test_completion_call(testCase)
            % Verify completion (current == total) prints newline
            text_progress_bar(5, 5, 'Done');
        end

        function test_zero_total(testCase)
            % Verify zero total does not error (no-op)
            text_progress_bar(0, 0, 'Empty');
        end

        function test_single_iteration(testCase)
            % Verify single iteration works
            text_progress_bar(1, 1, 'Single');
        end

        function test_default_label(testCase)
            % Verify default label when not provided
            text_progress_bar(1, 5);
            fprintf('\n');
        end

        function test_full_loop(testCase)
            % Run a complete loop and verify no errors
            for i = 1:10
                text_progress_bar(i, 10, 'Loop');
            end
        end

        function test_parfor_progress_factory(testCase)
            % Verify parfor_progress returns a callable function handle
            if exist('OCTAVE_VERSION', 'builtin')
                return;  % Skip in Octave (no DataQueue)
            end
            cb = parfor_progress(5, 'Test');
            testCase.verifyClass(cb, 'function_handle');
            % Call it manually to simulate afterEach callbacks
            cb(1);
            cb(2);
            fprintf('\n');
        end

    end
end
