classdef test_text_progress_bar < matlab.unittest.TestCase
    % TEST_TEXT_PROGRESS_BAR Unit tests for the text progress bar utilities.
    %
    % Validates text_progress_bar.m and parfor_progress.m (both in utils/).
    % These utilities provide console-based progress reporting for pipeline
    % loops and parallel (parfor) processing. Tests verify that various
    % calling patterns do not error, including edge cases like zero-length
    % loops and single iterations.

    methods (Test)

        function test_basic_call_no_error(testCase)
            % Smoke test: a single mid-loop call (iteration 1 of 10) with a
            % custom label should execute without error.
            text_progress_bar(1, 10, 'Test');
            fprintf('\n');  % Clean up the carriage-return line
        end

        function test_completion_call(testCase)
            % Verify that reaching completion (current == total) prints a
            % newline and does not error. The progress bar should display 100%.
            text_progress_bar(5, 5, 'Done');
        end

        function test_zero_total(testCase)
            % Edge case: zero-length loop (total=0). The function should
            % handle this gracefully as a no-op without division by zero.
            text_progress_bar(0, 0, 'Empty');
        end

        function test_single_iteration(testCase)
            % Edge case: a loop with exactly one iteration. Both the start
            % and completion states are the same call (1 of 1 = 100%).
            text_progress_bar(1, 1, 'Single');
        end

        function test_default_label(testCase)
            % Verify the function works when the optional label argument is
            % omitted, falling back to a default label internally.
            text_progress_bar(1, 5);
            fprintf('\n');
        end

        function test_full_loop(testCase)
            % Integration test: run a complete 10-iteration loop to ensure
            % the progress bar updates correctly from 10% to 100% without
            % accumulating display artifacts.
            for i = 1:10
                text_progress_bar(i, 10, 'Loop');
            end
        end

        function test_parfor_progress_factory(testCase)
            % Verify that parfor_progress() returns a callable function handle
            % suitable for use as a DataQueue afterEach callback. The returned
            % handle is called manually here to simulate two parallel workers
            % reporting completion.
            if exist('OCTAVE_VERSION', 'builtin')
                return;  % Skip in Octave (no DataQueue support)
            end
            cb = parfor_progress(5, 'Test');
            testCase.verifyClass(cb, 'function_handle');
            % Simulate two afterEach callbacks (as if 2 parfor iterations completed)
            cb(1);
            cb(2);
            fprintf('\n');
        end

    end
end
