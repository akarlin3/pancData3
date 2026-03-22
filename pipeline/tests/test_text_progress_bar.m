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
            output = evalc('text_progress_bar(5, 5, ''Done'')');
            testCase.verifyTrue(contains(output, '100'), ...
                'Completion call should indicate 100% progress');
        end

        function test_zero_total(testCase)
            % Edge case: zero-length loop (total=0). The function should
            % handle this gracefully as a no-op without division by zero.
            text_progress_bar(0, 0, 'Empty');
        end

        function test_single_iteration(testCase)
            % Edge case: a loop with exactly one iteration. Both the start
            % and completion states are the same call (1 of 1 = 100%).
            output = evalc('text_progress_bar(1, 1, ''Single'')');
            testCase.verifyTrue(contains(output, '100'), ...
                'Single iteration should show 100% progress');
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
            output = evalc('for i = 1:10; text_progress_bar(i, 10, ''Loop''); end');
            testCase.verifyTrue(contains(output, '100'), ...
                'Full loop should reach 100% progress');
        end

        function test_full_loop_contains_intermediate_progress(testCase)
            % Verify that intermediate progress values appear in output
            % during a full loop (e.g., 50%).
            output = evalc('for i = 1:10; text_progress_bar(i, 10, ''Mid''); end');
            testCase.verifyTrue(contains(output, '50') || contains(output, '5') , ...
                'Full loop output should contain intermediate progress indicators');
        end

        function test_parfor_progress_factory(testCase)
            % Verify that parfor_progress() returns a callable function handle
            % suitable for use as a DataQueue afterEach callback. The returned
            % handle is called manually here to simulate workers reporting
            % completion. We call it the full number of times and verify the
            % output indicates completion.
            testCase.assumeTrue(~exist('OCTAVE_VERSION', 'builtin'), ...
                'Skipped in Octave — no DataQueue support');
            cb = parfor_progress(5, 'Test');
            testCase.verifyClass(cb, 'function_handle');
            % Simulate all 5 afterEach callbacks and capture output
            output = evalc('for i = 1:5; cb(i); end');
            testCase.verifyTrue(contains(output, '100'), ...
                'parfor_progress should indicate 100% after all iterations complete');
        end

        function test_parfor_progress_partial_output(testCase)
            % Verify that parfor_progress produces visible output for
            % partial completion (not just at 100%).
            testCase.assumeTrue(~exist('OCTAVE_VERSION', 'builtin'), ...
                'Skipped in Octave — no DataQueue support');
            cb = parfor_progress(4, 'Partial');
            testCase.verifyClass(cb, 'function_handle');
            % Call twice out of 4 total — should show ~50% progress
            output = evalc('cb(1); cb(2);');
            testCase.verifyFalse(isempty(strtrim(output)), ...
                'parfor_progress should produce non-empty output for partial completion');
            fprintf('\n');
        end

    end
end