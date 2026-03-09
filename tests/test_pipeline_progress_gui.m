classdef test_pipeline_progress_gui < matlab.unittest.TestCase
    % TEST_PIPELINE_PROGRESS_GUI Unit tests for PipelineProgressGUI wrapper.
    %
    % PipelineProgressGUI maps internal pipeline step keys (e.g., 'load',
    % 'sanity') to human-readable display names and wraps the ProgressGUI
    % figure-based progress bar. These tests validate:
    %   - Constructor behavior with valid, unknown, and all known steps
    %   - Step lifecycle (start, complete with success/warning/skipped)
    %   - Graceful handling of unknown steps and double-close
    %   - isValid() state tracking after close

    properties
        OldVisible  % Saved DefaultFigureVisible setting for restoration
    end

    methods(TestMethodSetup)
        function suppressFigures(testCase)
            % Hide figure windows during tests to prevent GUI pop-ups
            testCase.OldVisible = get(0, 'DefaultFigureVisible');
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function restoreFigures(testCase)
            % Restore original figure visibility setting
            set(0, 'DefaultFigureVisible', testCase.OldVisible);
        end
    end

    methods(Test)
        function test_constructor_with_valid_steps(testCase)
            % Known step keys ('load', 'sanity') should be recognized
            % by STEP_MAP, producing a valid GUI instance.
            gui = PipelineProgressGUI({'load', 'sanity'}, 'Standard');
            testCase.verifyTrue(gui.isValid());
            gui.close();
        end

        function test_constructor_with_empty_steps(testCase)
            % Steps not in STEP_MAP should be filtered out, resulting in
            % zero tracked steps and an invalid (no-op) GUI.
            gui = PipelineProgressGUI({'nonexistent_step'}, 'Standard');
            testCase.verifyFalse(gui.isValid());
            gui.close();
        end

        function test_constructor_with_all_known_steps(testCase)
            % All 10 recognized pipeline step keys should be accepted,
            % and the DWI type label ('dnCNN') should be set correctly.
            allSteps = {'load', 'sanity', 'metrics_baseline', 'compare_cores', ...
                'metrics_longitudinal', 'metrics_dosimetry', ...
                'metrics_stats_comparisons', 'metrics_stats_predictive', ...
                'visualize', 'metrics_survival'};
            gui = PipelineProgressGUI(allSteps, 'dnCNN');
            testCase.verifyTrue(gui.isValid());
            gui.close();
        end

        function test_start_and_complete_step(testCase)
            % Full lifecycle: start a step, mark it complete with 'success',
            % then start and complete the next step. No errors expected.
            gui = PipelineProgressGUI({'load', 'sanity'}, 'Standard');
            gui.startStep('load');
            gui.completeStep('load', 'success');
            gui.startStep('sanity');
            gui.completeStep('sanity', 'success');
            gui.close();
        end

        function test_complete_with_warning_status(testCase)
            % A step that completes with 'warning' status (non-fatal
            % module failure) should be handled without error.
            gui = PipelineProgressGUI({'load', 'sanity'}, 'Standard');
            gui.startStep('load');
            gui.completeStep('load', 'warning');
            gui.close();
        end

        function test_complete_with_skipped_status(testCase)
            % A step that is skipped (e.g., data already loaded) should
            % be handled without error.
            gui = PipelineProgressGUI({'load', 'sanity'}, 'Standard');
            gui.startStep('load');
            gui.completeStep('load', 'skipped');
            gui.close();
        end

        function test_start_unknown_step_does_not_error(testCase)
            % Starting a step key not in STEP_MAP should be a no-op
            % (graceful degradation), not an error.
            gui = PipelineProgressGUI({'load'}, 'Standard');
            gui.startStep('nonexistent');
            gui.close();
        end

        function test_close_twice_does_not_error(testCase)
            % Calling close() on an already-closed GUI should be
            % idempotent (safe to call multiple times).
            gui = PipelineProgressGUI({'load'}, 'Standard');
            gui.close();
            gui.close();
        end

        function test_isValid_after_close(testCase)
            % After close(), isValid() should return false to indicate
            % the GUI is no longer usable.
            gui = PipelineProgressGUI({'load'}, 'Standard');
            gui.close();
            testCase.verifyFalse(gui.isValid());
        end
    end
end
