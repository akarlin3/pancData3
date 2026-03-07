classdef test_pipeline_progress_gui < matlab.unittest.TestCase
    % TEST_PIPELINE_PROGRESS_GUI Unit tests for PipelineProgressGUI wrapper.
    % Validates step mapping, progress tracking, and edge cases.

    properties
        OldVisible
    end

    methods(TestMethodSetup)
        function suppressFigures(testCase)
            testCase.OldVisible = get(0, 'DefaultFigureVisible');
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function restoreFigures(testCase)
            set(0, 'DefaultFigureVisible', testCase.OldVisible);
        end
    end

    methods(Test)
        function test_constructor_with_valid_steps(testCase)
            gui = PipelineProgressGUI({'load', 'sanity'}, 'Standard');
            testCase.verifyTrue(gui.isValid());
            gui.close();
        end

        function test_constructor_with_empty_steps(testCase)
            % Steps not in STEP_MAP should result in no tracked steps
            gui = PipelineProgressGUI({'nonexistent_step'}, 'Standard');
            testCase.verifyFalse(gui.isValid());
            gui.close();
        end

        function test_constructor_with_all_known_steps(testCase)
            allSteps = {'load', 'sanity', 'metrics_baseline', 'compare_cores', ...
                'metrics_longitudinal', 'metrics_dosimetry', ...
                'metrics_stats_comparisons', 'metrics_stats_predictive', ...
                'visualize', 'metrics_survival'};
            gui = PipelineProgressGUI(allSteps, 'dnCNN');
            testCase.verifyTrue(gui.isValid());
            gui.close();
        end

        function test_start_and_complete_step(testCase)
            gui = PipelineProgressGUI({'load', 'sanity'}, 'Standard');
            gui.startStep('load');
            gui.completeStep('load', 'success');
            gui.startStep('sanity');
            gui.completeStep('sanity', 'success');
            % Should not error
            gui.close();
        end

        function test_complete_with_warning_status(testCase)
            gui = PipelineProgressGUI({'load', 'sanity'}, 'Standard');
            gui.startStep('load');
            gui.completeStep('load', 'warning');
            gui.close();
        end

        function test_complete_with_skipped_status(testCase)
            gui = PipelineProgressGUI({'load', 'sanity'}, 'Standard');
            gui.startStep('load');
            gui.completeStep('load', 'skipped');
            gui.close();
        end

        function test_start_unknown_step_does_not_error(testCase)
            gui = PipelineProgressGUI({'load'}, 'Standard');
            gui.startStep('nonexistent');
            gui.close();
        end

        function test_close_twice_does_not_error(testCase)
            gui = PipelineProgressGUI({'load'}, 'Standard');
            gui.close();
            gui.close();  % Should not error
        end

        function test_isValid_after_close(testCase)
            gui = PipelineProgressGUI({'load'}, 'Standard');
            gui.close();
            testCase.verifyFalse(gui.isValid());
        end
    end
end
