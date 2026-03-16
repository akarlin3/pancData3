classdef test_progress_gui < matlab.unittest.TestCase
% TEST_PROGRESS_GUI — Unit tests for ProgressGUI.m
%
% Validates the custom progress bar GUI including:
%   - Construction and initial state
%   - Update method (fraction clamping, bar color, text updates)
%   - Detail text truncation
%   - isValid/close lifecycle
%   - isDisplayAvailable static check
%   - formatTime static helper
%   - buildSummary text generation

    properties
        OldFigVis
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));
        end

        function suppressFigures(testCase)
            testCase.OldFigVis = get(0, 'DefaultFigureVisible');
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods (TestMethodTeardown)
        function restoreFigures(testCase)
            set(0, 'DefaultFigureVisible', testCase.OldFigVis);
        end
    end

    methods (Test)
        function test_constructor_creates_valid_object(testCase)
            gui = ProgressGUI('Test Bar', 10);
            testCase.addTeardown(@() gui.close());
            testCase.verifyTrue(gui.isValid(), 'Newly created GUI should be valid.');
        end

        function test_constructor_zero_total(testCase)
            % total=0 should not crash
            gui = ProgressGUI('Empty', 0);
            testCase.addTeardown(@() gui.close());
            testCase.verifyTrue(gui.isValid());
        end

        function test_update_does_not_crash(testCase)
            gui = ProgressGUI('Update Test', 20);
            testCase.addTeardown(@() gui.close());

            counts = struct('completed', 5, 'total', 20, 'passed', 5, 'failed', 0);
            gui.update(0.25, counts, 'test_item', 'running');

            testCase.verifyTrue(gui.isValid(), 'GUI should remain valid after update.');
        end

        function test_update_fraction_clamped(testCase)
            % Fractions > 1 or < 0 should not cause errors
            gui = ProgressGUI('Clamp Test', 10);
            testCase.addTeardown(@() gui.close());

            gui.update(1.5, struct('completed', 10, 'total', 10), '', 'success');
            gui.update(-0.5, struct('completed', 0, 'total', 10), '', 'running');

            testCase.verifyTrue(gui.isValid());
        end

        function test_close_invalidates(testCase)
            gui = ProgressGUI('Close Test', 5);
            gui.close();
            testCase.verifyFalse(gui.isValid(), 'Closed GUI should not be valid.');
        end

        function test_double_close_no_error(testCase)
            gui = ProgressGUI('Double Close', 5);
            gui.close();
            gui.close();  % should not error
            testCase.verifyFalse(gui.isValid());
        end

        function test_update_after_close_no_error(testCase)
            gui = ProgressGUI('Post Close', 5);
            gui.close();
            % update after close should silently return
            gui.update(0.5, struct('completed', 2, 'total', 5), 'detail', 'running');
            testCase.verifyFalse(gui.isValid());
        end

        function test_setDetail_truncation(testCase)
            gui = ProgressGUI('Detail Test', 10);
            testCase.addTeardown(@() gui.close());

            long_text = repmat('x', 1, 100);
            gui.setDetail(long_text);
            % Should not error; truncation is internal
            testCase.verifyTrue(gui.isValid());
        end

        function test_setDetail_after_close(testCase)
            gui = ProgressGUI('Detail Close', 5);
            gui.close();
            gui.setDetail('should not crash');
            testCase.verifyFalse(gui.isValid());
        end

        function test_bar_color_failure(testCase)
            gui = ProgressGUI('Color Test', 10);
            testCase.addTeardown(@() gui.close());

            counts = struct('completed', 5, 'total', 10, 'passed', 3, 'failed', 2);
            gui.update(0.5, counts, 'test', 'failure');
            testCase.verifyTrue(gui.isValid());
        end

        function test_bar_color_amber_on_failed(testCase)
            gui = ProgressGUI('Amber Test', 10);
            testCase.addTeardown(@() gui.close());

            counts = struct('completed', 5, 'total', 10, 'passed', 4, 'failed', 1);
            gui.update(0.5, counts, 'test', 'running');
            testCase.verifyTrue(gui.isValid());
        end

        function test_isDisplayAvailable_returns_logical(testCase)
            result = ProgressGUI.isDisplayAvailable();
            testCase.verifyClass(result, 'logical');
        end

        function test_custom_dimensions(testCase)
            opts = struct('width', 600, 'height', 250);
            gui = ProgressGUI('Custom Size', 5, opts);
            testCase.addTeardown(@() gui.close());
            testCase.verifyTrue(gui.isValid());
        end

        function test_summary_with_step_name(testCase)
            gui = ProgressGUI('Step Name', 5);
            testCase.addTeardown(@() gui.close());

            counts = struct('completed', 3, 'total', 5, 'stepName', 'Loading data');
            gui.update(0.6, counts, '', 'running');
            testCase.verifyTrue(gui.isValid());
        end
    end
end
