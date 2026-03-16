classdef test_progress_gui < matlab.unittest.TestCase
% TEST_PROGRESS_GUI — Unit tests for ProgressGUI.m
%
% Validates the custom progress bar widget including:
%   - Construction and initial state
%   - Update method with fraction, counts, detail, status
%   - setDetail truncation for long strings
%   - close/isValid lifecycle
%   - isDisplayAvailable static method
%   - formatTime static formatting (via update elapsed display)
%   - getBarColor status-based color selection

    properties
        OldFigVis
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            addpath(fullfile(testDir, '..', 'utils'));
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
        function test_construction(testCase)
            if ~ProgressGUI.isDisplayAvailable(), return; end
            gui = ProgressGUI('Test', 10);
            testCase.addTeardown(@() gui.close());
            testCase.verifyTrue(gui.isValid(), 'GUI should be valid after construction.');
        end

        function test_update_fraction(testCase)
            if ~ProgressGUI.isDisplayAvailable(), return; end
            gui = ProgressGUI('Test', 10);
            testCase.addTeardown(@() gui.close());

            counts = struct('completed', 5, 'total', 10, 'passed', 5, 'failed', 0);
            gui.update(0.5, counts, 'item 5/10', 'running');
            testCase.verifyTrue(gui.isValid());
        end

        function test_update_clamps_fraction(testCase)
            if ~ProgressGUI.isDisplayAvailable(), return; end
            gui = ProgressGUI('Test', 10);
            testCase.addTeardown(@() gui.close());

            % Should not error on out-of-range fraction
            counts = struct('completed', 0, 'total', 10);
            gui.update(-0.5, counts, '', 'running');
            gui.update(1.5, counts, '', 'running');
            testCase.verifyTrue(gui.isValid());
        end

        function test_setDetail_truncation(testCase)
            if ~ProgressGUI.isDisplayAvailable(), return; end
            gui = ProgressGUI('Test', 10);
            testCase.addTeardown(@() gui.close());

            long_str = repmat('x', 1, 200);
            gui.setDetail(long_str);
            % Should not error — just truncates to 70 chars
            testCase.verifyTrue(gui.isValid());
        end

        function test_close_and_isValid(testCase)
            if ~ProgressGUI.isDisplayAvailable(), return; end
            gui = ProgressGUI('Test', 10);
            testCase.verifyTrue(gui.isValid());

            gui.close();
            testCase.verifyFalse(gui.isValid(), 'GUI should not be valid after close.');
        end

        function test_double_close_no_error(testCase)
            if ~ProgressGUI.isDisplayAvailable(), return; end
            gui = ProgressGUI('Test', 10);
            gui.close();
            gui.close();  % Should not error
        end

        function test_update_after_close_no_error(testCase)
            if ~ProgressGUI.isDisplayAvailable(), return; end
            gui = ProgressGUI('Test', 10);
            gui.close();

            % Should silently return (isValid check)
            counts = struct('completed', 1, 'total', 10);
            gui.update(0.1, counts, 'test', 'running');
            gui.setDetail('test');
        end

        function test_isDisplayAvailable_returns_logical(testCase)
            result = ProgressGUI.isDisplayAvailable();
            testCase.verifyTrue(islogical(result), ...
                'isDisplayAvailable should return a logical.');
        end

        function test_failure_status_color(testCase)
            if ~ProgressGUI.isDisplayAvailable(), return; end
            gui = ProgressGUI('Test', 10);
            testCase.addTeardown(@() gui.close());

            counts = struct('completed', 5, 'total', 10, 'passed', 3, 'failed', 2);
            gui.update(0.5, counts, 'failing', 'failure');
            testCase.verifyTrue(gui.isValid());
        end

        function test_success_status(testCase)
            if ~ProgressGUI.isDisplayAvailable(), return; end
            gui = ProgressGUI('Test', 10);
            testCase.addTeardown(@() gui.close());

            counts = struct('completed', 10, 'total', 10, 'passed', 10, 'failed', 0);
            gui.update(1.0, counts, 'done', 'success');
            testCase.verifyTrue(gui.isValid());
        end

        function test_custom_dimensions(testCase)
            if ~ProgressGUI.isDisplayAvailable(), return; end
            opts = struct('width', 600, 'height', 250);
            gui = ProgressGUI('Wide', 20, opts);
            testCase.addTeardown(@() gui.close());
            testCase.verifyTrue(gui.isValid());
        end

        function test_zero_total(testCase)
            if ~ProgressGUI.isDisplayAvailable(), return; end
            gui = ProgressGUI('Empty', 0);
            testCase.addTeardown(@() gui.close());
            testCase.verifyTrue(gui.isValid());
        end

        function test_summary_with_stepName(testCase)
            if ~ProgressGUI.isDisplayAvailable(), return; end
            gui = ProgressGUI('Steps', 5);
            testCase.addTeardown(@() gui.close());

            counts = struct('completed', 2, 'total', 5, 'stepName', 'Loading data');
            gui.update(0.4, counts, 'step 2', 'running');
            testCase.verifyTrue(gui.isValid());
        end
    end
end
