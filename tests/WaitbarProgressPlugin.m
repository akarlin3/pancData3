classdef WaitbarProgressPlugin < matlab.unittest.plugins.TestRunnerPlugin
%WAITBARPROGRESSPLUGIN Displays a ProgressGUI dialog during test execution.
%   Updates a ProgressGUI figure as each test completes. Receives an
%   existing ProgressGUI instance and an offset count so that progress
%   spans both parallel and serial test phases.
%
%   Shows precise percentage, pass/fail counts, elapsed time, and the
%   current test name in a professional styled progress window.
%
%   Usage:
%       gui = ProgressGUI('Running Tests', totalTests);
%       runner.addPlugin(WaitbarProgressPlugin(gui, totalTests, offset));

    properties (Access = private)
        GUI                 % ProgressGUI handle for the graphical progress window
        TotalTests double   % Total test count across all phases (parallel + serial)
        Offset double       % Number of tests already completed before this plugin starts
        CompletedTests double = 0  % Tests completed by this plugin instance
        PassedTests double = 0     % Tests passed by this plugin instance
        FailedTests double = 0     % Tests failed by this plugin instance
    end

    methods
        function plugin = WaitbarProgressPlugin(gui, totalTests, offset)
        %WAITBARPROGRESSPLUGIN Constructor.
        %   gui        - An existing ProgressGUI instance to update.
        %   totalTests - Total test count (parallel + serial combined).
        %   offset     - Number of tests already completed (e.g., from a
        %                prior parallel phase), so progress starts where
        %                the parallel phase left off.
            plugin.GUI = gui;
            plugin.TotalTests = totalTests;
            plugin.Offset = offset;
            plugin.CompletedTests = 0;
        end
    end

    methods (Access = protected)
        function runTestSuite(plugin, pluginData)
        %RUNTESTSUITE Called once at the start of the serial phase.
        %   Initializes the progress bar at the offset fraction, delegates
        %   to the superclass to run all tests, then sets the final status
        %   to 'success' or 'failure' depending on results.
            if ~isempty(plugin.GUI) && plugin.GUI.isValid()
                frac = plugin.Offset / max(plugin.TotalTests, 1);
                counts = struct('completed', plugin.Offset, 'total', plugin.TotalTests, ...
                                'passed', plugin.PassedTests, 'failed', plugin.FailedTests);
                plugin.GUI.update(frac, counts, 'Running serial tests...', 'running');
            end

            runTestSuite@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);

            % Update GUI with final status after all serial tests complete
            if ~isempty(plugin.GUI) && plugin.GUI.isValid()
                total_done = plugin.Offset + plugin.CompletedTests;
                counts = struct('completed', total_done, 'total', plugin.TotalTests, ...
                                'passed', plugin.PassedTests, 'failed', plugin.FailedTests);
                if plugin.FailedTests > 0
                    plugin.GUI.update(1, counts, 'Done', 'failure');
                else
                    plugin.GUI.update(1, counts, 'Done', 'success');
                end
            end
        end

        function runTest(plugin, pluginData)
        %RUNTEST Called for each individual test. Shows the test name in
        %   the GUI detail line, runs the test, increments pass/fail
        %   counters, and updates the progress fraction.
            % Show current test name before running
            testName = pluginData.Name;
            if ~isempty(plugin.GUI) && plugin.GUI.isValid()
                plugin.GUI.setDetail(testName);
            end

            runTest@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);

            % Update counters after the test completes
            plugin.CompletedTests = plugin.CompletedTests + 1;
            if pluginData.TestResult.Failed
                plugin.FailedTests = plugin.FailedTests + 1;
            else
                plugin.PassedTests = plugin.PassedTests + 1;
            end

            % Compute overall progress (offset from parallel phase + serial progress)
            total_done = plugin.Offset + plugin.CompletedTests;
            if ~isempty(plugin.GUI) && plugin.GUI.isValid()
                frac = total_done / max(plugin.TotalTests, 1);
                counts = struct('completed', total_done, 'total', plugin.TotalTests, ...
                                'passed', plugin.PassedTests, 'failed', plugin.FailedTests);
                if plugin.FailedTests > 0
                    status = 'failure';
                else
                    status = 'running';
                end
                plugin.GUI.update(frac, counts, testName, status);
            end
        end
    end
end
