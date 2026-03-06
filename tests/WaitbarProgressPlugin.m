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
        GUI                 % ProgressGUI instance
        TotalTests double   % Total test count (parallel + serial)
        Offset double       % Tests already completed (from parallel phase)
        CompletedTests double = 0
        PassedTests double = 0
        FailedTests double = 0
    end

    methods
        function plugin = WaitbarProgressPlugin(gui, totalTests, offset)
            plugin.GUI = gui;
            plugin.TotalTests = totalTests;
            plugin.Offset = offset;
            plugin.CompletedTests = 0;
        end
    end

    methods (Access = protected)
        function runTestSuite(plugin, pluginData)
            if ~isempty(plugin.GUI) && plugin.GUI.isValid()
                frac = plugin.Offset / max(plugin.TotalTests, 1);
                counts = struct('completed', plugin.Offset, 'total', plugin.TotalTests, ...
                                'passed', plugin.PassedTests, 'failed', plugin.FailedTests);
                plugin.GUI.update(frac, counts, 'Running serial tests...', 'running');
            end

            runTestSuite@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);

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
            % Show current test name before running
            testName = pluginData.Name;
            if ~isempty(plugin.GUI) && plugin.GUI.isValid()
                plugin.GUI.setDetail(testName);
            end

            runTest@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);

            plugin.CompletedTests = plugin.CompletedTests + 1;
            if pluginData.TestResult.Failed
                plugin.FailedTests = plugin.FailedTests + 1;
            else
                plugin.PassedTests = plugin.PassedTests + 1;
            end

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
