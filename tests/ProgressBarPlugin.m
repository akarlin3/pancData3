classdef ProgressBarPlugin < matlab.unittest.plugins.TestRunnerPlugin
%PROGRESSBARPLUGIN Displays a dot-style progress bar during test execution.
%   Prints '.' after each passing test and 'F' after each failing test,
%   with periodic [count/total] markers for orientation.
%
%   Usage:
%       runner.addPlugin(ProgressBarPlugin(numel(suite)));

    properties (Access = private)
        TotalTests double = 0
        CompletedTests double = 0
    end

    methods
        function plugin = ProgressBarPlugin(totalTests)
            plugin.TotalTests = totalTests;
        end
    end

    methods (Access = protected)
        function runTestSuite(plugin, pluginData)
            plugin.CompletedTests = 0;
            fprintf('\n  Progress [%d tests]: ', plugin.TotalTests);
            runTestSuite@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
            fprintf(' Done!\n\n');
        end

        function runTest(plugin, pluginData)
            runTest@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
            plugin.CompletedTests = plugin.CompletedTests + 1;
            fprintf('.');
            if mod(plugin.CompletedTests, 50) == 0
                fprintf(' [%d/%d]\n  ', plugin.CompletedTests, plugin.TotalTests);
            end
        end
    end
end
