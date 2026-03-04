classdef WaitbarProgressPlugin < matlab.unittest.plugins.TestRunnerPlugin
%WAITBARPROGRESSPLUGIN Displays a waitbar dialog during test execution.
%   Updates a MATLAB waitbar figure as each test completes. Receives an
%   existing waitbar handle and an offset count so that progress spans
%   both parallel and serial test phases.
%
%   Usage:
%       h = waitbar(0, 'Running tests...', 'Visible', 'on');
%       runner.addPlugin(WaitbarProgressPlugin(h, totalTests, offset));

    properties (Access = private)
        WaitbarHandle       % Handle to the waitbar figure
        TotalTests double   % Total test count (parallel + serial)
        Offset double       % Tests already completed (from parallel phase)
        CompletedTests double = 0
    end

    methods
        function plugin = WaitbarProgressPlugin(hWaitbar, totalTests, offset)
            plugin.WaitbarHandle = hWaitbar;
            plugin.TotalTests = totalTests;
            plugin.Offset = offset;
            plugin.CompletedTests = 0;
        end
    end

    methods (Access = protected)
        function runTestSuite(plugin, pluginData)
            if isvalid(plugin.WaitbarHandle)
                frac = plugin.Offset / max(plugin.TotalTests, 1);
                waitbar(frac, plugin.WaitbarHandle, ...
                    sprintf('Running serial tests... (%d/%d)', ...
                    plugin.Offset, plugin.TotalTests));
            end

            runTestSuite@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);

            if isvalid(plugin.WaitbarHandle)
                waitbar(1, plugin.WaitbarHandle, 'Tests complete!');
            end
        end

        function runTest(plugin, pluginData)
            runTest@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
            plugin.CompletedTests = plugin.CompletedTests + 1;
            total_done = plugin.Offset + plugin.CompletedTests;
            if isvalid(plugin.WaitbarHandle)
                frac = total_done / max(plugin.TotalTests, 1);
                waitbar(frac, plugin.WaitbarHandle, ...
                    sprintf('Test %d / %d', total_done, plugin.TotalTests));
            end
        end
    end
end
