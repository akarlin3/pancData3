classdef WaitbarProgressPlugin < matlab.unittest.plugins.TestRunnerPlugin
%WAITBARPROGRESSPLUGIN Displays a waitbar dialog during test execution.
%   Updates a MATLAB waitbar figure as each test completes. Receives an
%   existing waitbar handle and an offset count so that progress spans
%   both parallel and serial test phases.
%
%   Shows precise percentage, pass/fail counts, and the current test name.
%
%   Usage:
%       h = waitbar(0, 'Running tests...', 'Visible', 'on');
%       runner.addPlugin(WaitbarProgressPlugin(h, totalTests, offset));

    properties (Access = private)
        WaitbarHandle       % Handle to the waitbar figure
        TotalTests double   % Total test count (parallel + serial)
        Offset double       % Tests already completed (from parallel phase)
        CompletedTests double = 0
        PassedTests double = 0
        FailedTests double = 0
    end

    methods (Access = private, Static)
        function centerWaitbarText(hWaitbar)
        %CENTERWAITBARTEXT Re-center text objects after waitbar updates.
            if ~isvalid(hWaitbar); return; end
            hText = findobj(hWaitbar, 'Type', 'text');
            for ti = 1:numel(hText)
                set(hText(ti), 'Units', 'normalized', ...
                    'Position', [0.5, 0.5, 0], ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle');
            end
        end
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
                pct = frac * 100;
                waitbar(frac, plugin.WaitbarHandle, ...
                    sprintf('%.1f%% — Running serial tests... (%d/%d)', ...
                    pct, plugin.Offset, plugin.TotalTests));
                WaitbarProgressPlugin.centerWaitbarText(plugin.WaitbarHandle);
            end

            runTestSuite@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);

            if isvalid(plugin.WaitbarHandle)
                if plugin.FailedTests > 0
                    waitbar(1, plugin.WaitbarHandle, ...
                        sprintf('100.0%% — Done: %d passed, %d failed', ...
                        plugin.PassedTests, plugin.FailedTests));
                else
                    waitbar(1, plugin.WaitbarHandle, ...
                        sprintf('100.0%% — Done: %d passed', plugin.PassedTests));
                end
                WaitbarProgressPlugin.centerWaitbarText(plugin.WaitbarHandle);
            end
        end

        function runTest(plugin, pluginData)
            % Show current test name before running
            testName = pluginData.Name;
            total_done = plugin.Offset + plugin.CompletedTests;
            if isvalid(plugin.WaitbarHandle)
                frac = total_done / max(plugin.TotalTests, 1);
                pct = frac * 100;
                % Truncate long test names to fit the waitbar
                displayName = testName;
                if length(displayName) > 80
                    displayName = ['...' displayName(end-76:end)];
                end
                waitbar(frac, plugin.WaitbarHandle, ...
                    sprintf('%.1f%% (%d/%d) — %s', ...
                    pct, total_done, plugin.TotalTests, displayName));
                WaitbarProgressPlugin.centerWaitbarText(plugin.WaitbarHandle);
            end

            runTest@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);

            plugin.CompletedTests = plugin.CompletedTests + 1;
            if pluginData.TestResult.Failed
                plugin.FailedTests = plugin.FailedTests + 1;
            else
                plugin.PassedTests = plugin.PassedTests + 1;
            end

            total_done = plugin.Offset + plugin.CompletedTests;
            if isvalid(plugin.WaitbarHandle)
                frac = total_done / max(plugin.TotalTests, 1);
                pct = frac * 100;
                if plugin.FailedTests > 0
                    waitbar(frac, plugin.WaitbarHandle, ...
                        sprintf('%.1f%% (%d/%d) — %d passed, %d failed', ...
                        pct, total_done, plugin.TotalTests, ...
                        plugin.PassedTests, plugin.FailedTests));
                else
                    waitbar(frac, plugin.WaitbarHandle, ...
                        sprintf('%.1f%% (%d/%d) — %d passed', ...
                        pct, total_done, plugin.TotalTests, plugin.PassedTests));
                end
                WaitbarProgressPlugin.centerWaitbarText(plugin.WaitbarHandle);
            end
        end
    end
end
