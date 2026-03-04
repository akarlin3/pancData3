classdef ProgressBarPlugin < matlab.unittest.plugins.TestRunnerPlugin
%PROGRESSBARPLUGIN Displays a visual progress bar during test execution.
%   Shows an in-place progress bar with percentage, pass/fail counts,
%   and the name of the currently running test.
%
%   Usage:
%       runner.addPlugin(ProgressBarPlugin(numel(suite)));

    properties (Access = private)
        TotalTests double = 0
        CompletedTests double = 0
        PassedTests double = 0
        FailedTests double = 0
    end

    methods
        function plugin = ProgressBarPlugin(totalTests)
            plugin.TotalTests = totalTests;
        end
    end

    methods (Access = protected)
        function runTestSuite(plugin, pluginData)
            plugin.CompletedTests = 0;
            plugin.PassedTests = 0;
            plugin.FailedTests = 0;
            fprintf('\n');
            runTestSuite@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
            % Print final complete bar
            plugin.printBar('');
            if plugin.FailedTests > 0
                fprintf('  Result: %d passed, %d FAILED\n\n', ...
                    plugin.PassedTests, plugin.FailedTests);
            else
                fprintf('  Result: %d passed\n\n', plugin.PassedTests);
            end
        end

        function runTest(plugin, pluginData)
            % Show which test is about to run
            testName = pluginData.Name;
            plugin.printBar(testName);

            runTest@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);

            plugin.CompletedTests = plugin.CompletedTests + 1;
            if pluginData.TestResult.Failed
                plugin.FailedTests = plugin.FailedTests + 1;
            else
                plugin.PassedTests = plugin.PassedTests + 1;
            end
            plugin.printBar('');
        end
    end

    methods (Access = private)
        function printBar(plugin, currentTest)
            bar_width = 30;
            total = max(plugin.TotalTests, 1);
            fraction = plugin.CompletedTests / total;
            filled = round(bar_width * fraction);
            empty_count = bar_width - filled;

            fill_char = char(9608);   % U+2588 FULL BLOCK
            empty_char = char(9617);  % U+2591 LIGHT SHADE

            bar_str = [repmat(fill_char, 1, filled), ...
                       repmat(empty_char, 1, empty_count)];
            pct = fraction * 100;

            if isempty(currentTest)
                status = '';
            else
                % Truncate long test names
                maxLen = 40;
                if length(currentTest) > maxLen
                    currentTest = ['...' currentTest(end-maxLen+4:end)];
                end
                status = sprintf(' | %s', currentTest);
            end

            if plugin.FailedTests > 0
                counts = sprintf('%d/%d (%d failed)', ...
                    plugin.CompletedTests, plugin.TotalTests, plugin.FailedTests);
            else
                counts = sprintf('%d/%d', ...
                    plugin.CompletedTests, plugin.TotalTests);
            end

            fprintf('\r  Tests: |%s| %5.1f%% %s%s', ...
                bar_str, pct, counts, status);

            if plugin.CompletedTests >= plugin.TotalTests && isempty(currentTest)
                fprintf('\n');
            end
        end
    end
end
