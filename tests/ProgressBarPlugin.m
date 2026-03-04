classdef ProgressBarPlugin < matlab.unittest.plugins.TestRunnerPlugin
%PROGRESSBARPLUGIN Displays a clear one-line-per-test progress summary.
%   Prints each test result on its own line with pass/fail status, the
%   fully-qualified test name, and elapsed time.  At the end of the suite,
%   prints a summary block with counts and a list of any failures.
%
%   Usage:
%       runner.addPlugin(ProgressBarPlugin(numel(suite)));

    properties (Access = private)
        TotalTests double = 0
        CompletedTests double = 0
        PassedCount double = 0
        FailedCount double = 0
        FailedNames cell = {}
        TestTimer
        SuiteTimer
    end

    methods
        function plugin = ProgressBarPlugin(totalTests)
            plugin.TotalTests = totalTests;
        end
    end

    methods (Access = protected)
        function runTestSuite(plugin, pluginData)
            plugin.CompletedTests = 0;
            plugin.PassedCount = 0;
            plugin.FailedCount = 0;
            plugin.FailedNames = {};
            plugin.SuiteTimer = tic;

            fprintf('\n');
            runTestSuite@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);

            elapsed = toc(plugin.SuiteTimer);
            plugin.printSummary(elapsed);
        end

        function runTest(plugin, pluginData)
            plugin.TestTimer = tic;
            runTest@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
            testTime = toc(plugin.TestTimer);

            plugin.CompletedTests = plugin.CompletedTests + 1;

            testName = pluginData.Name;
            passed = ~pluginData.TestResult.Failed;

            if passed
                plugin.PassedCount = plugin.PassedCount + 1;
                fprintf('  ✅ PASS  [%3d/%d] %s (%.2fs)\n', ...
                    plugin.CompletedTests, plugin.TotalTests, testName, testTime);
            else
                plugin.FailedCount = plugin.FailedCount + 1;
                plugin.FailedNames{end+1} = testName;
                fprintf('  ❌ FAIL  [%3d/%d] %s (%.2fs)\n', ...
                    plugin.CompletedTests, plugin.TotalTests, testName, testTime);
            end
        end
    end

    methods (Access = private)
        function printSummary(plugin, elapsed)
            fprintf('\n');
            fprintf('───────────────────────────────────────────────────\n');
            fprintf('  Test Suite Summary\n');
            fprintf('───────────────────────────────────────────────────\n');
            fprintf('  Total:   %d tests in %.1fs\n', plugin.TotalTests, elapsed);
            fprintf('  Passed:  %d\n', plugin.PassedCount);
            if plugin.FailedCount > 0
                fprintf('  Failed:  %d\n', plugin.FailedCount);
                fprintf('\n  Failed tests:\n');
                for k = 1:numel(plugin.FailedNames)
                    fprintf('    ❌ %s\n', plugin.FailedNames{k});
                end
            else
                fprintf('  Failed:  0\n');
            end
            fprintf('───────────────────────────────────────────────────\n\n');
        end
    end
end
