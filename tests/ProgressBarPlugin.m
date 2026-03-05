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

            % Truncate test name to fit terminal width.
            % Line format: "  XX PASS  [NNN/NNN] <name> (X.XXs)\n"
            % Overhead is ~30 chars + counter digits; cap name to the remainder.
            termWidth = ProgressBarPlugin.getTerminalWidth();
            counterStr = sprintf('%3d/%d', plugin.CompletedTests, plugin.TotalTests);
            timeStr = sprintf('%.2fs', testTime);
            % 2 (indent) + 3 (emoji) + 6 (" PASS ") + 2 (" [") + counter + 2 ("] ") + 2 (" (") + time + 1 (")")
            overhead = 2 + 3 + 6 + 2 + length(counterStr) + 2 + 2 + length(timeStr) + 1;
            maxName = termWidth - overhead;
            displayName = ProgressBarPlugin.truncateName(testName, maxName);

            if passed
                plugin.PassedCount = plugin.PassedCount + 1;
                fprintf('  ✅ PASS  [%s] %s (%s)\n', ...
                    counterStr, displayName, timeStr);
            else
                plugin.FailedCount = plugin.FailedCount + 1;
                plugin.FailedNames{end+1} = testName;
                fprintf('  ❌ FAIL  [%s] %s (%s)\n', ...
                    counterStr, displayName, timeStr);
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

    methods (Static, Access = private)
        function w = getTerminalWidth()
        %GETTERMINALWIDTH Return the terminal column width (default 80).
            w = 80;
            try
                if isunix
                    [status, result] = system('tput cols 2>/dev/null');
                    if status == 0
                        val = str2double(strtrim(result));
                        if ~isnan(val) && val > 0
                            w = val;
                        end
                    end
                elseif ispc
                    [status, result] = system('mode con 2>nul');
                    if status == 0
                        tok = regexp(result, 'Columns:\s*(\d+)', 'tokens');
                        if ~isempty(tok)
                            val = str2double(tok{1}{1});
                            if ~isnan(val) && val > 0
                                w = val;
                            end
                        end
                    end
                end
            catch
                % Keep default
            end
        end

        function name = truncateName(name, maxLen)
        %TRUNCATENAME Shorten a test name to maxLen, preserving class/method.
            if maxLen < 10
                maxLen = 10;
            end
            if length(name) <= maxLen
                return;
            end
            % Try to keep the class/method boundary visible by splitting on '/'
            sep = strfind(name, '/');
            if ~isempty(sep)
                className = name(1:sep(end)-1);
                methodName = name(sep(end):end);  % includes the '/'
                % Prefer truncating the class path (package-qualified names)
                % while keeping the full method name visible.
                avail = maxLen - length(methodName) - 3;  % 3 for '...'
                if avail >= 8
                    name = ['...' className(end-avail+1:end) methodName];
                else
                    % Both parts are long; just truncate from the left
                    name = ['...' name(end-maxLen+4:end)];
                end
            else
                name = ['...' name(end-maxLen+4:end)];
            end
        end
    end
end
