classdef ProgressBarPlugin < matlab.unittest.plugins.TestRunnerPlugin
%PROGRESSBARPLUGIN Displays a clear one-line-per-test progress summary.
%   Prints each test result on its own line with pass/fail status, the
%   fully-qualified test name, and elapsed time.  At the end of the suite,
%   prints a summary block with counts and a list of any failures.
%
%   Usage:
%       runner.addPlugin(ProgressBarPlugin(numel(suite)));

    properties (Access = private)
        TotalTests double = 0       % Expected number of tests in the suite
        CompletedTests double = 0   % Running count of completed tests
        PassedCount double = 0      % Running count of passed tests
        FailedCount double = 0      % Running count of failed tests
        FailedNames cell = {}       % Cell array of fully-qualified names of failed tests
        TestTimer                   % tic handle for measuring individual test duration
        SuiteTimer                  % tic handle for measuring total suite duration
        DiaryFile char = ''         % Diary file to restart after each test (core modules hijack diary)
        CounterOffset double = 0    % Offset for display counter (e.g., after parallel phase)
    end

    methods
        function plugin = ProgressBarPlugin(totalTests, diaryFile, counterOffset)
        %PROGRESSBARPLUGIN Constructor. Accepts the total test count for
        %   progress fraction display (e.g., [3/42]).
        %   Optional diaryFile: path to restart diary after each test
        %   (core modules open their own diary, overriding the test log).
        %   Optional counterOffset: starting offset for display numbering
        %   (e.g., after parallel phase completes, so numbering continues).
            plugin.TotalTests = totalTests;
            if nargin >= 2 && ~isempty(diaryFile)
                plugin.DiaryFile = diaryFile;
            end
            if nargin >= 3 && ~isempty(counterOffset)
                plugin.CounterOffset = counterOffset;
            end
        end
    end

    methods (Access = protected)
        function runTestSuite(plugin, pluginData)
        %RUNTESTSUITE Called once when the suite starts. Resets counters,
        %   delegates to the superclass to execute all tests, then prints
        %   the final summary block.
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
        %RUNTEST Called for each individual test. Times the test, prints a
        %   single-line pass/fail result, and updates running counters.
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
            counterStr = sprintf('%3d/%d', plugin.CompletedTests + plugin.CounterOffset, plugin.TotalTests);
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

            % Core modules (visualize_results, sanity_checks, etc.) open
            % their own diary, overriding the test runner's diary. When
            % they call 'diary off', the test log is permanently lost.
            % Restart the diary after each test to recapture output.
            if ~isempty(plugin.DiaryFile)
                try
                    diary(plugin.DiaryFile);
                catch
                    % Diary restart failed — continue without logging
                end
            end
        end
    end

    methods (Access = private)
        function printSummary(plugin, elapsed)
        %PRINTSUMMARY Prints the final bordered summary block with total
        %   counts, elapsed time, and a list of any failed test names.
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
