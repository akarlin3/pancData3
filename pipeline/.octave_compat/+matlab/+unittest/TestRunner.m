% TESTRUNNER  Octave-compatible shim for MATLAB's matlab.unittest.TestRunner.
%
%   MATLAB's TestRunner orchestrates test execution, plugin management,
%   and result collection. This shim provides a minimal sequential runner
%   that iterates over the struct array produced by TestSuite.fromFolder(),
%   instantiates each test class, calls setup/test/teardown, and collects
%   pass/fail results.
%
%   Behavioral differences from MATLAB's TestRunner:
%   - Plugins are accepted (addPlugin) but ignored -- no coverage or
%     diagnostic output beyond simple text.
%   - withTextOutput() returns a plain TestRunner (text output is the only mode).
%   - Setup/teardown methods are discovered by naming convention ('setup' /
%     'teardown', case-insensitive) rather than by methods block attributes,
%     because Octave cannot introspect method attributes at runtime.
%   - Results are returned as a struct array, not matlab.unittest.TestResult objects.
classdef TestRunner < handle
    % TestRunner  Minimal shim for matlab.unittest.TestRunner under Octave.

    properties
        Plugins = {};  % Accepted but not used; present for API compatibility.
    end

    methods (Static)
        function runner = withTextOutput()
            % WITHTEXTOUTPUT  Create a runner with console text output.
            %   In MATLAB, this configures verbose text diagnostics. Here it
            %   simply returns a default TestRunner since text output is the
            %   only mode available.
            runner = matlab.unittest.TestRunner();
        end
    end

    methods
        function addPlugin(runner, plugin)
            % ADDPLUGIN  Register a plugin (no-op in this shim).
            %   Plugins are stored but never invoked. This allows test runner
            %   code that adds CodeCoveragePlugin to work without errors.
            runner.Plugins{end+1} = plugin;
        end

        function results = run(runner, suite)
            % RUN  Execute all tests in the suite sequentially.
            %   Iterates over each entry in the suite struct array, creates
            %   a test class instance, runs setup -> test -> teardown, and
            %   records the outcome. Teardown runs even if the test fails.
            % Run test suite entries and collect results
            nTests = numel(suite);
            results = struct('Name', {}, 'Passed', {}, 'Failed', {}, ...
                'Incomplete', {}, 'Duration', {}, 'Details', {});

            for i = 1:nTests
                entry = suite(i);
                testName = entry.Name;
                fprintf('  Running %s ... ', testName);
                t0 = tic;

                passed = false;
                failMsg = '';
                try
                    % Ensure the test file's folder is on path
                    if ~isempty(entry.Folder)
                        addpath(entry.Folder);
                    end

                    % Create instance of the test class
                    obj = feval(entry.TestClass);

                    % Run TestMethodSetup if it exists
                    runner.runSetup(obj);

                    % Run the test method
                    feval(entry.TestMethod, obj);

                    % Run TestMethodTeardown if it exists
                    runner.runTeardown(obj);

                    passed = true;
                    fprintf('PASSED (%.2fs)\n', toc(t0));
                catch e
                    % Run teardown even on failure
                    try
                        runner.runTeardown(obj);
                    catch
                    end
                    failMsg = sprintf('%s: %s', e.identifier, e.message);
                    fprintf('FAILED (%.2fs)\n', toc(t0));
                    fprintf('    %s\n', failMsg);
                end

                r.Name = testName;
                r.Passed = passed;
                r.Failed = ~passed;
                r.Incomplete = false;
                r.Duration = toc(t0);
                r.Details = failMsg;
                if isempty(fieldnames(results)) && i == 1
                    results = r;
                else
                    results(i) = r;
                end
            end
        end
    end

    methods (Access = private)
        function runSetup(runner, obj)
            % RUNSETUP  Find and call the test class's setup method if it exists.
            %   Uses a naming convention (case-insensitive 'setup') because
            %   Octave's classdef does not expose method block attributes
            %   (TestMethodSetup) for runtime introspection.
            % Call the setup method if the class defines one.
            % In Octave classdef, we look for a method named after
            % the TestMethodSetup convention. Since we cannot introspect
            % methods blocks by attribute, we use a naming convention:
            % look for a method called 'setup' (case-insensitive).
            meths = methods(obj);
            for k = 1:numel(meths)
                if strcmpi(meths{k}, 'setup')
                    feval(meths{k}, obj);
                    return;
                end
            end
        end

        function runTeardown(runner, obj)
            % RUNTEARDOWN  Find and call the test class's teardown method.
            %   Mirror of runSetup; looks for a method named 'teardown'.
            meths = methods(obj);
            for k = 1:numel(meths)
                if strcmpi(meths{k}, 'teardown')
                    feval(meths{k}, obj);
                    return;
                end
            end
        end
    end
end
