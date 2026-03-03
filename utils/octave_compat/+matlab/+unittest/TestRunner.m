classdef TestRunner < handle
    % TestRunner  Minimal shim for matlab.unittest.TestRunner under Octave.

    properties
        Plugins = {};
    end

    methods (Static)
        function runner = withTextOutput()
            runner = matlab.unittest.TestRunner();
        end
    end

    methods
        function addPlugin(runner, plugin)
            runner.Plugins{end+1} = plugin;
        end

        function results = run(runner, suite)
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
