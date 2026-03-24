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
%   - Setup/teardown methods are discovered by parsing source files for
%     methods(TestMethodSetup) and methods(TestMethodTeardown) blocks.
%   - Results are returned as a struct array, not matlab.unittest.TestResult objects.
classdef TestRunner < handle
    % TestRunner  Minimal shim for matlab.unittest.TestRunner under Octave.

    properties
        Plugins = {};  % Accepted but not used; present for API compatibility.
    end

    methods (Static)
        function runner = withTextOutput()
            % WITHTEXTOUTPUT  Create a runner with console text output.
            runner = matlab.unittest.TestRunner();
        end
    end

    methods
        function addPlugin(runner, plugin)
            % ADDPLUGIN  Register a plugin (no-op in this shim).
            runner.Plugins{end+1} = plugin;
        end

        function results = run(runner, suite)
            % RUN  Execute all tests in the suite sequentially.
            nTests = numel(suite);
            results = struct('Name', {}, 'Passed', {}, 'Failed', {}, ...
                'Incomplete', {}, 'Duration', {}, 'Details', {});

            % Cache parsed setup/teardown methods per file
            setup_cache = struct();
            teardown_cache = struct();

            for i = 1:nTests
                entry = suite(i);
                testName = entry.Name;
                fprintf('  Running %s ... ', testName);
                t0 = tic;

                passed = false;
                failMsg = '';
                obj = [];
                try
                    % Ensure the test file's folder is on path
                    if ~isempty(entry.Folder)
                        addpath(entry.Folder);
                    end

                    % Create instance of the test class
                    obj = feval(entry.TestClass);

                    % Reset verification failures before each test
                    if ismethod(obj, 'resetVerifications')
                        obj.resetVerifications();
                    end

                    % Run TestMethodSetup methods (parsed from source)
                    runner.runSetup(obj, entry);

                    % Run the test method
                    feval(entry.TestMethod, obj);

                    % Run TestMethodTeardown methods
                    runner.runTeardown(obj, entry);

                    % Check for accumulated verification failures
                    if ismethod(obj, 'hasVerificationFailures') && obj.hasVerificationFailures()
                        obj.assertNoFailures();
                    end

                    passed = true;
                    fprintf('PASSED (%.2fs)\n', toc(t0));
                catch e
                    % Run teardown even on failure
                    try
                        if ~isempty(obj)
                            runner.runTeardown(obj, entry);
                        end
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
        function runSetup(runner, obj, entry)
            % RUNSETUP  Find and call setup methods.
            %   First tries to parse methods(TestMethodSetup) blocks from the
            %   source file. Falls back to looking for a method named 'setup'.
            setup_meths = runner.parseBlockMethods(entry, 'TestMethodSetup');
            if ~isempty(setup_meths)
                for k = 1:numel(setup_meths)
                    if ismethod(obj, setup_meths{k})
                        feval(setup_meths{k}, obj);
                    end
                end
                return;
            end
            % Fallback: look for method named 'setup'
            meths = methods(obj);
            for k = 1:numel(meths)
                if strcmpi(meths{k}, 'setup')
                    feval(meths{k}, obj);
                    return;
                end
            end
        end

        function runTeardown(runner, obj, entry)
            % RUNTEARDOWN  Find and call teardown methods.
            teardown_meths = runner.parseBlockMethods(entry, 'TestMethodTeardown');
            if ~isempty(teardown_meths)
                for k = 1:numel(teardown_meths)
                    if ismethod(obj, teardown_meths{k})
                        feval(teardown_meths{k}, obj);
                    end
                end
                return;
            end
            % Fallback: look for method named 'teardown'
            meths = methods(obj);
            for k = 1:numel(meths)
                if strcmpi(meths{k}, 'teardown')
                    feval(meths{k}, obj);
                    return;
                end
            end
        end

        function method_names = parseBlockMethods(runner, entry, blockType)
            % PARSEBLOCKMETHODS  Parse method names from a specific methods block.
            %   Reads the test source file and extracts function names from
            %   methods(TestMethodSetup) or methods(TestMethodTeardown) blocks.
            method_names = {};
            if ~isfield(entry, 'FilePath') || isempty(entry.FilePath)
                return;
            end
            fid = fopen(entry.FilePath, 'r');
            if fid == -1; return; end
            txt = fread(fid, '*char')';
            fclose(fid);

            in_block = false;
            brace_depth = 0;
            lines = strsplit(txt, '\n');
            % Build regex pattern for the block type
            pat = ['^\s*methods\s*\(\s*' blockType '\s*\)'];
            for li = 1:numel(lines)
                ln = strtrim(lines{li});
                if ~isempty(regexp(ln, pat, 'once'))
                    in_block = true;
                    brace_depth = 0;
                    continue;
                end
                if in_block
                    if ~isempty(regexp(ln, '^\s*end\s*$', 'once'))
                        if brace_depth <= 0
                            in_block = false;
                            continue;
                        else
                            brace_depth = brace_depth - 1;
                            continue;
                        end
                    end
                    tokens = regexp(ln, '^\s*function\s+(?:\w+\s*=\s*)?(\w+)\s*\(', 'tokens');
                    if ~isempty(tokens)
                        method_names{end+1} = tokens{1}{1};
                        brace_depth = brace_depth + 1;
                    end
                end
            end
        end
    end
end
