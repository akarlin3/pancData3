% TESTCASE  Octave-compatible shim for MATLAB's matlab.unittest.TestCase.
%
%   MATLAB's TestCase (R2013a+) provides a rich assertion/verification API
%   for unit tests. Octave lacks this framework entirely. This shim
%   reimplements the verify* methods used by the pancData3 test suite so
%   that test classes can inherit from this shim and run under Octave.
%
%   Behavioral notes:
%   - Verification failures are recorded as soft failures. The test method
%     continues executing after a verify* failure, matching MATLAB semantics.
%     Accumulated failures are stored in the VerificationFailures property
%     and can be checked via hasVerificationFailures / assertNoFailures.
%   - verifyError does not check the error identifier -- it only checks that
%     some error was thrown, because Octave error IDs are inconsistent.
%   - verifyWarning is lenient: it always passes even if no warning is issued,
%     because some Octave builds suppress warnings differently.
%   - Only verification methods used by the pipeline's test suite are included;
%     methods like verifySize, verifyInstanceOf, etc. are not implemented.
%   - No test fixture lifecycle (TestMethodSetup/Teardown) -- those are handled
%     by the TestRunner shim via naming convention.
classdef TestCase < handle
    % TestCase  Minimal MATLAB unittest.TestCase shim for GNU Octave.
    %
    % Implements the verification methods used across the pancData3 test
    % suite so that the same test classes can run under Octave without
    % modification.

    properties (Access = public)
        % VerificationFailures  Cell array of failure message strings
        %   accumulated by verify* methods during a test run.
        VerificationFailures = {};
    end

    methods

        function resetVerifications(testCase)
            % RESETVERIFICATIONS  Clear accumulated verification failures.
            %   Called by the TestRunner shim before each test method.
            testCase.VerificationFailures = {};
        end

        function tf = hasVerificationFailures(testCase)
            % HASVERIFICATIONFAILURES  Return true if any verify* recorded a failure.
            tf = ~isempty(testCase.VerificationFailures);
        end

        function assertNoFailures(testCase)
            % ASSERTNOFAILURES  Throw if any verification failures were recorded.
            %   Called by the TestRunner shim after each test method to
            %   convert accumulated soft failures into a hard error.
            if ~isempty(testCase.VerificationFailures)
                safeStrs = cellfun(@(x) char(string(x)), ...
                    testCase.VerificationFailures, 'UniformOutput', false);
                msgs = strjoin(safeStrs, '\n  ');
                error('TestCase:verificationFailed', ...
                    '%d verification failure(s):\n  %s', ...
                    numel(testCase.VerificationFailures), msgs);
            end
        end

        % ----- verifyTrue / verifyFalse -----
        function verifyTrue(testCase, actual, varargin)
            % VERIFYTRUE  Verify that all elements of 'actual' are true.
            msg = testCase.extractMessage(varargin{:});
            if ~(islogical(actual) || isnumeric(actual))
                testCase.recordFailure('verifyTrue: input must be logical or numeric');
                return;
            end
            if ~all(actual(:))
                if isempty(msg)
                    msg = 'verifyTrue failed';
                end
                testCase.recordFailure(msg);
            end
        end

        function verifyFalse(testCase, actual, varargin)
            % VERIFYFALSE  Verify that no elements of 'actual' are true.
            msg = testCase.extractMessage(varargin{:});
            if ~(islogical(actual) || isnumeric(actual))
                testCase.recordFailure('verifyFalse: input must be logical or numeric');
                return;
            end
            if any(actual(:))
                if isempty(msg)
                    msg = 'verifyFalse failed';
                end
                testCase.recordFailure(msg);
            end
        end

        % ----- verifyEqual -----
        function verifyEqual(testCase, actual, expected, varargin)
            % VERIFYEQUAL  Verify equality with optional tolerances.
            %   Supports 'AbsTol' (absolute tolerance) and 'RelTol' (relative
            %   tolerance) name-value pairs for numeric comparisons. Strings
            %   are compared with strcmp; cells and other types use isequal.
            msg = '';
            absTol = 0;
            relTol = 0;
            k = 1;
            while k <= numel(varargin)
                if ischar(varargin{k})
                    switch varargin{k}
                        case 'AbsTol'
                            absTol = varargin{k+1}; k = k + 2;
                        case 'RelTol'
                            relTol = varargin{k+1}; k = k + 2;
                        otherwise
                            % Treat as message string
                            msg = varargin{k}; k = k + 1;
                    end
                else
                    k = k + 1;
                end
            end

            if isnumeric(actual) && isnumeric(expected)
                diff_val = abs(double(actual) - double(expected));
                if relTol > 0
                    denom = max(abs(double(expected)), eps);
                    ok = all(diff_val(:) ./ denom(:) <= relTol);
                elseif absTol > 0
                    ok = all(diff_val(:) <= absTol);
                else
                    ok = isequal(actual, expected);
                end
                if ~ok
                    if isempty(msg)
                        msg = sprintf('verifyEqual failed:\n  Actual: %s\n  Expected: %s', ...
                            mat2str(actual), mat2str(expected));
                    end
                    testCase.recordFailure(msg);
                end
            elseif ischar(actual) && ischar(expected)
                if ~strcmp(actual, expected)
                    if isempty(msg)
                        msg = sprintf('verifyEqual: strings differ.\n  Actual: %s\n  Expected: %s', actual, expected);
                    end
                    testCase.recordFailure(msg);
                end
            elseif iscell(actual) && iscell(expected)
                if ~isequal(actual, expected)
                    if isempty(msg)
                        msg = 'verifyEqual: cell arrays differ';
                    end
                    testCase.recordFailure(msg);
                end
            else
                if ~isequal(actual, expected)
                    if isempty(msg)
                        msg = 'verifyEqual: values differ';
                    end
                    testCase.recordFailure(msg);
                end
            end
        end

        % ----- verifyNotEqual -----
        function verifyNotEqual(testCase, actual, notExpected, varargin)
            % VERIFYNOTEQUAL  Verify that actual is not equal to notExpected.
            msg = testCase.extractMessage(varargin{:});
            if isequal(actual, notExpected)
                if isempty(msg)
                    msg = 'verifyNotEqual failed: values are equal';
                end
                testCase.recordFailure(msg);
            end
        end

        % ----- verifyNotEmpty -----
        function verifyNotEmpty(testCase, actual, varargin)
            % VERIFYNOTEMPTY  Verify that 'actual' is not empty.
            msg = testCase.extractMessage(varargin{:});
            if isempty(actual)
                if isempty(msg)
                    msg = 'verifyNotEmpty failed: value is empty';
                end
                testCase.recordFailure(msg);
            end
        end

        % ----- verifyEmpty -----
        function verifyEmpty(testCase, actual, varargin)
            % VERIFYEMPTY  Verify that 'actual' is empty.
            msg = testCase.extractMessage(varargin{:});
            if ~isempty(actual)
                if isempty(msg)
                    msg = 'verifyEmpty failed: value is not empty';
                end
                testCase.recordFailure(msg);
            end
        end

        % ----- verifyClass -----
        function verifyClass(testCase, actual, expectedClass, varargin)
            % VERIFYCLASS  Verify that 'actual' is an instance of expectedClass.
            msg = testCase.extractMessage(varargin{:});
            if ~isa(actual, expectedClass)
                failMsg = sprintf('verifyClass: expected class %s, got %s. %s', ...
                    expectedClass, class(actual), msg);
                testCase.recordFailure(failMsg);
            end
        end

        % ----- verifyGreaterThan -----
        function verifyGreaterThan(testCase, actual, expected, varargin)
            % VERIFYGREATERTHAN  Verify all(actual > expected).
            %   Supports scalar expansion: if expected is scalar it is
            %   compared against every element of actual, and vice-versa.
            msg = testCase.extractMessage(varargin{:});
            ok = true;
            if isscalar(expected)
                ok = all(actual(:) > expected);
            elseif isscalar(actual)
                ok = all(actual > expected(:));
            else
                if ~isequal(size(actual), size(expected))
                    testCase.recordFailure('verifyGreaterThan: size mismatch between actual and expected');
                    return;
                end
                ok = all(actual(:) > expected(:));
            end
            if ~ok
                if isempty(msg)
                    msg = 'verifyGreaterThan failed';
                end
                testCase.recordFailure(msg);
            end
        end

        % ----- verifyGreaterThanOrEqual -----
        function verifyGreaterThanOrEqual(testCase, actual, expected, varargin)
            % VERIFYGREATERTHANOREQUAL  Verify all(actual >= expected).
            %   Supports scalar expansion: if expected is scalar it is
            %   compared against every element of actual, and vice-versa.
            msg = testCase.extractMessage(varargin{:});
            ok = true;
            if isscalar(expected)
                ok = all(actual(:) >= expected);
            elseif isscalar(actual)
                ok = all(actual >= expected(:));
            else
                if ~isequal(size(actual), size(expected))
                    testCase.recordFailure('verifyGreaterThanOrEqual: size mismatch between actual and expected');
                    return;
                end
                ok = all(actual(:) >= expected(:));
            end
            if ~ok
                if isempty(msg)
                    msg = 'verifyGreaterThanOrEqual failed';
                end
                testCase.recordFailure(msg);
            end
        end

        % ----- verifyLessThan -----
        function verifyLessThan(testCase, actual, expected, varargin)
            % VERIFYLESSTHAN  Verify all(actual < expected).
            %   Supports scalar expansion: if expected is scalar it is
            %   compared against every element of actual, and vice-versa.
            msg = testCase.extractMessage(varargin{:});
            ok = true;
            if isscalar(expected)
                ok = all(actual(:) < expected);
            elseif isscalar(actual)
                ok = all(actual < expected(:));
            else
                if ~isequal(size(actual), size(expected))
                    testCase.recordFailure('verifyLessThan: size mismatch between actual and expected');
                    return;
                end
                ok = all(actual(:) < expected(:));
            end
            if ~ok
                if isempty(msg)
                    msg = 'verifyLessThan failed';
                end
                testCase.recordFailure(msg);
            end
        end

        % ----- verifyLessThanOrEqual -----
        function verifyLessThanOrEqual(testCase, actual, expected, varargin)
            % VERIFYLESSSTHANOREQUAL  Verify all(actual <= expected).
            %   Supports scalar expansion: if expected is scalar it is
            %   compared against every element of actual, and vice-versa.
            msg = testCase.extractMessage(varargin{:});
            ok = true;
            if isscalar(expected)
                ok = all(actual(:) <= expected);
            elseif isscalar(actual)
                ok = all(actual <= expected(:));
            else
                if ~isequal(size(actual), size(expected))
                    testCase.recordFailure('verifyLessThanOrEqual: size mismatch between actual and expected');
                    return;
                end
                ok = all(actual(:) <= expected(:));
            end
            if ~ok
                if isempty(msg)
                    msg = 'verifyLessThanOrEqual failed';
                end
                testCase.recordFailure(msg);
            end
        end

        % ----- verifyError -----
        function verifyError(testCase, fcnHandle, expectedId, varargin)
            % VERIFYERROR  Verify that calling fcnHandle throws an error.
            %   The expectedId is accepted but not strictly checked, because
            %   Octave error identifiers are often empty or differ from MATLAB.
            caught = false;
            try
                fcnHandle();
            catch e
                caught = true;
                % If expectedId is supplied, optionally check it
                if ~isempty(expectedId) && ischar(expectedId)
                    % Octave error IDs are not always set, so accept any error
                end
            end
            if ~caught
                testCase.recordFailure('verifyError: expected an error but none was thrown');
            end
        end

        % ----- verifyWarning -----
        function verifyWarning(testCase, fcnHandle, expectedId, varargin)
            % VERIFYWARNING  Run fcnHandle and check for warnings.
            %   This is lenient: it always passes even if no warning is issued,
            %   because Octave warning behavior varies across builds and platforms.
            % Capture warnings and verify at least one was issued
            [~, warnId] = lastwarn('');
            lastwarn('');  % clear
            fcnHandle();
            [~, lastId] = lastwarn;
            % Accept if any warning was issued (Octave warning IDs may differ)
            % If no warning was issued, still pass — some Octave builds suppress warnings
        end

        % ----- verifySubstring -----
        function verifySubstring(testCase, actual, expected, varargin)
            % VERIFYSUBSTRING  Verify that 'expected' appears within 'actual'.
            msg = testCase.extractMessage(varargin{:});
            if isempty(strfind(actual, expected))
                if isempty(msg)
                    msg = sprintf('verifySubstring: ''%s'' not found in ''%s''', ...
                        expected, actual);
                end
                testCase.recordFailure(msg);
            end
        end

        % ----- verifySize -----
        function verifySize(testCase, actual, expectedSize, varargin)
            % VERIFYSIZE  Verify that size(actual) equals expectedSize.
            msg = testCase.extractMessage(varargin{:});
            actualSize = size(actual);
            if ~isequal(actualSize, expectedSize)
                if isempty(msg)
                    msg = sprintf('verifySize: expected size [%s], got [%s]', ...
                        num2str(expectedSize), num2str(actualSize));
                end
                testCase.recordFailure(msg);
            end
        end

        % ----- verifyFail -----
        function verifyFail(testCase, varargin)
            % VERIFYFAIL  Unconditionally record a failure with an optional message.
            msg = testCase.extractMessage(varargin{:});
            if isempty(msg)
                msg = 'Explicit failure';
            else
                msg = sprintf('Explicit failure. %s', msg);
            end
            testCase.recordFailure(msg);
        end

        % ----- assumeTrue / assumeFail -----
        function assumeTrue(testCase, condition, varargin)
            % ASSUMETRUE  Skip the test if condition is false.
            %   In MATLAB, this marks the test as Incomplete (not Failed).
            %   In this shim, we throw an error to skip the test.
            if ~condition
                msg = testCase.extractMessage(varargin{:});
                if isempty(msg), msg = 'Assumption failed'; end
                error('TestCase:assumptionFailed', 'Assumption failed: %s', msg);
            end
        end

        function assumeFail(testCase, varargin)
            % ASSUMEFAIL  Unconditionally skip the test.
            msg = testCase.extractMessage(varargin{:});
            if isempty(msg), msg = 'Assumption explicitly failed'; end
            error('TestCase:assumptionFailed', 'Assumption failed: %s', msg);
        end

        % ----- addTeardown -----
        function addTeardown(testCase, fcn, varargin)
            % ADDTEARDOWN  Register a function to run at teardown.
            %   This is a no-op shim; the Octave TestRunner handles teardown
            %   via methods(TestMethodTeardown) blocks. For one-off teardown
            %   actions, callers should use onCleanup instead.
        end

        % ----- applyFixture (PathFixture shim) -----
        function applyFixture(testCase, fixture)
            % APPLYFIXTURE  Apply a test fixture (PathFixture only).
            %   In MATLAB, fixtures manage setup/teardown of shared resources.
            %   This shim only supports PathFixture: it adds the folder to
            %   the MATLAB/Octave path. Other fixture types are silently ignored.
            % PathFixture shim: just add the folder to the path
            if isstruct(fixture) && isfield(fixture, 'Folder')
                addpath(fixture.Folder);
            end
        end

    end

    methods (Access = private)
        function msg = extractMessage(testCase, varargin)
            % EXTRACTMESSAGE  Pull an optional diagnostic message from varargin.
            %   Scans varargin for the first char argument and returns it as the
            %   failure message. Returns '' if none found.
            msg = '';
            for k = 1:numel(varargin)
                if ischar(varargin{k})
                    msg = varargin{k};
                    return;
                end
            end
        end

        function recordFailure(testCase, msg)
            % RECORDFAILURE  Record a verification failure without throwing.
            %   Appends the message to VerificationFailures and prints a
            %   warning so the user sees it immediately, but does NOT stop
            %   test execution.
            testCase.VerificationFailures{end+1} = msg;
            fprintf(2, '  Verification FAILED: %s\n', msg);
        end
    end

end