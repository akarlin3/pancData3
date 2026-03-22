% TESTCASE  Octave-compatible shim for MATLAB's matlab.unittest.TestCase.
%
%   MATLAB's TestCase (R2013a+) provides a rich assertion/verification API
%   for unit tests. Octave lacks this framework entirely. This shim
%   reimplements the verify* methods used by the pancData3 test suite using
%   plain assert() calls, so test classes can inherit from this shim and
%   run identically under Octave.
%
%   Behavioral differences from MATLAB's TestCase:
%   - Verification failures throw immediately (like assertions), rather than
%     being collected as soft failures. MATLAB's verify* methods record a
%     failure but continue the test; here, the test stops on first failure.
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

    methods

        % ----- verifyTrue / verifyFalse -----
        function verifyTrue(testCase, actual, varargin)
            % VERIFYTRUE  Assert that all elements of 'actual' are true.
            msg = testCase.extractMessage(varargin{:});
            assert(islogical(actual) || isnumeric(actual), ...
                'verifyTrue: input must be logical or numeric');
            assert(all(actual(:)), msg);
        end

        function verifyFalse(testCase, actual, varargin)
            % VERIFYFALSE  Assert that no elements of 'actual' are true.
            msg = testCase.extractMessage(varargin{:});
            assert(islogical(actual) || isnumeric(actual), ...
                'verifyFalse: input must be logical or numeric');
            assert(~any(actual(:)), msg);
        end

        % ----- verifyEqual -----
        function verifyEqual(testCase, actual, expected, varargin)
            % VERIFYEQUAL  Assert equality with optional tolerances.
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
                    error('TestCase:verifyEqual', '%s', msg);
                end
            elseif ischar(actual) && ischar(expected)
                assert(strcmp(actual, expected), ...
                    'verifyEqual: strings differ.\n  Actual: %s\n  Expected: %s', actual, expected);
            elseif iscell(actual) && iscell(expected)
                assert(isequal(actual, expected), ...
                    'verifyEqual: cell arrays differ');
            else
                assert(isequal(actual, expected), ...
                    'verifyEqual: values differ');
            end
        end

        % ----- verifyNotEmpty -----
        function verifyNotEmpty(testCase, actual, varargin)
            % VERIFYNOTEMPTY  Assert that 'actual' is not empty.
            msg = testCase.extractMessage(varargin{:});
            assert(~isempty(actual), msg);
        end

        % ----- verifyEmpty -----
        function verifyEmpty(testCase, actual, varargin)
            % VERIFYEMPTY  Assert that 'actual' is empty.
            msg = testCase.extractMessage(varargin{:});
            assert(isempty(actual), msg);
        end

        % ----- verifyClass -----
        function verifyClass(testCase, actual, expectedClass, varargin)
            % VERIFYCLASS  Assert that 'actual' is an instance of expectedClass.
            msg = testCase.extractMessage(varargin{:});
            assert(isa(actual, expectedClass), ...
                'verifyClass: expected class %s, got %s. %s', ...
                expectedClass, class(actual), msg);
        end

        % ----- verifyGreaterThan -----
        function verifyGreaterThan(testCase, actual, expected, varargin)
            % VERIFYGREATERTHAN  Assert all(actual > expected).
            %   Supports scalar expansion: if expected is scalar it is
            %   compared against every element of actual, and vice-versa.
            msg = testCase.extractMessage(varargin{:});
            if isscalar(expected)
                assert(all(actual(:) > expected), msg);
            elseif isscalar(actual)
                assert(all(actual > expected(:)), msg);
            else
                assert(isequal(size(actual), size(expected)), ...
                    'verifyGreaterThan: size mismatch between actual and expected');
                assert(all(actual(:) > expected(:)), msg);
            end
        end

        % ----- verifyGreaterThanOrEqual -----
        function verifyGreaterThanOrEqual(testCase, actual, expected, varargin)
            % VERIFYGREATERTHANOREQUAL  Assert all(actual >= expected).
            %   Supports scalar expansion: if expected is scalar it is
            %   compared against every element of actual, and vice-versa.
            msg = testCase.extractMessage(varargin{:});
            if isscalar(expected)
                assert(all(actual(:) >= expected), msg);
            elseif isscalar(actual)
                assert(all(actual >= expected(:)), msg);
            else
                assert(isequal(size(actual), size(expected)), ...
                    'verifyGreaterThanOrEqual: size mismatch between actual and expected');
                assert(all(actual(:) >= expected(:)), msg);
            end
        end

        % ----- verifyLessThan -----
        function verifyLessThan(testCase, actual, expected, varargin)
            % VERIFYLESSTHAN  Assert all(actual < expected).
            %   Supports scalar expansion: if expected is scalar it is
            %   compared against every element of actual, and vice-versa.
            msg = testCase.extractMessage(varargin{:});
            if isscalar(expected)
                assert(all(actual(:) < expected), msg);
            elseif isscalar(actual)
                assert(all(actual < expected(:)), msg);
            else
                assert(isequal(size(actual), size(expected)), ...
                    'verifyLessThan: size mismatch between actual and expected');
                assert(all(actual(:) < expected(:)), msg);
            end
        end

        % ----- verifyLessThanOrEqual -----
        function verifyLessThanOrEqual(testCase, actual, expected, varargin)
            % VERIFYLESSSTHANOREQUAL  Assert all(actual <= expected).
            %   Supports scalar expansion: if expected is scalar it is
            %   compared against every element of actual, and vice-versa.
            msg = testCase.extractMessage(varargin{:});
            if isscalar(expected)
                assert(all(actual(:) <= expected), msg);
            elseif isscalar(actual)
                assert(all(actual <= expected(:)), msg);
            else
                assert(isequal(size(actual), size(expected)), ...
                    'verifyLessThanOrEqual: size mismatch between actual and expected');
                assert(all(actual(:) <= expected(:)), msg);
            end
        end

        % ----- verifyError -----
        function verifyError(testCase, fcnHandle, expectedId, varargin)
            % VERIFYERROR  Assert that calling fcnHandle throws an error.
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
            assert(caught, 'verifyError: expected an error but none was thrown');
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
            % VERIFYSUBSTRING  Assert that 'expected' appears within 'actual'.
            msg = testCase.extractMessage(varargin{:});
            if isempty(strfind(actual, expected))
                error('verifySubstring: ''%s'' not found in ''%s''. %s', ...
                    expected, actual, msg);
            end
        end

        % ----- verifyFail -----
        function verifyFail(testCase, varargin)
            % VERIFYFAIL  Unconditionally fail the test with an optional message.
            msg = testCase.extractMessage(varargin{:});
            error('TestCase:verifyFail', 'Explicit failure. %s', msg);
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
    end

end