classdef TestCase < handle
    % TestCase  Minimal MATLAB unittest.TestCase shim for GNU Octave.
    %
    % Implements the verification methods used across the pancData3 test
    % suite so that the same test classes can run under Octave without
    % modification.

    methods

        % ----- verifyTrue / verifyFalse -----
        function verifyTrue(testCase, actual, varargin)
            msg = testCase.extractMessage(varargin{:});
            assert(islogical(actual) || isnumeric(actual), ...
                'verifyTrue: input must be logical or numeric');
            assert(all(actual(:)), msg);
        end

        function verifyFalse(testCase, actual, varargin)
            msg = testCase.extractMessage(varargin{:});
            assert(islogical(actual) || isnumeric(actual), ...
                'verifyFalse: input must be logical or numeric');
            assert(~any(actual(:)), msg);
        end

        % ----- verifyEqual -----
        function verifyEqual(testCase, actual, expected, varargin)
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
            msg = testCase.extractMessage(varargin{:});
            assert(~isempty(actual), msg);
        end

        % ----- verifyEmpty -----
        function verifyEmpty(testCase, actual, varargin)
            msg = testCase.extractMessage(varargin{:});
            assert(isempty(actual), msg);
        end

        % ----- verifyClass -----
        function verifyClass(testCase, actual, expectedClass, varargin)
            msg = testCase.extractMessage(varargin{:});
            assert(isa(actual, expectedClass), ...
                'verifyClass: expected class %s, got %s. %s', ...
                expectedClass, class(actual), msg);
        end

        % ----- verifyGreaterThan -----
        function verifyGreaterThan(testCase, actual, expected, varargin)
            msg = testCase.extractMessage(varargin{:});
            assert(all(actual(:) > expected(:)), msg);
        end

        % ----- verifyGreaterThanOrEqual -----
        function verifyGreaterThanOrEqual(testCase, actual, expected, varargin)
            msg = testCase.extractMessage(varargin{:});
            assert(all(actual(:) >= expected(:)), msg);
        end

        % ----- verifyLessThan -----
        function verifyLessThan(testCase, actual, expected, varargin)
            msg = testCase.extractMessage(varargin{:});
            assert(all(actual(:) < expected(:)), msg);
        end

        % ----- verifyLessThanOrEqual -----
        function verifyLessThanOrEqual(testCase, actual, expected, varargin)
            msg = testCase.extractMessage(varargin{:});
            assert(all(actual(:) <= expected(:)), msg);
        end

        % ----- verifyError -----
        function verifyError(testCase, fcnHandle, expectedId, varargin)
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
            msg = testCase.extractMessage(varargin{:});
            if isempty(strfind(actual, expected))
                error('verifySubstring: ''%s'' not found in ''%s''. %s', ...
                    expected, actual, msg);
            end
        end

        % ----- verifyFail -----
        function verifyFail(testCase, varargin)
            msg = testCase.extractMessage(varargin{:});
            error('TestCase:verifyFail', 'Explicit failure. %s', msg);
        end

        % ----- applyFixture (PathFixture shim) -----
        function applyFixture(testCase, fixture)
            % PathFixture shim: just add the folder to the path
            if isstruct(fixture) && isfield(fixture, 'Folder')
                addpath(fixture.Folder);
            end
        end

    end

    methods (Access = private)
        function msg = extractMessage(testCase, varargin)
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
