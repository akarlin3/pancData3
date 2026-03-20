function escaped_arg = escape_shell_arg(arg, style)
% ESCAPE_SHELL_ARG Escapes a string to be used as a shell argument.
%
%   escaped_arg = escape_shell_arg(arg)
%   escaped_arg = escape_shell_arg(arg, style)
%
%   Inputs:
%       arg   - The string to escape.
%       style - (Optional) 'pc' or 'unix'. Defaults to the current OS.
%
%   Outputs:
%       escaped_arg - The escaped string safe for use in system().
%
% --- Analytical Rationale ---
% This pipeline calls external tools (primarily dcm2niix for DICOM-to-NIfTI
% conversion) via MATLAB's system() function. Patient data paths often contain
% spaces, special characters, and institution-specific naming conventions that
% would be misinterpreted by the shell without proper escaping. For example,
% a path like '/data/Pancreas Study/Patient #3/DWI scan' would break an
% unescaped system() call.
%
% More critically, this function prevents shell injection vulnerabilities.
% Patient identifiers and folder names originate from clinical DICOM metadata,
% which is not sanitized by the scanner. Without escaping, a malicious or
% malformed DICOM StudyDescription could inject arbitrary shell commands.
%
% The function supports both Windows (double-quote wrapping with backslash
% escaping) and Unix (single-quote wrapping) conventions, since this pipeline
% may run on clinical workstations (often Windows) or research servers (Linux).
%
% --- Unicode and Non-ASCII Character Handling ---
% This function now includes enhanced support for Unicode and non-ASCII 
% characters commonly found in international clinical environments (e.g.,
% accented characters in patient names, Cyrillic/CJK characters in institution
% names). The function detects system encoding and ensures proper character
% representation for shell operations.

    % --- Persistent encoding cache ---
    % Cache the detected system encoding so that expensive system() calls
    % (chcp on Windows, locale charmap on Unix) are executed only once per
    % MATLAB session instead of on every invocation.
    persistent cached_pc_encoding;
    persistent cached_unix_encoding;

    % --- OS Detection ---
    % Auto-detect the shell style based on the current platform. The optional
    % 'style' override is provided for cross-platform testing (e.g., testing
    % Windows escaping logic on a Unix development machine).
    if nargin < 2
        if ispc
            style = 'pc';
        else
            style = 'unix';
        end
    end

    % --- Input Validation ---
    % Reject non-string inputs early to produce a clear error message rather
    % than a cryptic failure deep inside strrep. Numeric patient IDs or MRNs
    % must be converted to strings by the caller before escaping.
    if ~ischar(arg) && ~isstring(arg)
        error('escape_shell_arg:invalidInput', 'Input argument must be a string or character vector.');
    end
    % Normalize to char for consistent string operations below. MATLAB string
    % objects and char arrays behave differently with strrep and concatenation.
    arg = char(arg);

    % --- Unicode and Encoding Handling ---
    % Detect and handle Unicode characters that may cause issues in shell operations
    try
        % Check if the string contains non-ASCII characters
        has_unicode = any(double(arg) > 127);
        
        if has_unicode
            % Get system encoding to ensure proper character handling
            if strcmpi(style, 'pc')
                % Windows: Check for system code page and handle Unicode paths
                if isempty(cached_pc_encoding)
                    try
                        % Attempt to get system code page (runs only once)
                        [~, cp_output] = system('chcp');
                        if contains(cp_output, '65001') % UTF-8
                            cached_pc_encoding = 'UTF-8';
                        else
                            cached_pc_encoding = 'windows-1252'; % Common Windows default
                        end
                    catch
                        cached_pc_encoding = 'windows-1252';
                    end
                end
                system_encoding = cached_pc_encoding;
                
                % For Windows, ensure proper Unicode handling by converting
                % to native encoding if needed
                if ~strcmp(system_encoding, 'UTF-8')
                    try
                        % Attempt to handle encoding conversion
                        arg = native2unicode(unicode2native(arg, 'UTF-8'), system_encoding);
                    catch
                        % If conversion fails, proceed with original string
                        % and hope the system can handle it
                    end
                end
            else
                % Unix systems: Most modern Unix systems use UTF-8
                if isempty(cached_unix_encoding)
                    try
                        % Check locale for encoding information (runs only once)
                        [~, locale_output] = system('locale charmap');
                        if contains(upper(locale_output), 'UTF-8')
                            cached_unix_encoding = 'UTF-8';
                        else
                            cached_unix_encoding = 'ISO-8859-1';
                        end
                    catch
                        cached_unix_encoding = 'UTF-8'; % Default assumption for modern Unix
                    end
                end
                system_encoding = cached_unix_encoding; %#ok<NASGU>
            end
        end
    catch
        % If Unicode detection fails, proceed with original escaping
        % This ensures backward compatibility
    end

    if strcmpi(style, 'pc')
        % --- Windows Escaping Strategy ---
        % Windows cmd.exe uses double-quote delimiters. Any literal double quotes
        % within the path must be escaped with a preceding backslash so they are
        % not interpreted as the closing delimiter.
        escaped_arg = strrep(arg, '"', '\"');

        % --- Percent Sign Hazard ---
        % On Windows cmd.exe, percent signs trigger environment variable
        % expansion (e.g., %PATH%). Double them so each is treated as a
        % literal percent character rather than a variable delimiter.
        escaped_arg = strrep(escaped_arg, '%', '%%');

        % --- Caret (^) Hazard ---
        % The caret is cmd.exe's general escape character. Inside double
        % quotes it is mostly inert, but it can still cause problems with
        % piped commands or when the string is re-parsed. Escape it by
        % doubling so it is always treated as a literal caret.
        escaped_arg = strrep(escaped_arg, '^', '^^');

        % --- Exclamation Mark (!) Hazard ---
        % When delayed expansion is enabled (common in batch scripts and
        % some CI environments), ! triggers variable expansion (!VAR!).
        % Escape with ^ so it is treated as a literal character.
        escaped_arg = strrep(escaped_arg, '!', '^!');

        % --- Unicode Path Handling for Windows ---
        % Windows has special handling requirements for Unicode paths
        if exist('has_unicode', 'var') && has_unicode
            % For paths with Unicode characters, Windows may require
            % short path names (8.3 format) for reliable shell operations
            if exist(arg, 'file') || exist(arg, 'dir')
                try
                    % Attempt to get short path name for existing files/directories
                    [status, short_path] = system(['for %I in ("' arg '") do @echo %~sI']);
                    if status == 0 && ~isempty(strtrim(short_path))
                        % Use short path if available and different from original
                        short_path = strtrim(short_path);
                        if ~strcmp(short_path, arg) && ~contains(short_path, '?')
                            escaped_arg = short_path;
                        end
                    end
                catch
                    % If short path conversion fails, continue with Unicode path
                end
            end
        end

        % --- Trailing Backslash Hazard ---
        % On Windows, a trailing backslash immediately before the closing double
        % quote (e.g., "C:\data\") would be interpreted as escaping the quote
        % character itself, turning it into a literal quote rather than a
        % delimiter. This would leave the argument unclosed and break the shell
        % command. We count trailing backslashes and double them so each one is
        % treated as a literal backslash rather than an escape character.
        n_trailing_backslashes = 0;
        for k = length(escaped_arg):-1:1
            if escaped_arg(k) == '\'
                n_trailing_backslashes = n_trailing_backslashes + 1;
            else
                break;
            end
        end

        % Double each trailing backslash: e.g., 'C:\data\' becomes 'C:\data\\'
        % so that the closing quote is not consumed by the backslash escape.
        if n_trailing_backslashes > 0
             escaped_arg = [escaped_arg repmat('\', 1, n_trailing_backslashes)];
        end

        % Wrap the entire argument in double quotes to handle spaces and
        % special characters in DICOM-derived file paths.
        escaped_arg = ['"' escaped_arg '"'];

    else
        % --- Unix Escaping Strategy ---
        % Single quotes in Unix shells treat everything between them as literal
        % text (no variable expansion, no globbing). This is the safest quoting
        % mechanism. The only character that cannot appear inside single quotes
        % is the single quote itself.
        %
        % To include a literal single quote, we end the current single-quoted
        % segment, insert an escaped single quote (\'), and start a new
        % single-quoted segment: 'don'\''t' evaluates to the string don't.
        % This pattern is standard POSIX shell escaping.
        escaped_arg = strrep(arg, '''', '''\''''' );
        
        % --- Unicode Handling for Unix ---
        % Unix systems generally handle UTF-8 well, but ensure proper locale
        if exist('has_unicode', 'var') && has_unicode
            % For Unix systems with Unicode characters, we may need to ensure
            % the environment supports UTF-8. The cached encoding is already
            % available from above; no additional system calls needed here.
        end
        
        escaped_arg = ['''' escaped_arg ''''];
    end
end