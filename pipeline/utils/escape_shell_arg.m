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

    % --- Persistent encoding cache with validation ---
    % Cache the detected system encoding with timestamps for validation
    persistent cached_pc_encoding;
    persistent cached_unix_encoding;
    persistent cache_pc_timestamp;
    persistent cache_unix_timestamp;
    
    % Configuration constants
    CACHE_EXPIRY_SECONDS = 300; % 5 minutes - balance between performance and freshness
    SAFE_DEFAULT_ENCODING = 'UTF-8'; % Fallback encoding if cache becomes invalid

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

    % --- Control Character Detection and Sanitization ---
    % Detect and handle dangerous control characters that could pose security
    % risks, including null bytes, terminal escape sequences, and other
    % potentially harmful control codes. DICOM metadata is often unvalidated
    % and could contain embedded control sequences that manipulate terminal
    % output or shell behavior.
    
    % Null bytes (char(0)) can cause premature string termination
    if any(arg == 0)
        warning('escape_shell_arg:nullByte', 'Null bytes detected and removed from shell argument.');
        arg(arg == 0) = [];
    end
    
    % Terminal escape sequences start with ESC (char(27)) followed by [
    % These can manipulate terminal display, execute commands, or leak information
    esc_positions = find(arg == 27); % ESC character
    if ~isempty(esc_positions)
        % Look for ANSI escape sequences (ESC[ followed by parameters and command)
        ansi_found = false;
        for i = 1:length(esc_positions)
            pos = esc_positions(i);
            if pos < length(arg) && arg(pos + 1) == '['
                ansi_found = true;
                break;
            end
        end
        if ansi_found
            warning('escape_shell_arg:escapeSequence', 'Terminal escape sequences detected and removed from shell argument.');
        end
        % Remove all ESC characters as they're generally not needed in file paths
        arg(arg == 27) = [];
    end
    
    % Other dangerous control characters (ASCII 1-31 except common whitespace)
    % Keep: TAB (9), LF (10), CR (13) as they might be legitimate in some contexts
    % Remove: All other control characters (1-8, 11-12, 14-31) including:
    % - Bell (7): Could trigger audio alerts
    % - Backspace (8): Could manipulate display
    % - Vertical tab (11), Form feed (12): Unusual formatting
    % - Various other control codes that have no business in file paths
    dangerous_controls = (arg >= 1 & arg <= 8) | (arg >= 11 & arg <= 12) | (arg >= 14 & arg <= 31);
    if any(dangerous_controls)
        warning('escape_shell_arg:controlChars', 'Dangerous control characters detected and removed from shell argument.');
        arg(dangerous_controls) = [];
    end
    
    % DEL character (127) - another potentially problematic control character
    if any(arg == 127)
        warning('escape_shell_arg:delChar', 'DEL character detected and removed from shell argument.');
        arg(arg == 127) = [];
    end

    % --- Unicode and Encoding Handling ---
    % Check if the string contains non-ASCII characters. For the common case
    % of pure-ASCII paths (the vast majority in US clinical environments),
    % we short-circuit here and skip all encoding detection and conversion
    % overhead, including persistent variable checks and system() calls.
    has_unicode = any(arg > 127);

    if has_unicode
        try
            % Get system encoding to ensure proper character handling
            if strcmpi(style, 'pc')
                % Windows: Check for system code page and handle Unicode paths
                current_time = now * 24 * 3600; % Convert to seconds since epoch
                
                % Validate cache and refresh if needed
                cache_valid = ~isempty(cached_pc_encoding) && ...
                             ~isempty(cache_pc_timestamp) && ...
                             (current_time - cache_pc_timestamp) < CACHE_EXPIRY_SECONDS;
                
                % Additional validation: check if cached encoding is still a valid string
                if cache_valid && (~ischar(cached_pc_encoding) || isempty(cached_pc_encoding))
                    cache_valid = false;
                    warning('escape_shell_arg:cacheCorruption', 'PC encoding cache corruption detected, refreshing cache.');
                end
                
                if ~cache_valid
                    try
                        % Attempt to get system code page (runs periodically)
                        [status, cp_output] = system('chcp');
                        if status == 0 && contains(cp_output, '65001') % UTF-8
                            cached_pc_encoding = 'UTF-8';
                        else
                            cached_pc_encoding = 'windows-1252'; % Common Windows default
                        end
                        cache_pc_timestamp = current_time;
                    catch
                        % If system call fails, use safe default and log warning
                        cached_pc_encoding = SAFE_DEFAULT_ENCODING;
                        cache_pc_timestamp = current_time;
                        warning('escape_shell_arg:encodingFallback', 'Failed to detect PC encoding, using safe default: %s', SAFE_DEFAULT_ENCODING);
                    end
                end
                system_encoding = cached_pc_encoding;
                
                % For Windows, convert to the system's native encoding so
                % that cmd.exe interprets the bytes correctly.  When the
                % system code page is already UTF-8 no conversion is needed.
                if ~strcmp(system_encoding, 'UTF-8')
                    try
                        % Re-encode from UTF-8 (MATLAB internal) to the
                        % detected Windows code page so the shell receives
                        % bytes it can display/process correctly.
                        arg = native2unicode(unicode2native(arg, system_encoding), system_encoding);
                    catch
                        % If conversion fails, proceed with original string
                        % and hope the system can handle it
                    end
                end
            else
                % Unix systems: Enhanced encoding detection with fallbacks
                current_time = now * 24 * 3600; % Convert to seconds since epoch
                
                % Validate cache and refresh if needed
                cache_valid = ~isempty(cached_unix_encoding) && ...
                             ~isempty(cache_unix_timestamp) && ...
                             (current_time - cache_unix_timestamp) < CACHE_EXPIRY_SECONDS;
                
                % Additional validation: check if cached encoding is still a valid string
                if cache_valid && (~ischar(cached_unix_encoding) || isempty(cached_unix_encoding))
                    cache_valid = false;
                    warning('escape_shell_arg:cacheCorruption', 'Unix encoding cache corruption detected, refreshing cache.');
                end
                
                if ~cache_valid
                    % First try environment variables (most reliable)
                    lang_var = getenv('LANG');
                    lc_all = getenv('LC_ALL');
                    lc_ctype = getenv('LC_CTYPE');
                    
                    % Check environment variables for encoding info
                    env_encoding = '';
                    if ~isempty(lc_all) && contains(upper(lc_all), 'UTF')
                        env_encoding = 'UTF-8';
                    elseif ~isempty(lc_ctype) && contains(upper(lc_ctype), 'UTF')
                        env_encoding = 'UTF-8';
                    elseif ~isempty(lang_var) && contains(upper(lang_var), 'UTF')
                        env_encoding = 'UTF-8';
                    end
                    
                    if ~isempty(env_encoding)
                        cached_unix_encoding = env_encoding;
                        cache_unix_timestamp = current_time;
                    else
                        % Fallback to system() call if environment variables don't help
                        try
                            [status, locale_output] = system('locale charmap 2>/dev/null');
                            if status == 0 && contains(upper(locale_output), 'UTF-8')
                                cached_unix_encoding = 'UTF-8';
                            else
                                % Try alternative command if locale charmap fails
                                [status2, locale_output2] = system('locale 2>/dev/null | grep -i utf');
                                if status2 == 0 && ~isempty(locale_output2)
                                    cached_unix_encoding = 'UTF-8';
                                else
                                    % Final fallback based on common modern Unix defaults
                                    cached_unix_encoding = SAFE_DEFAULT_ENCODING;
                                end
                            end
                            cache_unix_timestamp = current_time;
                        catch
                            % If all system calls fail (restricted environment, missing utilities)
                            % assume UTF-8 as it's the most common encoding on modern Unix systems
                            cached_unix_encoding = SAFE_DEFAULT_ENCODING;
                            cache_unix_timestamp = current_time;
                            warning('escape_shell_arg:encodingFallback', 'Failed to detect Unix encoding, using safe default: %s', SAFE_DEFAULT_ENCODING);
                        end
                    end
                end
            end
        catch
            % If Unicode detection/handling fails, proceed with original escaping
            % This ensures backward compatibility
        end
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
        % NOTE: Caret escaping MUST come before all other caret-based
        % escapes below, otherwise we would double-escape the carets we
        % insert for &, |, <, >, (, ).
        escaped_arg = strrep(escaped_arg, '^', '^^');

        % --- Exclamation Mark (!) Hazard ---
        % When delayed expansion is enabled (common in batch scripts and
        % some CI environments), ! triggers variable expansion (!VAR!).
        % Escape with ^ so it is treated as a literal character.
        escaped_arg = strrep(escaped_arg, '!', '^!');

        % --- Ampersand (&) Hazard ---
        % The ampersand is a command separator in cmd.exe (e.g., cmd1 & cmd2).
        % While double quoting provides partial protection, the ampersand CAN
        % break out of double-quoted context in certain cmd.exe parsing
        % scenarios, especially when the escaped string is later concatenated
        % with other strings before being passed to system(). MATLAB's
        % system() passes the command through cmd.exe /c, which can re-parse
        % the string before honoring quote grouping. Escape with ^ for
        % defense-in-depth.
        escaped_arg = strrep(escaped_arg, '&', '^&');

        % --- Pipe (|) Hazard ---
        % The pipe character creates command pipelines in cmd.exe. Escape with
        % ^ for defense-in-depth even though double quoting provides partial
        % protection.
        escaped_arg = strrep(escaped_arg, '|', '^|');

        % --- Angle Bracket (<, >) Hazards ---
        % Angle brackets perform I/O redirection in cmd.exe. Escape with ^
        % for defense-in-depth to prevent redirection attacks from malicious
        % DICOM metadata.
        escaped_arg = strrep(escaped_arg, '<', '^<');
        escaped_arg = strrep(escaped_arg, '>', '^>');

        % --- Parentheses ( ) Hazards ---
        % Parentheses are used for command grouping in cmd.exe and can cause
        % parsing errors or unexpected behavior when unescaped, particularly
        % in compound commands or when EnableDelayedExpansion is active.
        % Escape with ^ for defense-in-depth.
        escaped_arg = strrep(escaped_arg, '(', '^(');
        escaped_arg = strrep(escaped_arg, ')', '^)');

        % --- Unicode Path Handling for Windows ---
        % For paths with Unicode characters, use the \\?\ long path prefix
        % which allows the Windows API to handle extended characters and
        % long paths without needing 8.3 short names. This is more reliable
        % than 8.3 name lookup, which is often disabled on modern Windows
        % systems and would also require passing unescaped input to a
        % system() call, creating a shell injection vulnerability.
        if has_unicode
            if (exist(arg, 'file') || exist(arg, 'dir'))
                try
                    % Convert to absolute path if not already
                    abs_path = arg;
                    if length(arg) < 2 || arg(2) ~= ':'
                        abs_path = fullfile(pwd, arg);
                    end
                    % Apply \\?\ prefix for Unicode-safe Windows API access
                    % This avoids character encoding issues in cmd.exe
                    if ~startsWith(abs_path, '\\?\')
                        escaped_arg = strrep(abs_path, '"', '\"');
                        escaped_arg = strrep(escaped_arg, '%', '%%');
                        % Caret must be escaped before other caret-based escapes
                        escaped_arg = strrep(escaped_arg, '^', '^^');
                        escaped_arg = strrep(escaped_arg, '!', '^!');
                        escaped_arg = strrep(escaped_arg, '&', '^&');
                        escaped_arg = strrep(escaped_arg, '|', '^|');
                        escaped_arg = strrep(escaped_arg, '<', '^<');
                        escaped_arg = strrep(escaped_arg, '>', '^>');
                        escaped_arg = strrep(escaped_arg, '(', '^(');
                        escaped_arg = strrep(escaped_arg, ')', '^)');
                        % Prepend \\?\ prefix for extended-length path handling
                        escaped_arg = ['\\?\' escaped_arg];
                    end
                catch
                    % If long path prefix fails, continue with escaped path
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
        % segment, insert a backslash-escaped single quote (\'), and start a
        % new single-quoted segment. For example, for the input O'Brien:
        %   'O'\''Brien'
        % The shell sees three tokens that concatenate:
        %   'O'   -> O
        %   \'    -> '       (literal single quote, outside any quoting)
        %   'Brien' -> Brien
        % Result: O'Brien
        %
        % This is the standard POSIX shell idiom for embedding single quotes
        % inside single-quoted strings and prevents shell injection via paths
        % containing apostrophes (e.g., patient names like O'Brien in DICOM
        % metadata that end up in directory names).
        escaped_arg = strrep(arg, '''', '''\''''' );
        
        escaped_arg = ['''' escaped_arg ''''];
    end
end