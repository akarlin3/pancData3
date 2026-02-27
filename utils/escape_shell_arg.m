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

    if nargin < 2
        if ispc
            style = 'pc';
        else
            style = 'unix';
        end
    end

    if ~ischar(arg) && ~isstring(arg)
        error('escape_shell_arg:invalidInput', 'Input argument must be a string or character vector.');
    end
    arg = char(arg);

    if strcmpi(style, 'pc')
        % Windows escaping
        % Escape double quotes with backslash
        escaped_arg = strrep(arg, '"', '\"');

        % Handle trailing backslash checking:
        % If the string ends in a backslash, we need to escape it so it doesn't escape the closing quote.
        % We iterate from the end to count trailing backslashes.
        n_trailing_backslashes = 0;
        for k = length(escaped_arg):-1:1
            if escaped_arg(k) == '\'
                n_trailing_backslashes = n_trailing_backslashes + 1;
            else
                break;
            end
        end

        % If there are trailing backslashes, double them only if they precede the closing quote
        if n_trailing_backslashes > 0
             escaped_arg = [escaped_arg repmat('\', 1, n_trailing_backslashes)];
        end

        escaped_arg = ['"' escaped_arg '"'];

    else
        % Unix escaping
        % Wrap in single quotes. Replace ' with '\''
        escaped_arg = strrep(arg, "'", "'\''");
        escaped_arg = ['''' escaped_arg ''''];
    end
end
