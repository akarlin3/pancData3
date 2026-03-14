function json_str = json_set_field(json_str, field_name, new_value)
%JSON_SET_FIELD Replace a field's value in a raw JSON string.
%   out = json_set_field(json_str, field_name, new_value) performs a
%   targeted regex replacement of the named field's value, preserving all
%   surrounding formatting (indentation, field order, comments, etc.).
%
%   Supported value types:
%     - character vector / string  -> written as a JSON string
%     - logical / numeric scalar   -> written as true/false or number
%     - empty array []             -> written as []
%
%   Example:
%     s = '{"dwi_type": "Standard", "skip_to_reload": false}';
%     s = json_set_field(s, 'dwi_type', 'dnCNN');
%     s = json_set_field(s, 'skip_to_reload', true);

    % Build the JSON literal string for the replacement value.
    % Each MATLAB type maps to its JSON equivalent:
    %   char/string -> "value"
    %   logical     -> true/false (not 0/1)
    %   empty []    -> []
    %   scalar num  -> full-precision number (%.17g preserves double precision)
    if ischar(new_value) || isstring(new_value)
        val_literal = ['"' char(new_value) '"'];
    elseif islogical(new_value)
        if new_value
            val_literal = 'true';
        else
            val_literal = 'false';
        end
    elseif isnumeric(new_value) && isempty(new_value)
        val_literal = '[]';
    elseif isnumeric(new_value) && isscalar(new_value)
        val_literal = num2str(new_value, '%.17g');
    else
        error('json_set_field:unsupportedType', ...
            'Only string, logical, scalar numeric, and empty values are supported.');
    end

    % Regex pattern matches: "field_name" : <old_value>
    % Capture group $1 = everything up to and including the colon + whitespace,
    % so the replacement preserves the original formatting/indentation.
    % The value alternation covers all JSON value types:
    %   "..."           - quoted string (with escaped chars)
    %   true|false|null - JSON keywords
    %   -?[\d.eE+\-]+  - integers, floats, scientific notation
    %   \[\]            - empty array
    pattern = ['("' field_name '"\s*:\s*)("(?:[^"\\]|\\.)*"|true|false|null|-?[\d.eE+\-]+|\[\])'];
    replacement = ['$1' val_literal];

    % 'once' ensures only the first occurrence is replaced (safe for
    % non-nested JSON; nested objects with duplicate keys are not supported)
    result = regexprep(json_str, pattern, replacement, 'once');

    if strcmp(result, json_str)
        % No replacement occurred. Distinguish between "field exists but
        % value pattern didn't match" vs "field not found at all".
        if isempty(regexp(json_str, ['"' field_name '"'], 'once'))
            error('json_set_field:fieldNotFound', ...
                'Field "%s" not found in JSON string.', field_name);
        end
        % Field exists but value didn't match any pattern — silently
        % return unchanged (may indicate an unsupported value format)
    end

    json_str = result;
end
