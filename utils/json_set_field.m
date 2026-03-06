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

    % Build the replacement value literal
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

    % Pattern: "field_name" <colon> <optional whitespace> <old value>
    % Old value is one of: a quoted string, true/false/null, a number, or []
    pattern = ['("' field_name '"\s*:\s*)("(?:[^"\\]|\\.)*"|true|false|null|-?[\d.eE+\-]+|\[\])'];
    replacement = ['$1' val_literal];

    result = regexprep(json_str, pattern, replacement, 'once');

    if strcmp(result, json_str)
        % Check if the field exists at all
        if isempty(regexp(json_str, ['"' field_name '"'], 'once'))
            error('json_set_field:fieldNotFound', ...
                'Field "%s" not found in JSON string.', field_name);
        end
    end

    json_str = result;
end
