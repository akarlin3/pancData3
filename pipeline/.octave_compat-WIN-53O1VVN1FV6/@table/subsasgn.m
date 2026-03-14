% SUBSASGN  Subscripted assignment for the Octave table shim.
%
%   Replaces MATLAB's table subscripted assignment. Only dot-reference
%   assignment is supported (e.g., t.NewCol = values or t.Col(idx) = val).
%   Parenthetical and brace assignment (e.g., t{:,1} = ...) are not
%   implemented because the pipeline does not use them.
%
%   Behavioral differences from MATLAB's table:
%   - No {:,...} or (...) assignment support.
%   - No type checking or size validation on assigned values.
%   - Assigning to a new field name implicitly adds a column (like MATLAB).
function obj = subsasgn(obj, S, val)
    % Unwrap the table object back to its underlying struct for manipulation.
    s = struct(obj);

    if strcmp(S(1).type, '.')
        field_name = S(1).subs;

        if length(S) == 1
            % Simple dot assignment: t.ColName = values
            % e.g. obj.ADC_z = val
            s.(field_name) = val;
        else
            % Chained indexing: t.ColName(idx) = val
            % Retrieve existing column (or initialize empty), apply the
            % remaining subscript chain via built-in subsasgn, then store back.
            % e.g. obj.ADC_z(idx) = val
            if isfield(s, field_name)
                temp = s.(field_name);
            else
                temp = [];
            end
            temp = subsasgn(temp, S(2:end), val);
            s.(field_name) = temp;
        end
    else
        % Parenthetical or brace assignment not supported in this shim.
        error('Unsupported subsasgn for table');
    end

    % Rewrap the modified struct as a table object.
    obj = class(s, 'table');
end
