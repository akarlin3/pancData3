% SUBSREF  Subscripted reference for the Octave table shim.
%
%   Replaces MATLAB's table subscripted reference. Supports two access
%   patterns used by the pipeline:
%     1. Dot indexing:  t.ColName  or  t.ColName(idx)
%     2. Row slicing:   t(row_idx, :)
%
%   Behavioral differences from MATLAB's table:
%   - Brace indexing (t{:, 'col'}) is not supported.
%   - Column selection by name via parentheses is not supported; only
%     row slicing with all columns is implemented.
%   - No VariableNames property access (use fieldnames(struct(t)) instead).
function out = subsref(obj, S)
    % Handle property indexing like table.VarName and slicing like table(idx, :)

    % Unwrap to the underlying struct for field access.
    s = struct(obj);

    if strcmp(S(1).type, '.')
        % Dot-reference: retrieve a single column by variable name.
        field_name = S(1).subs;
        out = s.(field_name);

        % If there are more indices, apply them (e.g. table.VarName(1:5))
        % This delegates to MATLAB/Octave's built-in subsref for the column data.
        if length(S) > 1
            out = subsref(out, S(2:end));
        end
    elseif strcmp(S(1).type, '()')
        % Parenthetical row slicing: t(rows, :)
        % Extracts the specified rows from every column, returning a new table.
        row_idx = S(1).subs{1};
        % Col idx is S(1).subs{2} but we assume we want all columns if it's ':'

        fields = fieldnames(s);
        new_s = struct();

        for i = 1:length(fields)
            f = fields{i};
            % Slice each column by row index, preserving column orientation.
            new_s.(f) = s.(f)(row_idx, :); % keep it a column vector sliced
        end

        out = class(new_s, 'table');
    else
        error('Unsupported table indexing type: %s', S(1).type);
    end
end
