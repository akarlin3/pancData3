function out = subsref(obj, S)
    % Handle property indexing like table.VarName and slicing like table(idx, :)

    s = struct(obj);

    if strcmp(S(1).type, '.')
        % Property access
        field_name = S(1).subs;
        out = s.(field_name);

        % If there are more indices, apply them (e.g. table.VarName(1:5))
        if length(S) > 1
            out = subsref(out, S(2:end));
        end
    elseif strcmp(S(1).type, '()')
        % Slicing table(row_idx, col_idx)
        row_idx = S(1).subs{1};
        % Col idx is S(1).subs{2} but we assume we want all columns if it's ':'

        fields = fieldnames(s);
        new_s = struct();

        for i = 1:length(fields)
            f = fields{i};
            new_s.(f) = s.(f)(row_idx, :); % keep it a column vector sliced
        end

        out = class(new_s, 'table');
    else
        error('Unsupported table indexing type: %s', S(1).type);
    end
end
