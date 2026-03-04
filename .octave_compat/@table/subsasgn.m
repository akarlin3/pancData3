function obj = subsasgn(obj, S, val)
    s = struct(obj);

    if strcmp(S(1).type, '.')
        field_name = S(1).subs;

        if length(S) == 1
            % e.g. obj.ADC_z = val
            s.(field_name) = val;
        else
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
        error('Unsupported subsasgn for table');
    end

    obj = class(s, 'table');
end
