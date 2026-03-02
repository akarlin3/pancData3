function obj = table(varargin)
    % TABLE class constructor for Octave compatibility
    s = struct();
    varNames = {};
    for i = 1:length(varargin)
        if ischar(varargin{i}) && strcmp(varargin{i}, 'VariableNames')
            if i < length(varargin)
                varNames = varargin{i+1};
            end
            break;
        end
    end

    idx = 1;
    for i = 1:length(varargin)
        if ischar(varargin{i}) && strcmp(varargin{i}, 'VariableNames')
            break;
        end
        if idx <= length(varNames)
            name = varNames{idx};
        else
            name = sprintf('Var%d', idx);
        end

        val = varargin{i};
        if isnumeric(val) || islogical(val) || iscell(val)
            val = val(:); % Ensure column vector
        end

        s.(name) = val;
        idx = idx + 1;
    end

    obj = class(s, 'table');
end
