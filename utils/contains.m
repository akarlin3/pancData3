function tf = contains(str, pattern, varargin)
    % CONTAINS True if pattern is found in text.
    % Custom implementation for Octave compatibility.

    ignoreCase = false;
    for i = 1:2:length(varargin)
        if strcmpi(varargin{i}, 'IgnoreCase')
            ignoreCase = varargin{i+1};
        end
    end

    if ignoreCase
        if iscell(str)
            str = cellfun(@lower, str, 'UniformOutput', false);
        else
            str = lower(str);
        end
        if iscell(pattern)
            pattern = cellfun(@lower, pattern, 'UniformOutput', false);
        else
            pattern = lower(pattern);
        end
    end

    if iscell(str)
        tf = ~cellfun(@isempty, strfind(str, pattern));
    else
        tf = ~isempty(strfind(str, pattern));
    end
end
