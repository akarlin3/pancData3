% CONTAINS  Octave-compatible shim for MATLAB's contains (R2016b+).
%
%   tf = contains(str, pattern) returns true if pattern is found within str.
%   tf = contains(str, pattern, 'IgnoreCase', true) performs case-insensitive
%       matching.
%
%   Supports both scalar strings and cell arrays of strings for str.
%   When str is a cell array, returns a logical array of the same size.
%
%   Behavioral differences from MATLAB's contains:
%   - Does not support string arrays (only char vectors and cell arrays of
%     char vectors), because Octave's string type support is limited.
%   - Pattern matching uses strfind (exact substring), matching MATLAB's
%     behavior. Multiple patterns (cell array of patterns) are not fully
%     supported -- only a single pattern string is matched.
function tf = contains(str, pattern, varargin)
    % Parse optional 'IgnoreCase' name-value pair.
    ignoreCase = false;
    for i = 1:2:length(varargin)
        if strcmpi(varargin{i}, 'IgnoreCase')
            ignoreCase = varargin{i+1};
        end
    end

    % If case-insensitive, convert both str and pattern to lowercase.
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

    % Use strfind for substring matching; convert empty results to false.
    if iscell(str)
        % Cell array input: return a logical array, one element per cell.
        tf = ~cellfun(@isempty, strfind(str, pattern));
    else
        % Scalar string input: return a scalar logical.
        tf = ~isempty(strfind(str, pattern));
    end
end
