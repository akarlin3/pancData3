% TABLE  Octave-compatible shim for MATLAB's built-in table class.
%
%   MATLAB's table data type (introduced in R2013b) is not available in
%   GNU Octave. This shim implements a minimal table class backed by a
%   plain struct, where each field corresponds to a table variable (column).
%
%   Supported syntax:
%       t = table(col1, col2, ..., 'VariableNames', {'A', 'B', ...})
%
%   Behavioral differences from MATLAB's table:
%   - No RowNames support.
%   - No variable type enforcement or heterogeneous column handling.
%   - Indexing is limited to dot-access and row slicing (see subsref/subsasgn).
%   - height(), width(), and other table methods are not implemented.
%   - Columns without names are auto-named 'Var1', 'Var2', etc.
function obj = table(varargin)
    % TABLE class constructor for Octave compatibility
    % Creates a table object from column data and optional variable names.
    s = struct();
    % First pass: scan varargin for the 'VariableNames' name-value pair.
    % Everything before it is treated as column data.
    varNames = {};
    for i = 1:length(varargin)
        if ischar(varargin{i}) && strcmp(varargin{i}, 'VariableNames')
            if i < length(varargin)
                varNames = varargin{i+1};
            end
            break;
        end
    end

    % Second pass: assign each data argument to a struct field.
    % Stop when we hit the 'VariableNames' key (it and its value are not data).
    idx = 1;
    for i = 1:length(varargin)
        if ischar(varargin{i}) && strcmp(varargin{i}, 'VariableNames')
            break;
        end
        % Use the user-supplied name if available; otherwise auto-generate.
        if idx <= length(varNames)
            name = varNames{idx};
        else
            name = sprintf('Var%d', idx);
        end

        val = varargin{i};
        % Force column orientation to match MATLAB table convention where
        % each variable is stored as a column vector.
        if isnumeric(val) || islogical(val) || iscell(val)
            val = val(:); % Ensure column vector
        end

        s.(name) = val;
        idx = idx + 1;
    end

    % Octave's old-style class() constructor wraps the struct as a 'table' object.
    obj = class(s, 'table');
end
