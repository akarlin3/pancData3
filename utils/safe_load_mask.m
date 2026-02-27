function mask = safe_load_mask(filepath, varname)
% SAFE_LOAD_MASK Safely loads a mask variable from a .mat file.
%
%   mask = SAFE_LOAD_MASK(filepath, varname)
%
%   INPUTS:
%       filepath - Path to the .mat file.
%       varname  - (Optional) Name of the variable to load. Defaults to 'Stvol3d'.
%
%   OUTPUTS:
%       mask     - The loaded variable content if safe, or [] if unsafe/missing.
%
%   SECURITY:
%       This function inspects the .mat file content using 'whos' without
%       loading it into the workspace. It verifies that the target variable:
%       1. Exists in the file.
%       2. Is of a safe class (numeric, logical).
%       3. Is not a complex object or script that could execute code.
%
%   EXAMPLE:
%       mask = safe_load_mask('patient_data.mat', 'Stvol3d');

    if nargin < 2
        varname = 'Stvol3d';
    end

    mask = [];

    if ~exist(filepath, 'file')
        warning('safe_load_mask:FileNotFound', 'File not found: %s', filepath);
        return;
    end

    % Inspect file contents without loading
    try
        file_info = whos('-file', filepath);
    catch ME
        warning('safe_load_mask:FileReadError', 'Could not read file info: %s', ME.message);
        return;
    end

    % Find the target variable
    idx = find(strcmp({file_info.name}, varname));

    if isempty(idx)
        warning('safe_load_mask:VariableNotFound', ...
            'Variable ''%s'' not found in %s', varname, filepath);
        return;
    end

    target_info = file_info(idx);

    % Define safe classes
    safe_classes = {'double', 'single', 'logical', ...
                    'int8', 'uint8', 'int16', 'uint16', ...
                    'int32', 'uint32', 'int64', 'uint64'};

    if ~ismember(target_info.class, safe_classes)
        warning('safe_load_mask:SecurityRisk', ...
            'Security Risk: Variable ''%s'' is of unsafe class ''%s''. skipping.', ...
            varname, target_info.class);
        return;
    end

    % Load only the specific variable
    try
        tmp = load(filepath, varname);
        if isfield(tmp, varname)
            mask = tmp.(varname);
        end
    catch ME
        warning('safe_load_mask:LoadError', 'Error loading variable: %s', ME.message);
    end

end
