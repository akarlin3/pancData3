function mask = safe_load_mask(filepath, varname, max_file_size_mb)
% SAFE_LOAD_MASK Safely loads a mask variable from a .mat file.
%
%   mask = SAFE_LOAD_MASK(filepath, varname, max_file_size_mb)
%
%   INPUTS:
%       filepath        - Path to the .mat file.
%       varname         - (Optional) Name of the variable to load. Defaults to 'Stvol3d'.
%       max_file_size_mb- (Optional) Maximum file size in MB. Defaults to 100 MB.
%
%   OUTPUTS:
%       mask            - The loaded variable content if safe, or [] if unsafe/missing.
%
%   SECURITY:
%       This function inspects the .mat file content using 'whos' without
%       loading it into the workspace. It verifies that the target variable:
%       1. Exists in the file.
%       2. Is of a safe class (numeric, logical).
%       3. Is not a complex object or script that could execute code.
%       4. The file size is within acceptable limits.
%       5. The file has a valid .mat structure.
%
%   Analytical Rationale — Why Safe Loading is Critical:
%   ----------------------------------------------------
%   GTV mask files (.mat) are often produced by external contouring
%   software or shared across institutions.  MATLAB's `load` function can
%   deserialize arbitrary objects, including classes with custom `loadobj`
%   methods that execute code upon loading.  A maliciously crafted .mat
%   file could exploit this to run arbitrary code in the MATLAB session
%   that has access to patient data.  By inspecting the variable class via
%   `whos('-file',...)` BEFORE loading, we reject any non-numeric class
%   (e.g., function_handle, java objects, custom classes) without ever
%   executing their deserialization logic.
%
%   The default variable name 'Stvol3d' corresponds to the standard GTV
%   structure volume variable produced by the contouring pipeline.
%
%   EXAMPLE:
%       mask = safe_load_mask('patient_data.mat', 'Stvol3d', 200);

    if nargin < 2
        varname = 'Stvol3d';
    end
    
    if nargin < 3
        max_file_size_mb = 100;  % Default 100 MB limit
    end

    mask = [];

    if ~exist(filepath, 'file')
        warning('safe_load_mask:FileNotFound', 'File not found: %s', filepath);
        return;
    end

    % Check file size first to prevent loading oversized files
    file_info_sys = dir(filepath);
    if isempty(file_info_sys)
        warning('safe_load_mask:FileAccessError', 'Cannot access file: %s', filepath);
        return;
    end
    
    file_size_mb = file_info_sys.bytes / (1024 * 1024);
    if file_size_mb > max_file_size_mb
        warning('safe_load_mask:FileTooLarge', ...
            'File size (%.1f MB) exceeds maximum allowed size (%.1f MB)', ...
            file_size_mb, max_file_size_mb);
        return;
    end

    % Inspect file contents without loading into the workspace.  whos
    % reads only the variable metadata (name, size, class) from the .mat
    % file header — it does NOT deserialize the data or execute any
    % loadobj methods, making it safe for untrusted files.
    try
        file_info = whos('-file', filepath);
    catch ME
        warning('safe_load_mask:FileReadError', ...
            'Could not read file info (possibly corrupted .mat file): %s', ME.message);
        return;
    end
    
    % Validate .mat file structure - should have at least one variable
    if isempty(file_info)
        warning('safe_load_mask:InvalidStructure', ...
            'File appears to be empty or not a valid .mat file: %s', filepath);
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
    
    % Additional validation: check if variable size is reasonable
    var_bytes = target_info.bytes;
    var_size_mb = var_bytes / (1024 * 1024);
    
    % Variable size should not exceed 80% of the maximum file size limit
    max_var_size_mb = max_file_size_mb * 0.8;
    if var_size_mb > max_var_size_mb
        warning('safe_load_mask:VariableTooLarge', ...
            'Variable ''%s'' size (%.1f MB) exceeds safety limit (%.1f MB)', ...
            varname, var_size_mb, max_var_size_mb);
        return;
    end

    % Define safe classes — only primitive numeric and logical types are
    % permitted.  These types have no custom deserialization behavior and
    % cannot execute arbitrary code.  GTV masks are always representable
    % as logical (binary mask) or numeric (probability maps, label maps)
    % arrays, so no legitimate mask file should contain other classes.
    % Notably excluded: 'struct' (could contain nested unsafe types),
    % 'cell' (could contain function_handles), 'function_handle', and
    % any user-defined class.
    safe_classes = {'double', 'single', 'logical', ...
                    'int8', 'uint8', 'int16', 'uint16', ...
                    'int32', 'uint32', 'int64', 'uint64'};

    if ~ismember(target_info.class, safe_classes)
        warning('safe_load_mask:SecurityRisk', ...
            'Security Risk: Variable ''%s'' is of unsafe class ''%s''. skipping.', ...
            varname, target_info.class);
        return;
    end
    
    % Additional structure validation: check for suspicious variables
    % that might indicate a malicious file
    suspicious_classes = {'function_handle', 'onCleanup', 'timer', 'java'};
    for i = 1:length(file_info)
        if any(contains(file_info(i).class, suspicious_classes))
            warning('safe_load_mask:SuspiciousContent', ...
                'File contains suspicious variable of class ''%s''. Aborting load for security.', ...
                file_info(i).class);
            return;
        end
    end

    % Load only the specific variable by name to minimize memory usage
    % and avoid deserializing any other (potentially unsafe) variables
    % in the file.  The whos check above has already verified this
    % variable's class is safe, so the load is secure.
    try
        tmp = load(filepath, varname);
        if isfield(tmp, varname)
            mask = tmp.(varname);
        end
    catch ME
        warning('safe_load_mask:LoadError', 'Error loading variable: %s', ME.message);
    end

end