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

    % Check for duplicate variable name entries. In well-formed .mat files
    % each variable name appears exactly once in the whos listing. However,
    % .mat v7.3 (HDF5-based) files that are manually constructed or
    % corrupted via HDF5 tools or matfile operations can contain multiple
    % entries with the same name but different sizes or classes. If we only
    % validated the first entry, an unsafe duplicate could slip through
    % undetected. We therefore validate ALL duplicate entries or reject the
    % file outright.
    if numel(idx) > 1
        warning('safe_load_mask:DuplicateVariable', ...
            ['Variable ''%s'' appears %d times in %s. ' ...
             'This may indicate a corrupted or tampered .mat file. ' ...
             'All duplicate entries will be validated before proceeding.'], ...
            varname, numel(idx), filepath);

        % Define safe classes here as well so we can validate all entries
        safe_classes_dup = {'double', 'single', 'logical', ...
                        'int8', 'uint8', 'int16', 'uint16', ...
                        'int32', 'uint32', 'int64', 'uint64'};

        % Validate every duplicate entry for class safety and size limits
        for dup_i = 1:numel(idx)
            dup_info = file_info(idx(dup_i));

            % Check class safety for this duplicate
            if ~ismember(dup_info.class, safe_classes_dup)
                error('safe_load_mask:SecurityRisk', ...
                    ['Security Risk: Duplicate entry %d of variable ''%s'' has unsafe class ''%s''. ' ...
                     'Aborting load due to potentially corrupted/tampered file.'], ...
                    dup_i, varname, dup_info.class);
            end

            % Check size limit for this duplicate
            dup_size_mb = dup_info.bytes / (1024 * 1024);
            if dup_size_mb > max_file_size_mb
                error('safe_load_mask:VariableTooLarge', ...
                    ['Duplicate entry %d of variable ''%s'' in-memory size (%.1f MB) ' ...
                     'exceeds safety limit (%.1f MB). Aborting load.'], ...
                    dup_i, varname, dup_size_mb, max_file_size_mb);
            end
        end

        % Check that all duplicates agree on class and size. If they
        % differ, the file is likely corrupted or tampered with and we
        % should not trust any single entry.
        dup_classes = {file_info(idx).class};
        dup_sizes = {file_info(idx).size};
        classes_match = all(strcmp(dup_classes, dup_classes{1}));
        sizes_match = true;
        for dup_i = 2:numel(idx)
            if ~isequal(dup_sizes{dup_i}, dup_sizes{1})
                sizes_match = false;
                break;
            end
        end

        if ~classes_match || ~sizes_match
            error('safe_load_mask:DuplicateMismatch', ...
                ['Duplicate entries for variable ''%s'' have inconsistent class or size. ' ...
                 'File is likely corrupted or tampered. Aborting load for security.'], ...
                varname);
        end
    end

    % Use the first matching entry for subsequent checks (all duplicates
    % have been validated above if there were multiple).
    target_info = file_info(idx(1));
    
    % Additional validation: check if the reported in-memory variable
    % size is reasonable.  We use the same threshold as the file-level
    % check (max_file_size_mb).
    %
    % NOTE: For .mat v7.3 (HDF5-based) files, the 'bytes' field
    % returned by whos('-file',...) reports the *uncompressed* in-memory
    % size, NOT the on-disk size.  A compressed .mat file that passes the
    % on-disk file size check above could therefore contain a variable
    % whose reported bytes exceed the file's physical size on disk.  Using
    % the same max_file_size_mb limit (rather than a fraction of it)
    % avoids false rejections of legitimate high-resolution mask data
    % stored in compressed .mat files.
    var_bytes = target_info.bytes;
    var_size_mb = var_bytes / (1024 * 1024);
    
    if var_size_mb > max_file_size_mb
        warning('safe_load_mask:VariableTooLarge', ...
            'Variable ''%s'' in-memory size (%.1f MB) exceeds safety limit (%.1f MB)', ...
            varname, var_size_mb, max_file_size_mb);
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
    % that might indicate a malicious file.
    %
    % We use exact matching (ismember) for known dangerous MATLAB class
    % names and startsWith for Java classes (which always begin with a
    % Java package prefix such as 'java.', 'javax.', 'com.', 'org.', 'net.').
    % This avoids false positives from substring matching (e.g. a
    % legitimate class named 'javascript_parser' would NOT be flagged).
    %
    % Defense-in-depth: 'struct' and 'cell' variables are also checked
    % for non-target variables. whos('-file') only reports the top-level
    % class; a struct or cell could contain nested unsafe types (e.g.,
    % function_handle) that are invisible to whos inspection. Since we
    % only load the target variable (which has already been validated as
    % a safe primitive type), the presence of struct/cell in OTHER
    % variables does not block loading but is logged as a warning so
    % that analysts are aware the file may contain potentially unsafe
    % nested content.
    exact_suspicious_classes = {'function_handle', 'onCleanup', 'timer'};
    java_prefixes = {'java.', 'javax.', 'com.', 'org.', 'net.'};
    opaque_container_classes = {'struct', 'cell'};
    for i = 1:length(file_info)
        var_class = file_info(i).class;
        is_exact_suspicious = ismember(var_class, exact_suspicious_classes);
        is_java_class = any(startsWith(var_class, java_prefixes));
        if is_exact_suspicious || is_java_class
            warning('safe_load_mask:SuspiciousContent', ...
                'File contains suspicious variable of class ''%s''. Aborting load for security.', ...
                var_class);
            return;
        end
        % For non-target variables, warn about opaque container types
        % (struct, cell) that could hide nested unsafe types such as
        % function_handle. whos('-file') cannot inspect their contents,
        % so we log a defense-in-depth warning. Loading proceeds because
        % only the individually validated target variable is loaded.
        is_non_target = ~strcmp(file_info(i).name, varname);
        is_opaque_container = ismember(var_class, opaque_container_classes);
        if is_non_target && is_opaque_container
            warning('safe_load_mask:OpaqueContainerInFile', ...
                ['Non-target variable ''%s'' has class ''%s'' which may contain ' ...
                 'nested unsafe types (e.g., function_handle) invisible to whos inspection. ' ...
                 'Only the validated target variable ''%s'' will be loaded.'], ...
                file_info(i).name, var_class, varname);
        end
    end

    % Validate that the variable name is a valid MATLAB identifier.
    % whos('-file') can find variables with names that are not valid
    % identifiers (e.g., names starting with a digit like '3d_mask', or
    % containing spaces). Such variables can be created via matfile
    % objects or HDF5 tools.  load() returns a struct, and dynamic field
    % access via loaded_struct.(varname) would fail for non-identifier
    % names.  We check isvarname() to choose the appropriate access
    % strategy.
    varname_is_valid_identifier = isvarname(varname);

    % Load only the specific variable by name to minimize memory usage
    % and avoid deserializing any other (potentially unsafe) variables
    % in the file.  The whos check above has already verified this
    % variable's class is safe, so the load is secure.
    try
        if varname_is_valid_identifier
            % Standard path: load returns a struct and we access the
            % field via dynamic field reference.
            tmp = load(filepath, varname);
            if isfield(tmp, varname)
                mask = tmp.(varname);
            else
                % Unexpected: load succeeded but field not present.
                % Fall back to matfile access.
                warning('safe_load_mask:FieldMissing', ...
                    'load() succeeded but struct field ''%s'' not found. Trying matfile fallback.', ...
                    varname);
                mask = safe_load_mask_via_matfile(filepath, varname);
            end
        else
            % Non-standard variable name (not a valid MATLAB identifier).
            % Dynamic struct field access would fail, so use matfile()
            % which can handle arbitrary variable names.
            warning('safe_load_mask:NonStandardVarName', ...
                ['Variable name ''%s'' is not a valid MATLAB identifier. ' ...
                 'Using matfile() for access.'], varname);
            mask = safe_load_mask_via_matfile(filepath, varname);
        end
    catch ME
        % If the primary method failed, attempt matfile fallback (only
        % if we haven't already tried it).
        if varname_is_valid_identifier
            try
                warning('safe_load_mask:LoadErrorFallback', ...
                    'Primary load failed (%s). Attempting matfile() fallback.', ME.message);
                mask = safe_load_mask_via_matfile(filepath, varname);
            catch ME2
                warning('safe_load_mask:LoadError', ...
                    'Error loading variable via both methods: primary [%s], fallback [%s]', ...
                    ME.message, ME2.message);
                mask = [];
            end
        else
            warning('safe_load_mask:LoadError', 'Error loading variable: %s', ME.message);
            mask = [];
        end
    end

end


function data = safe_load_mask_via_matfile(filepath, varname)
% SAFE_LOAD_MASK_VIA_MATFILE Loads a variable using matfile() object.
%   This handles variable names that are not valid MATLAB identifiers,
%   which cannot be accessed via dynamic struct field reference after
%   a standard load() call.
%
%   matfile() uses partial loading (does not load the entire file) and
%   supports arbitrary variable names through its parenthetical access
%   syntax.

    mf = matfile(filepath);
    try
        % matfile objects support dynamic property access for valid
        % identifiers and subsref-based access for others.
        % Use the who() method to confirm the variable is accessible.
        vars_in_file = who(mf);
        if ~ismember(varname, vars_in_file)
            warning('safe_load_mask:MatfileVarNotFound', ...
                'Variable ''%s'' not found via matfile access.', varname);
            data = [];
            return;
        end

        % Access the variable. For valid identifiers, dynamic property
        % access works. For non-identifiers, we use subsref directly.
        if isvarname(varname)
            data = mf.(varname);
        else
            % Use subsref with '.' type referencing, which matfile
            % supports even for non-standard variable names.
            s.type = '.';
            s.subs = varname;
            data = subsref(mf, s);
        end
    catch ME
        warning('safe_load_mask:MatfileAccessError', ...
            'matfile() access failed for variable ''%s'': %s', varname, ME.message);
        data = [];
    end

end