function global_struct = align_and_assign_struct(global_struct, new_struct, index)
    % ALIGN_AND_ASSIGN_STRUCT Helper to assign struct arrays with potentially missing fields
    %   Extracted local helper from load_dwi_data.m.
    %   Different patients may have different sets of available data (e.g.,
    %   one patient has DnCNN results while another does not), resulting in
    %   struct arrays with different field sets. MATLAB requires all elements
    %   of a struct array to have identical fields. This helper reconciles
    %   field differences by adding empty placeholders for missing fields in
    %   both the global and new structs before performing the assignment.

    % Per-patient checkpoint data is stored as nFx x nRp (no patient dim).
    % Reshape to 1 x nFx x nRp so it can be slotted into the global
    % nPatients x nFx x nRp array at (index, :, :).
    if ndims(new_struct) <= 2 && size(new_struct, 1) > 1
        new_struct = reshape(new_struct, [1, size(new_struct)]);
    end

    if isempty(fieldnames(global_struct))
        % Initialise global struct: create a matching-fields template so
        % MATLAB can perform subscripted assignment at any index.
        fields = fieldnames(new_struct);
        empty_vals = repmat({[]}, numel(fields), 1);
        template = cell2struct(empty_vals, fields, 1);
        sz_new = size(new_struct);
        dims = [index, sz_new(2:end)];
        global_struct = repmat(template, dims);
        global_struct(index, :, :) = new_struct;
        return;
    end

    fields_global = fieldnames(global_struct);
    fields_new = fieldnames(new_struct);

    % Add any fields that exist in global_struct but are missing in new_struct
    missing_in_new = setdiff(fields_global, fields_new);
    for i = 1:length(missing_in_new)
        [new_struct.(missing_in_new{i})] = deal([]);
    end

    % Add any fields that exist in new_struct but are missing in global_struct
    missing_in_global = setdiff(fields_new, fields_global);
    for i = 1:length(missing_in_global)
        [global_struct.(missing_in_global{i})] = deal([]);
    end

    % Order the new struct fields to match global_struct
    new_struct = orderfields(new_struct, global_struct);

    % Perform the assignment safely
    global_struct(index, :, :) = new_struct;
end
