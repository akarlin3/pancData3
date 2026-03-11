function [pat_normalized, id_list_normalized] = normalize_patient_ids(T_Pat, id_list)
% NORMALIZE_PATIENT_IDS — Normalize patient IDs for cross-source matching.
%
%   Converts clinical spreadsheet patient IDs (T.Pat) and file-system folder
%   names to a common normalized form by:
%     1. Converting categorical/char arrays to cellstr (Octave-compatible)
%     2. Stripping Excel-embedded single quotes
%     3. Replacing underscores with hyphens (so 'P_01' matches 'P-01')
%
%   Parameters
%   ----------
%   T_Pat : categorical, cellstr, char, or numeric
%       The T.Pat column from the clinical spreadsheet.
%   id_list : cellstr
%       Patient folder names from the file system.
%
%   Returns
%   -------
%   pat_normalized : cellstr
%       Normalized spreadsheet patient IDs.
%   id_list_normalized : cellstr
%       Normalized folder patient IDs.
%
%   See also: load_dwi_data, metrics_baseline

    % --- Convert T_Pat to cellstr (Octave compatibility) ---
    if exist('OCTAVE_VERSION', 'builtin')
        % iscategorical may be missing or mocked; T.Pat might be char array
        if ischar(T_Pat)
            if size(T_Pat, 1) > 1
                T_Pat_cell = {};
                for i_pat_row = 1:size(T_Pat, 1)
                    T_Pat_cell{i_pat_row} = strtrim(T_Pat(i_pat_row, :));
                end
            else
                T_Pat_cell = {T_Pat};
            end
        elseif isnumeric(T_Pat)
            T_Pat_cell = {};
        elseif iscell(T_Pat)
            T_Pat_cell = T_Pat;
        else
            T_Pat_cell = {T_Pat};
        end
    else
        if iscategorical(T_Pat)
            T_Pat_cell = cellstr(T_Pat);
        else
            T_Pat_cell = T_Pat;
        end
    end

    % --- Normalize spreadsheet patient IDs ---
    if isempty(T_Pat_cell)
        pat_normalized = {};
    else
        try
            % Strip leading/trailing single quotes that Excel may embed
            T_Pat_cell = strrep(T_Pat_cell, '''', '');
            pat_normalized = strrep(T_Pat_cell, '_', '-');
        catch
            pat_normalized = {};
        end
    end

    % --- Normalize folder names ---
    if isempty(id_list)
        id_list_normalized = {};
    else
        try
            id_list_normalized = strrep(id_list, '_', '-');
        catch
            id_list_normalized = {};
        end
    end
end
