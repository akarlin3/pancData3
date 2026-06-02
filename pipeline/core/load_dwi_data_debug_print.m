function load_dwi_data_debug_print(actual_var_names, pat_col_found, T_Pat_normalized, T, id_list_normalized)
    % LOAD_DWI_DATA_DEBUG_PRINT Diagnostic dump of spreadsheet vs folder patient IDs.
    %   Extracted verbatim from load_dwi_data.m Section 2. Side-effect only (fprintf);
    %   prints the matched patient-ID column, the first few normalized spreadsheet and
    %   folder IDs, and the runtime class of the binary-outcome columns so cell-format
    %   mismatches (e.g. Excel "Text" cells silently parsed as 0) are visible at startup.

% --- DEBUG: print spreadsheet vs folder patient IDs for matching diagnosis ---
fprintf('\n--- DEBUG: Spreadsheet columns: %s ---\n', strjoin(actual_var_names, ', '));
if ~isempty(pat_col_found)
    fprintf('--- DEBUG: Using column ''%s'' for patient matching ---\n', pat_col_found);
end
fprintf('--- DEBUG: Clinical spreadsheet Pat column (first 5) ---\n');
for dbg_i = 1:min(5, numel(T_Pat_normalized))
    fprintf('  Spreadsheet[%d]: "%s"\n', dbg_i, T_Pat_normalized{dbg_i});
end
% Surface the runtime type readtable inferred for the binary-outcome
% columns. When values look numeric in Excel but the cell format is
% "Text", readtable returns a cell-of-char and naive index access
% silently produces 0 for every row. This print makes the failure mode
% visible the first time the pipeline starts up.
for dbg_col = {'LF', 'LocalFailure', 'Immuno', 'IO'}
    cn = dbg_col{1};
    if has_table_var(T, cn)
        col = T.(cn);
        fprintf('--- DEBUG: T.%s class=%s, first 10 raw values ---\n', cn, class(col));
        for dbg_k = 1:min(10, numel(col))
            try
                v = col(dbg_k);
                if iscell(v); v = v{1}; end
                if isnumeric(v) || islogical(v)
                    fprintf('  T.%s(%d) = %g\n', cn, dbg_k, double(v));
                else
                    fprintf('  T.%s(%d) = "%s"\n', cn, dbg_k, char(string(v)));
                end
            catch dbg_err  %#ok<NASGU>
                fprintf('  T.%s(%d) = <unreadable>\n', cn, dbg_k);
            end
        end
    end
end

fprintf('--- DEBUG: Folder id_list (first 5) ---\n');
for dbg_i = 1:min(5, numel(id_list_normalized))
    fprintf('  Folder[%d]: "%s"\n', dbg_i, id_list_normalized{dbg_i});
end
fprintf('--- DEBUG: Total spreadsheet entries: %d, Total folder entries: %d ---\n\n', numel(T_Pat_normalized), numel(id_list_normalized));
end
