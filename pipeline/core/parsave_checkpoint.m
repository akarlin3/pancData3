function parsave_checkpoint(fname, data, lock_file)
    % PARSAVE_CHECKPOINT Parallel-safe checkpoint writer with lock-file protocol.
    %   Extracted local helper from load_dwi_data.m.
    % In a parfor loop, multiple workers may finish near-simultaneously.
    % The .lock file acts as a write-ahead indicator: if the pipeline
    % crashes between lock creation and .mat completion, the recovery
    % logic in the main function detects the orphaned .lock and re-processes
    % the patient. This ensures data integrity even after unclean shutdowns.
    % Create lock sentinel BEFORE writing to prevent race conditions.
    % The lock is removed only after the .mat write completes successfully.
    if nargin >= 3 && ~isempty(lock_file)
        fid = fopen(lock_file, 'w');
        if fid > 0, fclose(fid); end
    end
    save(fname, '-struct', 'data');
    if nargin >= 3 && ~isempty(lock_file) && exist(lock_file, 'file')
        delete(lock_file);
    end
end
