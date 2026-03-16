function cleared = clear_pipeline_cache(config_struct)
% CLEAR_PIPELINE_CACHE  Remove pipeline-generated .mat cache files.
%
%   cleared = clear_pipeline_cache(config_struct)
%
%   When config_struct.clear_cache is true and the cache has not already
%   been cleared in this MATLAB session, deletes all pipeline-generated
%   .mat files from config_struct.dataloc so the pipeline recomputes
%   everything from scratch.  Also removes the per-patient checkpoint
%   directory if it was created by the pipeline (verified via sentinel).
%
%   A persistent variable ensures this runs at most once per session,
%   avoiding redundant clearing when execute_all_workflows calls
%   run_dwi_pipeline multiple times (once per DWI type).
%
%   Manually curated files (e.g., dwi_vectors_ea.mat) are protected.
%
%   Inputs:
%     config_struct - Pipeline configuration struct.  Must contain:
%                     .clear_cache (logical) - whether to clear
%                     .dataloc     (string)  - data directory path
%
%   Outputs:
%     cleared - true if cache files were actually deleted in this call,
%               false if skipped (already cleared or clear_cache is false)

    cleared = false;

    if ~isfield(config_struct, 'clear_cache') || ~config_struct.clear_cache
        return;
    end

    persistent cache_cleared_this_session;
    if ~isempty(cache_cleared_this_session) && cache_cleared_this_session
        return;
    end

    dataloc = config_struct.dataloc;
    cache_patterns = {
        fullfile(dataloc, 'dwi_vectors*.mat'), ...
        fullfile(dataloc, 'summary_metrics*.mat'), ...
        fullfile(dataloc, 'adc_vectors.mat')
    };
    % Protect manually curated files from cache clearing
    protected_files = {'dwi_vectors_ea.mat'};
    n_deleted = 0;
    for cp = 1:numel(cache_patterns)
        cached = dir(cache_patterns{cp});
        for cf = 1:numel(cached)
            if any(strcmpi(cached(cf).name, protected_files))
                fprintf('  🛡️ Skipping protected file: %s\n', cached(cf).name);
                continue;
            end
            delete(fullfile(cached(cf).folder, cached(cf).name));
            n_deleted = n_deleted + 1;
        end
    end
    % Remove per-patient checkpoint directory — only if it was
    % created by the pipeline (contains .pipeline_created sentinel).
    checkpoint_dir = fullfile(dataloc, 'processed_patients');
    if isfolder(checkpoint_dir) && exist(fullfile(checkpoint_dir, '.pipeline_created'), 'file')
        rmdir(checkpoint_dir, 's');
        fprintf('  🗑️ Removed per-patient checkpoint directory.\n');
    elseif isfolder(checkpoint_dir)
        fprintf('  🛡️ Skipping checkpoint directory (no pipeline sentinel): %s\n', checkpoint_dir);
    end
    if n_deleted > 0
        fprintf('  🗑️ Cleared %d cached .mat file(s) from %s\n', n_deleted, dataloc);
    else
        fprintf('  💡 No cached files found to clear.\n');
    end
    cache_cleared_this_session = true;
    cleared = true;
end
