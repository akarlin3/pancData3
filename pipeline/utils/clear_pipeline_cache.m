function cleared = clear_pipeline_cache(config_struct)
% CLEAR_PIPELINE_CACHE  Remove pipeline-generated derived .mat cache files.
%
%   cleared = clear_pipeline_cache(config_struct)
%
%   When config_struct.clear_cache is true and the cache has not already
%   been cleared in this MATLAB session, deletes derived .mat files from
%   config_struct.dataloc so the pipeline recomputes downstream metrics
%   from scratch.  Also removes the per-patient checkpoint directory if
%   it was created by the pipeline (verified via sentinel).
%
%   A persistent variable ensures this runs at most once per session,
%   avoiding redundant clearing when execute_all_workflows calls
%   run_dwi_pipeline multiple times (once per DWI type).
%
%   The pipeline's own voxel cache lives in pipeline_voxels*.mat.  That
%   is the expensive voxel-extraction checkpoint (hours of DICOM
%   conversion + model fitting) and is cleared here when clear_cache is
%   true so the next run regenerates it.
%
%   dwi_vectors*.mat is reserved for the user's own curated files and is
%   NEVER touched: never read, never written, never deleted by the
%   pipeline.  The defensive guard below skips any file whose name
%   starts with "dwi_vectors" as a belt-and-braces safety in case a
%   future glob edit accidentally matches one.
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
        fprintf(['  💡 clear_pipeline_cache: already ran this session — ' ...
                 'skipping (persistent short-circuit). Call ' ...
                 '`clear clear_pipeline_cache` to force a re-run.\n']);
        return;
    end

    dataloc = config_struct.dataloc;
    % Pipeline-owned caches: pipeline_voxels*.mat (upstream voxel
    % checkpoint), summary_metrics*.mat (derived), adc_vectors.mat.
    % dwi_vectors*.mat is intentionally excluded — see header docstring.
    cache_patterns = {
        fullfile(dataloc, 'pipeline_voxels*.mat'), ...
        fullfile(dataloc, 'summary_metrics*.mat'), ...
        fullfile(dataloc, 'adc_vectors.mat')
    };
    n_deleted = 0;
    for cp = 1:numel(cache_patterns)
        cached = dir(cache_patterns{cp});
        for cf = 1:numel(cached)
            % Defensive guard: if the glob ever accidentally matches a
            % dwi_vectors*.mat file, skip it rather than delete.
            if startsWith(lower(cached(cf).name), 'dwi_vectors')
                fprintf('  🛡️ Skipping protected file: %s\n', cached(cf).name);
                continue;
            end
            delete(fullfile(cached(cf).folder, cached(cf).name));
            n_deleted = n_deleted + 1;
        end
    end
    % Remove per-patient checkpoint directory — only if it was
    % created by the pipeline (contains .pipeline_created sentinel).
    % Legacy dirs that predate the sentinel convention get healed first
    % via backfill_checkpoint_sentinel, which writes the sentinel iff
    % the contents look pipeline-owned (patient_NNN_*.mat with no
    % foreign files).  This lets a single clear_cache:true run sweep
    % a legacy processed_patients/ instead of forcing the user to
    % rm -rf by hand.
    checkpoint_dir = fullfile(dataloc, 'processed_patients');
    if isfolder(checkpoint_dir)
        backfill_checkpoint_sentinel(checkpoint_dir);
    end
    if isfolder(checkpoint_dir) && exist(fullfile(checkpoint_dir, '.pipeline_created'), 'file')
        rmdir(checkpoint_dir, 's');
        fprintf('  🗑️ Removed per-patient checkpoint directory.\n');
    elseif isfolder(checkpoint_dir)
        fprintf('  🛡️ Skipping checkpoint directory (no pipeline sentinel; contents do not match pipeline naming): %s\n', checkpoint_dir);
    end
    if n_deleted > 0
        fprintf('  🗑️ Cleared %d cached .mat file(s) from %s\n', n_deleted, dataloc);
    else
        fprintf('  💡 No cached files found to clear.\n');
    end
    cache_cleared_this_session = true;
    cleared = true;
end
