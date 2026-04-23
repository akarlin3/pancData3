function backfilled = backfill_checkpoint_sentinel(checkpoint_dir)
% BACKFILL_CHECKPOINT_SENTINEL  Add .pipeline_created to legacy checkpoint dirs.
%
%   backfilled = backfill_checkpoint_sentinel(checkpoint_dir)
%
%   Older pipeline versions created processed_patients/ without the
%   .pipeline_created sentinel that clear_pipeline_cache.m looks for
%   before sweeping the directory.  Without the sentinel, those legacy
%   folders survive every clear_cache:true call, leaving stale per-
%   patient checkpoints in place.
%
%   This helper backfills the sentinel for an existing checkpoint
%   directory IFF its contents look pipeline-owned: contains at least
%   one patient_NNN_<id>.mat file AND every entry matches the
%   pipeline's own naming convention (patient_*.mat, patient_*.lock,
%   or .pipeline_created itself).  This conservative content check
%   prevents the pipeline from retroactively claiming an unrelated
%   user folder that happens to share the name.
%
%   No-ops cleanly when:
%     - The directory does not exist.
%     - The sentinel already exists.
%     - The directory is empty.
%     - Any file does not match the pipeline naming convention.
%
%   Inputs:
%     checkpoint_dir - Absolute path to processed_patients (or equivalent)
%
%   Outputs:
%     backfilled - true if a sentinel was newly written, false otherwise.

    backfilled = false;

    if ~isfolder(checkpoint_dir)
        return;
    end
    sentinel_path = fullfile(checkpoint_dir, '.pipeline_created');
    if exist(sentinel_path, 'file')
        return;
    end

    cp_entries = dir(checkpoint_dir);
    cp_entries = cp_entries(~[cp_entries.isdir]);
    if isempty(cp_entries)
        return;
    end

    names_lower = lower({cp_entries.name});
    has_patient_mat = any(~cellfun('isempty', ...
        regexp(names_lower, '^patient_\d+_.*\.mat$', 'once')));
    all_recognized = all(cellfun(@(nm) ...
        ~isempty(regexp(nm, '^patient_\d+_.*\.(mat|lock)$', 'once')) || ...
        strcmp(nm, '.pipeline_created'), names_lower));

    if ~(has_patient_mat && all_recognized)
        return;
    end

    sent_fid = fopen(sentinel_path, 'w');
    if sent_fid > 0
        fprintf(sent_fid, 'Backfilled by backfill_checkpoint_sentinel\n');
        fclose(sent_fid);
        backfilled = true;
        fprintf('  💡 Backfilled .pipeline_created sentinel for existing %s\n', ...
            checkpoint_dir);
    end
end
