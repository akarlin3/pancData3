function master_folder = setup_output_folders(pipeline_dir, master_output_folder)
% SETUP_OUTPUT_FOLDERS  Create or reuse the master pipeline output folder.
%
%   master_folder = setup_output_folders(pipeline_dir, master_output_folder)
%
%   Handles the two output folder modes:
%     1. Explicit folder (master_output_folder non-empty): reuses the given
%        path, creating it if needed.  Used by execute_all_workflows to share
%        a single timestamped folder across Standard/dnCNN/IVIMnet runs.
%     2. Standalone run (master_output_folder empty): creates a fresh
%        timestamped folder (saved_files_YYYYMMDD_HHMMSS) with a
%        .pipeline_created sentinel for safe cleanup.
%
%   Inputs:
%     pipeline_dir         - Absolute path to the pipeline/ directory
%     master_output_folder - Explicit master folder path (string), or ''
%                            to auto-create a timestamped folder
%
%   Outputs:
%     master_folder - Resolved master output folder path

    if ~isempty(master_output_folder)
        % User/execute_all_workflows explicitly provided a folder
        if ~exist(master_output_folder, 'dir'), mkdir(master_output_folder); end
        master_folder = master_output_folder;
        fprintf('      📁 Using explicitly provided master output folder: %s\n', master_output_folder);
    else
        % Standalone run — always create a fresh folder
        timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
        master_folder = fullfile(pipeline_dir, '..', sprintf('saved_files_%s', timestamp_str));
        if ~exist(master_folder, 'dir'), mkdir(master_folder); end
        % Write provenance sentinel so cleanup tools know this directory
        % was created by the pipeline and is safe to delete.
        sent_fid = fopen(fullfile(master_folder, '.pipeline_created'), 'w');
        if sent_fid > 0
            fprintf(sent_fid, 'Created by run_dwi_pipeline at %s\n', timestamp_str);
            fclose(sent_fid);
        end
        fprintf('      📁 Created NEW master output folder: %s\n', master_folder);
    end
end
