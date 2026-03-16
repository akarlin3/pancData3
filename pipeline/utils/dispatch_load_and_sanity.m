function [validated_data_gtvp, validated_data_gtvn, summary_metrics, abort] = dispatch_load_and_sanity(session)
% DISPATCH_LOAD_AND_SANITY  Execute or skip the load and sanity pipeline steps.
%
%   [validated_data_gtvp, validated_data_gtvn, summary_metrics, abort] = dispatch_load_and_sanity(session)
%
%   Runs the two "fatal" pipeline steps: data loading (DICOM conversion +
%   model fitting) and sanity checks (parameter validation).  Failure in
%   either step sets abort=true, signaling the caller to halt the pipeline.
%
%   When a step is not in session.steps_to_run, previously computed results
%   are loaded from disk via load_data_from_disk.
%
%   Inputs:
%     session - Session struct from prepare_pipeline_session containing:
%       steps_to_run, config_struct, pipeGUI, log_fid, master_diary_file,
%       current_name, current_dtype, dwi_vectors_file,
%       fallback_dwi_vectors_file, summary_metrics_file
%
%   Outputs:
%     validated_data_gtvp  - Validated GTVp data vectors
%     validated_data_gtvn  - Validated GTVn data vectors
%     summary_metrics      - Summary metrics struct
%     abort                - true if the pipeline should halt

    steps_to_run = session.steps_to_run;
    config_struct = session.config_struct;
    pipeGUI = session.pipeGUI;
    log_fid = session.log_fid;
    master_diary_file = session.master_diary_file;
    current_name = session.current_name;
    current_dtype = session.current_dtype;
    dwi_vectors_file = session.dwi_vectors_file;
    fallback_dwi_vectors_file = session.fallback_dwi_vectors_file;
    summary_metrics_file = session.summary_metrics_file;

    validated_data_gtvp = [];
    validated_data_gtvn = [];
    summary_metrics = [];
    abort = false;

    % =====================================================================
    % Step 2: Load DWI Data
    % =====================================================================
    % [ANALYTICAL RATIONALE]:
    % This is the computationally dominant step.  For each patient and
    % treatment fraction, load_dwi_data converts DICOM to NIFTI, loads
    % GTV masks, optionally applies DnCNN denoising, and fits IVIM/ADC
    % models voxel-by-voxel.
    if ismember('load', steps_to_run)
        if ~isempty(pipeGUI), pipeGUI.startStep('load'); end
        try
            fprintf('\xe2\x9a\x99\xef\xb8\x8f [2/5] [%s] Loading DWI data... \n', current_name);
            [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = load_dwi_data(config_struct);

            save(summary_metrics_file, 'summary_metrics');
            fprintf('      \xf0\x9f\x92\xbe Saved summary_metrics to %s\n', summary_metrics_file);

            fprintf('      \xe2\x9c\x85 Done: Successfully loaded data.\n');
            if ~isempty(pipeGUI), pipeGUI.completeStep('load', 'success'); end
            [warn_msg, warn_id] = lastwarn;
            lastwarn('');
            if ~isempty(warn_msg) && log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] During load_dwi_data: %s (id: %s)\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), warn_msg, warn_id);
            end
        catch ME
            if ~isempty(ME.stack)
                stack_line = ME.stack(1).line;
                stack_file = ME.stack(1).name;
                fprintf('\xe2\x9d\x8c FAILED at %s:%d.\n', stack_file, stack_line);
            else
                fprintf('\xe2\x9d\x8c FAILED.\n');
            end
            fprintf('\xe2\x9d\x8c Error during data loading: %s\n', ME.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [ERROR] Data loading failed: %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
                if ~isempty(ME.stack)
                    fprintf(log_fid, '         at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
                end
            end
            abort = true;
            return;
        end
    else
        % Load from disk (checkpoint recovery)
        if ~isempty(pipeGUI), pipeGUI.completeStep('load', 'skipped'); end
        fprintf('\xe2\x8f\xad\xef\xb8\x8f [2/5] [%s] Skipping Load Step. Loading from disk...\n', current_name);
        try
            [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = ...
                load_data_from_disk(dwi_vectors_file, fallback_dwi_vectors_file, ...
                                    summary_metrics_file, current_dtype, current_name);
            if isempty(summary_metrics)
                error('Data files not found. Please run "load" step first.');
            end
        catch ME
             fprintf('\xe2\x9d\x8c FAILED.\n');
             fprintf('\xe2\x9d\x8c Error loading data from disk: %s\n', ME.message);
             if log_fid > 0
                 fprintf(log_fid, '[%s] [ERROR] Loading data from disk failed: %s\n', ...
                     datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
             end
             abort = true;
             return;
        end
    end
    diary(master_diary_file);

    % =====================================================================
    % Step 3: Sanity Checks
    % =====================================================================
    % [ANALYTICAL RATIONALE]:
    % Validates that fitted diffusion parameters are physically plausible.
    % Failure halts the pipeline to prevent scientifically invalid results.
    if ismember('sanity', steps_to_run)
        if ~isempty(pipeGUI), pipeGUI.startStep('sanity'); end
        try
            fprintf('\xe2\x9a\x99\xef\xb8\x8f [3/5] [%s] Running sanity checks...\n', current_name);
            [is_valid, validation_msg, validated_data_gtvp, validated_data_gtvn] = sanity_checks(data_vectors_gtvp, data_vectors_gtvn, summary_metrics, config_struct);

            if ~is_valid
                error('Sanity checks failed: %s', validation_msg);
            end
            fprintf('      \xe2\x9c\x85 Passed.\n');
            if ~isempty(pipeGUI), pipeGUI.completeStep('sanity', 'success'); end

            sanity_results_file = fullfile(config_struct.output_folder, sprintf('sanity_checks_results_%s.txt', current_name));
            fid = fopen(sanity_results_file, 'w');
            if fid < 0
                warning('run_dwi_pipeline:fileWriteFailed', 'Cannot write %s', sanity_results_file);
            else
                fprintf(fid, 'is_valid: %d\nvalidation_msg: %s\n', is_valid, validation_msg);
                fclose(fid);
            end
            fprintf('      \xf0\x9f\x92\xbe Saved sanity check results to %s\n', sanity_results_file);
            [warn_msg, warn_id] = lastwarn;
            lastwarn('');
            if ~isempty(warn_msg) && log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] During sanity_checks: %s (id: %s)\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), warn_msg, warn_id);
            end
        catch ME
            fprintf('\xe2\x9d\x8c FAILED.\n');
            fprintf('\xe2\x9d\x8c Pipeline halted due to sanity check failure: %s\n', ME.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [ERROR] Sanity checks failed: %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
            end
            abort = true;
            return;
        end
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('sanity', 'skipped'); end
        fprintf('\xe2\x8f\xad\xef\xb8\x8f [3/5] [%s] Skipping Sanity Checks.\n', current_name);
        if exist('data_vectors_gtvp', 'var')
            validated_data_gtvp = data_vectors_gtvp;
            validated_data_gtvn = data_vectors_gtvn;
        else
            try
                [validated_data_gtvp, validated_data_gtvn, sm_tmp] = ...
                    load_data_from_disk(dwi_vectors_file, fallback_dwi_vectors_file, ...
                                        summary_metrics_file, current_dtype, current_name);
                if ~exist('summary_metrics', 'var') || isempty(summary_metrics)
                    if ~isempty(sm_tmp)
                        summary_metrics = sm_tmp;
                    end
                end
            catch ME_load
                fprintf('\xe2\x9d\x8c Error loading data vectors from disk: %s\n', ME_load.message);
                if log_fid > 0
                    fprintf(log_fid, '[%s] [ERROR] Loading data vectors for visualize failed: %s\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME_load.message);
                end
                abort = true;
                return;
            end
        end
    end
    diary(master_diary_file);
    lastwarn('');
end
