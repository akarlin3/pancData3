function [success, err_msg] = execute_pipeline_step(step_name, step_fn, pipeGUI, log_fid, master_diary_file, output_folder, current_name)
% EXECUTE_PIPELINE_STEP  Generic executor for non-fatal pipeline steps.
%
%   [success, err_msg] = execute_pipeline_step(step_name, step_fn, pipeGUI, log_fid, master_diary_file, output_folder, current_name)
%
%   Wraps the boilerplate try-catch-GUI-warning-diary-logging pattern used
%   by non-fatal pipeline steps in run_dwi_pipeline.  The caller provides a
%   function handle (step_fn) that contains ALL step-specific logic
%   (console output, function calls, saves).  This wrapper handles:
%     - GUI startStep / completeStep
%     - try-catch with non-fatal error handling
%     - Warning capture to error log
%     - Diary restart after the step completes (success or failure)
%     - lastwarn('') reset between steps
%
%   Inputs:
%     step_name         - Pipeline step key for GUI (e.g., 'metrics_longitudinal')
%     step_fn           - Function handle (@() ...) containing the step logic.
%                         Called inside a try block; any error triggers non-fatal
%                         handling.
%     pipeGUI           - PipelineProgressGUI object (may be empty)
%     log_fid           - File handle for error.log (may be -1 if not open)
%     master_diary_file - Path to the master diary file for diary restart
%     output_folder     - Type output folder path (currently unused but
%                         available for future step-specific logging)
%     current_name      - Human-readable DWI type name (currently unused but
%                         available for future step-specific logging)
%
%   Outputs:
%     success   - true if step_fn completed without error, false otherwise
%     err_msg   - Empty string on success, error message string on failure

    err_msg = '';

    if ~isempty(pipeGUI), pipeGUI.startStep(step_name); end
    try
        step_fn();
        if ~isempty(pipeGUI), pipeGUI.completeStep(step_name, 'success'); end
        [warn_msg, warn_id] = lastwarn;  % read current warning
        lastwarn('');                       % then clear
        if ~isempty(warn_msg) && log_fid > 0
            fprintf(log_fid, '[%s] [WARNING] During %s: %s (id: %s)\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), step_name, warn_msg, warn_id);
        end
        success = true;
    catch ME
        fprintf('\u26a0\ufe0f FAILED (Non-Fatal).\n');
        fprintf('\u26a0\ufe0f Error during %s: %s\n', step_name, ME.message);
        if ~isempty(pipeGUI), pipeGUI.completeStep(step_name, 'warning'); end
        if log_fid > 0
            fprintf(log_fid, '[%s] [ERROR] %s failed: %s\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), step_name, ME.message);
            if ~isempty(ME.stack)
                fprintf(log_fid, '         at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
            end
        end
        success = false;
        err_msg = ME.message;
    end
    % Restart master diary after step completes (each core module may
    % override the diary with its own file)
    diary(master_diary_file);
    lastwarn('');  % reset warning tracker between steps
end
