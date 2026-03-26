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

    % Notify the progress GUI that this step has started (updates visual indicator)
    if ~isempty(pipeGUI), pipeGUI.startStep(step_name); end
    try
        % Execute the step-specific logic provided by the caller.
        % step_fn is a zero-argument function handle that encapsulates
        % all module calls, console output, and file saves for this step.
        step_fn();

        % Mark step as successful in the GUI
        if ~isempty(pipeGUI), pipeGUI.completeStep(step_name, 'success'); end

        % Capture any non-fatal warnings emitted during the step and
        % write them to the error log for post-run debugging. Then clear
        % the warning state so the next step starts clean.
        [warn_msg, warn_id] = lastwarn;
        lastwarn('');
        if ~isempty(warn_msg) && log_fid > 0
            fprintf(log_fid, '[%s] [WARNING] During %s: %s (id: %s)\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), step_name, warn_msg, warn_id);
        end
        success = true;
    catch ME
        % Non-fatal error handling: log the error but do NOT rethrow.
        % This allows the pipeline to continue with subsequent steps
        % (e.g., survival analysis can still run even if dosimetry fails).
        fprintf('\xe2\x9a\xa0\xef\xb8\x8f FAILED (Non-Fatal).\n');
        fprintf('\xe2\x9a\xa0\xef\xb8\x8f Error during %s: %s\n', step_name, ME.message);
        if ~isempty(pipeGUI), pipeGUI.completeStep(step_name, 'warning'); end
        if log_fid > 0
            fprintf(log_fid, '[%s] [ERROR] %s failed: %s\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), step_name, ME.message);
            % Log the first stack frame for quick debugging without a full traceback
            if ~isempty(ME.stack)
                fprintf(log_fid, '         at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
            end
        end
        success = false;
        err_msg = ME.message;
    end
    % Restart master diary after step completes. Each core module opens its
    % own diary file (e.g., sanity_checks_output.txt), which overrides the
    % orchestrator's diary. Re-opening the master diary here resumes capture
    % of inter-step console output in the orchestrator's log file.
    diary(master_diary_file);
    lastwarn('');  % reset warning tracker between steps
end
