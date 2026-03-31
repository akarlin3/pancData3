function run_dwi_pipeline(config_path, steps_to_run, master_output_folder)
% RUN_DWI_PIPELINE Master orchestrator function for the DWI analysis pipeline
%
% [ORCHESTRATOR PATTERN]:
% This function acts as the central controller for the DWI analysis workflow.
% Instead of monolithic scripts, the pipeline is broken down into modular steps
% (Load -> Sanity -> Visualize -> Metrics). The orchestrator manages:
%   1. Data Flow: Passing outputs from one module as inputs to the next.
%   2. Error Handling: Catching exceptions and halting execution gracefully.
%   3. Environment: Setting up paths and verifying dependencies dynamically.
% Usage:
%   run_dwi_pipeline('config.json');
%   run_dwi_pipeline('config.json', {'test', 'load', 'sanity', 'visualize', 'metrics_baseline', 'compare_cores', 'cross_pipeline_dice', 'core_failure_rates', 'metrics_longitudinal', 'metrics_dosimetry', 'core_method_outcomes', 'metrics_stats_comparisons', 'metrics_stats_predictive', 'metrics_survival'});
%   run_dwi_pipeline('config.json', {'load'}, 'path/to/my_output_folder');
%
% This function sequentially calls the module scripts:
%   1. load_dwi_data
%   2. sanity_checks
%   3. visualize_results
%   4. metrics (baseline, longitudinal, dosimetry, stats_comparisons, stats_predictive, survival)
%
% It explicitly passes data between modules to avoid workspace pollution
% and includes error handling to halt execution if checks fail.
%
% [ANALYTICAL DESIGN RATIONALE]:
% The pipeline follows a strict sequential dependency chain that mirrors the
% scientific workflow for DWI-based treatment response assessment in pancreatic
% cancer radiotherapy:
%
%   Load: Convert raw DICOM DWI images -> fit voxel-wise IVIM/ADC models.
%         This is the most computationally expensive step (hours for large cohorts)
%         because each voxel in each b-value image must be fit to either the
%         mono-exponential ADC model (S = S0*exp(-b*ADC)) or the bi-exponential
%         IVIM model (S = S0*((1-f)*exp(-b*D) + f*exp(-b*D*))). Checkpointing
%         ensures that a crash after fitting patient 40/60 does not lose results.
%
%   Sanity: Validate that fitted parameters are physically meaningful. IVIM D
%           values should be ~1e-3 mm^2/s (tissue diffusion), f should be 0-1
%           (perfusion fraction), D* should be ~10x D (capillary pseudo-diffusion).
%           Invalid fits indicate convergence failure, corrupt DICOM data, or
%           misregistered GTV masks.
%
%   Metrics: The analytical core. Baseline metrics establish pre-treatment
%            diffusion characteristics of the tumor. Longitudinal metrics track
%            how these change across fractions (Fx1-Fx5) -- early diffusion
%            increases may indicate tumor cell kill. Dosimetry correlates RT
%            dose heterogeneity with diffusion sub-volumes. Statistical and
%            survival analyses test whether diffusion changes predict outcomes.
%
%   Visualize: Generates parameter maps and distributions AFTER predictive
%              modeling so that risk stratification overlays are available.
%
% The step ordering matters: metrics_baseline must precede longitudinal/dosimetry
% because percent-change calculations require baseline values. Dosimetry must
% precede stats because dose-response sub-volume features (D95, V50 within
% diffusion-defined tumor sub-regions) are candidate predictors. Survival
% analysis runs last because it uses risk scores from predictive modeling.
    % --- Initialization Block ---
    pipeline_dir = fileparts(mfilename('fullpath'));

    if nargin < 1
        config_path = fullfile(pipeline_dir, '..', 'config.json');
    end

    if nargin < 2
        steps_to_run = {'test', 'load', 'sanity', 'visualize', 'metrics_baseline', 'metrics_longitudinal', 'metrics_dosimetry', 'metrics_stats_comparisons', 'metrics_stats_predictive', 'metrics_survival'};
    end

    if nargin < 3
        master_output_folder = '';
    end

    % Pre-flight tests (once per session via dedicated caching function).
    has_tests_passed = check_tests_cached(pipeline_dir, config_path, steps_to_run, master_output_folder);

    if ~has_tests_passed
        error('run_dwi_pipeline:PreFlightFailed', ...
            '%s Pre-flight initialization failed. Cannot proceed with pipeline execution. Review test output above for details.', ...
            safe_icon('fail'));
    end

    config_path = resolve_config_path(pipeline_dir, config_path);

    % If 'test' was the only requested step, stop here.
    other_steps = setdiff(steps_to_run, {'test'});
    if isempty(other_steps)
        fprintf('%s Test-only run complete. No pipeline steps to execute.\n', safe_icon('ok'));
        return;
    end

    % ----------------------------
    % Session Setup
    % ----------------------------
    % Delegate config parsing, DWI type resolution, output folder setup,
    % diary/error log initialization, GUI, file paths, and cache clearing
    % to prepare_pipeline_session.  The session struct includes a cleanup
    % function handle that encapsulates all resource teardown logic, so
    % the orchestrator does not depend on internal session field names.
    session = prepare_pipeline_session(pipeline_dir, config_path, master_output_folder, steps_to_run);

    % ----------------------------
    % Session Cleanup Setup
    % ----------------------------
    % The session cleanup logic uses a structured approach with two fallback levels:
    %
    % 1. PREFERRED: session.cleanup function handle (modern structured approach)
    %    - prepare_pipeline_session() provides a dedicated cleanup function handle that
    %      encapsulates ALL resource teardown logic (diary, GUI, log files, figure settings)
    %    - This ensures the orchestrator doesn't need to know about internal session structure
    %    - Provides better error isolation and cleaner separation of concerns
    %
    % 2. FALLBACK: fallback_cleanup() local function (backwards compatibility)
    %    - Used only when prepare_pipeline_session() doesn't provide session.cleanup
    %    - Maintains compatibility with older session preparation code
    %    - Directly accesses session fields (log_fid, prev_fig_vis, pipeGUI) for cleanup
    %    - Defined as a local function within this file to guarantee availability
    %
    % onCleanup guard: MUST be created before the abort check so that early returns
    % still trigger all resource cleanup. The cleanup logic handles failures gracefully.
    
    if isfield(session, 'cleanup') && isa(session.cleanup, 'function_handle')
        % Modern structured cleanup: use the session-provided cleanup handle
        cleanup_guard = onCleanup(session.cleanup); %#ok<NASGU>
    else
        % Backwards compatibility: fallback to field-based cleanup using local function.
        % The fallback_cleanup local function is defined at the bottom of this file,
        % so it is always available regardless of the MATLAB path.
        warning('run_dwi_pipeline:UsingFallbackCleanup', ...
            'Session does not provide structured cleanup handle. Using fallback cleanup method. Consider updating prepare_pipeline_session() to provide session.cleanup function handle.');
        % Capture session by value for the onCleanup callback; wrap in try-catch
        % so that errors in cleanup are logged to stderr rather than silently
        % swallowed by onCleanup's error suppression.
        session_snapshot = session;
        cleanup_guard = onCleanup(@() fallback_cleanup_safe(session_snapshot)); %#ok<NASGU>
    end

    if session.abort
        return;
    end

    % ----------------------------
    % Load & Sanity Steps
    % ----------------------------
    % These are "fatal" steps: failure halts the pipeline entirely.
    [validated_data_gtvp, validated_data_gtvn, summary_metrics, abort] = ...
        dispatch_load_and_sanity(session);
    if abort, return; end

    % ----------------------------
    % Analysis Steps
    % ----------------------------
    % Metrics baseline through survival, plus visualization.
    dispatch_pipeline_steps(session, validated_data_gtvp, validated_data_gtvn, summary_metrics);

    % ----------------------------
    % Completion
    % ----------------------------
    fprintf('=======================================================\n');
    fprintf('%s Pipeline Execution Complete for parameter %s\n', safe_icon('done'), session.current_name);
    fprintf('=======================================================\n');

    if isfield(session, 'log_fid') && session.log_fid > 0
        fprintf(session.log_fid, '[%s] ===== Pipeline run completed (type: %s) =====\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), session.current_name);
    end

    diary off;
end

%% ===== Utility functions (remain in orchestrator) =====

function fallback_cleanup_safe(session)
%FALLBACK_CLEANUP_SAFE  Wrapper that calls fallback_cleanup with stderr logging.
%
%   onCleanup suppresses errors in callbacks, so this wrapper catches any
%   exception from fallback_cleanup and writes the error details to stderr.
%   This ensures resource-leak failures are visible to the user rather than
%   being silently swallowed.
    try
        fallback_cleanup(session);
    catch ME
        % Log to stderr so the user is aware cleanup failed, even though
        % onCleanup would have suppressed the error.
        try
            fprintf(2, '[ERROR] fallback_cleanup failed during onCleanup:\n');
            fprintf(2, '  Identifier: %s\n', ME.identifier);
            fprintf(2, '  Message:    %s\n', ME.message);
            for k = 1:numel(ME.stack)
                fprintf(2, '  In %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
        catch
            % Last resort: if even stderr logging fails, there is nothing
            % more we can do.
        end
    end
end

function fallback_cleanup(session)
%FALLBACK_CLEANUP  Perform resource cleanup using session fields (backwards compatibility).
%   
%   This function is used only when prepare_pipeline_session() does not provide 
%   a structured session.cleanup function handle. It directly accesses known 
%   session fields to perform cleanup operations.
%
%   BACKWARDS COMPATIBILITY NOTICE:
%   This fallback approach is maintained for compatibility with older session 
%   preparation code. New implementations should migrate to providing a structured
%   session.cleanup function handle from prepare_pipeline_session().
%
%   Each cleanup operation is wrapped in try-catch to ensure partial cleanup 
%   succeeds even if individual operations fail. This prevents cascading failures
%   during resource teardown.
%
%   Cleaned resources:
%   - MATLAB diary logging
%   - Figure visibility settings (session.prev_fig_vis)
%   - Log file descriptor (session.log_fid) 
%   - GUI window handle (session.pipeGUI)

    try 
        diary('off'); 
    catch ME
        warning('fallback_cleanup:DiaryFailed', 'Failed to turn off diary: %s', ME.message);
    end
    
    if isfield(session, 'prev_fig_vis')
        try 
            set(0, 'DefaultFigureVisible', session.prev_fig_vis); 
        catch ME
            warning('fallback_cleanup:FigureVisibilityFailed', 'Failed to restore figure visibility: %s', ME.message);
        end
    end
    
    if isfield(session, 'log_fid')
        try
            safe_fclose_log(session.log_fid);
        catch ME
            warning('fallback_cleanup:LogCloseFailed', 'Failed to close log file: %s', ME.message);
        end
    end
    
    if isfield(session, 'pipeGUI')
        try
            closeIfValid(session.pipeGUI);
        catch ME
            warning('fallback_cleanup:GUICloseFailed', 'Failed to close GUI: %s', ME.message);
        end
    end
end

function closeIfValid(gui)
%CLOSEIFVALID  Close a PipelineProgressGUI if it is still valid.
    if ~isempty(gui)
        try gui.close(); catch, end
    end
end

function safe_fclose_log(fid)
%SAFE_FCLOSE_LOG  Close a file handle only if it is valid.
    if fid > 0
        try fclose(fid); catch, end
    end
end

function resolved = resolve_config_path(pipeline_dir, config_path)
%RESOLVE_CONFIG_PATH  Return an absolute config path, resolving relative paths.
    if ~isempty(config_path) && exist(config_path, 'file')
        resolved = config_path;
    else
        candidate = fullfile(pipeline_dir, '..', config_path);
        if exist(candidate, 'file')
            resolved = candidate;
        else
            resolved = config_path;  % pass through; downstream will error
        end
    end
end

function has_passed = check_tests_cached(pipeline_dir, config_path, steps_to_run, master_output_folder)
%CHECK_TESTS_CACHED  Run pre-flight tests once per session, caching the result.
%
%   Encapsulates the test-caching logic with its own persistent variables so
%   that the orchestrator does not manage persistent state directly. Tests are
%   re-run only if the cache is empty or older than 1 hour.
    persistent cached_passed;
    persistent cached_timestamp;

    cache_valid = false;
    if ~isempty(cached_passed) && ~isempty(cached_timestamp)
        elapsed_hours = (now - cached_timestamp) * 24;
        if cached_passed && elapsed_hours < 1
            cache_valid = true;
        end
    end

    if cache_valid
        fprintf('%s Pre-flight tests already passed this session (%.1f min ago). Skipping.\n', ...
            safe_icon('ok'), (now - cached_timestamp) * 24 * 60);
        has_passed = true;
        return;
    end

    % Delegate to initialize_pipeline without passing persistent state.
    % initialize_pipeline handles path setup, running tests if 'test' is in
    % steps_to_run, and returns the (possibly updated) config path.
    try
        config_path = initialize_pipeline(pipeline_dir, config_path, steps_to_run, master_output_folder, cached_passed, cached_timestamp);
        cached_passed = true;
        cached_timestamp = now;
        has_passed = true;
    catch ME
        cached_passed = false;
        cached_timestamp = [];
        fprintf('%s Pre-flight initialization failed: %s\n', safe_icon('fail'), ME.message);
        has_passed = false;
    end
end

function icon = safe_icon(name)
%SAFE_ICON  Return a Unicode icon string if the console supports UTF-8,
%   otherwise return an ASCII-safe fallback.
    is_utf8 = false;
    try
        charset = feature('DefaultCharacterSet');
        if strcmpi(charset, 'UTF-8')
            is_utf8 = true;
        end
    catch
    end

    switch lower(name)
        case 'ok'
            if is_utf8, icon = char([10004]); else, icon = '[OK]'; end
        case 'done'
            if is_utf8, icon = char([9989]); else, icon = '[DONE]'; end
        case 'fail'
            if is_utf8, icon = char([10060]); else, icon = '[FAIL]'; end
        otherwise
            icon = '';
    end
end