function session = prepare_pipeline_session(pipeline_dir, config_path, master_output_folder, steps_to_run)
% PREPARE_PIPELINE_SESSION  Set up all session state for a pipeline run.
%
%   session = prepare_pipeline_session(pipeline_dir, config_path, master_output_folder, steps_to_run)
%
%   Performs all setup between initialization and step dispatch:
%     1. Parses config.json via parse_config
%     2. Resolves DWI type to a scalar index and human-readable name
%     3. Sets up master and type-specific output folders
%     4. Opens the master diary and error log
%     5. Conditionally injects compare_cores into steps_to_run
%     6. Initializes the pipeline progress GUI (if display available)
%     7. Builds type-specific file paths for pipeline artifacts
%     8. Clears cached files if requested by config
%
%   Inputs:
%     pipeline_dir         - Absolute path to the pipeline/ directory
%     config_path          - Resolved path to config.json
%     master_output_folder - Parent output folder (empty for auto-creation)
%     steps_to_run         - Cell array of pipeline step names
%
%   Outputs:
%     session - Struct with fields:
%       config_struct             - Parsed configuration struct
%       steps_to_run              - Possibly modified step list
%       current_dtype             - Scalar DWI type index (1/2/3)
%       current_name              - Human-readable DWI type name
%       master_output_folder      - Resolved master output folder path
%       type_output_folder        - DWI-type-specific output subfolder
%       master_diary_file         - Path to the master diary log
%       log_fid                   - File handle for error.log (-1 if failed)
%       pipeGUI                   - PipelineProgressGUI object (may be [])
%       prev_fig_vis              - Previous DefaultFigureVisible setting
%       voxel_cache_file          - Type-specific voxel-cache .mat path
%       voxel_cache_fallback_file - Un-typed voxel-cache .mat fallback path
%       summary_metrics_file      - Type-specific summary_metrics .mat path
%       results_file              - Calculated results .mat path
%       baseline_results_file     - Baseline results .mat path
%       dosimetry_results_file    - Dosimetry results .mat path
%       predictive_results_file   - Predictive results .mat path
%
%   Errors:
%     Returns session with session.abort = true and prints error message
%     if config parsing fails. Caller should check session.abort and return.

    session = struct();
    session.abort = false;
    session.log_fid = -1;
    session.pipeGUI = [];

    % Guard ownership pattern: This function creates the error-log file
    % handle (log_fid) and modifies DefaultFigureVisible.  Under normal
    % operation the *caller* owns cleanup via onCleanup tied to the
    % returned session struct.  However, if this function throws before
    % returning, the caller never receives the struct, so its onCleanup
    % guards never fire.  The try-catch below ensures that on any error
    % after the global-state change we (a) close log_fid if it was opened
    % and (b) restore DefaultFigureVisible, then re-throw so the caller
    % still sees the failure.

    % Suppress figure windows — all plots are saved to disk via saveas()
    session.prev_fig_vis = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');

    fprintf('=======================================================\n');
    fprintf('🚀 Starting Master DWI Pipeline Orchestrator\n');
    fprintf('=======================================================\n');

    % Step 1: Parse configuration
    try
        fprintf('⚙️ [1/5] Parsing configuration from %s...\n', config_path);
        config_struct = parse_config(config_path);
        fprintf('      ✅ Done.\n');

        master_output_folder = setup_output_folders(pipeline_dir, master_output_folder);
        config_struct.master_output_folder = master_output_folder;
    catch ME
        fprintf('❌ FAILED.\n');
        fprintf('❌ Error parsing configuration: %s\n', ME.message);
        fb_fid = fopen(fullfile(pipeline_dir, '..', 'error.log'), 'a');
        if fb_fid > 0
            fprintf(fb_fid, '[%s] [ERROR] Error parsing configuration: %s\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
            fclose(fb_fid);
        end
        % Restore figure visibility — caller will never see session.prev_fig_vis
        set(0, 'DefaultFigureVisible', session.prev_fig_vis);
        session.abort = true;
        return;
    end

    % Everything below can throw (mkdir, diary, fopen, clear_pipeline_cache,
    % etc.).  Wrap in try-catch so we clean up global state before
    % re-throwing — the caller's onCleanup cannot help because it never
    % receives the session struct on an error path.
    log_fid = -1;
    try
        % Resolve DWI type to scalar
        current_dtype = config_struct.dwi_types_to_run;
        if ~isscalar(current_dtype)
            current_dtype = current_dtype(1);
            config_struct.dwi_types_to_run = current_dtype;
        end
        dwi_type_names = {'Standard', 'dnCNN', 'IVIMnet'};
        current_name = dwi_type_names{current_dtype};

        fprintf('\n=======================================================\n');
        fprintf('\xf0\x9f\x8e\xaf EXECUTING PIPELINE FOR TARGET: %s\n', upper(current_name));
        fprintf('=======================================================\n');

        config_struct.dwi_type_name = current_name;

        % Type-specific output folder
        type_output_folder = fullfile(master_output_folder, current_name);
        if ~exist(type_output_folder, 'dir'), mkdir(type_output_folder); end
        config_struct.output_folder = type_output_folder;

        % Master diary
        master_diary_file = fullfile(type_output_folder, sprintf('pipeline_log_%s.txt', current_name));
        if exist(master_diary_file, 'file'), delete(master_diary_file); end
        diary(master_diary_file);

        % Error log
        error_log_file = fullfile(master_output_folder, 'error.log');
        log_fid = fopen(error_log_file, 'a');
        if log_fid > 0
            fprintf(log_fid, '\n[%s] ===== Pipeline run started (type: %s) =====\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), current_name);
        end
        lastwarn('');
        fprintf('      \xf0\x9f\x93\x8b Logging errors/warnings to: %s\n', error_log_file);

        % Conditionally inject config-gated steps into steps_to_run.
        % Each step is inserted at its correct position in the dependency
        % chain when its config flag is true and it is not already present.

        % compare_cores: after metrics_baseline (needs validated voxel data)
        if config_struct.run_compare_cores && ~ismember('compare_cores', steps_to_run)
            idx = find(strcmp(steps_to_run, 'metrics_baseline'));
            if ~isempty(idx)
                steps_to_run = [steps_to_run(1:idx), {'compare_cores'}, steps_to_run(idx+1:end)];
            else
                steps_to_run{end+1} = 'compare_cores';
            end
        end

        % cross_pipeline_dice: after compare_cores / metrics_baseline
        if config_struct.run_cross_pipeline_dice && ~ismember('cross_pipeline_dice', steps_to_run)
            idx = find(strcmp(steps_to_run, 'compare_cores'));
            if isempty(idx), idx = find(strcmp(steps_to_run, 'metrics_baseline')); end
            if ~isempty(idx)
                steps_to_run = [steps_to_run(1:idx), {'cross_pipeline_dice'}, steps_to_run(idx+1:end)];
            else
                steps_to_run{end+1} = 'cross_pipeline_dice';
            end
        end

        % core_failure_rates: after cross_pipeline_dice (or metrics_baseline)
        if config_struct.run_core_failure_rates && ~ismember('core_failure_rates', steps_to_run)
            idx = find(strcmp(steps_to_run, 'cross_pipeline_dice'));
            if isempty(idx), idx = find(strcmp(steps_to_run, 'compare_cores')); end
            if isempty(idx), idx = find(strcmp(steps_to_run, 'metrics_baseline')); end
            if ~isempty(idx)
                steps_to_run = [steps_to_run(1:idx), {'core_failure_rates'}, steps_to_run(idx+1:end)];
            else
                steps_to_run{end+1} = 'core_failure_rates';
            end
        end

        % core_method_outcomes: after metrics_dosimetry (needs per_method_dosimetry)
        if config_struct.run_core_method_outcomes && ~ismember('core_method_outcomes', steps_to_run)
            idx = find(strcmp(steps_to_run, 'metrics_dosimetry'));
            if ~isempty(idx)
                steps_to_run = [steps_to_run(1:idx), {'core_method_outcomes'}, steps_to_run(idx+1:end)];
            else
                steps_to_run{end+1} = 'core_method_outcomes';
            end
        end

        % per_method_cor: after cross_pipeline_dice (needs Fx1 repeat data)
        if isfield(config_struct, 'run_per_method_cor') && config_struct.run_per_method_cor && ...
                ~ismember('per_method_cor', steps_to_run)
            idx = find(strcmp(steps_to_run, 'cross_pipeline_dice'));
            if isempty(idx), idx = find(strcmp(steps_to_run, 'core_failure_rates')); end
            if isempty(idx), idx = find(strcmp(steps_to_run, 'metrics_baseline')); end
            if ~isempty(idx)
                steps_to_run = [steps_to_run(1:idx), {'per_method_cor'}, steps_to_run(idx+1:end)];
            else
                steps_to_run{end+1} = 'per_method_cor';
            end
        end

        % subvolume_stability: after core_failure_rates (needs validated voxel data)
        if isfield(config_struct, 'run_subvolume_stability') && config_struct.run_subvolume_stability && ...
                ~ismember('subvolume_stability', steps_to_run)
            idx = find(strcmp(steps_to_run, 'core_failure_rates'));
            if isempty(idx), idx = find(strcmp(steps_to_run, 'metrics_baseline')); end
            if ~isempty(idx)
                steps_to_run = [steps_to_run(1:idx), {'subvolume_stability'}, steps_to_run(idx+1:end)];
            else
                steps_to_run{end+1} = 'subvolume_stability';
            end
        end

        % dose_response_roc: after core_method_outcomes (needs per_method_dosimetry)
        if isfield(config_struct, 'run_dose_response_roc') && config_struct.run_dose_response_roc && ...
                ~ismember('dose_response_roc', steps_to_run)
            idx = find(strcmp(steps_to_run, 'core_method_outcomes'));
            if isempty(idx), idx = find(strcmp(steps_to_run, 'metrics_dosimetry')); end
            if ~isempty(idx)
                steps_to_run = [steps_to_run(1:idx), {'dose_response_roc'}, steps_to_run(idx+1:end)];
            else
                steps_to_run{end+1} = 'dose_response_roc';
            end
        end

        % gtv_confounding: after core_method_outcomes (needs per_method_dosimetry + baseline)
        if isfield(config_struct, 'run_gtv_confounding') && config_struct.run_gtv_confounding && ...
                ~ismember('gtv_confounding', steps_to_run)
            idx = find(strcmp(steps_to_run, 'dose_response_roc'));
            if isempty(idx), idx = find(strcmp(steps_to_run, 'core_method_outcomes')); end
            if isempty(idx), idx = find(strcmp(steps_to_run, 'metrics_dosimetry')); end
            if ~isempty(idx)
                steps_to_run = [steps_to_run(1:idx), {'gtv_confounding'}, steps_to_run(idx+1:end)];
            else
                steps_to_run{end+1} = 'gtv_confounding';
            end
        end

        % risk_dose_concordance: after core_method_outcomes (needs predictive + dosimetry)
        if isfield(config_struct, 'run_risk_dose_concordance') && config_struct.run_risk_dose_concordance && ...
                ~ismember('risk_dose_concordance', steps_to_run)
            idx = find(strcmp(steps_to_run, 'gtv_confounding'));
            if isempty(idx), idx = find(strcmp(steps_to_run, 'dose_response_roc')); end
            if isempty(idx), idx = find(strcmp(steps_to_run, 'core_method_outcomes')); end
            if ~isempty(idx)
                steps_to_run = [steps_to_run(1:idx), {'risk_dose_concordance'}, steps_to_run(idx+1:end)];
            else
                steps_to_run{end+1} = 'risk_dose_concordance';
            end
        end

        % Pipeline progress GUI (optional — failure here must not abort the pipeline)
        pipeGUI = [];
        if ProgressGUI.isDisplayAvailable()
            try
                pipeGUI = PipelineProgressGUI(steps_to_run, current_name);
            catch ME_gui
                fprintf('⚠️ PipelineProgressGUI creation failed: %s\n', ME_gui.message);
                pipeGUI = [];
            end
        end

        % Build type-specific file paths
        voxel_cache_file = fullfile(config_struct.dataloc, sprintf('pipeline_voxels_%s.mat', current_name));
        voxel_cache_fallback_file = fullfile(config_struct.dataloc, 'pipeline_voxels.mat');
        summary_metrics_file = fullfile(config_struct.output_folder, sprintf('summary_metrics_%s.mat', current_name));
        results_file = fullfile(config_struct.output_folder, sprintf('calculated_results_%s.mat', current_name));
        baseline_results_file = fullfile(config_struct.output_folder, sprintf('metrics_baseline_results_%s.mat', current_name));
        dosimetry_results_file = fullfile(config_struct.output_folder, sprintf('metrics_dosimetry_results_%s.mat', current_name));
        predictive_results_file = fullfile(config_struct.output_folder, sprintf('metrics_stats_predictive_results_%s.mat', current_name));

        % Clear cached files if requested
        clear_pipeline_cache(config_struct);

    catch ME
        % Clean up resources this function acquired — the caller's
        % onCleanup will never fire because session was not returned.
        if log_fid > 0
            fclose(log_fid);
        end
        set(0, 'DefaultFigureVisible', session.prev_fig_vis);
        rethrow(ME);
    end

    % Pack session struct
    session.config_struct = config_struct;
    session.steps_to_run = steps_to_run;
    session.current_dtype = current_dtype;
    session.current_name = current_name;
    session.master_output_folder = master_output_folder;
    session.type_output_folder = type_output_folder;
    session.master_diary_file = master_diary_file;
    session.log_fid = log_fid;
    session.pipeGUI = pipeGUI;
    session.voxel_cache_file = voxel_cache_file;
    session.voxel_cache_fallback_file = voxel_cache_fallback_file;
    session.summary_metrics_file = summary_metrics_file;
    session.results_file = results_file;
    session.baseline_results_file = baseline_results_file;
    session.dosimetry_results_file = dosimetry_results_file;
    session.predictive_results_file = predictive_results_file;
end
