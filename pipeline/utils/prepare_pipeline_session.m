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
%       dwi_vectors_file          - Type-specific dwi_vectors .mat path
%       fallback_dwi_vectors_file - Legacy dwi_vectors.mat path
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

    % Suppress figure windows — all plots are saved to disk via saveas()
    session.prev_fig_vis = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');

    fprintf('=======================================================\n');
    fprintf('\xf0\x9f\x9a\x80 Starting Master DWI Pipeline Orchestrator\n');
    fprintf('=======================================================\n');

    % Step 1: Parse configuration
    try
        fprintf('\xe2\x9a\x99\xef\xb8\x8f [1/5] Parsing configuration from %s...\n', config_path);
        config_struct = parse_config(config_path);
        fprintf('      \xe2\x9c\x85 Done.\n');

        master_output_folder = setup_output_folders(pipeline_dir, master_output_folder);
        config_struct.master_output_folder = master_output_folder;
    catch ME
        fprintf('\xe2\x9d\x8c FAILED.\n');
        fprintf('\xe2\x9d\x8c Error parsing configuration: %s\n', ME.message);
        fb_fid = fopen(fullfile(pipeline_dir, '..', 'error.log'), 'a');
        if fb_fid > 0
            fprintf(fb_fid, '[%s] [ERROR] Error parsing configuration: %s\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
            fclose(fb_fid);
        end
        session.abort = true;
        return;
    end

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

    % Conditionally inject compare_cores
    if config_struct.run_compare_cores && ~ismember('compare_cores', steps_to_run)
        idx = find(strcmp(steps_to_run, 'metrics_baseline'));
        if ~isempty(idx)
            steps_to_run = [steps_to_run(1:idx), {'compare_cores'}, steps_to_run(idx+1:end)];
        else
            steps_to_run{end+1} = 'compare_cores';
        end
    end

    % Pipeline progress GUI
    pipeGUI = [];
    if ProgressGUI.isDisplayAvailable()
        pipeGUI = PipelineProgressGUI(steps_to_run, current_name);
    end

    % Build type-specific file paths
    dwi_vectors_file = fullfile(config_struct.dataloc, sprintf('dwi_vectors_%s.mat', current_name));
    fallback_dwi_vectors_file = fullfile(config_struct.dataloc, 'dwi_vectors.mat');
    summary_metrics_file = fullfile(config_struct.output_folder, sprintf('summary_metrics_%s.mat', current_name));
    results_file = fullfile(config_struct.output_folder, sprintf('calculated_results_%s.mat', current_name));
    baseline_results_file = fullfile(config_struct.output_folder, sprintf('metrics_baseline_results_%s.mat', current_name));
    dosimetry_results_file = fullfile(config_struct.output_folder, sprintf('metrics_dosimetry_results_%s.mat', current_name));
    predictive_results_file = fullfile(config_struct.output_folder, sprintf('metrics_stats_predictive_results_%s.mat', current_name));

    % Clear cached files if requested
    clear_pipeline_cache(config_struct);

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
    session.dwi_vectors_file = dwi_vectors_file;
    session.fallback_dwi_vectors_file = fallback_dwi_vectors_file;
    session.summary_metrics_file = summary_metrics_file;
    session.results_file = results_file;
    session.baseline_results_file = baseline_results_file;
    session.dosimetry_results_file = dosimetry_results_file;
    session.predictive_results_file = predictive_results_file;
end
