function run_dwi_pipeline(config_path, steps_to_run)
% RUN_DWI_PIPELINE Master orchestrator function for the DWI analysis pipeline
%
% [ORCHESTRATOR PATTERN]:
% This function acts as the central controller for the DWI analysis workflow.
% Instead of monolithic scripts, the pipeline is broken down into modular steps
% (Load -> Sanity -> Metrics -> Visualize). The orchestrator manages:
%   1. Data Flow: Passing outputs from one module as inputs to the next.
%   2. Error Handling: Catching exceptions and halting execution gracefully.
%   3. Environment: Setting up paths and verifying dependencies dynamically.
%
% Usage:
%   run_dwi_pipeline('config.json');
%   run_dwi_pipeline('config.json', {'load', 'sanity', 'metrics', 'visualize'});
%
% This function sequentially calls the module scripts:
%   1. load_dwi_data
%   2. sanity_checks
%   3. metrics
%   4. visualize_results
% 
% It explicitly passes data between modules to avoid workspace pollution
% and includes error handling to halt execution if checks fail.
    % --- Initialization Block ---
    % 1) Dynamically add folders to the MATLAB path
    % This ensures that helper functions in 'core', 'utils', and 'dependencies'
    % are accessible regardless of the current working directory.
    pipeline_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(pipeline_dir, 'core'));
    addpath(fullfile(pipeline_dir, 'utils'));
    addpath(fullfile(pipeline_dir, 'dependencies'));

    % 2) Programmatically check for required toolboxes
    % The pipeline relies on specific toolboxes (Stats, Image).
    % Verification prevents obscure runtime errors deep within the code.
    if ~license('test', 'Statistics_Toolbox')
        error('InitializationError:MissingToolbox', ...
            'The "Statistics and Machine Learning Toolbox" is required but not installed or licensed.');
    end
    
    if ~license('test', 'Image_Toolbox')
        error('InitializationError:MissingToolbox', ...
            'The "Image Processing Toolbox" is required but not installed or licensed.');
    end
    % ----------------------------

    if nargin < 1
        config_path = 'config.json';
    end

    if nargin < 2
        steps_to_run = {'load', 'sanity', 'metrics', 'visualize'};
    end

    fprintf('=======================================================\n');
    fprintf('Starting Master DWI Pipeline Orchestrator\n');
    fprintf('=======================================================\n');

    % Step 1: Parse configuration
    try
        fprintf('[1/5] Parsing configuration from %s... ', config_path);
        config_struct = parse_config(config_path);
        
        % Generate timestamp and create output folder
        timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
        master_output_folder = fullfile(pwd, sprintf('saved_figures_%s', timestamp_str));
        if ~exist(master_output_folder, 'dir'), mkdir(master_output_folder); end
        config_struct.master_output_folder = master_output_folder;
        fprintf('Done. Outputs will be saved to %s\n', master_output_folder);
    catch ME
        fprintf('FAILED.\n');
        fprintf('Error parsing configuration: %s\n', ME.message);
        return; % Halt pipeline
    end

    current_dtype = config_struct.dwi_types_to_run;
    % Ensure it's a scalar for single-pass mode
    if ~isscalar(current_dtype)
        current_dtype = current_dtype(1);
        config_struct.dwi_types_to_run = current_dtype;
    end
    
    dwi_type_names = {'Standard', 'dnCNN', 'IVIMnet'};
    current_name = dwi_type_names{current_dtype};
    
    fprintf('\n=======================================================\n');
    fprintf('  EXECUTING PIPELINE FOR TARGET: %s\n', upper(current_name));
    fprintf('=======================================================\n');
    
    config_struct.dwi_type_name = current_name;
    
    type_output_folder = fullfile(master_output_folder, current_name);
    if ~exist(type_output_folder, 'dir'), mkdir(type_output_folder); end
    config_struct.output_folder = type_output_folder;

    % Set isolated file names
    dwi_vectors_file = fullfile(config_struct.dataloc, sprintf('dwi_vectors_%s.mat', current_name));
    summary_metrics_file = fullfile(config_struct.dataloc, sprintf('summary_metrics_%s.mat', current_name));
    results_file = fullfile(config_struct.dataloc, sprintf('calculated_results_%s.mat', current_name));

    % Step 2: Load DWI Data
    if ismember('load', steps_to_run)
        try
            fprintf('[2/5] [%s] Loading DWI data... \n', current_name);
            [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = load_dwi_data(config_struct);

            % Save intermediate results for subsequent steps
            save(summary_metrics_file, 'summary_metrics');
            fprintf('      Saved summary_metrics to %s\n', summary_metrics_file);

            fprintf('      Done: Successfully loaded data.\n');
        catch ME
            fprintf('FAILED.\n');
            fprintf('Error during data loading: %s\n', ME.message);
            return; % Halt pipeline
        end
    else
        fprintf('[2/5] [%s] Skipping Load Step. Loading from disk...\n', current_name);
        try
            if exist(dwi_vectors_file, 'file') && exist(summary_metrics_file, 'file')
                tmp_vectors = load(dwi_vectors_file, 'data_vectors_gtvp', 'data_vectors_gtvn');
                data_vectors_gtvp = tmp_vectors.data_vectors_gtvp;
                data_vectors_gtvn = tmp_vectors.data_vectors_gtvn;

                tmp_metrics = load(summary_metrics_file, 'summary_metrics');
                summary_metrics = tmp_metrics.summary_metrics;
                fprintf('      Loaded data from disk.\n');
            else
                error('Data files not found. Please run "load" step first.');
            end
        catch ME
             fprintf('FAILED.\n');
             fprintf('Error loading data from disk: %s\n', ME.message);
             return;
        end
    end

    % Step 3: Sanity Checks
    if ismember('sanity', steps_to_run)
        try
            fprintf('[3/5] [%s] Running sanity checks... ', current_name);
            [is_valid, validation_msg, validated_data_gtvp, validated_data_gtvn] = sanity_checks(data_vectors_gtvp, data_vectors_gtvn, summary_metrics, config_struct);

            if ~is_valid
                error('Sanity checks failed: %s', validation_msg);
            end
            fprintf('Passed.\n');
        catch ME
            fprintf('FAILED.\n');
            fprintf('Pipeline halted due to sanity check failure: %s\n', ME.message);
            return; % Halt pipeline MUST STOP HERE IF DATA IS CORRUPT
        end
    else
        fprintf('[3/5] [%s] Skipping Sanity Checks.\n', current_name);
        validated_data_gtvp = data_vectors_gtvp;
        validated_data_gtvn = data_vectors_gtvn;
    end

    % Step 4: Calculate Metrics
    if ismember('metrics', steps_to_run)
        try
            fprintf('[4/5] [%s] Calculating metrics... ', current_name);
            calculated_results = metrics(validated_data_gtvp, validated_data_gtvn, summary_metrics, config_struct);

            % Save calculated results
            save(results_file, 'calculated_results');
            fprintf('      Saved calculated_results to %s\n', results_file);

            fprintf('Done.\n');
        catch ME
            fprintf('FAILED.\n');
            fprintf('Error during metrics calculation: %s\n', ME.message);
            return; % Halt pipeline
        end
    else
        fprintf('[4/5] [%s] Skipping Metrics Calculation.\n', current_name);
        if ismember('visualize', steps_to_run)
             if exist(results_file, 'file')
                 tmp_results = load(results_file, 'calculated_results');
                 calculated_results = tmp_results.calculated_results;
                 fprintf('      Loaded calculated_results from disk.\n');
             else
                 fprintf('      Warning: calculated_results file not found. Visualizations may fail.\n');
                 calculated_results = struct(); % Empty struct fallback
             end
        end
    end

    % Step 5: Visualize Results
    if ismember('visualize', steps_to_run)
        try
            fprintf('\n[5/5] [%s] Visualizing results...\n', current_name);
            visualize_results(data_vectors_gtvp, summary_metrics, calculated_results, config_struct);
            fprintf('      Done: Visualizations generated.\n');
        catch ME
            fprintf('FAILED (Non-Fatal).\n');
            fprintf('Error generating visualizations: %s\n', ME.message);
        end
    else
        fprintf('[5/5] [%s] Skipping Visualization.\n', current_name);
    end

    fprintf('=======================================================\n');
    fprintf('Pipeline Execution Complete for parameter %s\n', current_name);
    fprintf('=======================================================\n');

end
