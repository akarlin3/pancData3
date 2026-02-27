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
    % Reads the JSON configuration file to determine input/output paths
    % and analysis parameters.
    try
        fprintf('[1/5] Parsing configuration from %s... ', config_path);
        config_struct = parse_config(config_path);
        fprintf('Done.\n');
    catch ME
        fprintf('FAILED.\n');
        fprintf('Error parsing configuration: %s\n', ME.message);
        return; % Halt pipeline
    end

    % Step 2: Load DWI Data
    % Calls core/load_dwi_data.m to process raw DICOMs or reload pre-processed data.
    % If successful, returns structured data arrays for downstream analysis.
    if ismember('load', steps_to_run)
        try
            fprintf('[2/5] Loading DWI data... \n');
            [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = load_dwi_data(config_struct);

            % Save intermediate results for subsequent steps
            summary_metrics_file = fullfile(config_struct.dataloc, 'summary_metrics.mat');
            save(summary_metrics_file, 'summary_metrics');
            fprintf('      Saved summary_metrics to %s\n', summary_metrics_file);

            fprintf('      Done: Successfully loaded data.\n');
        catch ME
            fprintf('FAILED.\n');
            fprintf('Error during data loading: %s\n', ME.message);
            return; % Halt pipeline
        end
    else
        fprintf('[2/5] Skipping Load Step. Loading from disk...\n');
        try
            dwi_vectors_file = fullfile(config_struct.dataloc, 'dwi_vectors.mat');
            summary_metrics_file = fullfile(config_struct.dataloc, 'summary_metrics.mat');

            if exist(dwi_vectors_file, 'file') && exist(summary_metrics_file, 'file')
                % Load data using structure assignment to avoid unsafe deserialization
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
    % Validates data integrity (checking for NaNs, Infs, spatial alignment issues).
    % Prevents garbage-in-garbage-out by halting if critical checks fail.
    if ismember('sanity', steps_to_run)
        try
            fprintf('[3/5] Running sanity checks... ');
            % NOTE: Assumes sanity_checks.m accepts the loaded data as inputs
            % and returns a flag (is_valid) and/or sanitized data.
            % Adjust the function signature below if sanity_checks.m requires
            % different input arguments based on its actual implementation.
            [is_valid, validation_msg, validated_data_gtvp, validated_data_gtvn] = sanity_checks(data_vectors_gtvp, data_vectors_gtvn, summary_metrics);

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
        fprintf('[3/5] Skipping Sanity Checks.\n');
        % If skipped, assume data is valid and use loaded data
        validated_data_gtvp = data_vectors_gtvp;
        validated_data_gtvn = data_vectors_gtvn;
    end

    % Step 4: Calculate Metrics
    % Computes statistical metrics (survival analysis, feature selection,
    % longitudinal trending) on the validated data.
    if ismember('metrics', steps_to_run)
        try
            fprintf('[4/5] Calculating metrics... ');
            % NOTE: Assumes metrics.m accepts the validated data and original
            % summary metrics, returning a struct of computed results.
            % Adjust as needed based on the actual metrics.m signature.
            calculated_results = metrics(validated_data_gtvp, validated_data_gtvn, summary_metrics, config_struct);

            % Save calculated results
            results_file = fullfile(config_struct.dataloc, 'calculated_results.mat');
            save(results_file, 'calculated_results');
            fprintf('      Saved calculated_results to %s\n', results_file);

            fprintf('Done.\n');
        catch ME
            fprintf('FAILED.\n');
            fprintf('Error during metrics calculation: %s\n', ME.message);
            return; % Halt pipeline
        end
    else
        fprintf('[4/5] Skipping Metrics Calculation.\n');
        if ismember('visualize', steps_to_run)
             % Load results if visualize is requested
             results_file = fullfile(config_struct.dataloc, 'calculated_results.mat');
             if exist(results_file, 'file')
                 tmp_results = load(results_file, 'calculated_results');
                 calculated_results = tmp_results.calculated_results;
                 fprintf('      Loaded calculated_results from disk.\n');
             else
                 fprintf('      Warning: calculated_results.mat not found. Visualizations may fail.\n');
                 calculated_results = struct(); % Empty struct fallback
             end
        end
    end

    % Step 5: Visualize Results
    % Generates plots, reports, and visual summaries of the analysis.
    % Errors here are non-fatal to preserve the calculated results.
    if ismember('visualize', steps_to_run)
        try
            fprintf('\n[5/5] Visualizing results...\n');
            % NOTE: Assumes visualize_results.m accepts the raw data, summary metrics,
            % calculated results, and configuration to produce plots.
            % Adjust as needed based on the actual visualize_results.m signature.
            visualize_results(data_vectors_gtvp, summary_metrics, calculated_results, config_struct);
            fprintf('      Done: Visualizations generated.\n');
        catch ME
            fprintf('FAILED (Non-Fatal).\n');
            fprintf('Error generating visualizations: %s\n', ME.message);
            % We don't necessarily need to return here since it's the last step.
        end
    else
        fprintf('[5/5] Skipping Visualization.\n');
    end

    fprintf('=======================================================\n');
    fprintf('Pipeline Execution Complete.\n');
    fprintf('=======================================================\n');

end
