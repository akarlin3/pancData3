function run_dwi_pipeline(config_path)
% RUN_DWI_PIPELINE Master orchestrator function for the DWI analysis pipeline
%
% Usage:
%   run_dwi_pipeline('config.json');
%
% This function sequentially calls the module scripts:
%   1. load_dwi_data
%   2. sanity_checks
%   3. metrics
%   4. visualize_results
% 
% It explicitly passes data between modules to avoid workspace pollution
% and includes error handling to halt execution if checks fail.
% Dynamically set up the MATLAB path based on the orchestrator's location
    pipeline_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(pipeline_dir, 'core'));
    addpath(fullfile(pipeline_dir, 'utils'));
    addpath(fullfile(pipeline_dir, 'dependencies'));
    % (Do not add the 'tests' directory to the main path to keep production clean)

    if nargin < 1
        config_path = 'config.json';
    end

    % Ensure all subdirectories (core, utils, dependencies) are on the path
    addpath(genpath(pwd));

    fprintf('=======================================================\n');
    fprintf('Starting Master DWI Pipeline Orchestrator\n');
    fprintf('=======================================================\n');

    % Step 1: Parse configuration
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
    try
        fprintf('[2/5] Loading DWI data... \n');
        [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = load_dwi_data(config_struct);
        fprintf('      Done: Successfully loaded data.\n');
    catch ME
        fprintf('FAILED.\n');
        fprintf('Error during data loading: %s\n', ME.message);
        return; % Halt pipeline
    end

    % Step 3: Sanity Checks
    try
        fprintf('[3/5] Running sanity checks... ');
        % NOTE: Assumes sanity_checks.m accepts the loaded data as inputs
        % and returns a flag (is_valid) and/or sanitized data.
        % Adjust the function signature below if sanity_checks.m requires
        % different input arguments based on its actual implementation.
        [is_valid, validation_msg, validated_data_gtvp, validated_data_gtvn] = sanity_checks(data_vectors_gtvp, data_vectors_gtvn);
        
        if ~is_valid
            error('Sanity checks failed: %s', validation_msg);
        end
        fprintf('Passed.\n');
    catch ME
        fprintf('FAILED.\n');
        fprintf('Pipeline halted due to sanity check failure: %s\n', ME.message);
        return; % Halt pipeline MUST STOP HERE IF DATA IS CORRUPT
    end

    % Step 4: Calculate Metrics
    try
        fprintf('[4/5] Calculating metrics... ');
        % NOTE: Assumes metrics.m accepts the validated data and original
        % summary metrics, returning a struct of computed results.
        % Adjust as needed based on the actual metrics.m signature.
        calculated_results = metrics(validated_data_gtvp, validated_data_gtvn, summary_metrics, config_struct);
        fprintf('Done.\n');
    catch ME
        fprintf('FAILED.\n');
        fprintf('Error during metrics calculation: %s\n', ME.message);
        return; % Halt pipeline
    end

    % Step 5: Visualize Results
    try
        fprintf('[5/5] Generating visualizations... ');
        % NOTE: Assumes visualize_results.m accepts the calculated results
        % and configuration to produce plots.
        % Adjust as needed based on the actual visualize_results.m signature.
        visualize_results(calculated_results, config_struct);
        fprintf('Done.\n');
    catch ME
        fprintf('FAILED (Non-Fatal).\n');
        fprintf('Error generating visualizations: %s\n', ME.message);
        % We don't necessarily need to return here since it's the last step.
    end

    fprintf('=======================================================\n');
    fprintf('Pipeline Execution Complete.\n');
    fprintf('=======================================================\n');

end
