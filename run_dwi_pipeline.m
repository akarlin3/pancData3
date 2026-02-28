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
%   run_dwi_pipeline('config.json', {'load', 'sanity', 'visualize', 'metrics_baseline', 'metrics_longitudinal', 'metrics_dosimetry', 'metrics_stats_comparisons', 'metrics_stats_predictive', 'metrics_survival'});
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
        steps_to_run = {'load', 'sanity', 'visualize', 'metrics_baseline', 'metrics_longitudinal', 'metrics_dosimetry', 'metrics_stats_comparisons', 'metrics_stats_predictive', 'metrics_survival'};
    end

    if nargin < 3
        master_output_folder = '';
    end

    fprintf('=======================================================\n');
    fprintf('🚀 Starting Master DWI Pipeline Orchestrator\n');
    fprintf('=======================================================\n');

    % Step 1: Parse configuration
    try
        fprintf('⚙️ [1/5] Parsing configuration from %s... ', config_path);
        config_struct = parse_config(config_path);
        fprintf('✅ Done.\n');
        
        % --- Master Output Folder Logic ---
        global MASTER_OUTPUT_FOLDER;
        
        if isempty(master_output_folder) && isempty(MASTER_OUTPUT_FOLDER)
            % First run in this MATLAB session without a specific folder
            timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
            master_output_folder = fullfile(pwd, sprintf('saved_figures_%s', timestamp_str));
            if ~exist(master_output_folder, 'dir'), mkdir(master_output_folder); end
            MASTER_OUTPUT_FOLDER = master_output_folder;
            fprintf('      📁 Created NEW master output folder: %s\n', master_output_folder);
            fprintf('      💡 (To start a new session folder, run: clear global MASTER_OUTPUT_FOLDER)\n');
        elseif isempty(master_output_folder) && ~isempty(MASTER_OUTPUT_FOLDER)
            % Subsequent run in the same session
            master_output_folder = MASTER_OUTPUT_FOLDER;
            fprintf('      📁 Reusing EXISTING master output folder: %s\n', master_output_folder);
        else
            % User explicitly provided a folder
            if ~exist(master_output_folder, 'dir'), mkdir(master_output_folder); end
            MASTER_OUTPUT_FOLDER = master_output_folder;
            fprintf('      📁 Using explicitly provided master output folder: %s\n', master_output_folder);
        end
        
        config_struct.master_output_folder = master_output_folder;
    catch ME
        fprintf('❌ FAILED.\n');
        fprintf('❌ Error parsing configuration: %s\n', ME.message);
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
    fprintf('🎯 EXECUTING PIPELINE FOR TARGET: %s\n', upper(current_name));
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
            fprintf('⚙️ [2/5] [%s] Loading DWI data... \n', current_name);
            [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = load_dwi_data(config_struct);

            % Save intermediate results for subsequent steps
            save(summary_metrics_file, 'summary_metrics');
            fprintf('      💾 Saved summary_metrics to %s\n', summary_metrics_file);

            fprintf('      ✅ Done: Successfully loaded data.\n');
        catch ME
            fprintf('❌ FAILED.\n');
            fprintf('❌ Error during data loading: %s\n', ME.message);
            return; % Halt pipeline
        end
    else
        fprintf('⏭️ [2/5] [%s] Skipping Load Step. Loading from disk...\n', current_name);
        try
            fallback_dwi_vectors_file = fullfile(config_struct.dataloc, 'dwi_vectors.mat');
            if exist(dwi_vectors_file, 'file')
                target_dwi_file = dwi_vectors_file;
            elseif exist(fallback_dwi_vectors_file, 'file')
                target_dwi_file = fallback_dwi_vectors_file;
            else
                target_dwi_file = '';
            end
            
            if ~isempty(target_dwi_file) && exist(summary_metrics_file, 'file')
                tmp_vectors = load(target_dwi_file, 'data_vectors_gtvp', 'data_vectors_gtvn');
                data_vectors_gtvp = tmp_vectors.data_vectors_gtvp;
                data_vectors_gtvn = tmp_vectors.data_vectors_gtvn;

                tmp_metrics = load(summary_metrics_file, 'summary_metrics');
                summary_metrics = tmp_metrics.summary_metrics;
                fprintf('      💾 Loaded data from disk (%s).\n', target_dwi_file);
            else
                error('Data files not found. Please run "load" step first.');
            end
        catch ME
             fprintf('❌ FAILED.\n');
             fprintf('❌ Error loading data from disk: %s\n', ME.message);
             return;
        end
    end

    % Step 3: Sanity Checks
    if ismember('sanity', steps_to_run)
        try
            fprintf('⚙️ [3/5] [%s] Running sanity checks... ', current_name);
            [is_valid, validation_msg, validated_data_gtvp, validated_data_gtvn] = sanity_checks(data_vectors_gtvp, data_vectors_gtvn, summary_metrics, config_struct);

            if ~is_valid
                error('Sanity checks failed: %s', validation_msg);
            end
            fprintf('✅ Passed.\n');
        catch ME
            fprintf('❌ FAILED.\n');
            fprintf('❌ Pipeline halted due to sanity check failure: %s\n', ME.message);
            return; % Halt pipeline MUST STOP HERE IF DATA IS CORRUPT
        end
    else
        fprintf('⏭️ [3/5] [%s] Skipping Sanity Checks.\n', current_name);
        validated_data_gtvp = data_vectors_gtvp;
        validated_data_gtvn = data_vectors_gtvn;
    end

    % Step 4: Visualize Results
    if ismember('visualize', steps_to_run)
        try
            fprintf('\n⚙️ [4/5] [%s] Visualizing results...\n', current_name);
            % Load existing results if visualize is run independently
            if ~exist('calculated_results', 'var')
                if exist(results_file, 'file')
                    tmp_results = load(results_file, 'calculated_results');
                    calculated_results = tmp_results.calculated_results;
                    fprintf('      💾 Loaded calculated_results from disk for visualization.\n');
                else
                    fprintf('      ⚠️ Warning: calculated_results file not found. Visualizations may be incomplete.\n');
                    calculated_results = struct(); % Empty struct fallback
                end
            end
            visualize_results(data_vectors_gtvp, summary_metrics, calculated_results, config_struct);
            fprintf('      ✅ Done: Visualizations generated.\n');
        catch ME
            fprintf('⚠️ FAILED (Non-Fatal).\n');
            fprintf('⚠️ Error generating visualizations: %s\n', ME.message);
        end
    else
        fprintf('⏭️ [4/5] [%s] Skipping Visualization.\n', current_name);
    end

    % Step 5: Calculate Metrics
    baseline_results_file = fullfile(config_struct.dataloc, sprintf('metrics_baseline_results_%s.mat', current_name));

    if ismember('metrics_baseline', steps_to_run)
        try
            fprintf('\n⚙️ [5.1/5] [%s] Running metrics_baseline... ', current_name);
            [m_lf, m_total_time, m_total_follow_up_time, m_gtv_vol, m_adc_mean, m_d_mean, m_f_mean, m_dstar_mean, ...
             m_id_list, m_mrn_list, m_d95_gtvp, m_v50gy_gtvp, m_data_vectors_gtvp, lf_group, valid_pts, ...
             ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_pct, Dstar_pct, ...
             nTp, metric_sets, set_names, time_labels, dtype_label, dl_provenance] = ...
             metrics_baseline(validated_data_gtvp, validated_data_gtvn, summary_metrics, config_struct);
             
            save(baseline_results_file, 'm_lf', 'm_total_time', 'm_total_follow_up_time', 'm_gtv_vol', 'm_adc_mean', 'm_d_mean', 'm_f_mean', 'm_dstar_mean', ...
             'm_id_list', 'm_mrn_list', 'm_d95_gtvp', 'm_v50gy_gtvp', 'm_data_vectors_gtvp', 'lf_group', 'valid_pts', ...
             'ADC_abs', 'D_abs', 'f_abs', 'Dstar_abs', 'ADC_pct', 'D_pct', 'f_pct', 'Dstar_pct', ...
             'nTp', 'metric_sets', 'set_names', 'time_labels', 'dtype_label', 'dl_provenance');
            fprintf('✅ Done.\n');
        catch ME
            fprintf('❌ FAILED.\n');
            fprintf('❌ Error during metrics_baseline: %s\n', ME.message);
            return;
        end
    else
        metrics_steps = {'metrics_longitudinal', 'metrics_dosimetry', 'metrics_stats_comparisons', 'metrics_stats_predictive', 'metrics_survival'};
        if any(ismember(metrics_steps, steps_to_run))
            fprintf('\n⏭️ [5.1/5] [%s] Skipping metrics_baseline. Loading from disk...\n', current_name);
            if exist(baseline_results_file, 'file')
                load(baseline_results_file);
            else
                fprintf('      ⚠️ Warning: metrics_baseline results not found. Subsequent metrics steps will fail.\n');
            end
        end
    end

    if ismember('metrics_longitudinal', steps_to_run)
        try
            fprintf('⚙️ [5.2/5] [%s] Running metrics_longitudinal... ', current_name);
            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_pct, Dstar_pct, ...
                                 nTp, dtype_label, config_struct.output_folder);
            fprintf('✅ Done.\n');
        catch ME
            fprintf('❌ FAILED.\n');
            fprintf('❌ Error during metrics_longitudinal: %s\n', ME.message);
        end
    else
        fprintf('⏭️ [5.2/5] [%s] Skipping metrics_longitudinal.\n', current_name);
    end
    
    dosimetry_results_file = fullfile(config_struct.dataloc, sprintf('metrics_dosimetry_results_%s.mat', current_name));
    if ismember('metrics_dosimetry', steps_to_run)
        try
            fprintf('⚙️ [5.3/5] [%s] Running metrics_dosimetry... ', current_name);
            [d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub] = ...
                metrics_dosimetry(m_id_list, summary_metrics.id_list, nTp, config_struct, ...
                                  m_data_vectors_gtvp, summary_metrics.gtv_locations);
            
            save(dosimetry_results_file, 'd95_adc_sub', 'v50_adc_sub', 'd95_d_sub', 'v50_d_sub', 'd95_f_sub', 'v50_f_sub', 'd95_dstar_sub', 'v50_dstar_sub');
            fprintf('✅ Done.\n');
        catch ME
            fprintf('❌ FAILED.\n');
            fprintf('❌ Error during metrics_dosimetry: %s\n', ME.message);
        end
    else
        if any(ismember({'metrics_stats_comparisons', 'metrics_stats_predictive'}, steps_to_run))
            fprintf('⏭️ [5.3/5] [%s] Skipping metrics_dosimetry. Loading from disk...\n', current_name);
            if exist(dosimetry_results_file, 'file')
                load(dosimetry_results_file);
            else
                fprintf('      ⚠️ Warning: metrics_dosimetry results not found. metrics_stats may fail.\n');
                % define defaults just to prevent hard crash occasionally
                d95_adc_sub=[]; v50_adc_sub=[]; d95_d_sub=[]; v50_d_sub=[]; d95_f_sub=[]; v50_f_sub=[]; d95_dstar_sub=[]; v50_dstar_sub=[];
            end
        else
            fprintf('⏭️ [5.3/5] [%s] Skipping metrics_dosimetry.\n', current_name);
        end
    end
    
    predictive_results_file = fullfile(config_struct.dataloc, sprintf('metrics_stats_predictive_results_%s.mat', current_name));
    
    if ismember('metrics_stats_comparisons', steps_to_run)
        try
            fprintf('⚙️ [5.4a/5] [%s] Running metrics_stats_comparisons... ', current_name);
            metrics_stats_comparisons(valid_pts, lf_group, ...
                metric_sets, set_names, time_labels, dtype_label, config_struct.output_folder, config_struct.dataloc, nTp, ...
                ADC_abs, D_abs, f_abs, Dstar_abs);
            
            fprintf('✅ Done.\n');
        catch ME
            fprintf('❌ FAILED.\n');
            fprintf('❌ Error during metrics_stats_comparisons: %s\n', ME.message);
        end
    else
        fprintf('⏭️ [5.4a/5] [%s] Skipping metrics_stats_comparisons.\n', current_name);
    end

    if ismember('metrics_stats_predictive', steps_to_run)
        try
            fprintf('⚙️ [5.4b/5] [%s] Running metrics_stats_predictive... ', current_name);
            [risk_scores_all, is_high_risk, times_km, events_km] = metrics_stats_predictive(valid_pts, lf_group, ...
                dtype_label, config_struct.output_folder, config_struct.dataloc, nTp, ...
                m_gtv_vol, summary_metrics.adc_sd, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_pct, Dstar_pct, ...
                m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, ...
                d95_dstar_sub, v50_dstar_sub, summary_metrics.id_list, config_struct.dwi_types_to_run, ...
                dl_provenance, time_labels, m_lf, m_total_time, m_total_follow_up_time);
            
            save(predictive_results_file, 'risk_scores_all', 'is_high_risk', 'times_km', 'events_km');

            % Build a calculated_results struct for any downstream visualize_results steps
            calculated_results = struct();
            calculated_results.risk_scores_all = risk_scores_all;
            calculated_results.is_high_risk = is_high_risk;
            calculated_results.times_km = times_km;
            calculated_results.events_km = events_km;
            calculated_results.m_lf = m_lf;
            calculated_results.m_id_list = m_id_list;

            % Save calculated results
            save(results_file, 'calculated_results');
            fprintf('      💾 Saved calculated_results to %s\n', results_file);
            
            fprintf('✅ Done.\n');
        catch ME
            fprintf('❌ FAILED.\n');
            fprintf('❌ Error during metrics_stats_predictive: %s\n', ME.message);
        end
    else
        if ismember('metrics_survival', steps_to_run)
            fprintf('⏭️ [5.4b/5] [%s] Skipping metrics_stats_predictive. Loading from disk...\n', current_name);
            if exist(predictive_results_file, 'file')
                load(predictive_results_file);
            else
                fprintf('      ⚠️ Warning: metrics_stats_predictive results not found. metrics_survival will fail.\n');
            end
        else
            fprintf('⏭️ [5.4b/5] [%s] Skipping metrics_stats_predictive.\n', current_name);
        end
    end
    
    if ismember('metrics_survival', steps_to_run)
        try
            fprintf('⚙️ [5.5/5] [%s] Running metrics_survival... ', current_name);
            metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, ...
                             m_total_follow_up_time, nTp, 'Survival', dtype_label);
            fprintf('✅ Done.\n');
        catch ME
            fprintf('❌ FAILED.\n');
            fprintf('❌ Error during metrics_survival: %s\n', ME.message);
        end
    else
        fprintf('⏭️ [5.5/5] [%s] Skipping metrics_survival.\n', current_name);
    end

    fprintf('=======================================================\n');
    fprintf('🎉 Pipeline Execution Complete for parameter %s\n', current_name);
    fprintf('=======================================================\n');

end
