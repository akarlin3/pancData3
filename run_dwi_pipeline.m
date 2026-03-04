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
%   run_dwi_pipeline('config.json', {'test', 'load', 'sanity', 'visualize', 'metrics_baseline', 'metrics_longitudinal', 'metrics_dosimetry', 'metrics_stats_comparisons', 'metrics_stats_predictive', 'metrics_survival'});
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
    if exist('OCTAVE_VERSION', 'builtin')
        addpath(fullfile(pipeline_dir, '.octave_compat'));
    end
    addpath(fullfile(pipeline_dir, 'dependencies'));

    if nargin < 1
        config_path = fullfile(pipeline_dir, 'config.json');
    end

    % Resolve relative config path against pipeline directory
    if ~isfile(config_path) && isfile(fullfile(pipeline_dir, config_path))
        config_path = fullfile(pipeline_dir, config_path);
    end

    if nargin < 2
        steps_to_run = {'test', 'load', 'sanity', 'visualize', 'metrics_baseline', 'metrics_longitudinal', 'metrics_dosimetry', 'metrics_stats_comparisons', 'metrics_stats_predictive', 'metrics_survival'};
    end

    if nargin < 3
        master_output_folder = '';
    end

    % 1.5) Run test suite before pipeline execution (once per session)
    % Uses a persistent variable so tests only run on the first call,
    % avoiding redundant re-runs when execute_all_workflows calls this
    % function multiple times.
    persistent tests_passed_this_session;
    skip_preflight = strcmp(getenv('SKIP_PIPELINE_PREFLIGHT'), '1');
    if ~skip_preflight
        try
            pf_raw = fileread(config_path);
            pf_cfg = jsondecode(pf_raw);
            if isfield(pf_cfg, 'skip_tests') && pf_cfg.skip_tests
                skip_preflight = true;
            end
        catch
            % If config can't be read here, let parse_config handle the error later
        end
    end
    if ismember('test', steps_to_run)
        if ~skip_preflight && (isempty(tests_passed_this_session) || ~tests_passed_this_session)
            try
                fprintf('⚙️ [Pre-flight] Running unit tests before pipeline...\n');
                % Run tests directly instead of via run_all_tests.m to:
                %  - Exclude integration tests (test_dwi_pipeline, test_modularity)
                %    that call run_dwi_pipeline themselves, which would be circular.
                %  - Skip the CodeCoveragePlugin (pre-flight only needs pass/fail).
                %  - Keep pre-flight fast.
                tests_dir = fullfile(pipeline_dir, 'tests');
                addpath(tests_dir);
                pf_suite = matlab.unittest.TestSuite.fromFolder(tests_dir, 'IncludingSubfolders', true);
                % Exclude integration tests that invoke the pipeline
                integration_names = {'test_dwi_pipeline', 'test_modularity'};
                keep = true(size(pf_suite));
                for ii = 1:numel(pf_suite)
                    for jj = 1:numel(integration_names)
                        if startsWith(pf_suite(ii).Name, integration_names{jj})
                            keep(ii) = false;
                        end
                    end
                end
                pf_suite = pf_suite(keep);
                % Capture test output to its own diary file
                if ~isempty(master_output_folder) && exist(master_output_folder, 'dir')
                    pf_diary = fullfile(master_output_folder, 'preflight_tests_output.log');
                    diary(pf_diary);
                end
                old_fig_vis = get(0, 'DefaultFigureVisible');
                set(0, 'DefaultFigureVisible', 'off');
                pf_results = run(pf_suite);
                set(0, 'DefaultFigureVisible', old_fig_vis);
                diary off;
                if any([pf_results.Failed])
                    error('PreFlight:TestFailure', '%d test(s) failed.', sum([pf_results.Failed]));
                end
                tests_passed_this_session = true;
                fprintf('      ✅ %d unit tests passed.\n', numel(pf_results));
            catch ME
                diary off;
                tests_passed_this_session = false;
                fprintf('❌ Test suite failed: %s\n', ME.message);
                error('PipelineAborted:TestFailure', ...
                    'Pipeline aborted because the test suite did not pass.');
            end
        else
            fprintf('⏭️ [Pre-flight] Test suite already passed this session or configured to skip. Skipping.\n');
        end
    else
        fprintf('⏭️ [Pre-flight] Skipping test step.\n');
    end

    % If 'test' was the only requested step, stop here.
    other_steps = setdiff(steps_to_run, {'test'});
    if isempty(other_steps)
        fprintf('✅ Test-only run complete. No pipeline steps to execute.\n');
        return;
    end

    % 2) Programmatically check for required toolboxes
    % The pipeline relies on specific toolboxes (Stats, Image).
    % Verification prevents obscure runtime errors deep within the code.
    if ~exist('OCTAVE_VERSION', 'builtin')
        if ~license('test', 'Statistics_Toolbox')
            error('InitializationError:MissingToolbox', ...
                'The "Statistics and Machine Learning Toolbox" is required but not installed or licensed.');
        end

        if ~license('test', 'Image_Toolbox')
            error('InitializationError:MissingToolbox', ...
                'The "Image Processing Toolbox" is required but not installed or licensed.');
        end
    end
    % ----------------------------

    log_fid = -1; % Error log file handle (opened after output folder is determined)

    % Suppress figure windows for the entire pipeline run — all plots are
    % saved to disk via saveas(), so visible windows are unnecessary.
    prev_fig_vis = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    cleanup_fig_vis = onCleanup(@() set(0, 'DefaultFigureVisible', prev_fig_vis));

    fprintf('=======================================================\n');
    fprintf('🚀 Starting Master DWI Pipeline Orchestrator\n');
    fprintf('=======================================================\n');

    % Step 1: Parse configuration
    try
        fprintf('⚙️ [1/5] Parsing configuration from %s...\n', config_path);
        config_struct = parse_config(config_path);
        fprintf('      ✅ Done.\n');

        % --- Master Output Folder Logic ---
        % Use a persistent variable scoped to this function instead of a
        % global, which can leak across unrelated MATLAB sessions.
        % Persistent variable is ONLY used when execute_all_workflows
        % passes an explicit folder (3rd arg) and subsequent calls within
        % the same multi-DWI-type session need to reuse it.  Direct
        % manual calls (no 3rd arg) always create a fresh timestamped
        % folder, avoiding silent reuse of a previous run's directory.
        persistent MASTER_OUTPUT_FOLDER;

        if ~isempty(master_output_folder)
            % User/execute_all_workflows explicitly provided a folder
            if ~exist(master_output_folder, 'dir'), mkdir(master_output_folder); end
            MASTER_OUTPUT_FOLDER = master_output_folder;
            fprintf('      📁 Using explicitly provided master output folder: %s\n', master_output_folder);
        elseif ~isempty(MASTER_OUTPUT_FOLDER) && exist(MASTER_OUTPUT_FOLDER, 'dir')
            % Subsequent call within the same execute_all_workflows session
            master_output_folder = MASTER_OUTPUT_FOLDER;
            fprintf('      📁 Reusing EXISTING master output folder: %s\n', master_output_folder);
        else
            % Standalone run — always create a fresh folder
            timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
            master_output_folder = fullfile(pipeline_dir, sprintf('saved_files_%s', timestamp_str));
            if ~exist(master_output_folder, 'dir'), mkdir(master_output_folder); end
            MASTER_OUTPUT_FOLDER = master_output_folder;
            fprintf('      📁 Created NEW master output folder: %s\n', master_output_folder);
        end

        config_struct.master_output_folder = master_output_folder;
    catch ME
        fprintf('❌ FAILED.\n');
        fprintf('❌ Error parsing configuration: %s\n', ME.message);
        fb_fid = fopen(fullfile(pipeline_dir, 'error.log'), 'a');
        if fb_fid > 0
            fprintf(fb_fid, '[%s] [ERROR] Error parsing configuration: %s\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
            fclose(fb_fid);
        end
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

    % Master diary: capture orchestrator-level console output.
    % Each core module runs its own diary (overriding this one), then calls
    % diary off.  We restart this master diary after every module returns.
    master_diary_file = fullfile(type_output_folder, sprintf('pipeline_log_%s.txt', current_name));
    if exist(master_diary_file, 'file'), delete(master_diary_file); end
    diary(master_diary_file);

    % Open error log file in master output folder
    error_log_file = fullfile(master_output_folder, 'error.log');
    log_fid = fopen(error_log_file, 'a');
    if log_fid > 0
        fprintf(log_fid, '\n[%s] ===== Pipeline run started (type: %s) =====\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), current_name);
    end
    cleanup_log = onCleanup(@() fclose(log_fid));
    lastwarn(''); % Reset MATLAB warning tracker
    fprintf('      📋 Logging errors/warnings to: %s\n', error_log_file);

    % Set isolated file names
    dwi_vectors_file = fullfile(config_struct.dataloc, sprintf('dwi_vectors_%s.mat', current_name));
    summary_metrics_file = fullfile(config_struct.output_folder, sprintf('summary_metrics_%s.mat', current_name));
    results_file = fullfile(config_struct.output_folder, sprintf('calculated_results_%s.mat', current_name));

    % Step 2: Load DWI Data
    if ismember('load', steps_to_run)
        try
            fprintf('⚙️ [2/5] [%s] Loading DWI data... \n', current_name);
            [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = load_dwi_data(config_struct);

            % Save intermediate results for subsequent steps
            save(summary_metrics_file, 'summary_metrics');
            fprintf('      💾 Saved summary_metrics to %s\n', summary_metrics_file);

            fprintf('      ✅ Done: Successfully loaded data.\n');
            [warn_msg, warn_id] = lastwarn('');
            if ~isempty(warn_msg) && log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] During load_dwi_data: %s (id: %s)\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), warn_msg, warn_id);
            end
        catch ME
            if ~isempty(ME.stack)
                stack_line = ME.stack(1).line;
                stack_file = ME.stack(1).name;
                fprintf('❌ FAILED at %s:%d.\n', stack_file, stack_line);
            else
                fprintf('❌ FAILED.\n');
            end
            fprintf('❌ Error during data loading: %s\n', ME.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [ERROR] Data loading failed: %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
                if ~isempty(ME.stack)
                    fprintf(log_fid, '         at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
                end
            end
            return; % Halt pipeline
        end
    else
        fprintf('⏭️ [2/5] [%s] Skipping Load Step. Loading from disk...\n', current_name);
        try
            fallback_dwi_vectors_file = fullfile(config_struct.dataloc, 'dwi_vectors.mat');
            if exist(dwi_vectors_file, 'file')
                target_dwi_file = dwi_vectors_file;
            elseif exist(fallback_dwi_vectors_file, 'file')
                % Legacy file without DWI-type suffix exists.  Only allow
                % fallback for Standard (type 1) to prevent silently loading
                % data from a different DWI processing method.
                if current_dtype == 1
                    fprintf('  💡 Using legacy dwi_vectors.mat (no type suffix) for Standard.\n');
                    target_dwi_file = fallback_dwi_vectors_file;
                else
                    fprintf('  ❌ Type-specific file %s not found and legacy dwi_vectors.mat cannot be used for %s (risk of cross-contamination).\n', ...
                        dwi_vectors_file, current_name);
                    target_dwi_file = '';
                end
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
             if log_fid > 0
                 fprintf(log_fid, '[%s] [ERROR] Loading data from disk failed: %s\n', ...
                     datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
             end
             return;
        end
    end

    % Step 3: Sanity Checks
    if ismember('sanity', steps_to_run)
        try
            fprintf('⚙️ [3/5] [%s] Running sanity checks...\n', current_name);
            [is_valid, validation_msg, validated_data_gtvp, validated_data_gtvn] = sanity_checks(data_vectors_gtvp, data_vectors_gtvn, summary_metrics, config_struct);

            if ~is_valid
                error('Sanity checks failed: %s', validation_msg);
            end
            fprintf('      ✅ Passed.\n');

            sanity_results_file = fullfile(config_struct.output_folder, sprintf('sanity_checks_results_%s.txt', current_name));
            fid = fopen(sanity_results_file, 'w');
            fprintf(fid, 'is_valid: %d\nvalidation_msg: %s\n', is_valid, validation_msg);
            fclose(fid);
            fprintf('      💾 Saved sanity check results to %s\n', sanity_results_file);
            [warn_msg, warn_id] = lastwarn('');
            if ~isempty(warn_msg) && log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] During sanity_checks: %s (id: %s)\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), warn_msg, warn_id);
            end
        catch ME
            fprintf('❌ FAILED.\n');
            fprintf('❌ Pipeline halted due to sanity check failure: %s\n', ME.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [ERROR] Sanity checks failed: %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
            end
            return; % Halt pipeline MUST STOP HERE IF DATA IS CORRUPT
        end
    else
        fprintf('⏭️ [3/5] [%s] Skipping Sanity Checks.\n', current_name);
        if exist('data_vectors_gtvp', 'var')
            validated_data_gtvp = data_vectors_gtvp;
            validated_data_gtvn = data_vectors_gtvn;
        else
            % When running visualize alone without load or sanity, try loading
            % from the DWI vectors file on disk.
            try
                fallback_dwi_vectors_file = fullfile(config_struct.dataloc, 'dwi_vectors.mat');
                if exist(dwi_vectors_file, 'file')
                    target_dwi_file = dwi_vectors_file;
                elseif exist(fallback_dwi_vectors_file, 'file')
                    target_dwi_file = fallback_dwi_vectors_file;
                else
                    target_dwi_file = '';
                end
                if ~isempty(target_dwi_file)
                    tmp_vectors = load(target_dwi_file, 'data_vectors_gtvp', 'data_vectors_gtvn');
                    validated_data_gtvp = tmp_vectors.data_vectors_gtvp;
                    validated_data_gtvn = tmp_vectors.data_vectors_gtvn;
                    fprintf('      💾 Loaded validated_data_gtvp from disk (%s).\n', target_dwi_file);
                else
                    error('Data vectors not found. Please run "load" step first.');
                end
            catch ME_load
                fprintf('❌ Error loading data vectors from disk: %s\n', ME_load.message);
                if log_fid > 0
                    fprintf(log_fid, '[%s] [ERROR] Loading data vectors for visualize failed: %s\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME_load.message);
                end
                return;
            end
        end
    end
    diary(master_diary_file);  % restart master diary after sanity_checks
    lastwarn('');  % reset warning tracker between steps

    % Step 5: Calculate Metrics
    baseline_results_file = fullfile(config_struct.output_folder, sprintf('metrics_baseline_results_%s.mat', current_name));

    if ismember('metrics_baseline', steps_to_run)
        try
            fprintf('\n⚙️ [5.1/5] [%s] Running metrics_baseline...\n', current_name);
            [m_lf, m_total_time, m_total_follow_up_time, m_gtv_vol, m_adc_mean, m_d_mean, m_f_mean, m_dstar_mean, ...
             m_id_list, m_mrn_list, m_d95_gtvp, m_v50gy_gtvp, m_data_vectors_gtvp, lf_group, valid_pts, ...
             ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
             nTp, metric_sets, set_names, time_labels, dtype_label, dl_provenance] = ...
             metrics_baseline(validated_data_gtvp, validated_data_gtvn, summary_metrics, config_struct);

            save(baseline_results_file, 'm_lf', 'm_total_time', 'm_total_follow_up_time', 'm_gtv_vol', 'm_adc_mean', 'm_d_mean', 'm_f_mean', 'm_dstar_mean', ...
             'm_id_list', 'm_mrn_list', 'm_d95_gtvp', 'm_v50gy_gtvp', 'm_data_vectors_gtvp', 'lf_group', 'valid_pts', ...
             'ADC_abs', 'D_abs', 'f_abs', 'Dstar_abs', 'ADC_pct', 'D_pct', 'f_delta', 'Dstar_pct', ...
             'nTp', 'metric_sets', 'set_names', 'time_labels', 'dtype_label', 'dl_provenance');
            fprintf('      ✅ Done.\n');
            [warn_msg, warn_id] = lastwarn('');
            if ~isempty(warn_msg) && log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] During metrics_baseline: %s (id: %s)\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), warn_msg, warn_id);
            end
        catch ME
            fprintf('❌ FAILED.\n');
            fprintf('❌ Error during metrics_baseline: %s\n', ME.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [ERROR] metrics_baseline failed: %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
                if ~isempty(ME.stack)
                    fprintf(log_fid, '         at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
                end
            end
            return;
        end
    else
        metrics_steps = {'metrics_longitudinal', 'metrics_dosimetry', 'metrics_stats_comparisons', 'metrics_stats_predictive', 'metrics_survival'};
        if any(ismember(metrics_steps, steps_to_run))
            fprintf('\n⏭️ [5.1/5] [%s] Skipping metrics_baseline. Loading from disk...\n', current_name);
            if exist(baseline_results_file, 'file')
                tmp_base = load(baseline_results_file);
                m_lf = tmp_base.m_lf; m_total_time = tmp_base.m_total_time; m_total_follow_up_time = tmp_base.m_total_follow_up_time; m_gtv_vol = tmp_base.m_gtv_vol; m_adc_mean = tmp_base.m_adc_mean; m_d_mean = tmp_base.m_d_mean; m_f_mean = tmp_base.m_f_mean; m_dstar_mean = tmp_base.m_dstar_mean;
                m_id_list = tmp_base.m_id_list; m_mrn_list = tmp_base.m_mrn_list; m_d95_gtvp = tmp_base.m_d95_gtvp; m_v50gy_gtvp = tmp_base.m_v50gy_gtvp; m_data_vectors_gtvp = tmp_base.m_data_vectors_gtvp; lf_group = tmp_base.lf_group; valid_pts = tmp_base.valid_pts;
                ADC_abs = tmp_base.ADC_abs; D_abs = tmp_base.D_abs; f_abs = tmp_base.f_abs; Dstar_abs = tmp_base.Dstar_abs; ADC_pct = tmp_base.ADC_pct; D_pct = tmp_base.D_pct; f_delta = tmp_base.f_delta; Dstar_pct = tmp_base.Dstar_pct;
                nTp = tmp_base.nTp; metric_sets = tmp_base.metric_sets; set_names = tmp_base.set_names; time_labels = tmp_base.time_labels; dtype_label = tmp_base.dtype_label; dl_provenance = tmp_base.dl_provenance;
            else
                fprintf('❌ metrics_baseline results not found at: %s\n', baseline_results_file);
                fprintf('❌ Downstream metrics steps require baseline results. Halting pipeline.\n');
                if log_fid > 0
                    fprintf(log_fid, '[%s] [ERROR] metrics_baseline results not found at: %s. Halting pipeline.\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), baseline_results_file);
                end
                return;
            end
        end
    end
    diary(master_diary_file);  % restart master diary after metrics_baseline
    lastwarn('');  % reset warning tracker between steps

    if ismember('metrics_longitudinal', steps_to_run)
        try
            fprintf('⚙️ [5.2/5] [%s] Running metrics_longitudinal...\n', current_name);
            metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
                                 nTp, dtype_label, config_struct.output_folder);
            longitudinal_results_file = fullfile(config_struct.output_folder, sprintf('metrics_longitudinal_results_%s.txt', current_name));
            fid = fopen(longitudinal_results_file, 'w');
            fprintf(fid, 'Longitudinal metrics generated successfully.\n');
            fclose(fid);
            fprintf('      💾 Saved longitudinal results log to %s\n', longitudinal_results_file);
            fprintf('      ✅ Done.\n');
            [warn_msg, warn_id] = lastwarn('');
            if ~isempty(warn_msg) && log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] During metrics_longitudinal: %s (id: %s)\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), warn_msg, warn_id);
            end
        catch ME
            fprintf('⚠️ FAILED (Non-Fatal).\n');
            fprintf('⚠️ Error during metrics_longitudinal: %s\n', ME.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [ERROR] metrics_longitudinal failed: %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
                if ~isempty(ME.stack)
                    fprintf(log_fid, '         at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
                end
            end
        end
    else
        fprintf('⏭️ [5.2/5] [%s] Skipping metrics_longitudinal.\n', current_name);
    end
    diary(master_diary_file);  % restart master diary after metrics_longitudinal
    lastwarn('');  % reset warning tracker between steps

    dosimetry_results_file = fullfile(config_struct.output_folder, sprintf('metrics_dosimetry_results_%s.mat', current_name));
    if ismember('metrics_dosimetry', steps_to_run)
        try
            fprintf('⚙️ [5.3/5] [%s] Running metrics_dosimetry...\n', current_name);
            [d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub] = ...
                metrics_dosimetry(m_id_list, summary_metrics.id_list, nTp, config_struct, ...
                                  m_data_vectors_gtvp, summary_metrics.gtv_locations);

            save(dosimetry_results_file, 'd95_adc_sub', 'v50_adc_sub', 'd95_d_sub', 'v50_d_sub', 'd95_f_sub', 'v50_f_sub', 'd95_dstar_sub', 'v50_dstar_sub');
            fprintf('      ✅ Done.\n');
            [warn_msg, warn_id] = lastwarn('');
            if ~isempty(warn_msg) && log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] During metrics_dosimetry: %s (id: %s)\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), warn_msg, warn_id);
            end
        catch ME
            fprintf('⚠️ FAILED (Non-Fatal).\n');
            fprintf('⚠️ Error during metrics_dosimetry: %s\n', ME.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [ERROR] metrics_dosimetry failed: %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
                if ~isempty(ME.stack)
                    fprintf(log_fid, '         at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
                end
            end
        end
    else
        if any(ismember({'metrics_stats_comparisons', 'metrics_stats_predictive'}, steps_to_run))
            fprintf('⏭️ [5.3/5] [%s] Skipping metrics_dosimetry. Loading from disk...\n', current_name);
            if exist(dosimetry_results_file, 'file')
                tmp_dosimetry = load(dosimetry_results_file);
                d95_adc_sub = tmp_dosimetry.d95_adc_sub; v50_adc_sub = tmp_dosimetry.v50_adc_sub; d95_d_sub = tmp_dosimetry.d95_d_sub; v50_d_sub = tmp_dosimetry.v50_d_sub; d95_f_sub = tmp_dosimetry.d95_f_sub; v50_f_sub = tmp_dosimetry.v50_f_sub; d95_dstar_sub = tmp_dosimetry.d95_dstar_sub; v50_dstar_sub = tmp_dosimetry.v50_dstar_sub;
            else
                fprintf('      ⚠️ Warning: metrics_dosimetry results not found. metrics_stats may fail.\n');
                if log_fid > 0
                    fprintf(log_fid, '[%s] [WARNING] metrics_dosimetry results not found at: %s\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), dosimetry_results_file);
                end
                % define defaults just to prevent hard crash
                nTp_val = nTp;
                d95_adc_sub = nan(length(m_id_list), nTp_val); v50_adc_sub = nan(length(m_id_list), nTp_val);
                d95_d_sub = nan(length(m_id_list), nTp_val); v50_d_sub = nan(length(m_id_list), nTp_val);
                d95_f_sub = nan(length(m_id_list), nTp_val); v50_f_sub = nan(length(m_id_list), nTp_val);
                d95_dstar_sub = nan(length(m_id_list), nTp_val); v50_dstar_sub = nan(length(m_id_list), nTp_val);
            end
        else
            fprintf('⏭️ [5.3/5] [%s] Skipping metrics_dosimetry.\n', current_name);
        end
    end
    diary(master_diary_file);  % restart master diary after metrics_dosimetry
    lastwarn('');  % reset warning tracker between steps

    predictive_results_file = fullfile(config_struct.output_folder, sprintf('metrics_stats_predictive_results_%s.mat', current_name));

    if ismember('metrics_stats_comparisons', steps_to_run)
        try
            fprintf('⚙️ [5.4a/5] [%s] Running metrics_stats_comparisons...\n', current_name);
            metrics_stats_comparisons(valid_pts, lf_group, ...
                metric_sets, set_names, time_labels, dtype_label, config_struct.output_folder, config_struct.dataloc, nTp, ...
                ADC_abs, D_abs, f_abs, Dstar_abs);

            comparisons_results_file = fullfile(config_struct.output_folder, sprintf('metrics_stats_comparisons_results_%s.txt', current_name));
            fid = fopen(comparisons_results_file, 'w');
            fprintf(fid, 'Stats Comparisons generated successfully.\n');
            fclose(fid);
            fprintf('      💾 Saved comparisons results log to %s\n', comparisons_results_file);
            fprintf('      ✅ Done.\n');
            [warn_msg, warn_id] = lastwarn('');
            if ~isempty(warn_msg) && log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] During metrics_stats_comparisons: %s (id: %s)\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), warn_msg, warn_id);
            end
        catch ME
            fprintf('⚠️ FAILED (Non-Fatal).\n');
            fprintf('⚠️ Error during metrics_stats_comparisons: %s\n', ME.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [ERROR] metrics_stats_comparisons failed: %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
                if ~isempty(ME.stack)
                    fprintf(log_fid, '         at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
                end
            end
        end
    else
        fprintf('⏭️ [5.4a/5] [%s] Skipping metrics_stats_comparisons.\n', current_name);
    end
    diary(master_diary_file);  % restart master diary after metrics_stats_comparisons
    lastwarn('');  % reset warning tracker between steps

    if ismember('metrics_stats_predictive', steps_to_run)
        try
            fprintf('⚙️ [5.4b/5] [%s] Running metrics_stats_predictive...\n', current_name);
            [risk_scores_all, is_high_risk, times_km, events_km] = metrics_stats_predictive(valid_pts, lf_group, ...
                dtype_label, config_struct.output_folder, config_struct.dataloc, nTp, ...
                m_gtv_vol, summary_metrics.adc_sd, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
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

            fprintf('      ✅ Done.\n');
            [warn_msg, warn_id] = lastwarn('');
            if ~isempty(warn_msg) && log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] During metrics_stats_predictive: %s (id: %s)\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), warn_msg, warn_id);
            end
        catch ME
            fprintf('⚠️ FAILED (Non-Fatal).\n');
            fprintf('⚠️ Error during metrics_stats_predictive: %s\n', ME.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [ERROR] metrics_stats_predictive failed: %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
                if ~isempty(ME.stack)
                    fprintf(log_fid, '         at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
                end
            end
        end
    else
        if ismember('metrics_survival', steps_to_run)
            fprintf('⏭️ [5.4b/5] [%s] Skipping metrics_stats_predictive. Loading from disk...\n', current_name);
            if exist(predictive_results_file, 'file')
                tmp_pred = load(predictive_results_file);
                risk_scores_all = tmp_pred.risk_scores_all; is_high_risk = tmp_pred.is_high_risk; times_km = tmp_pred.times_km; events_km = tmp_pred.events_km;
            else
                fprintf('      ⚠️ Warning: metrics_stats_predictive results not found. metrics_survival will fail.\n');
                if log_fid > 0
                    fprintf(log_fid, '[%s] [WARNING] metrics_stats_predictive results not found at: %s\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), predictive_results_file);
                end
            end
        else
            fprintf('⏭️ [5.4b/5] [%s] Skipping metrics_stats_predictive.\n', current_name);
        end
    end
    diary(master_diary_file);  % restart master diary after metrics_stats_predictive
    lastwarn('');  % reset warning tracker between steps

    % Moved Step: Visualize Results (MUST run after calculated_results is prepared)
    if ismember('visualize', steps_to_run)
        try
            fprintf('\n⚙️ [5.4c/5] [%s] Visualizing results...\n', current_name);
            % Load existing results if visualize is run independently
            if ~exist('calculated_results', 'var')
                if exist(results_file, 'file')
                    tmp_results = load(results_file, 'calculated_results');
                    calculated_results = tmp_results.calculated_results;
                    fprintf('      💾 Loaded calculated_results from disk for visualization.\n');
                else
                    calculated_results = struct(); % Empty struct fallback
                end
            end
            visualize_results(validated_data_gtvp, summary_metrics, calculated_results, config_struct);
            visualize_results_file = fullfile(config_struct.output_folder, sprintf('visualize_results_state_%s.txt', current_name));
            fid = fopen(visualize_results_file, 'w');
            fprintf(fid, 'Visualizations generated successfully for: %s\n', current_name);
            fclose(fid);
            fprintf('      ✅ Done: Visualizations generated and state saved to %s.\n', visualize_results_file);
            [warn_msg, warn_id] = lastwarn('');
            if ~isempty(warn_msg) && log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] During visualize_results: %s (id: %s)\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), warn_msg, warn_id);
            end
        catch ME
            fprintf('⚠️ FAILED (Non-Fatal).\n');
            fprintf('⚠️ Error generating visualizations: %s\n', ME.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] visualize_results failed (non-fatal): %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
                if ~isempty(ME.stack)
                    fprintf(log_fid, '         at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
                end
            end
        end
    else
        fprintf('⏭️ [5.4c/5] [%s] Skipping Visualization.\n', current_name);
    end
    diary(master_diary_file);  % restart master diary after visualize_results
    lastwarn('');  % reset warning tracker between steps

    if ismember('metrics_survival', steps_to_run)
        try
            fprintf('⚙️ [5.5/5] [%s] Running metrics_survival...\n', current_name);
            % Pass actual scan days from config if available; otherwise
            % metrics_survival uses defaults and emits a warning.
            td_scan_days_cfg = [];
            if isfield(config_struct, 'td_scan_days') && ~isempty(config_struct.td_scan_days)
                td_scan_days_cfg = config_struct.td_scan_days;
            end
            metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, ...
                             m_total_follow_up_time, nTp, 'Survival', dtype_label, m_gtv_vol, config_struct.output_folder, td_scan_days_cfg);

            survival_results_file = fullfile(config_struct.output_folder, sprintf('metrics_survival_results_%s.txt', current_name));
            fid = fopen(survival_results_file, 'w');
            fprintf(fid, 'Survival metrics generated successfully.\n');
            fclose(fid);
            fprintf('      💾 Saved survival results log to %s\n', survival_results_file);
            fprintf('      ✅ Done.\n');
            [warn_msg, warn_id] = lastwarn('');
            if ~isempty(warn_msg) && log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] During metrics_survival: %s (id: %s)\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), warn_msg, warn_id);
            end
        catch ME
            fprintf('⚠️ FAILED (Non-Fatal).\n');
            fprintf('⚠️ Error during metrics_survival: %s\n', ME.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [ERROR] metrics_survival failed: %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME.message);
                if ~isempty(ME.stack)
                    fprintf(log_fid, '         at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
                end
            end
        end
    else
        fprintf('⏭️ [5.5/5] [%s] Skipping metrics_survival.\n', current_name);
    end
    diary(master_diary_file);  % restart master diary after metrics_survival
    lastwarn('');  % reset warning tracker between steps

    fprintf('=======================================================\n');
    fprintf('🎉 Pipeline Execution Complete for parameter %s\n', current_name);
    fprintf('=======================================================\n');

    if log_fid > 0
        fprintf(log_fid, '[%s] ===== Pipeline run completed (type: %s) =====\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), current_name);
    end

    diary off;  % close master diary at end of pipeline run
end
