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
    % 1) Dynamically add folders to the MATLAB path
    % This ensures that helper functions in 'core', 'utils', and 'dependencies'
    % are accessible regardless of the current working directory.
    % The pipeline uses fileparts(mfilename('fullpath')) rather than pwd() so
    % that it works correctly when invoked from any working directory -- a
    % common scenario in cluster/batch environments at MSK.
    pipeline_dir = fileparts(mfilename('fullpath'));

    if nargin < 1
        config_path = fullfile(pipeline_dir, 'config.json');
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
    % function multiple times (once per DWI type: Standard, dnCNN, IVIMnet).
    % The timestamp tracks WHEN tests last passed so we can invalidate
    % if test files have been modified since (interactive development).
    %
    % [ANALYTICAL RATIONALE]: Pre-flight testing is essential because
    % incorrect IVIM fitting, leaky cross-validation, or broken imputation
    % could produce scientifically invalid results that would only be caught
    % much later during manuscript review. Running unit tests before each
    % pipeline execution ensures that core statistical and physical constraints
    % (e.g., patient-stratified CV folds, temporal leakage bounds in KNN
    % imputation, IVIM parameter range validation) are intact before
    % committing to a multi-hour computation on the full patient cohort.
    persistent tests_passed_this_session;
    persistent tests_passed_timestamp;

    % Delegate path setup, pre-flight tests, and toolbox checks to
    % initialize_pipeline.  Persistent variables cannot live in the helper
    % (they are function-scoped), so they are passed as I/O arguments.
    [config_path, tests_passed_this_session, tests_passed_timestamp] = ...
        initialize_pipeline(pipeline_dir, config_path, steps_to_run, ...
                            master_output_folder, tests_passed_this_session, ...
                            tests_passed_timestamp);

    % If 'test' was the only requested step, stop here.
    other_steps = setdiff(steps_to_run, {'test'});
    if isempty(other_steps)
        fprintf('✅ Test-only run complete. No pipeline steps to execute.\n');
        return;
    end

    % ----------------------------

    log_fid = -1; % Error log file handle (opened after output folder is determined)

    % Suppress figure windows for the entire pipeline run — all plots are
    % saved to disk via saveas(), so visible windows are unnecessary.
    % [RATIONALE]: On headless cluster nodes (common in hospital HPC
    % environments), creating visible figure windows causes MATLAB to
    % crash or hang. Even on workstations, a full pipeline run generates
    % hundreds of parameter maps, distribution plots, and survival curves
    % — spawning that many figure windows would overwhelm the desktop.
    % The onCleanup guard restores the original visibility setting even
    % if the pipeline errors out, preventing side effects on subsequent
    % interactive MATLAB use.
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

    % [DWI TYPE RESOLUTION]:
    % dwi_types_to_run is a numeric index (1=Standard, 2=dnCNN, 3=IVIMnet)
    % parsed from the config.json dwi_type string by parse_config.m.
    % When called from execute_all_workflows, this is always a scalar
    % because the outer script sets dwi_type explicitly for each run.
    % The scalar enforcement below is a safety check for direct calls
    % where a user might accidentally set dwi_type to an array — in that
    % case, only the first type is processed (the orchestrator pattern
    % expects one type per run_dwi_pipeline invocation).
    current_dtype = config_struct.dwi_types_to_run;
    % Ensure it's a scalar for single-pass mode
    if ~isscalar(current_dtype)
        current_dtype = current_dtype(1);
        config_struct.dwi_types_to_run = current_dtype;
    end

    % Map numeric index to human-readable name used in output folder names,
    % filenames, plot titles, and log messages.
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
    cleanup_log = onCleanup(@() safe_fclose_log(log_fid));
    lastwarn(''); % Reset MATLAB warning tracker
    fprintf('      📋 Logging errors/warnings to: %s\n', error_log_file);

    % --- Conditionally inject compare_cores into default steps ---
    % When run_compare_cores is true in config, automatically include the
    % compare_cores step (pairwise Dice/Hausdorff across all 11 core
    % methods) after metrics_baseline in the pipeline.
    if config_struct.run_compare_cores && ~ismember('compare_cores', steps_to_run)
        idx = find(strcmp(steps_to_run, 'metrics_baseline'));
        if ~isempty(idx)
            steps_to_run = [steps_to_run(1:idx), {'compare_cores'}, steps_to_run(idx+1:end)];
        else
            steps_to_run{end+1} = 'compare_cores';
        end
    end

    % --- Pipeline Progress GUI (GUI environments only) ---
    pipeGUI = [];
    if ProgressGUI.isDisplayAvailable()
        pipeGUI = PipelineProgressGUI(steps_to_run, current_name);
    end
    cleanup_gui = onCleanup(@() closeIfValid(pipeGUI));

    % [TYPE-ISOLATED FILE NAMING]:
    % Each DWI processing method produces its own set of output files, suffixed
    % with the method name (e.g., dwi_vectors_Standard.mat, dwi_vectors_dnCNN.mat).
    % This prevents cross-contamination: loading Standard diffusion parameters
    % when the pipeline is configured for dnCNN would produce meaningless results,
    % since the underlying voxel-wise IVIM fits are completely different.
    %
    % dwi_vectors_*.mat contains per-voxel fitted IVIM/ADC parameters within the
    % GTV mask for every patient and fraction. This is the primary data artifact.
    % summary_metrics_*.mat contains per-patient aggregated metrics (mean ADC,
    % mean D, etc.) derived from the voxel distributions.
    % calculated_results_*.mat contains predictive model outputs (risk scores,
    % high-risk classification) for downstream survival analysis and visualization.
    dwi_vectors_file = fullfile(config_struct.dataloc, sprintf('dwi_vectors_%s.mat', current_name));
    fallback_dwi_vectors_file = fullfile(config_struct.dataloc, 'dwi_vectors.mat');
    summary_metrics_file = fullfile(config_struct.output_folder, sprintf('summary_metrics_%s.mat', current_name));
    results_file = fullfile(config_struct.output_folder, sprintf('calculated_results_%s.mat', current_name));

    % --- Clear cached files if requested ---
    % When clear_cache is true, remove all pipeline-generated .mat files
    % from the data directory so the pipeline recomputes everything from
    % scratch.  This is a one-time operation on the first DWI type run;
    % subsequent types in the same execute_all_workflows session reuse
    % freshly generated caches.
    persistent cache_cleared_this_session;
    if isfield(config_struct, 'clear_cache') && config_struct.clear_cache
        if isempty(cache_cleared_this_session) || ~cache_cleared_this_session
            dataloc = config_struct.dataloc;
            cache_patterns = {
                fullfile(dataloc, 'dwi_vectors*.mat'), ...
                fullfile(dataloc, 'summary_metrics*.mat'), ...
                fullfile(dataloc, 'adc_vectors.mat')
            };
            n_deleted = 0;
            for cp = 1:numel(cache_patterns)
                cached = dir(cache_patterns{cp});
                for cf = 1:numel(cached)
                    delete(fullfile(cached(cf).folder, cached(cf).name));
                    n_deleted = n_deleted + 1;
                end
            end
            % Remove per-patient checkpoint directory
            checkpoint_dir = fullfile(dataloc, 'processed_patients');
            if isfolder(checkpoint_dir)
                rmdir(checkpoint_dir, 's');
                fprintf('  🗑️ Removed per-patient checkpoint directory.\n');
            end
            if n_deleted > 0
                fprintf('  🗑️ Cleared %d cached .mat file(s) from %s\n', n_deleted, dataloc);
            else
                fprintf('  💡 No cached files found to clear.\n');
            end
            cache_cleared_this_session = true;
        end
    end

    % Step 2: Load DWI Data
    % [ANALYTICAL RATIONALE — DATA LOADING AND MODEL FITTING]:
    % This is the computationally dominant step. For each patient and each
    % treatment fraction (Fx1-Fx5 + post-treatment), load_dwi_data:
    %   1. Converts DICOM to NIFTI via dcm2niix (standardizes orientation/headers)
    %   2. Loads GTV mask (GTVp = primary tumor, GTVn = nodal disease)
    %   3. Optionally applies DnCNN denoising (if dwi_type = 'dnCNN')
    %   4. Fits IVIM model: S(b) = S0*((1-f)*exp(-b*D) + f*exp(-b*D*))
    %      where D = true tissue diffusion (~1e-3 mm^2/s in pancreatic tissue),
    %      f = perfusion fraction (blood volume fraction, typically 5-30%),
    %      D* = pseudo-diffusion from capillary blood flow (~10x D)
    %   5. Fits ADC model: S(b) = S0*exp(-b*ADC) using only high b-values
    %      (above ivim_bthr, typically b > 100 s/mm^2) to minimize perfusion
    %      contamination in the diffusion estimate
    %   6. Extracts voxel-wise parameters within the GTV mask
    %
    % The output separates GTVp and GTVn because primary tumor and nodal
    % disease may show different diffusion characteristics and treatment
    % responses — collapsing them would obscure potentially distinct
    % biological behaviors. GTVp is the primary analysis target for
    % pancreatic cancer response assessment.
    if ismember('load', steps_to_run)
        if ~isempty(pipeGUI), pipeGUI.startStep('load'); end
        try
            fprintf('⚙️ [2/5] [%s] Loading DWI data... \n', current_name);
            [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = load_dwi_data(config_struct);

            % Save intermediate results for subsequent steps
            save(summary_metrics_file, 'summary_metrics');
            fprintf('      💾 Saved summary_metrics to %s\n', summary_metrics_file);

            fprintf('      ✅ Done: Successfully loaded data.\n');
            if ~isempty(pipeGUI), pipeGUI.completeStep('load', 'success'); end
            [warn_msg, warn_id] = lastwarn;  % read current warning
            lastwarn('');                       % then clear
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
        % [CHECKPOINT RECOVERY]:
        % When the load step is skipped (e.g., rerunning only metrics after
        % fixing a downstream bug), we must reload the previously computed
        % voxel-wise parameter vectors from disk. The type-specific filename
        % ensures we load parameters from the correct processing method.
        % The legacy fallback (dwi_vectors.mat without type suffix) exists
        % for backward compatibility with older pipeline runs that predated
        % the multi-type architecture, but is ONLY allowed for Standard to
        % prevent the dangerous scenario of loading Standard-fitted parameters
        % when running a dnCNN or IVIMnet analysis — this would silently
        % invalidate all downstream results.
        if ~isempty(pipeGUI), pipeGUI.completeStep('load', 'skipped'); end
        fprintf('⏭️ [2/5] [%s] Skipping Load Step. Loading from disk...\n', current_name);
        try
            [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = ...
                load_data_from_disk(dwi_vectors_file, fallback_dwi_vectors_file, ...
                                    summary_metrics_file, current_dtype, current_name);
            if isempty(summary_metrics)
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
    diary(master_diary_file);  % restart master diary after load step

    % Step 3: Sanity Checks
    % [ANALYTICAL RATIONALE — PARAMETER VALIDATION]:
    % After model fitting, sanity_checks validates that the fitted diffusion
    % parameters are physically and biologically plausible:
    %   - D (true diffusion) should be ~0.5-3.0 x 10^-3 mm^2/s for pancreatic
    %     tissue. Values outside this range indicate fitting failure (local
    %     minimum), corrupt input data, or misregistered GTV masks where
    %     voxels outside the tumor (e.g., bowel gas, ducts) are included.
    %   - f (perfusion fraction) must be in [0, 1]. Negative values or f > 1
    %     are nonphysical and indicate the bi-exponential IVIM fit diverged.
    %   - D* (pseudo-diffusion) should be ~5-100 x 10^-3 mm^2/s. Very large
    %     values suggest the optimizer hit the upper bound constraint.
    %   - ADC should be positive and typically ~1-2 x 10^-3 mm^2/s.
    %   - NaN/Inf values indicate voxels where the fit failed entirely.
    %
    % Sanity check FAILURE halts the pipeline entirely (return, not continue).
    % This is deliberate: if a significant fraction of voxels have implausible
    % parameters, all downstream metrics (mean ADC, percent change, survival
    % correlations) would be scientifically invalid. It is better to halt and
    % investigate than to produce misleading results.
    if ismember('sanity', steps_to_run)
        if ~isempty(pipeGUI), pipeGUI.startStep('sanity'); end
        try
            fprintf('⚙️ [3/5] [%s] Running sanity checks...\n', current_name);
            [is_valid, validation_msg, validated_data_gtvp, validated_data_gtvn] = sanity_checks(data_vectors_gtvp, data_vectors_gtvn, summary_metrics, config_struct);

            if ~is_valid
                error('Sanity checks failed: %s', validation_msg);
            end
            fprintf('      ✅ Passed.\n');
            if ~isempty(pipeGUI), pipeGUI.completeStep('sanity', 'success'); end

            sanity_results_file = fullfile(config_struct.output_folder, sprintf('sanity_checks_results_%s.txt', current_name));
            fid = fopen(sanity_results_file, 'w');
            if fid < 0
                warning('run_dwi_pipeline:fileWriteFailed', 'Cannot write %s', sanity_results_file);
            else
                fprintf(fid, 'is_valid: %d\nvalidation_msg: %s\n', is_valid, validation_msg);
                fclose(fid);
            end
            fprintf('      💾 Saved sanity check results to %s\n', sanity_results_file);
            [warn_msg, warn_id] = lastwarn;  % read current warning
            lastwarn('');                       % then clear
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
        if ~isempty(pipeGUI), pipeGUI.completeStep('sanity', 'skipped'); end
        fprintf('⏭️ [3/5] [%s] Skipping Sanity Checks.\n', current_name);
        if exist('data_vectors_gtvp', 'var')
            validated_data_gtvp = data_vectors_gtvp;
            validated_data_gtvn = data_vectors_gtvn;
        else
            % When running visualize alone without load or sanity, try loading
            % from the DWI vectors file on disk.
            try
                [validated_data_gtvp, validated_data_gtvn, sm_tmp] = ...
                    load_data_from_disk(dwi_vectors_file, fallback_dwi_vectors_file, ...
                                        summary_metrics_file, current_dtype, current_name);
                % Also load summary_metrics which is needed by
                % visualize_results and metrics_baseline downstream.
                if ~exist('summary_metrics', 'var') && ~isempty(sm_tmp)
                    summary_metrics = sm_tmp;
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
    % [ANALYTICAL RATIONALE — METRICS BASELINE]:
    % metrics_baseline is the analytical foundation. It:
    %   1. Reads the clinical spreadsheet to obtain outcome data (local failure,
    %      overall survival, follow-up time) and clinical covariates (GTV volume,
    %      dose coverage D95/V50) for each patient.
    %   2. Aggregates voxel-wise IVIM/ADC parameter distributions within each
    %      patient's GTV into summary statistics (mean, median, percentiles).
    %      This voxel-to-patient reduction is necessary because survival analysis
    %      operates at the patient level, not the voxel level.
    %   3. Computes ABSOLUTE values at each timepoint (Fx1-Fx5, post) and
    %      PERCENT CHANGE relative to baseline (Fx1). Percent change in ADC
    %      during the first week of RT is a candidate early-response biomarker:
    %      increases in ADC may reflect tumor cell kill (reduced cellularity
    %      increases water diffusivity).
    %   4. For the perfusion fraction f, ABSOLUTE DELTA (not percent change)
    %      is used because f values near zero make percent change numerically
    %      unstable and clinically uninterpretable.
    %   5. Identifies "valid_pts" — patients with complete baseline data needed
    %      for downstream analysis. Patients missing Fx1 data are excluded
    %      because percent change and longitudinal metrics require a baseline.
    %   6. Loads DL provenance to flag any overlap between IVIMnet training
    %      patients and analysis patients (data leakage check).
    %
    % The large number of output variables reflects the comprehensive set of
    % features extracted for subsequent statistical and predictive modeling.
    baseline_results_file = fullfile(config_struct.output_folder, sprintf('metrics_baseline_results_%s.mat', current_name));

    if ismember('metrics_baseline', steps_to_run)
        if ~isempty(pipeGUI), pipeGUI.startStep('metrics_baseline'); end
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
            if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_baseline', 'success'); end
            [warn_msg, warn_id] = lastwarn;  % read current warning
            lastwarn('');                       % then clear
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
            if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_baseline', 'skipped'); end
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

    % --- Compare Core Methods ---
    % Runs all 11 tumor core delineation methods on each patient/timepoint
    % and computes pairwise Dice/Hausdorff agreement metrics with summary
    % figures.  Can be auto-included by setting "run_compare_cores": true
    % in config.json, or invoked explicitly:
    %   run_dwi_pipeline('config.json', {'compare_cores'});
    if ismember('compare_cores', steps_to_run)
        execute_pipeline_step('compare_cores', @() run_compare_cores_step( ...
            validated_data_gtvp, summary_metrics, config_struct, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('compare_cores', 'skipped'); end
    end

    % [ANALYTICAL RATIONALE — LONGITUDINAL METRICS]:
    % metrics_longitudinal analyzes how diffusion parameters change across
    % treatment fractions (Fx1 through Fx5 and post-treatment). In pancreatic
    % cancer RT, the temporal trajectory of diffusion parameters carries
    % biological information:
    %   - Rising ADC/D during early fractions (Fx2-Fx3) suggests tumor cell
    %     death and reduced cellularity — a favorable treatment response signal.
    %   - Stable or decreasing ADC may indicate treatment resistance.
    %   - Changes in f (perfusion fraction) may reflect vascular damage or
    %     normalization from RT.
    %   - D* changes are the noisiest and least reliable but may capture
    %     microvascular flow alterations.
    %
    % This module generates temporal trajectory plots and computes paired
    % statistical tests (e.g., Fx2 vs Fx1) to identify significant early
    % changes that could serve as decision points for adaptive RT.
    %
    % Non-fatal: if longitudinal analysis fails (e.g., too few patients with
    % multiple fractions), the pipeline continues because downstream modules
    % (dosimetry, predictive modeling, survival) can still operate on
    % baseline and absolute metrics alone.
    if ismember('metrics_longitudinal', steps_to_run)
        execute_pipeline_step('metrics_longitudinal', @() run_metrics_longitudinal_step( ...
            ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
            nTp, dtype_label, config_struct, m_lf, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_longitudinal', 'skipped'); end
        fprintf('⏭️ [5.2/5] [%s] Skipping metrics_longitudinal.\n', current_name);
    end

    % [ANALYTICAL RATIONALE — DOSIMETRY METRICS]:
    % metrics_dosimetry bridges the gap between diffusion MRI and radiation
    % therapy dose distributions. It computes dose coverage metrics within
    % DIFFUSION-DEFINED tumor sub-volumes rather than the conventional
    % whole-GTV approach. The key innovation:
    %
    %   Traditional dosimetry: D95 and V50 are computed over the entire GTV
    %   contour drawn by the radiation oncologist.
    %
    %   Diffusion sub-volume dosimetry: The GTV is partitioned into sub-regions
    %   based on diffusion parameter thresholds (e.g., voxels with ADC below
    %   median = "restricted diffusion" sub-volume, potentially more cellular/
    %   aggressive tissue). D95 and V50 are then computed within each sub-volume.
    %
    % This tests the hypothesis that dose coverage of the most biologically
    % aggressive tumor sub-regions (identified by diffusion characteristics)
    % is a better predictor of local control than whole-GTV dose coverage.
    %
    % D95 = minimum dose covering 95% of the sub-volume (Gy)
    % V50 = fraction of sub-volume receiving >= 50 Gy
    %
    % Output: 8 matrices (patients x timepoints), one per combination of
    % {D95, V50} x {ADC, D, f, D*} sub-volume definitions.
    dosimetry_results_file = fullfile(config_struct.output_folder, sprintf('metrics_dosimetry_results_%s.mat', current_name));
    if ismember('metrics_dosimetry', steps_to_run)
        execute_pipeline_step('metrics_dosimetry', @() run_metrics_dosimetry_step( ...
            m_id_list, summary_metrics, nTp, config_struct, m_data_vectors_gtvp, ...
            dosimetry_results_file, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
        % Reload dosimetry results from disk (they were saved inside the step function)
        if exist(dosimetry_results_file, 'file')
            tmp_dosimetry = load(dosimetry_results_file);
            d95_adc_sub = tmp_dosimetry.d95_adc_sub; v50_adc_sub = tmp_dosimetry.v50_adc_sub;
            d95_d_sub = tmp_dosimetry.d95_d_sub; v50_d_sub = tmp_dosimetry.v50_d_sub;
            d95_f_sub = tmp_dosimetry.d95_f_sub; v50_f_sub = tmp_dosimetry.v50_f_sub;
            d95_dstar_sub = tmp_dosimetry.d95_dstar_sub; v50_dstar_sub = tmp_dosimetry.v50_dstar_sub;
        end
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_dosimetry', 'skipped'); end
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

    % [METRIC SET ASSEMBLY FOR STATISTICAL ANALYSIS]:
    % Append dosimetry sub-volume metrics to metric_sets for univariate
    % analysis.  metrics_baseline creates sets 1-2 (absolute diffusion values
    % + percent change from baseline); sets 3-4 (D95 + V50 dose coverage
    % within diffusion-defined sub-volumes) require dosimetry results which
    % are only available after metrics_dosimetry completes.
    %
    % The metric_sets cell array organizes all candidate features into
    % thematic groups for systematic univariate testing:
    %   Set 1: Absolute diffusion parameter values at each timepoint
    %   Set 2: Percent change in diffusion parameters from baseline
    %   Set 3: D95 dose coverage — whole GTV and each diffusion sub-volume
    %   Set 4: V50 dose coverage — whole GTV and each diffusion sub-volume
    %
    % Including both whole-GTV and sub-volume dose metrics in the same
    % analysis framework allows direct comparison of their predictive value.
    % If sub-volume D95 outperforms whole-GTV D95, it supports the hypothesis
    % that diffusion-guided tumor characterization adds value to RT planning.
    if exist('d95_adc_sub', 'var') && exist('metric_sets', 'var')
        metric_sets{3} = {m_d95_gtvp, d95_adc_sub, d95_d_sub, d95_f_sub, d95_dstar_sub};
        metric_sets{4} = {m_v50gy_gtvp, v50_adc_sub, v50_d_sub, v50_f_sub, v50_dstar_sub};
        set_names{3} = {'D95 GTVp (whole)', 'D95 Sub(ADC)', 'D95 Sub(D)', 'D95 Sub(f)', 'D95 Sub(D*)'};
        set_names{4} = {'V50 GTVp (whole)', 'V50 Sub(ADC)', 'V50 Sub(D)', 'V50 Sub(f)', 'V50 Sub(D*)'};
    end

    predictive_results_file = fullfile(config_struct.output_folder, sprintf('metrics_stats_predictive_results_%s.mat', current_name));

    % [ANALYTICAL RATIONALE — STATISTICAL GROUP COMPARISONS]:
    % metrics_stats_comparisons performs univariate hypothesis testing to
    % identify which diffusion-derived features differ significantly between
    % clinical outcome groups (e.g., local failure vs. no local failure).
    % Uses non-parametric Wilcoxon rank-sum (Mann-Whitney U) tests because:
    %   (a) Diffusion parameter distributions in pancreatic tumors are
    %       typically non-Gaussian (skewed by necrotic/cystic regions).
    %   (b) Sample sizes in pancreatic cancer studies are small (N ~ 30-60),
    %       making parametric assumptions unreliable.
    %   (c) The rank-sum test is robust to outliers from fitting artifacts.
    %
    % This step is exploratory/descriptive — it identifies promising features
    % for subsequent multivariate predictive modeling (metrics_stats_predictive).
    % Multiple comparison correction is applied to control false discovery rate.
    if ismember('metrics_stats_comparisons', steps_to_run)
        execute_pipeline_step('metrics_stats_comparisons', @() run_metrics_stats_comparisons_step( ...
            valid_pts, lf_group, metric_sets, set_names, time_labels, dtype_label, ...
            config_struct, nTp, ADC_abs, D_abs, f_abs, Dstar_abs, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_stats_comparisons', 'skipped'); end
        fprintf('⏭️ [5.4a/5] [%s] Skipping metrics_stats_comparisons.\n', current_name);
    end

    % [ANALYTICAL RATIONALE — PREDICTIVE MODELING]:
    % metrics_stats_predictive builds multivariate predictive models for
    % clinical outcomes (e.g., local failure prediction) using the full
    % feature set: diffusion parameters (absolute + change), dosimetric
    % features (D95/V50 for whole-GTV and sub-volumes), and clinical
    % covariates (GTV volume, ADC heterogeneity via adc_sd).
    %
    % Key methodological safeguards:
    %   - Patient-stratified cross-validation (make_grouped_folds): ensures
    %     no patient appears in both train and test sets, preventing
    %     intra-patient data leakage from longitudinal measurements.
    %   - KNN imputation with temporal bounds (knn_impute_train_test):
    %     missing timepoint data is imputed from other patients' same-
    %     timepoint data only, never from the same patient's other
    %     timepoints (which would leak temporal information).
    %   - Collinearity filtering (filter_collinear_features): removes
    %     highly correlated features before model fitting to stabilize
    %     coefficient estimates in the small-N, high-p regime typical of
    %     pancreatic cancer cohorts.
    %   - DL provenance checking: for IVIMnet runs, ensures analysis
    %     patients were not in the neural network's training set.
    %
    % The adc_sd (standard deviation of ADC within the GTV) is included as
    % a tumor heterogeneity feature — more heterogeneous tumors (higher
    % intra-tumoral ADC variance) may indicate mixed cellularity/necrosis
    % and different treatment response patterns.
    if ismember('metrics_stats_predictive', steps_to_run)
        execute_pipeline_step('metrics_stats_predictive', @() run_metrics_stats_predictive_step( ...
            valid_pts, lf_group, dtype_label, config_struct, nTp, m_gtv_vol, summary_metrics, ...
            m_id_list, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
            m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, ...
            d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub, dl_provenance, time_labels, ...
            m_lf, m_total_time, m_total_follow_up_time, predictive_results_file, results_file, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
        % Load calculated_results from disk if saved by the step
        if exist(results_file, 'file')
            tmp_results = load(results_file, 'calculated_results');
            calculated_results = tmp_results.calculated_results;
        end
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_stats_predictive', 'skipped'); end
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

    % [ANALYTICAL RATIONALE — VISUALIZATION AFTER PREDICTION]:
    % Visualization is intentionally placed AFTER predictive modeling rather
    % than in its traditional position after sanity checks. This ordering
    % enables the visualization module to incorporate risk stratification
    % overlays from the predictive model — for example, coloring patient
    % trajectories by predicted risk group or annotating parameter maps
    % with model confidence scores. Without the calculated_results struct,
    % visualize_results falls back to basic parameter distributions without
    % risk annotations (using an empty struct).
    if ismember('visualize', steps_to_run)
        execute_pipeline_step('visualize', @() run_visualize_step( ...
            validated_data_gtvp, summary_metrics, config_struct, results_file, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('visualize', 'skipped'); end
        fprintf('⏭️ [5.4c/5] [%s] Skipping Visualization.\n', current_name);
    end

    % [ANALYTICAL RATIONALE — SURVIVAL ANALYSIS]:
    % metrics_survival performs time-to-event analysis for pancreatic cancer
    % patients using diffusion parameters as candidate prognostic biomarkers.
    %
    % Key analytical components:
    %   - Cox proportional hazards regression: tests whether baseline or
    %     early-change diffusion parameters are independently associated with
    %     overall survival or local failure time after adjusting for known
    %     clinical prognostic factors (GTV volume, dose coverage).
    %   - Competing risks modeling via Cause-Specific Hazards (CSH): pancreatic
    %     cancer patients can experience local failure, distant metastasis, or
    %     death from other causes. Standard Cox models treat competing events
    %     as censored, which biases the hazard estimate. CSH models each
    %     event type separately while properly accounting for the competing
    %     events, yielding unbiased cause-specific hazard ratios.
    %   - Kaplan-Meier survival curves stratified by risk group (from
    %     metrics_stats_predictive) with log-rank tests.
    %
    % [SCAN DAY RESOLUTION — AVOIDING IMMORTAL TIME BIAS]:
    % The td_scan_days variable specifies the actual calendar days (relative to
    % treatment start) when each fraction's MRI was acquired. This is critical
    % for time-dependent Cox models: using assumed/default scan days (e.g.,
    % Fx1=day 0, Fx2=day 7, etc.) when actual scan dates differ introduces
    % "immortal time bias" — a well-known pitfall where the time window between
    % a covariate measurement and the event is misspecified, systematically
    % biasing hazard ratio estimates. The three-level fallback (DICOM dates ->
    % config -> defaults with warning) reflects decreasing data quality.
    if ismember('metrics_survival', steps_to_run)
        execute_pipeline_step('metrics_survival', @() run_metrics_survival_step( ...
            valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, ...
            m_total_follow_up_time, nTp, dtype_label, m_gtv_vol, config_struct, ...
            summary_metrics, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_survival', 'skipped'); end
        fprintf('⏭️ [5.5/5] [%s] Skipping metrics_survival.\n', current_name);
    end

    fprintf('=======================================================\n');
    fprintf('🎉 Pipeline Execution Complete for parameter %s\n', current_name);
    fprintf('=======================================================\n');

    if log_fid > 0
        fprintf(log_fid, '[%s] ===== Pipeline run completed (type: %s) =====\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'), current_name);
    end

    diary off;  % close master diary at end of pipeline run
end

%% ===== Local step wrapper functions =====
% Each wrapper encapsulates all the logic that was previously inline in the
% pipeline body.  They are called as function handles by execute_pipeline_step.

function run_compare_cores_step(validated_data_gtvp, summary_metrics, config_struct, current_name)
    fprintf('\n⚙️ [%s] Running core method comparison...\n', current_name);
    compare_results = compare_core_methods(validated_data_gtvp, summary_metrics, config_struct); %#ok<NASGU>
    fprintf('      ✅ Done.\n');
end

function run_metrics_longitudinal_step(ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
    nTp, dtype_label, config_struct, m_lf, current_name)
    fprintf('⚙️ [5.2/5] [%s] Running metrics_longitudinal...\n', current_name);
    metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
                         nTp, dtype_label, config_struct.output_folder, m_lf);
    longitudinal_results_file = fullfile(config_struct.output_folder, sprintf('metrics_longitudinal_results_%s.txt', current_name));
    fid = fopen(longitudinal_results_file, 'w');
    if fid < 0
        warning('run_dwi_pipeline:fileWriteFailed', 'Cannot write %s', longitudinal_results_file);
    else
        fprintf(fid, 'Longitudinal metrics generated successfully.\n');
        fclose(fid);
    end
    fprintf('      💾 Saved longitudinal results log to %s\n', longitudinal_results_file);
    fprintf('      ✅ Done.\n');
end

function run_metrics_dosimetry_step(m_id_list, summary_metrics, nTp, config_struct, ...
    m_data_vectors_gtvp, dosimetry_results_file, current_name)
    fprintf('⚙️ [5.3/5] [%s] Running metrics_dosimetry...\n', current_name);
    if config_struct.run_all_core_methods
        [d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub, per_method_dosimetry] = ...
            metrics_dosimetry(m_id_list, summary_metrics.id_list, nTp, config_struct, ...
                              m_data_vectors_gtvp, summary_metrics.gtv_locations);
    else
        [d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub] = ...
            metrics_dosimetry(m_id_list, summary_metrics.id_list, nTp, config_struct, ...
                              m_data_vectors_gtvp, summary_metrics.gtv_locations);
        per_method_dosimetry = struct();
    end

    save(dosimetry_results_file, 'd95_adc_sub', 'v50_adc_sub', 'd95_d_sub', 'v50_d_sub', 'd95_f_sub', 'v50_f_sub', 'd95_dstar_sub', 'v50_dstar_sub', 'per_method_dosimetry');
    fprintf('      ✅ Done.\n');
end

function run_metrics_stats_comparisons_step(valid_pts, lf_group, metric_sets, set_names, ...
    time_labels, dtype_label, config_struct, nTp, ADC_abs, D_abs, f_abs, Dstar_abs, current_name)
    fprintf('⚙️ [5.4a/5] [%s] Running metrics_stats_comparisons...\n', current_name);
    metrics_stats_comparisons(valid_pts, lf_group, ...
        metric_sets, set_names, time_labels, dtype_label, config_struct.output_folder, config_struct.dataloc, nTp, ...
        ADC_abs, D_abs, f_abs, Dstar_abs);

    comparisons_results_file = fullfile(config_struct.output_folder, sprintf('metrics_stats_comparisons_results_%s.txt', current_name));
    fid = fopen(comparisons_results_file, 'w');
    if fid < 0
        warning('run_dwi_pipeline:fileWriteFailed', 'Cannot write %s', comparisons_results_file);
    else
        fprintf(fid, 'Stats Comparisons generated successfully.\n');
        fclose(fid);
    end
    fprintf('      💾 Saved comparisons results log to %s\n', comparisons_results_file);
    fprintf('      ✅ Done.\n');
end

function run_metrics_stats_predictive_step(valid_pts, lf_group, dtype_label, config_struct, nTp, ...
    m_gtv_vol, summary_metrics, m_id_list, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, ...
    f_delta, Dstar_pct, m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, ...
    d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub, dl_provenance, time_labels, ...
    m_lf, m_total_time, m_total_follow_up_time, predictive_results_file, results_file, current_name)
    fprintf('⚙️ [5.4b/5] [%s] Running metrics_stats_predictive...\n', current_name);
    % Subset adc_sd to baseline-valid patients to match the
    % dimensions of m_id_list, ADC_abs, etc.  summary_metrics.adc_sd
    % has N rows (all patients) but valid_pts has M elements
    % (baseline-valid patients).  Without this mapping, logical
    % indexing by valid_pts silently selects wrong rows.
    [~, adc_sd_valid_idx] = ismember(m_id_list, summary_metrics.id_list);
    m_adc_sd = summary_metrics.adc_sd(adc_sd_valid_idx, :, :);

    [risk_scores_all, is_high_risk, times_km, events_km] = metrics_stats_predictive(valid_pts, lf_group, ...
        dtype_label, config_struct.output_folder, config_struct.dataloc, nTp, ...
        m_gtv_vol, m_adc_sd, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
        m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, ...
        d95_dstar_sub, v50_dstar_sub, m_id_list, config_struct.dwi_types_to_run, ...
        dl_provenance, time_labels, m_lf, m_total_time, m_total_follow_up_time, config_struct);

    save(predictive_results_file, 'risk_scores_all', 'is_high_risk', 'times_km', 'events_km');

    % Build a calculated_results struct for any downstream visualize_results
    % steps. This struct packages the predictive model outputs so that
    % visualization can overlay risk stratification on parameter maps —
    % e.g., highlighting high-risk patients' diffusion trajectories in
    % red and low-risk in blue on longitudinal plots, or annotating
    % Kaplan-Meier curves with model-predicted risk groups.
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
end

function run_visualize_step(validated_data_gtvp, summary_metrics, config_struct, results_file, current_name)
    fprintf('\n⚙️ [5.4c/5] [%s] Visualizing results...\n', current_name);
    % Load existing results if visualize is run independently
    if exist(results_file, 'file')
        tmp_results = load(results_file, 'calculated_results');
        calculated_results = tmp_results.calculated_results;
        fprintf('      💾 Loaded calculated_results from disk for visualization.\n');
    else
        calculated_results = struct(); % Empty struct fallback
    end
    visualize_results(validated_data_gtvp, summary_metrics, calculated_results, config_struct);
    visualize_results_file = fullfile(config_struct.output_folder, sprintf('visualize_results_state_%s.txt', current_name));
    fid = fopen(visualize_results_file, 'w');
    if fid < 0
        warning('run_dwi_pipeline:fileWriteFailed', 'Cannot write %s', visualize_results_file);
    else
        fprintf(fid, 'Visualizations generated successfully for: %s\n', current_name);
        fclose(fid);
    end
    fprintf('      ✅ Done: Visualizations generated and state saved to %s.\n', visualize_results_file);
end

function run_metrics_survival_step(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, ...
    m_total_follow_up_time, nTp, dtype_label, m_gtv_vol, config_struct, summary_metrics, current_name)
    fprintf('⚙️ [5.5/5] [%s] Running metrics_survival...\n', current_name);
    % Derive actual scan days from DICOM StudyDate headers stored
    % in summary_metrics.fx_dates (patients x fractions cell matrix
    % of 'YYYYMMDD' strings).  Fall back to config.json td_scan_days
    % if fx_dates are unavailable, then to built-in defaults in
    % metrics_survival (which emits an immortal-time-bias warning).
    td_scan_days_cfg = [];
    if isfield(summary_metrics, 'fx_dates') && ~isempty(summary_metrics.fx_dates)
        n_dates = sum(~cellfun('isempty', summary_metrics.fx_dates(:)));
        fprintf('      💡 Found %d DICOM StudyDate entries across cohort.\n', n_dates);
        td_scan_days_cfg = compute_scan_days_from_dates(summary_metrics.fx_dates);
        if isempty(td_scan_days_cfg)
            fprintf('      ⚠️  Could not derive scan days from DICOM dates (insufficient valid dates or non-monotonic).\n');
        end
    else
        fprintf('      💡 No DICOM StudyDate data available (fx_dates empty).\n');
    end
    if isempty(td_scan_days_cfg) && isfield(config_struct, 'td_scan_days') && ~isempty(config_struct.td_scan_days)
        td_scan_days_cfg = config_struct.td_scan_days;
        fprintf('      💡 Using td_scan_days from config.json.\n');
    end
    if isempty(td_scan_days_cfg)
        fprintf('      💡 Set "td_scan_days" in config.json to override defaults.\n');
    end
    metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, ...
                     m_total_follow_up_time, nTp, 'Survival', dtype_label, m_gtv_vol, config_struct.output_folder, td_scan_days_cfg);

    survival_results_file = fullfile(config_struct.output_folder, sprintf('metrics_survival_results_%s.txt', current_name));
    fid = fopen(survival_results_file, 'w');
    if fid < 0
        warning('run_dwi_pipeline:fileWriteFailed', 'Cannot write %s', survival_results_file);
    else
        fprintf(fid, 'Survival metrics generated successfully.\n');
        fclose(fid);
    end
    fprintf('      💾 Saved survival results log to %s\n', survival_results_file);
    fprintf('      ✅ Done.\n');
end

%% ===== Utility functions (remain in orchestrator) =====

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
