function dispatch_pipeline_steps(session, validated_data_gtvp, validated_data_gtvn, summary_metrics)
% DISPATCH_PIPELINE_STEPS  Execute all analysis pipeline steps in sequence.
%
%   dispatch_pipeline_steps(session, validated_data_gtvp, validated_data_gtvn, summary_metrics)
%
%   Runs the analysis portion of the DWI pipeline: metrics_baseline through
%   metrics_survival and visualize_results.  Each step is conditionally
%   executed based on session.steps_to_run.  Data flows between steps via
%   a baseline_results struct and a dosimetry_results struct rather than
%   individual variables.
%
%   Non-fatal steps use execute_pipeline_step for uniform error handling.
%   Fatal steps (metrics_baseline) use inline try-catch and halt the pipeline.
%
%   Inputs:
%     session                - Session struct from prepare_pipeline_session
%     validated_data_gtvp    - Validated GTVp data vectors
%     validated_data_gtvn    - Validated GTVn data vectors
%     summary_metrics        - Summary metrics struct

    steps_to_run = session.steps_to_run;
    config_struct = session.config_struct;
    pipeGUI = session.pipeGUI;
    log_fid = session.log_fid;
    master_diary_file = session.master_diary_file;
    type_output_folder = session.type_output_folder;
    current_name = session.current_name;

    % =====================================================================
    % Metrics Baseline
    % =====================================================================
    % metrics_baseline is the analytical foundation: aggregates voxel-wise
    % parameters into patient-level summary statistics, computes percent
    % change from baseline, and identifies patients with complete data.
    baseline_results = struct();
    if ismember('metrics_baseline', steps_to_run)
        if ~isempty(pipeGUI), pipeGUI.startStep('metrics_baseline'); end
        try
            fprintf('\n\xe2\x9a\x99\xef\xb8\x8f [5.1/5] [%s] Running metrics_baseline...\n', current_name);
            baseline = metrics_baseline(validated_data_gtvp, validated_data_gtvn, summary_metrics, config_struct);

            save(session.baseline_results_file, '-struct', 'baseline');
            fprintf('      \xe2\x9c\x85 Done.\n');
            if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_baseline', 'success'); end
            [warn_msg, warn_id] = lastwarn;
            lastwarn('');
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
            try
                baseline = load_baseline_from_disk(session.baseline_results_file);
            catch ME_base
                fprintf('❌ %s\n', ME_base.message);
                fprintf('❌ Downstream metrics steps require baseline results. Halting pipeline.\n');
                if log_fid > 0
                    fprintf(log_fid, '[%s] [ERROR] %s. Halting pipeline.\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME_base.message);
                end
                return;
            end
        end
    end

    % Populate baseline_results from the loaded/computed baseline struct so
    % that downstream helper functions can access fields via baseline_results.
    if exist('baseline', 'var')
        baseline_results = baseline;
    end

    diary(master_diary_file);
    lastwarn('');

    % =====================================================================
    % Compare Core Methods
    % =====================================================================
    if ismember('compare_cores', steps_to_run)
        execute_pipeline_step('compare_cores', @() run_compare_cores_step( ...
            validated_data_gtvp, summary_metrics, config_struct, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('compare_cores', 'skipped'); end
    end

    % =====================================================================
    % Cross-Pipeline Dice (Fraction 1)
    % =====================================================================
    if ismember('cross_pipeline_dice', steps_to_run)
        execute_pipeline_step('cross_pipeline_dice', @() run_cross_pipeline_dice_step( ...
            validated_data_gtvp, summary_metrics, config_struct, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('cross_pipeline_dice', 'skipped'); end
    end

    % =====================================================================
    % Core Method Failure Rates
    % =====================================================================
    if ismember('core_failure_rates', steps_to_run)
        execute_pipeline_step('core_failure_rates', @() run_core_failure_rates_step( ...
            validated_data_gtvp, summary_metrics, config_struct, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('core_failure_rates', 'skipped'); end
    end

    % =====================================================================
    % Core Method Pruning
    % =====================================================================
    % If failure rates were computed and a pruning threshold is set,
    % filter the method list and store in config_struct for downstream steps.
    ALL_CORE_METHODS_DEFAULT = {'adc_threshold', 'd_threshold', 'df_intersection', ...
        'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
        'percentile', 'spectral', 'fdm'};

    failure_rates_file = fullfile(type_output_folder, ...
        sprintf('core_failure_rates_%s.mat', config_struct.dwi_type_name));

    max_fail_rate = 1.0;
    if isfield(config_struct, 'max_core_failure_rate')
        max_fail_rate = config_struct.max_core_failure_rate;
    end
    excl_methods = {};
    if isfield(config_struct, 'excluded_core_methods')
        excl_methods = config_struct.excluded_core_methods;
    end
    if max_fail_rate < 1.0 || ~isempty(excl_methods)
        if exist(failure_rates_file, 'file')
            loaded = load(failure_rates_file, 'failure_table');
            [active_methods, pruned_info, retained_with_warning] = filter_core_methods( ...
                ALL_CORE_METHODS_DEFAULT, loaded.failure_table, config_struct);
            config_struct.active_core_methods = active_methods;

            % Log pruning results
            if ~isempty(pruned_info)
                fprintf('\n🪓 Core method pruning: %d/%d methods removed\n', ...
                    numel(pruned_info), numel(ALL_CORE_METHODS_DEFAULT));
                for pi = 1:numel(pruned_info)
                    fprintf('   ✂ %s — %s (%.1f%% failure rate)\n', ...
                        pruned_info(pi).name, pruned_info(pi).reason, ...
                        pruned_info(pi).failure_rate * 100);
                end
                fprintf('   ✅ %d methods retained: %s\n', ...
                    numel(active_methods), strjoin(active_methods, ', '));
            end

            % Log retained-with-warning methods
            if ~isempty(retained_with_warning)
                for wi = 1:numel(retained_with_warning)
                    fprintf('   ⚠️ %s — %s\n', ...
                        retained_with_warning(wi).name, retained_with_warning(wi).reason);
                end
            end

            % Save pruning results (including retained_with_warning)
            min_core_voxels_used = config_struct.min_core_voxels;
            save(fullfile(type_output_folder, ...
                sprintf('core_pruning_%s.mat', config_struct.dwi_type_name)), ...
                'active_methods', 'pruned_info', 'retained_with_warning', 'min_core_voxels_used');
        else
            fprintf('⚠️ Core failure rates not found at %s. Skipping pruning.\n', failure_rates_file);
            config_struct.active_core_methods = ALL_CORE_METHODS_DEFAULT;
        end
    else
        config_struct.active_core_methods = ALL_CORE_METHODS_DEFAULT;
    end

    % =====================================================================
    % Metrics Longitudinal
    % =====================================================================
    if ismember('metrics_longitudinal', steps_to_run)
        execute_pipeline_step('metrics_longitudinal', @() run_metrics_longitudinal_step( ...
            session, baseline_results), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_longitudinal', 'skipped'); end
        fprintf('⏭️ [5.2/5] [%s] Skipping metrics_longitudinal.\n', current_name);
    end

    % =====================================================================
    % Metrics Dosimetry
    % =====================================================================
    dosimetry_results = struct();
    if ismember('metrics_dosimetry', steps_to_run)
        execute_pipeline_step('metrics_dosimetry', @() run_metrics_dosimetry_step( ...
            session, baseline_results, summary_metrics), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
        if exist(session.dosimetry_results_file, 'file')
            dosimetry_results = load(session.dosimetry_results_file);
        end
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_dosimetry', 'skipped'); end
        if any(ismember({'metrics_stats_comparisons', 'metrics_stats_predictive'}, steps_to_run))
            fprintf('⏭️ [5.3/5] [%s] Skipping metrics_dosimetry. Loading from disk...\n', current_name);
            if exist(session.dosimetry_results_file, 'file')
                dosimetry_results = load(session.dosimetry_results_file);
            else
                fprintf('      ⚠️ Warning: metrics_dosimetry results not found. metrics_stats may fail.\n');
                if log_fid > 0
                    fprintf(log_fid, '[%s] [WARNING] metrics_dosimetry results not found at: %s\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), session.dosimetry_results_file);
                end
                dosimetry_results = make_empty_dosimetry(baseline_results);
            end
        else
            fprintf('⏭️ [5.3/5] [%s] Skipping metrics_dosimetry.\n', current_name);
        end
    end

    % Append dosimetry sub-volume metrics to metric_sets for univariate analysis
    if isfield(dosimetry_results, 'd95_adc_sub') && isfield(baseline_results, 'metric_sets')
        baseline_results.metric_sets{3} = {baseline_results.m_d95_gtvp, dosimetry_results.d95_adc_sub, dosimetry_results.d95_d_sub, dosimetry_results.d95_f_sub, dosimetry_results.d95_dstar_sub};
        baseline_results.metric_sets{4} = {baseline_results.m_v50gy_gtvp, dosimetry_results.v50_adc_sub, dosimetry_results.v50_d_sub, dosimetry_results.v50_f_sub, dosimetry_results.v50_dstar_sub};
        baseline_results.set_names{3} = {'D95 GTVp (whole)', 'D95 Sub(ADC)', 'D95 Sub(D)', 'D95 Sub(f)', 'D95 Sub(D*)'};
        baseline_results.set_names{4} = {'V50 GTVp (whole)', 'V50 Sub(ADC)', 'V50 Sub(D)', 'V50 Sub(f)', 'V50 Sub(D*)'};
    end

    % =====================================================================
    % Dose vs Delta D scatter plots (WU3)
    % =====================================================================
    if isfield(baseline_results, 'D_pct') && ~isempty(fieldnames(dosimetry_results))
        try
            plot_dose_vs_delta(baseline_results, dosimetry_results, config_struct);
        catch ME_delta
            fprintf('  \xe2\x9a\xa0\xef\xb8\x8f plot_dose_vs_delta failed: %s\n', ME_delta.message);
        end
    end

    % =====================================================================
    % ADC Threshold Optimization (WU4)
    % =====================================================================
    % Sweeps ADC threshold 0.8e-3..2.0e-3 and finds the value maximizing
    % Fx1 repeat Dice. Output is saved to type_output_folder (DWI subfolder)
    % as adc_threshold_optimization_{dwi}.mat (variable: opt_results).
    if isfield(config_struct, 'run_optimize_threshold') && config_struct.run_optimize_threshold
        try
            opt_results = optimize_adc_threshold(validated_data_gtvp, config_struct, ...
                summary_metrics.id_list, summary_metrics.gtv_locations); %#ok<NASGU>
            save(fullfile(type_output_folder, ...
                sprintf('adc_threshold_optimization_%s.mat', current_name)), 'opt_results');
        catch ME_opt
            fprintf('  \xe2\x9a\xa0\xef\xb8\x8f optimize_adc_threshold failed: %s\n', ME_opt.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] optimize_adc_threshold failed: %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME_opt.message);
            end
        end
    end

    % =====================================================================
    % Stats Comparisons
    % =====================================================================
    if ismember('metrics_stats_comparisons', steps_to_run)
        execute_pipeline_step('metrics_stats_comparisons', @() run_metrics_stats_comparisons_step( ...
            session, baseline_results), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_stats_comparisons', 'skipped'); end
        fprintf('⏭️ [5.4a/5] [%s] Skipping metrics_stats_comparisons.\n', current_name);
    end

    % =====================================================================
    % Stats Predictive
    % =====================================================================
    if ismember('metrics_stats_predictive', steps_to_run)
        execute_pipeline_step('metrics_stats_predictive', @() run_metrics_stats_predictive_step( ...
            session, baseline_results, dosimetry_results, summary_metrics), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
        if exist(session.results_file, 'file')
            tmp_results = load(session.results_file, 'calculated_results');
            calculated_results = tmp_results.calculated_results; %#ok<NASGU>
        end
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_stats_predictive', 'skipped'); end
        if ismember('metrics_survival', steps_to_run)
            fprintf('⏭️ [5.4b/5] [%s] Skipping metrics_stats_predictive. Loading from disk...\n', current_name);
            if exist(session.predictive_results_file, 'file')
                tmp_pred = load(session.predictive_results_file);
                risk_scores_all = tmp_pred.risk_scores_all; is_high_risk = tmp_pred.is_high_risk; times_km = tmp_pred.times_km; events_km = tmp_pred.events_km; %#ok<NASGU>
            else
                fprintf('      ⚠️ Warning: metrics_stats_predictive results not found. metrics_survival will fail.\n');
                if log_fid > 0
                    fprintf(log_fid, '[%s] [WARNING] metrics_stats_predictive results not found at: %s\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), session.predictive_results_file);
                end
            end
        else
            fprintf('⏭️ [5.4b/5] [%s] Skipping metrics_stats_predictive.\n', current_name);
        end
    end

    % =====================================================================
    % External Validation (export + apply)
    % =====================================================================
    % Runs after metrics_stats_predictive so the trained model is available.
    % Both sub-steps are non-fatal: failures are logged and the pipeline continues.
    if ismember('metrics_stats_predictive', steps_to_run)
        validation_model_path = fullfile(config_struct.output_folder, ...
            sprintf('validation_model_%s.mat', current_name));

        % --- Export validation model ---
        if config_struct.export_validation_model
            execute_pipeline_step('export_validation_model', ...
                @() run_export_validation_step(session, validation_model_path), ...
                pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
        end

        % --- Apply saved model to external dataset ---
        ext_data_path = config_struct.external_validation_data;
        if ~isempty(ext_data_path) && ischar(ext_data_path)
            if ~exist(ext_data_path, 'dir')
                fprintf('\xe2\x9a\xa0\xef\xb8\x8f External validation data file not found: %s\n', ext_data_path);
                if log_fid > 0
                    fprintf(log_fid, '[%s] [WARNING] External validation data not found: %s\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), ext_data_path);
                end
            else
                execute_pipeline_step('apply_external_validation', ...
                    @() run_apply_external_validation_step(session, ext_data_path, validation_model_path), ...
                    pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
            end
        end
    end

    % =====================================================================
    % Core Method Outcome Analysis
    % =====================================================================
    if ismember('core_method_outcomes', steps_to_run) && ...
            isfield(config_struct, 'active_core_methods') && ...
            ~isempty(fieldnames(dosimetry_results)) && ...
            isfield(dosimetry_results, 'per_method_dosimetry')
        execute_pipeline_step('core_method_outcomes', @() run_core_method_outcomes_step( ...
            dosimetry_results, baseline_results, config_struct, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('core_method_outcomes', 'skipped'); end
    end

    % =====================================================================
    % Per-Method Coefficient of Reproducibility (Prompt D)
    % =====================================================================
    if ismember('per_method_cor', steps_to_run)
        execute_pipeline_step('per_method_cor', @() run_per_method_cor_step( ...
            validated_data_gtvp, summary_metrics, config_struct, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('per_method_cor', 'skipped'); end
    end

    % =====================================================================
    % Sub-Volume Stability Over Fractions (Prompt A)
    % =====================================================================
    if ismember('subvolume_stability', steps_to_run)
        execute_pipeline_step('subvolume_stability', @() run_subvolume_stability_step( ...
            validated_data_gtvp, summary_metrics, config_struct, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('subvolume_stability', 'skipped'); end
    end

    % =====================================================================
    % Dose-Response ROC (Prompt B)
    % =====================================================================
    if ismember('dose_response_roc', steps_to_run) && ...
            ~isempty(fieldnames(dosimetry_results)) && ...
            isfield(dosimetry_results, 'per_method_dosimetry')
        execute_pipeline_step('dose_response_roc', @() run_dose_response_roc_step( ...
            dosimetry_results, baseline_results, config_struct, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('dose_response_roc', 'skipped'); end
    end

    % =====================================================================
    % GTV Volume Confounding Check (Prompt E)
    % =====================================================================
    if ismember('gtv_confounding', steps_to_run) && ...
            ~isempty(fieldnames(dosimetry_results)) && ...
            isfield(dosimetry_results, 'per_method_dosimetry')
        execute_pipeline_step('gtv_confounding', @() run_gtv_confounding_step( ...
            dosimetry_results, baseline_results, config_struct, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('gtv_confounding', 'skipped'); end
    end

    % =====================================================================
    % Risk-Dose Concordance (Prompt C)
    % =====================================================================
    if ismember('risk_dose_concordance', steps_to_run)
        predictive_results_loaded = struct();
        if exist(session.predictive_results_file, 'file')
            predictive_results_loaded = load(session.predictive_results_file);
        end
        if isfield(predictive_results_loaded, 'is_high_risk') && ...
                ~isempty(fieldnames(dosimetry_results)) && ...
                isfield(dosimetry_results, 'per_method_dosimetry')
            execute_pipeline_step('risk_dose_concordance', @() run_risk_dose_concordance_step( ...
                predictive_results_loaded, dosimetry_results, baseline_results, config_struct, current_name), ...
                pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
        else
            if ~isempty(pipeGUI), pipeGUI.completeStep('risk_dose_concordance', 'skipped'); end
        end
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('risk_dose_concordance', 'skipped'); end
    end

    % =====================================================================
    % Visualize Results
    % =====================================================================
    if ismember('visualize', steps_to_run)
        execute_pipeline_step('visualize', @() run_visualize_step( ...
            validated_data_gtvp, summary_metrics, config_struct, session.results_file, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('visualize', 'skipped'); end
        fprintf('⏭️ [5.4c/5] [%s] Skipping Visualization.\n', current_name);
    end

    % =====================================================================
    % Survival Analysis
    % =====================================================================
    if ismember('metrics_survival', steps_to_run)
        execute_pipeline_step('metrics_survival', @() run_metrics_survival_step( ...
            session, baseline_results, summary_metrics), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_survival', 'skipped'); end
        fprintf('⏭️ [5.5/5] [%s] Skipping metrics_survival.\n', current_name);
    end

    % =====================================================================
    % Baseline vs Delta Cox PH comparison (WU7)
    % =====================================================================
    if isfield(config_struct, 'run_baseline_vs_delta') && config_struct.run_baseline_vs_delta ...
            && isfield(baseline_results, 'D_pct')
        try
            comparison = compare_baseline_vs_delta(baseline_results, config_struct); %#ok<NASGU>
            save(fullfile(config_struct.output_folder, ...
                sprintf('baseline_vs_delta_%s.mat', current_name)), 'comparison');
        catch ME_bvd
            fprintf('  \xe2\x9a\xa0\xef\xb8\x8f compare_baseline_vs_delta failed: %s\n', ME_bvd.message);
            if log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] compare_baseline_vs_delta failed: %s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME_bvd.message);
            end
        end
    end
end
