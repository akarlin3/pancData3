function dispatch_pipeline_steps(session, validated_data_gtvp, validated_data_gtvn, summary_metrics)
% DISPATCH_PIPELINE_STEPS  Execute all analysis pipeline steps in sequence.
%
%   dispatch_pipeline_steps(session, validated_data_gtvp, validated_data_gtvn, summary_metrics)
%
%   Runs the analysis portion of the DWI pipeline: metrics_baseline through
%   metrics_survival and visualize_results.  Each step is conditionally
%   executed based on session.steps_to_run.  Data flows between steps via
%   local variables (metrics_baseline outputs feed into downstream steps).
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
    current_dtype = session.current_dtype;
    baseline_results_file = session.baseline_results_file;
    dosimetry_results_file = session.dosimetry_results_file;
    predictive_results_file = session.predictive_results_file;
    results_file = session.results_file;

    % =====================================================================
    % Metrics Baseline
    % =====================================================================
    % metrics_baseline is the analytical foundation: aggregates voxel-wise
    % parameters into patient-level summary statistics, computes percent
    % change from baseline, and identifies patients with complete data.
    if ismember('metrics_baseline', steps_to_run)
        if ~isempty(pipeGUI), pipeGUI.startStep('metrics_baseline'); end
        try
            fprintf('\n\xe2\x9a\x99\xef\xb8\x8f [5.1/5] [%s] Running metrics_baseline...\n', current_name);
            [m_lf, m_total_time, m_total_follow_up_time, m_gtv_vol, m_adc_mean, m_d_mean, m_f_mean, m_dstar_mean, ...
             m_id_list, m_mrn_list, m_d95_gtvp, m_v50gy_gtvp, m_data_vectors_gtvp, lf_group, valid_pts, ...
             ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
             nTp, metric_sets, set_names, time_labels, dtype_label, dl_provenance] = ...
             metrics_baseline(validated_data_gtvp, validated_data_gtvn, summary_metrics, config_struct);

            save(baseline_results_file, 'm_lf', 'm_total_time', 'm_total_follow_up_time', 'm_gtv_vol', 'm_adc_mean', 'm_d_mean', 'm_f_mean', 'm_dstar_mean', ...
             'm_id_list', 'm_mrn_list', 'm_d95_gtvp', 'm_v50gy_gtvp', 'm_data_vectors_gtvp', 'lf_group', 'valid_pts', ...
             'ADC_abs', 'D_abs', 'f_abs', 'Dstar_abs', 'ADC_pct', 'D_pct', 'f_delta', 'Dstar_pct', ...
             'nTp', 'metric_sets', 'set_names', 'time_labels', 'dtype_label', 'dl_provenance');
            fprintf('      \xe2\x9c\x85 Done.\n');
            if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_baseline', 'success'); end
            [warn_msg, warn_id] = lastwarn;
            lastwarn('');
            if ~isempty(warn_msg) && log_fid > 0
                fprintf(log_fid, '[%s] [WARNING] During metrics_baseline: %s (id: %s)\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), warn_msg, warn_id);
            end
        catch ME
            fprintf('\xe2\x9d\x8c FAILED.\n');
            fprintf('\xe2\x9d\x8c Error during metrics_baseline: %s\n', ME.message);
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
            fprintf('\n\xe2\x8f\xad\xef\xb8\x8f [5.1/5] [%s] Skipping metrics_baseline. Loading from disk...\n', current_name);
            try
                tmp_base = load_baseline_from_disk(baseline_results_file);
                m_lf = tmp_base.m_lf; m_total_time = tmp_base.m_total_time; m_total_follow_up_time = tmp_base.m_total_follow_up_time; m_gtv_vol = tmp_base.m_gtv_vol; m_adc_mean = tmp_base.m_adc_mean; m_d_mean = tmp_base.m_d_mean; m_f_mean = tmp_base.m_f_mean; m_dstar_mean = tmp_base.m_dstar_mean;
                m_id_list = tmp_base.m_id_list; m_mrn_list = tmp_base.m_mrn_list; m_d95_gtvp = tmp_base.m_d95_gtvp; m_v50gy_gtvp = tmp_base.m_v50gy_gtvp; m_data_vectors_gtvp = tmp_base.m_data_vectors_gtvp; lf_group = tmp_base.lf_group; valid_pts = tmp_base.valid_pts;
                ADC_abs = tmp_base.ADC_abs; D_abs = tmp_base.D_abs; f_abs = tmp_base.f_abs; Dstar_abs = tmp_base.Dstar_abs; ADC_pct = tmp_base.ADC_pct; D_pct = tmp_base.D_pct; f_delta = tmp_base.f_delta; Dstar_pct = tmp_base.Dstar_pct;
                nTp = tmp_base.nTp; metric_sets = tmp_base.metric_sets; set_names = tmp_base.set_names; time_labels = tmp_base.time_labels; dtype_label = tmp_base.dtype_label; dl_provenance = tmp_base.dl_provenance;
            catch ME_base
                fprintf('\xe2\x9d\x8c %s\n', ME_base.message);
                fprintf('\xe2\x9d\x8c Downstream metrics steps require baseline results. Halting pipeline.\n');
                if log_fid > 0
                    fprintf(log_fid, '[%s] [ERROR] %s. Halting pipeline.\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), ME_base.message);
                end
                return;
            end
        end
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
    % Metrics Longitudinal
    % =====================================================================
    if ismember('metrics_longitudinal', steps_to_run)
        execute_pipeline_step('metrics_longitudinal', @() run_metrics_longitudinal_step( ...
            ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
            nTp, dtype_label, config_struct, m_lf, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_longitudinal', 'skipped'); end
        fprintf('\xe2\x8f\xad\xef\xb8\x8f [5.2/5] [%s] Skipping metrics_longitudinal.\n', current_name);
    end

    % =====================================================================
    % Metrics Dosimetry
    % =====================================================================
    if ismember('metrics_dosimetry', steps_to_run)
        execute_pipeline_step('metrics_dosimetry', @() run_metrics_dosimetry_step( ...
            m_id_list, summary_metrics, nTp, config_struct, m_data_vectors_gtvp, ...
            dosimetry_results_file, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
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
            fprintf('\xe2\x8f\xad\xef\xb8\x8f [5.3/5] [%s] Skipping metrics_dosimetry. Loading from disk...\n', current_name);
            if exist(dosimetry_results_file, 'file')
                tmp_dosimetry = load(dosimetry_results_file);
                d95_adc_sub = tmp_dosimetry.d95_adc_sub; v50_adc_sub = tmp_dosimetry.v50_adc_sub; d95_d_sub = tmp_dosimetry.d95_d_sub; v50_d_sub = tmp_dosimetry.v50_d_sub; d95_f_sub = tmp_dosimetry.d95_f_sub; v50_f_sub = tmp_dosimetry.v50_f_sub; d95_dstar_sub = tmp_dosimetry.d95_dstar_sub; v50_dstar_sub = tmp_dosimetry.v50_dstar_sub;
            else
                fprintf('      \xe2\x9a\xa0\xef\xb8\x8f Warning: metrics_dosimetry results not found. metrics_stats may fail.\n');
                if log_fid > 0
                    fprintf(log_fid, '[%s] [WARNING] metrics_dosimetry results not found at: %s\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), dosimetry_results_file);
                end
                nTp_val = nTp;
                d95_adc_sub = nan(length(m_id_list), nTp_val); v50_adc_sub = nan(length(m_id_list), nTp_val);
                d95_d_sub = nan(length(m_id_list), nTp_val); v50_d_sub = nan(length(m_id_list), nTp_val);
                d95_f_sub = nan(length(m_id_list), nTp_val); v50_f_sub = nan(length(m_id_list), nTp_val);
                d95_dstar_sub = nan(length(m_id_list), nTp_val); v50_dstar_sub = nan(length(m_id_list), nTp_val);
            end
        else
            fprintf('\xe2\x8f\xad\xef\xb8\x8f [5.3/5] [%s] Skipping metrics_dosimetry.\n', current_name);
        end
    end

    % Append dosimetry sub-volume metrics to metric_sets for univariate analysis
    if exist('d95_adc_sub', 'var') && exist('metric_sets', 'var')
        metric_sets{3} = {m_d95_gtvp, d95_adc_sub, d95_d_sub, d95_f_sub, d95_dstar_sub};
        metric_sets{4} = {m_v50gy_gtvp, v50_adc_sub, v50_d_sub, v50_f_sub, v50_dstar_sub};
        set_names{3} = {'D95 GTVp (whole)', 'D95 Sub(ADC)', 'D95 Sub(D)', 'D95 Sub(f)', 'D95 Sub(D*)'};
        set_names{4} = {'V50 GTVp (whole)', 'V50 Sub(ADC)', 'V50 Sub(D)', 'V50 Sub(f)', 'V50 Sub(D*)'};
    end

    % =====================================================================
    % Stats Comparisons
    % =====================================================================
    if ismember('metrics_stats_comparisons', steps_to_run)
        execute_pipeline_step('metrics_stats_comparisons', @() run_metrics_stats_comparisons_step( ...
            valid_pts, lf_group, metric_sets, set_names, time_labels, dtype_label, ...
            config_struct, nTp, ADC_abs, D_abs, f_abs, Dstar_abs, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_stats_comparisons', 'skipped'); end
        fprintf('\xe2\x8f\xad\xef\xb8\x8f [5.4a/5] [%s] Skipping metrics_stats_comparisons.\n', current_name);
    end

    % =====================================================================
    % Stats Predictive
    % =====================================================================
    if ismember('metrics_stats_predictive', steps_to_run)
        execute_pipeline_step('metrics_stats_predictive', @() run_metrics_stats_predictive_step( ...
            valid_pts, lf_group, dtype_label, config_struct, nTp, m_gtv_vol, summary_metrics, ...
            m_id_list, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
            m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, ...
            d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub, dl_provenance, time_labels, ...
            m_lf, m_total_time, m_total_follow_up_time, predictive_results_file, results_file, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
        if exist(results_file, 'file')
            tmp_results = load(results_file, 'calculated_results');
            calculated_results = tmp_results.calculated_results; %#ok<NASGU>
        end
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_stats_predictive', 'skipped'); end
        if ismember('metrics_survival', steps_to_run)
            fprintf('\xe2\x8f\xad\xef\xb8\x8f [5.4b/5] [%s] Skipping metrics_stats_predictive. Loading from disk...\n', current_name);
            if exist(predictive_results_file, 'file')
                tmp_pred = load(predictive_results_file);
                risk_scores_all = tmp_pred.risk_scores_all; is_high_risk = tmp_pred.is_high_risk; times_km = tmp_pred.times_km; events_km = tmp_pred.events_km; %#ok<NASGU>
            else
                fprintf('      \xe2\x9a\xa0\xef\xb8\x8f Warning: metrics_stats_predictive results not found. metrics_survival will fail.\n');
                if log_fid > 0
                    fprintf(log_fid, '[%s] [WARNING] metrics_stats_predictive results not found at: %s\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS'), predictive_results_file);
                end
            end
        else
            fprintf('\xe2\x8f\xad\xef\xb8\x8f [5.4b/5] [%s] Skipping metrics_stats_predictive.\n', current_name);
        end
    end

    % =====================================================================
    % Visualize Results
    % =====================================================================
    if ismember('visualize', steps_to_run)
        execute_pipeline_step('visualize', @() run_visualize_step( ...
            validated_data_gtvp, summary_metrics, config_struct, results_file, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('visualize', 'skipped'); end
        fprintf('\xe2\x8f\xad\xef\xb8\x8f [5.4c/5] [%s] Skipping Visualization.\n', current_name);
    end

    % =====================================================================
    % Survival Analysis
    % =====================================================================
    if ismember('metrics_survival', steps_to_run)
        execute_pipeline_step('metrics_survival', @() run_metrics_survival_step( ...
            valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, ...
            m_total_follow_up_time, nTp, dtype_label, m_gtv_vol, config_struct, ...
            summary_metrics, current_name), ...
            pipeGUI, log_fid, master_diary_file, type_output_folder, current_name);
    else
        if ~isempty(pipeGUI), pipeGUI.completeStep('metrics_survival', 'skipped'); end
        fprintf('\xe2\x8f\xad\xef\xb8\x8f [5.5/5] [%s] Skipping metrics_survival.\n', current_name);
    end
end

%% ===== Local step wrapper functions =====

function run_compare_cores_step(validated_data_gtvp, summary_metrics, config_struct, current_name)
    fprintf('\n\xe2\x9a\x99\xef\xb8\x8f [%s] Running core method comparison...\n', current_name);
    compare_results = compare_core_methods(validated_data_gtvp, summary_metrics, config_struct); %#ok<NASGU>
    fprintf('      \xe2\x9c\x85 Done.\n');
end

function run_metrics_longitudinal_step(ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
    nTp, dtype_label, config_struct, m_lf, current_name)
    fprintf('\xe2\x9a\x99\xef\xb8\x8f [5.2/5] [%s] Running metrics_longitudinal...\n', current_name);
    metrics_longitudinal(ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
                         nTp, dtype_label, config_struct.output_folder, m_lf);
    write_sentinel_file(config_struct.output_folder, 'metrics_longitudinal_results', ...
        'Longitudinal metrics generated successfully.', current_name);
    fprintf('      \xe2\x9c\x85 Done.\n');
end

function run_metrics_dosimetry_step(m_id_list, summary_metrics, nTp, config_struct, ...
    m_data_vectors_gtvp, dosimetry_results_file, current_name)
    fprintf('\xe2\x9a\x99\xef\xb8\x8f [5.3/5] [%s] Running metrics_dosimetry...\n', current_name);
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
    fprintf('      \xe2\x9c\x85 Done.\n');
end

function run_metrics_stats_comparisons_step(valid_pts, lf_group, metric_sets, set_names, ...
    time_labels, dtype_label, config_struct, nTp, ADC_abs, D_abs, f_abs, Dstar_abs, current_name)
    fprintf('\xe2\x9a\x99\xef\xb8\x8f [5.4a/5] [%s] Running metrics_stats_comparisons...\n', current_name);
    metrics_stats_comparisons(valid_pts, lf_group, ...
        metric_sets, set_names, time_labels, dtype_label, config_struct.output_folder, config_struct.dataloc, nTp, ...
        ADC_abs, D_abs, f_abs, Dstar_abs);
    write_sentinel_file(config_struct.output_folder, 'metrics_stats_comparisons_results', ...
        'Stats Comparisons generated successfully.', current_name);
    fprintf('      \xe2\x9c\x85 Done.\n');
end

function run_metrics_stats_predictive_step(valid_pts, lf_group, dtype_label, config_struct, nTp, ...
    m_gtv_vol, summary_metrics, m_id_list, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, ...
    f_delta, Dstar_pct, m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, ...
    d95_f_sub, v50_f_sub, d95_dstar_sub, v50_dstar_sub, dl_provenance, time_labels, ...
    m_lf, m_total_time, m_total_follow_up_time, predictive_results_file, results_file, current_name)
    fprintf('\xe2\x9a\x99\xef\xb8\x8f [5.4b/5] [%s] Running metrics_stats_predictive...\n', current_name);
    [~, adc_sd_valid_idx] = ismember(m_id_list, summary_metrics.id_list);
    m_adc_sd = summary_metrics.adc_sd(adc_sd_valid_idx, :, :);

    [risk_scores_all, is_high_risk, times_km, events_km] = metrics_stats_predictive(valid_pts, lf_group, ...
        dtype_label, config_struct.output_folder, config_struct.dataloc, nTp, ...
        m_gtv_vol, m_adc_sd, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, ...
        m_d95_gtvp, m_v50gy_gtvp, d95_adc_sub, v50_adc_sub, d95_d_sub, v50_d_sub, d95_f_sub, v50_f_sub, ...
        d95_dstar_sub, v50_dstar_sub, m_id_list, config_struct.dwi_types_to_run, ...
        dl_provenance, time_labels, m_lf, m_total_time, m_total_follow_up_time, config_struct);

    save(predictive_results_file, 'risk_scores_all', 'is_high_risk', 'times_km', 'events_km');

    calculated_results = struct();
    calculated_results.risk_scores_all = risk_scores_all;
    calculated_results.is_high_risk = is_high_risk;
    calculated_results.times_km = times_km;
    calculated_results.events_km = events_km;
    calculated_results.m_lf = m_lf;
    calculated_results.m_id_list = m_id_list;
    save(results_file, 'calculated_results');
    fprintf('      \xf0\x9f\x92\xbe Saved calculated_results to %s\n', results_file);
    fprintf('      \xe2\x9c\x85 Done.\n');
end

function run_visualize_step(validated_data_gtvp, summary_metrics, config_struct, results_file, current_name)
    fprintf('\n\xe2\x9a\x99\xef\xb8\x8f [5.4c/5] [%s] Visualizing results...\n', current_name);
    if exist(results_file, 'file')
        tmp_results = load(results_file, 'calculated_results');
        calculated_results = tmp_results.calculated_results;
        fprintf('      \xf0\x9f\x92\xbe Loaded calculated_results from disk for visualization.\n');
    else
        calculated_results = struct();
    end
    visualize_results(validated_data_gtvp, summary_metrics, calculated_results, config_struct);
    write_sentinel_file(config_struct.output_folder, 'visualize_results_state', ...
        sprintf('Visualizations generated successfully for: %s', current_name), current_name);
    fprintf('      \xe2\x9c\x85 Done.\n');
end

function run_metrics_survival_step(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, ...
    m_total_follow_up_time, nTp, dtype_label, m_gtv_vol, config_struct, summary_metrics, current_name)
    fprintf('\xe2\x9a\x99\xef\xb8\x8f [5.5/5] [%s] Running metrics_survival...\n', current_name);
    td_scan_days_cfg = resolve_scan_days(summary_metrics, config_struct);
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
    fprintf('      \xf0\x9f\x92\xbe Saved survival results log to %s\n', survival_results_file);
    fprintf('      \xe2\x9c\x85 Done.\n');
end
