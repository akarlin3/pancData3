function run_metrics_survival_step(session, baseline_results, summary_metrics)
    fprintf('⚙️ [5.5/5] [%s] Running metrics_survival...\n', session.current_name);
    config_struct = session.config_struct;
    td_scan_days_cfg = resolve_scan_days(summary_metrics, config_struct);
    metrics_survival(baseline_results.valid_pts, baseline_results.ADC_abs, baseline_results.D_abs, ...
        baseline_results.f_abs, baseline_results.Dstar_abs, baseline_results.m_lf, ...
        baseline_results.m_total_time, baseline_results.m_total_follow_up_time, ...
        baseline_results.nTp, 'Survival', baseline_results.dtype_label, ...
        baseline_results.m_gtv_vol, config_struct.output_folder, td_scan_days_cfg);

    survival_results_file = fullfile(config_struct.output_folder, sprintf('metrics_survival_results_%s.txt', session.current_name));
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
