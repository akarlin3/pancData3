function run_core_failure_rates_step(validated_data_gtvp, summary_metrics, config_struct, current_name)
    fprintf('\n⚙️ [%s] Computing core method failure rates...\n', current_name);
    failure_table = compute_core_failure_rates(validated_data_gtvp, config_struct, ...
        summary_metrics.id_list, summary_metrics.gtv_locations);
    save(fullfile(config_struct.output_folder, ...
        sprintf('core_failure_rates_%s.mat', config_struct.dwi_type_name)), 'failure_table');
    fprintf('      ✅ Done.\n');
end
