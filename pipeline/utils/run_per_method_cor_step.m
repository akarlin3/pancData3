function run_per_method_cor_step(validated_data_gtvp, summary_metrics, config_struct, current_name)
    fprintf('\n⚙️ [%s] Computing per-method Coefficient of Reproducibility...\n', current_name);
    cor_results = compute_per_method_cor(validated_data_gtvp, config_struct, ...
        summary_metrics.id_list, summary_metrics.gtv_locations); %#ok<NASGU>
    save(fullfile(config_struct.output_folder, ...
        sprintf('per_method_cor_%s.mat', config_struct.dwi_type_name)), 'cor_results');
    fprintf('      ✅ Done.\n');
end
