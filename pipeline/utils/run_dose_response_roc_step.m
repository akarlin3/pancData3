function run_dose_response_roc_step(dosimetry_results, baseline_results, config_struct, current_name)
    fprintf('\n⚙️ [%s] Running dose-response ROC analysis...\n', current_name);
    active_methods = config_struct.active_core_methods;
    roc_results = compute_dose_response_roc( ...
        dosimetry_results.per_method_dosimetry, baseline_results, ...
        active_methods, config_struct); %#ok<NASGU>
    save(fullfile(config_struct.output_folder, ...
        sprintf('dose_response_roc_%s.mat', config_struct.dwi_type_name)), 'roc_results');
    fprintf('      ✅ Done.\n');
end
