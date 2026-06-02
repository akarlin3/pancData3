function run_gtv_confounding_step(dosimetry_results, baseline_results, config_struct, current_name)
    fprintf('\n⚙️ [%s] Running GTV confounding check...\n', current_name);
    active_methods = config_struct.active_core_methods;
    confound = compute_gtv_confounding( ...
        dosimetry_results.per_method_dosimetry, baseline_results, ...
        active_methods, config_struct); %#ok<NASGU>
    save(fullfile(config_struct.output_folder, ...
        sprintf('gtv_confounding_%s.mat', config_struct.dwi_type_name)), 'confound');
    fprintf('      ✅ Done.\n');
end
