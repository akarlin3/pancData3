function run_core_method_outcomes_step(dosimetry_results, baseline_results, config_struct, current_name)
    fprintf('\n⚙️ [%s] Running core method outcome analysis...\n', current_name);

    if ~isfield(dosimetry_results, 'per_method_dosimetry') || isempty(fieldnames(dosimetry_results.per_method_dosimetry))
        fprintf('      ⚠️ per_method_dosimetry not available. Set run_all_core_methods=true. Skipping.\n');
        return;
    end

    outcome_results = analyze_core_method_outcomes( ...
        dosimetry_results.per_method_dosimetry, baseline_results, ...
        config_struct.active_core_methods, config_struct);

    save(fullfile(config_struct.output_folder, ...
        sprintf('core_method_outcomes_%s.mat', config_struct.dwi_type_name)), 'outcome_results');
    fprintf('      ✅ Done.\n');
end
