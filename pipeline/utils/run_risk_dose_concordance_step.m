function run_risk_dose_concordance_step(predictive_results, dosimetry_results, baseline_results, config_struct, current_name)
    fprintf('\n⚙️ [%s] Running risk-dose concordance analysis...\n', current_name);
    active_methods = config_struct.active_core_methods;
    concordance = compute_risk_dose_concordance( ...
        predictive_results, dosimetry_results.per_method_dosimetry, ...
        baseline_results, active_methods, config_struct); %#ok<NASGU>
    save(fullfile(config_struct.output_folder, ...
        sprintf('risk_dose_concordance_%s.mat', config_struct.dwi_type_name)), 'concordance');
    fprintf('      ✅ Done.\n');
end
