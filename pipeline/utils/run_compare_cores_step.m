function run_compare_cores_step(validated_data_gtvp, summary_metrics, config_struct, current_name)
    fprintf('\n⚙️ [%s] Running core method comparison...\n', current_name);
    compare_results = compare_core_methods(validated_data_gtvp, summary_metrics, config_struct); %#ok<NASGU>
    fprintf('      ✅ Done.\n');
end
