function run_visualize_step(validated_data_gtvp, summary_metrics, config_struct, results_file, current_name)
    fprintf('\n⚙️ [5.4c/5] [%s] Visualizing results...\n', current_name);
    if exist(results_file, 'file')
        tmp_results = load(results_file, 'calculated_results');
        calculated_results = tmp_results.calculated_results;
        fprintf('      💾 Loaded calculated_results from disk for visualization.\n');
    else
        calculated_results = struct();
    end
    visualize_results(validated_data_gtvp, summary_metrics, calculated_results, config_struct);
    write_sentinel_file(config_struct.output_folder, 'visualize_results_state', ...
        sprintf('Visualizations generated successfully for: %s', current_name), current_name);
    fprintf('      ✅ Done.\n');
end
