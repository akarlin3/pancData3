function run_subvolume_stability_step(validated_data_gtvp, summary_metrics, config_struct, current_name)
    fprintf('\n⚙️ [%s] Computing sub-volume stability over fractions...\n', current_name);
    stability = compute_subvolume_stability(validated_data_gtvp, config_struct, ...
        summary_metrics.id_list, summary_metrics.gtv_locations); %#ok<NASGU>
    save(fullfile(config_struct.output_folder, ...
        sprintf('subvolume_stability_%s.mat', config_struct.dwi_type_name)), 'stability');
    fprintf('      ✅ Done.\n');
end
