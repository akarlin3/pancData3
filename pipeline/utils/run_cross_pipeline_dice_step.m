function run_cross_pipeline_dice_step(validated_data_gtvp, summary_metrics, config_struct, current_name)
    fprintf('\n⚙️ [%s] Running cross-pipeline Dice (Fx1)...\n', current_name);
    dice_results = compute_cross_pipeline_dice(validated_data_gtvp, config_struct, ...
        summary_metrics.id_list, summary_metrics.gtv_locations);

    % Generate summary figure: heatmap of mean Dice (methods x pipeline pairs)
    mean_dice = nanmean_safe(dice_results.dice, 3);  % average across patients
    fig = figure('Visible', 'off', 'Position', [100 100 700 500]);
    imagesc(mean_dice);
    colormap(parula); colorbar;
    caxis([0 1]);
    set(gca, 'XTick', 1:3, 'XTickLabel', dice_results.pipeline_pair_labels, ...
        'YTick', 1:11, 'YTickLabel', dice_results.method_names, ...
        'FontSize', 8, 'XTickLabelRotation', 30);
    title(sprintf('Cross-Pipeline Dice at Fx1 (%s)', config_struct.dwi_type_name));
    % Overlay text values
    for i = 1:11
        for j_idx = 1:3
            val = mean_dice(i, j_idx);
            if ~isnan(val)
                if val < 0.5, txt_color = [1 1 1]; else, txt_color = [0 0 0]; end
                text(j_idx, i, sprintf('%.2f', val), 'HorizontalAlignment', 'center', ...
                    'FontSize', 7, 'Color', txt_color);
            end
        end
    end
    saveas(fig, fullfile(config_struct.output_folder, ...
        sprintf('cross_pipeline_dice_heatmap_%s.png', config_struct.dwi_type_name)));
    close(fig);

    fprintf('      ✅ Done.\n');
end
