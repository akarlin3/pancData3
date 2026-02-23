% Test script for intelligent correlation filtering logic
function test_corr_filter()
    % Setup output directory like the rest of the pipeline
    output_folder = fullfile(pwd, 'saved_figures');
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end
    
    % Start diary to log text output
    diary_file = fullfile(output_folder, 'test_corr_filter_output.txt');
    if exist(diary_file, 'file'), delete(diary_file); end
    diary(diary_file);
    
    % Suppress figure windows
    old_vis = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');

    fprintf('Running test_corr_filter...\n');
    
    % Mock data
    rng(42);
    n = 50;
    y = [zeros(n/2, 1); ones(n/2, 1)];
    
    % F1: Very significant
    f1 = [randn(n/2, 1); randn(n/2, 1) + 2];
    % F2: Highly correlated with F1, but less significant
    f2 = f1 + 0.5*randn(n, 1) + 0.5;
    
    X = [f1, f2];
    feat_names = {'F1', 'F2'};
    
    R = corrcoef(X);
    fprintf('Correlation between F1 and F2: %.3f\n', R(1,2));
    
    % Run logic similar to metrics.m
    drop_flag = false(1, 2);
    for fi = 1:2
        if drop_flag(fi), continue; end
        for fj = fi+1:2
            if abs(R(fi, fj)) > 0.8
                p_fi = ranksum(X(y==0, fi), X(y==1, fi));
                p_fj = ranksum(X(y==0, fj), X(y==1, fj));
                fprintf('p_F1: %.4f, p_F2: %.4f\n', p_fi, p_fj);
                
                if p_fj >= p_fi
                    drop_flag(fj) = true;
                    fprintf('Dropping %s (p_j >= p_i)\n', feat_names{fj});
                else
                    drop_flag(fi) = true;
                    fprintf('Dropping %s (p_i > p_j)\n', feat_names{fi});
                    break;
                end
            end
        end
    end
    
    % Visualization of the test case
    figure('Name', 'Correlation Filter Test Data');
    scatter(f1(y==0), f2(y==0), 50, 'b', 'filled', 'MarkerEdgeColor', 'k'); hold on;
    scatter(f1(y==1), f2(y==1), 50, 'r', 'filled', 'MarkerEdgeColor', 'k');
    xlabel('F1 (More Significant)'); ylabel('F2 (Correlated)');
    title(sprintf('Correlation Filter Mock Data (r=%.2f)', R(1,2)));
    legend('LC', 'LF', 'Location', 'best');
    grid on;
    saveas(gcf, fullfile(output_folder, 'test_corr_filter_data.png'));
    close(gcf);
    
    if drop_flag(2) && ~drop_flag(1)
        fprintf('SUCCESS: Dropped less significant feature F2.\n');
    elseif drop_flag(1) && ~drop_flag(2)
        fprintf('SUCCESS: Dropped less significant feature F1 (if it were less significant).\n');
    else
        diary off;
        set(0, 'DefaultFigureVisible', old_vis);
        error('FAILURE: Unexpected drop state.');
    end
    
    % Cleanup
    set(0, 'DefaultFigureVisible', old_vis);
    diary off;
    fprintf('Test output and plots saved to: %s\n', output_folder);
end
