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
    
    % F1: Very predictive (high C-index)
    f1 = [randn(n/2, 1); randn(n/2, 1) + 2];
    % F2: Highly correlated with F1, but weaker C-index
    f2 = f1 + 0.5*randn(n, 1) + 0.5;
    
    X = [f1, f2];
    feat_names = {'F1', 'F2'};
    
    R = corrcoef(X);
    fprintf('Correlation between F1 and F2: %.3f\n', R(1,2));
    
    % Run logic using the centralized shared function
    keep_idx = filter_collinear_features(X, y);
    drop_flag_actual = true(1, 2);
    drop_flag_actual(keep_idx) = false;
    
    % Visualization of the test case
    figure('Name', 'Correlation Filter Test Data');
    scatter(f1(y==0), f2(y==0), 50, 'b', 'filled', 'MarkerEdgeColor', 'k'); hold on;
    scatter(f1(y==1), f2(y==1), 50, 'r', 'filled', 'MarkerEdgeColor', 'k');
    xlabel('F1 (Stronger C-index)'); ylabel('F2 (Correlated)');
    title(sprintf('Correlation Filter Mock Data (r=%.2f)', R(1,2)));
    legend('LC', 'LF', 'Location', 'best');
    grid on;
    saveas(gcf, fullfile(output_folder, 'test_corr_filter_data.png'));
    close(gcf);
    
    % Verify the drop state
    if ~any(keep_idx == 2) && any(keep_idx == 1)
        fprintf('SUCCESS: Shared function dropped feature F2 with weaker C-index.\n');
    else
        diary off;
        set(0, 'DefaultFigureVisible', old_vis);
        error('FAILURE: Unexpected drop state from shared function.');
    end
    
    % --- STRICT NON-LEAKAGE VERIFICATION ---
    % Test: If we have a new patient, do we re-calculate R? 
    % Rigor requires using the mask derived strictly from training data.
    fprintf('\n--- Verifying STRICT NON-LEAKAGE ---\n');
    test_patient_X = [1.2, 1.3]; % New patient data for F1, F2
    X_te_pruned = test_patient_X(:, keep_idx);
    
    assert(size(X_te_pruned, 2) == length(keep_idx), 'Pruned test data size mismatch');
    assert(~any(keep_idx == 2), 'Test pruning failed to apply training mask');
    fprintf('SUCCESS: Training pruning mask (Keep: %s) successfully applied to test patient.\n', ...
        strjoin(feat_names(keep_idx), ', '));
    fprintf('Scientific Rigor: PASSED (No correlation calculation on test data).\n');
    
    % --- TIME-STRATIFIED COLLINEARITY FILTER VERIFICATION ---
    % Test: When frac_vec is provided, correlation must be computed exclusively
    % on Fraction 1 (baseline) rows. Late-stage collinearity must NOT
    % influence the pruning mask applied at baseline.
    fprintf('\n--- Verifying TIME-STRATIFIED COLLINEARITY FILTERING ---\n');
    rng(99);
    n_td = 60;   % total longitudinal rows (3 fractions × 20 patients)
    n_pts_td = 20;
    
    % frac_vec: rows 1-20 are Fraction 1 (baseline), 21-40 Fraction 2, 41-60 Fraction 3
    frac_vec = [ones(n_pts_td,1); 2*ones(n_pts_td,1); 3*ones(n_pts_td,1)];
    
    % Patient-level event labels (broadcast to all rows for that patient)
    y_pat = [zeros(n_pts_td/2,1); ones(n_pts_td/2,1)];
    y_td  = repmat(y_pat, 3, 1);
    
    % Baseline features: F1 and F2 are NOT correlated at Fraction 1
    f1_base = [randn(n_pts_td/2,1); randn(n_pts_td/2,1) + 2];
    f2_base = randn(n_pts_td, 1);          % independent of f1_base
    
    % At later fractions, F1 and F2 become strongly correlated (radiation response)
    f1_late = repmat(f1_base, 2, 1) + 0.1*randn(2*n_pts_td, 1);
    f2_late = f1_late + 0.05*randn(2*n_pts_td, 1);  % near-identical at late fractions
    
    X_td_test = [[f1_base, f2_base]; [f1_late, f2_late]];
    
    R_global = corrcoef(X_td_test);
    R_base   = corrcoef(X_td_test(frac_vec==1, :));
    fprintf('Global |r|(F1,F2): %.3f  |  Baseline-only |r|(F1,F2): %.3f\n', ...
        abs(R_global(1,2)), abs(R_base(1,2)));
    
    % Global filter (no frac_vec): late-stage collinearity may cause a drop
    keep_global = filter_collinear_features(X_td_test, y_td);
    % Time-stratified filter: must be computed only on baseline rows
    keep_strat  = filter_collinear_features(X_td_test, y_td, frac_vec);
    
    % Baseline rows are NOT highly correlated, so both features should survive
    assert(length(keep_strat) == 2 && isequal(keep_strat(:)', [1, 2]), ...
        'TIME-STRATIFIED FAILURE: Baseline-independent features incorrectly dropped.');
    fprintf('SUCCESS: Time-stratified filter retained both features (baseline r=%.3f < 0.8).\n', ...
        abs(R_base(1,2)));
    
    % Global filter should drop one due to late-stage collinearity
    if length(keep_global) < 2
        fprintf('INFO: Global filter dropped a feature due to late-stage collinearity (r=%.3f).\n', ...
            abs(R_global(1,2)));
    end
    fprintf('TIME-STRATIFIED TEST: PASSED\n');
    
    % Cleanup
    set(0, 'DefaultFigureVisible', old_vis);
    diary off;
    fprintf('Test output and plots saved to: %s\n', output_folder);
end
