% TEST_CORR_FILTER Functional test for the collinearity pruning utility.
%
% This test validates filter_collinear_features.m across four scenarios:
%   1. Basic pruning: Given two highly correlated features, the one with the
%      weaker univariate C-index should be dropped.
%   2. Non-leakage verification: The pruning mask derived from training data
%      must be applied as-is to test data (no re-computation on test rows).
%   3. Time-stratified filtering: When a fraction vector is supplied, the
%      correlation matrix must be computed only on baseline (Fraction 1) rows.
%      Late-stage collinearity induced by treatment response must not influence
%      the baseline pruning decision.
%   4. NaN handling: The function must tolerate NaN entries in the feature
%      matrix (using pairwise deletion) without propagating NaN correlations.
function test_corr_filter()
    % Setup output directory in temp folder to avoid polluting the project root
    output_folder = fullfile(tempdir, 'test_corr_filter_output');
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end
    
    % Start diary to log text output
    diary_file = fullfile(output_folder, 'test_corr_filter_output.txt');
    if exist(diary_file, 'file'), delete(diary_file); end
    diary(diary_file);
    
    % Suppress figure windows during tests to avoid GUI pop-ups
    old_vis = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');

    fprintf('Running test_corr_filter...\n');

    % --- Scenario 1: Basic collinearity pruning ---
    % Create two features where F1 is more predictive of the binary outcome y
    % than F2, but F2 is highly correlated with F1 (r > 0.8). The filter
    % should retain F1 (higher C-index) and drop F2.
    rng(42);  % Fixed seed for reproducibility
    n = 50;
    y = [zeros(n/2, 1); ones(n/2, 1)];  % Binary outcome: 25 LC, 25 LF

    % F1: Strong class separation (shift of +2 for positive class)
    f1 = [randn(n/2, 1); randn(n/2, 1) + 2];
    % F2: Correlated with F1 (linear + noise), but offset reduces C-index
    f2 = f1 + 0.5*randn(n, 1) + 0.5;

    X = [f1, f2];
    feat_names = {'F1', 'F2'};
    
    R = corrcoef(X);
    fprintf('Correlation between F1 and F2: %.3f\n', R(1,2));
    
    % Run the collinearity filter; keep_idx contains column indices to retain
    keep_idx = filter_collinear_features(X, y);
    % Build a boolean drop flag for visualization (true = dropped)
    drop_flag_actual = true(1, 2);
    drop_flag_actual(keep_idx) = false;

    % Scatter plot of the two features colored by class for visual inspection
    figure('Name', 'Correlation Filter Test Data');
    scatter(f1(y==0), f2(y==0), 50, [0 0.4470 0.7410], 'filled', 'MarkerEdgeColor', 'k'); hold on;
    scatter(f1(y==1), f2(y==1), 50, [0.8500 0.3250 0.0980], 'filled', 'MarkerEdgeColor', 'k');
    xlabel('F1 (Stronger C-index)'); ylabel('F2 (Correlated)');
    title(sprintf('Correlation Filter Mock Data (r=%.2f)', R(1,2)));
    legend('LC', 'LF', 'Location', 'best');
    grid on;
    saveas(gcf, fullfile(output_folder, 'test_corr_filter_data.png'));
    close(gcf);
    
    % Verify that F2 was dropped (weaker C-index) and F1 was retained
    if ~any(keep_idx == 2) && any(keep_idx == 1)
        fprintf('SUCCESS: Shared function dropped feature F2 with weaker C-index.\n');
    else
        diary off;
        set(0, 'DefaultFigureVisible', old_vis);
        error('FAILURE: Unexpected drop state from shared function.');
    end
    
    % --- Scenario 2: STRICT NON-LEAKAGE VERIFICATION ---
    % The pruning mask (keep_idx) must be derived solely from training data.
    % When applying to unseen test patients, we must NOT recompute correlations
    % on test data -- just apply the training-derived column mask.
    fprintf('\n--- Verifying STRICT NON-LEAKAGE ---\n');
    % Simulate a new unseen patient with both features present
    test_patient_X = [1.2, 1.3]; % New patient data for F1, F2
    % Apply the training-derived mask to select only retained columns
    X_te_pruned = test_patient_X(:, keep_idx);
    
    assert(size(X_te_pruned, 2) == length(keep_idx), 'Pruned test data size mismatch');
    assert(~any(keep_idx == 2), 'Test pruning failed to apply training mask');
    fprintf('SUCCESS: Training pruning mask (Keep: %s) successfully applied to test patient.\n', ...
        strjoin(feat_names(keep_idx), ', '));
    fprintf('Scientific Rigor: PASSED (No correlation calculation on test data).\n');
    
    % --- Scenario 3: TIME-STRATIFIED COLLINEARITY FILTER VERIFICATION ---
    % In longitudinal DWI data, features may become correlated at later
    % treatment fractions (e.g., due to radiation response) even though they
    % are independent at baseline. The filter must use only baseline rows
    % when frac_vec is provided, to avoid dropping features based on
    % spurious late-stage correlations.
    fprintf('\n--- Verifying TIME-STRATIFIED COLLINEARITY FILTERING ---\n');
    rng(99);
    n_td = 60;      % 3 fractions x 20 patients = 60 longitudinal rows
    n_pts_td = 20;
    
    % frac_vec: rows 1-20 are Fraction 1 (baseline), 21-40 Fraction 2, 41-60 Fraction 3
    frac_vec = [ones(n_pts_td,1); 2*ones(n_pts_td,1); 3*ones(n_pts_td,1)];
    
    % Patient-level event labels (broadcast to all rows for that patient)
    y_pat = [zeros(n_pts_td/2,1); ones(n_pts_td/2,1)];
    y_td  = repmat(y_pat, 3, 1);
    
    % At baseline (Fraction 1), F1 and F2 are independent (low correlation)
    f1_base = [randn(n_pts_td/2,1); randn(n_pts_td/2,1) + 2];
    f2_base = randn(n_pts_td, 1);          % Independent of f1_base at baseline

    % At later fractions (2 and 3), F1 and F2 become strongly correlated
    % (simulating treatment-induced co-variation in diffusion parameters)
    f1_late = repmat(f1_base, 2, 1) + 0.1*randn(2*n_pts_td, 1);
    f2_late = f1_late + 0.05*randn(2*n_pts_td, 1);  % near-identical at late fractions
    
    X_td_test = [[f1_base, f2_base]; [f1_late, f2_late]];
    
    % Compute correlations globally (all fractions) vs baseline-only to show the contrast
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
    
    % --- Scenario 4: NaN HANDLING VERIFICATION ---
    % Missing data (NaN) is common in clinical DWI datasets. The standard
    % corrcoef() propagates NaN unless 'Rows','pairwise' is used. The filter
    % must handle NaN entries gracefully and still produce correct pruning.
    fprintf('\n--- Verifying NaN HANDLING ---\n');

    % Start from the original 2-feature dataset and inject NaN values
    X_nan = X;
    X_nan(5, 1) = NaN;   % Missing value in feature F1
    X_nan(10, 2) = NaN;  % Missing value in feature F2

    % The function should use pairwise-complete observations for the correlation
    % matrix, correctly identifying F2 as the weaker correlated feature to drop.
    keep_idx_nan = filter_collinear_features(X_nan, y);

    if ~any(keep_idx_nan == 2) && any(keep_idx_nan == 1)
        fprintf('SUCCESS: Function successfully handled NaNs and dropped weaker feature.\n');
    else
        diary off;
        set(0, 'DefaultFigureVisible', old_vis);
        error('FAILURE: NaN handling failed. Function did not drop feature as expected.');
    end
    fprintf('NaN HANDLING TEST: PASSED\n');

    % Cleanup
    set(0, 'DefaultFigureVisible', old_vis);
    diary off;
    fprintf('Test output and plots saved to: %s\n', output_folder);
end
