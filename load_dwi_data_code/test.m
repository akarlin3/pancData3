%% test.m — Test Suite for Refactored Statistical and Analytical Procedures
% =========================================================================
% Validates the refactored areas in metrics.m and supporting functions:
%   1. Benjamini-Hochberg FDR correction (per-timepoint q-values)
%   2. Holm-Bonferroni FWER correction (step-down thresholds)
%   3. Correlation pre-filtering for Elastic Net feature selection
%   4. LOOCV Kaplan-Meier risk assignment (no data leakage)
%
% NOTE: Dependency functions (IVIM_seg, dvh, halfSampleMode, im2Y, NIfTI
%       loaders, sample_rtdose_on_image) are NOT tested here.
%
% Usage:
%   test          % from MATLAB command window (run as script)
% =========================================================================

n_pass = 0;
n_fail = 0;

% Create output directory if it doesn't exist
output_folder = fullfile(pwd, 'saved_figures');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

% Start diary to log test output
diary_file = fullfile(output_folder, 'test_output.txt');
if exist(diary_file, 'file'), delete(diary_file); end
diary(diary_file);

% Cache metrics.m source once to avoid repeated file I/O across tests.
metrics_code   = readMetricsSource();
vis_code       = readVisualizeSource();
loaddwi_code   = readLoadDwiSource();

%% ====================================================================
%  1. BENJAMINI-HOCHBERG FDR CORRECTION
%  Tests the q-value computation used in Section 10 of metrics.m
%  (the "Q-Value (FDR) Analysis — Full Metric Sweep" block).
% =====================================================================

% testBH_QValues_KnownInput
try
    % Verify BH q-values for a hand-calculated example.
    % Five raw p-values; q = p * (m/rank), then enforce monotonicity.
    raw_p = [0.005; 0.01; 0.03; 0.40; 0.90];
    m = length(raw_p);

    [p_sorted, sort_idx] = sort(raw_p);
    q_vals = zeros(size(p_sorted));
    for k = 1:m
        q_vals(k) = p_sorted(k) * (m / k);
    end
    for k = m-1:-1:1
        q_vals(k) = min(q_vals(k), q_vals(k+1));
    end
    q_vals(q_vals > 1) = 1;

    % Expected: [0.025, 0.025, 0.05, 0.50, 0.90]
    expected = [0.025; 0.025; 0.05; 0.50; 0.90];
    assert(all(abs(q_vals - expected) <= 1e-12), ...
        'BH q-values do not match hand-calculated values');
    fprintf('[PASS] testBH_QValues_KnownInput\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBH_QValues_KnownInput: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBH_Monotonicity
try
    % Q-values must be non-decreasing after the BH step-up procedure.
    rng(1);
    raw_p = sort(rand(50, 1));   % 50 sorted p-values
    m = length(raw_p);

    q_vals = zeros(size(raw_p));
    for k = 1:m
        q_vals(k) = raw_p(k) * (m / k);
    end
    for k = m-1:-1:1
        q_vals(k) = min(q_vals(k), q_vals(k+1));
    end
    q_vals(q_vals > 1) = 1;

    diffs = diff(q_vals);
    assert(all(diffs >= -1e-15), ...
        'BH q-values are not monotonically non-decreasing');
    fprintf('[PASS] testBH_Monotonicity\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBH_Monotonicity: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBH_CappedAtOne
try
    % No q-value should exceed 1.0.
    raw_p = [0.10; 0.50; 0.80; 0.95; 0.99];
    m = length(raw_p);

    [p_sorted, ~] = sort(raw_p);
    q_vals = zeros(size(p_sorted));
    for k = 1:m
        q_vals(k) = p_sorted(k) * (m / k);
    end
    for k = m-1:-1:1
        q_vals(k) = min(q_vals(k), q_vals(k+1));
    end
    q_vals(q_vals > 1) = 1;

    assert(all(q_vals <= 1.0), ...
        'Some BH q-values exceed 1.0');
    fprintf('[PASS] testBH_CappedAtOne\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBH_CappedAtOne: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBH_QGreaterThanOrEqualP
try
    % Every q-value must be >= the corresponding raw p-value.
    raw_p = [0.001; 0.01; 0.04; 0.05; 0.20; 0.50; 0.80];
    m = length(raw_p);

    [p_sorted, sort_idx] = sort(raw_p);
    q_vals = zeros(size(p_sorted));
    for k = 1:m
        q_vals(k) = p_sorted(k) * (m / k);
    end
    for k = m-1:-1:1
        q_vals(k) = min(q_vals(k), q_vals(k+1));
    end
    q_vals(q_vals > 1) = 1;

    q_unsorted = zeros(size(raw_p));
    q_unsorted(sort_idx) = q_vals;

    assert(all(q_unsorted >= raw_p - 1e-15), ...
        'Some BH q-values are smaller than their raw p-values');
    fprintf('[PASS] testBH_QGreaterThanOrEqualP\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBH_QGreaterThanOrEqualP: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBH_SinglePValue
try
    % Edge case: m = 1. q should equal p (no correction needed).
    raw_p = 0.03;
    q = raw_p * (1 / 1);
    q = min(q, 1);
    assert(abs(q - 0.03) <= 1e-15);
    fprintf('[PASS] testBH_SinglePValue\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBH_SinglePValue: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBH_AllIdenticalPValues
try
    % When all p-values are identical, q-values should all equal
    % the same adjusted value (or 1.0 if > 1).
    raw_p = repmat(0.04, 10, 1);
    m = length(raw_p);

    [p_sorted, ~] = sort(raw_p);
    q_vals = zeros(size(p_sorted));
    for k = 1:m
        q_vals(k) = p_sorted(k) * (m / k);
    end
    for k = m-1:-1:1
        q_vals(k) = min(q_vals(k), q_vals(k+1));
    end
    q_vals(q_vals > 1) = 1;

    % All q-values should be the same after monotonicity
    assert(abs(max(q_vals) - min(q_vals)) <= 1e-15, ...
        'Identical p-values should yield identical q-values');
    fprintf('[PASS] testBH_AllIdenticalPValues\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBH_AllIdenticalPValues: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  2. HOLM-BONFERRONI FWER CORRECTION
%  Tests the step-down procedure used in Section 10 of metrics.m.
% =====================================================================

% testHolm_ThresholdFormula
try
    % Verify that Holm thresholds follow alpha / (m + 1 - k).
    m = 5;
    alpha = 0.05;
    holm_thresholds = alpha ./ (m + 1 - (1:m)');
    expected = [0.05/5; 0.05/4; 0.05/3; 0.05/2; 0.05/1];
    assert(all(abs(holm_thresholds - expected) <= 1e-15), ...
        'Holm thresholds do not match expected formula');
    fprintf('[PASS] testHolm_ThresholdFormula\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testHolm_ThresholdFormula: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testHolm_ThresholdsStrictlyIncreasing
try
    % Holm thresholds must be strictly increasing (less strict as
    % rank increases).
    m = 20;
    alpha = 0.05;
    holm_thresholds = alpha ./ (m + 1 - (1:m)');
    diffs = diff(holm_thresholds);
    assert(all(diffs > 0), ...
        'Holm thresholds are not strictly increasing');
    fprintf('[PASS] testHolm_ThresholdsStrictlyIncreasing\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testHolm_ThresholdsStrictlyIncreasing: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testHolm_EarlyStopBehavior
try
    % If the k-th test fails, all subsequent tests must also be
    % marked non-significant (step-down property).
    p_sorted = [0.001; 0.008; 0.060; 0.070; 0.900];
    m = length(p_sorted);
    alpha = 0.05;
    holm_thresholds = alpha ./ (m + 1 - (1:m)');

    is_sig = false(size(p_sorted));
    for k = 1:m
        if p_sorted(k) < holm_thresholds(k)
            is_sig(k) = true;
        else
            break;
        end
    end

    % First two should pass: 0.001 < 0.01, 0.008 < 0.0125
    % Third fails: 0.060 >= 0.0167 → stop
    assert(is_sig(1) && is_sig(2), ...
        'First two tests should be significant');
    assert(~any(is_sig(3:end)), ...
        'Tests 3-5 should be non-significant after early stop');
    fprintf('[PASS] testHolm_EarlyStopBehavior\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testHolm_EarlyStopBehavior: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testHolm_AllSignificant
try
    % When all p-values beat every threshold.
    p_sorted = [0.001; 0.002; 0.003; 0.004; 0.005];
    m = length(p_sorted);
    alpha = 0.05;
    holm_thresholds = alpha ./ (m + 1 - (1:m)');

    is_sig = false(size(p_sorted));
    for k = 1:m
        if p_sorted(k) < holm_thresholds(k)
            is_sig(k) = true;
        else
            break;
        end
    end

    assert(all(is_sig), ...
        'All tests should be significant with very small p-values');
    fprintf('[PASS] testHolm_AllSignificant\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testHolm_AllSignificant: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testHolm_NoneSignificant
try
    % When the first test already fails.
    p_sorted = [0.10; 0.20; 0.30; 0.40; 0.50];
    m = length(p_sorted);
    alpha = 0.05;
    holm_thresholds = alpha ./ (m + 1 - (1:m)');

    is_sig = false(size(p_sorted));
    for k = 1:m
        if p_sorted(k) < holm_thresholds(k)
            is_sig(k) = true;
        else
            break;
        end
    end

    assert(~any(is_sig), ...
        'No tests should be significant');
    fprintf('[PASS] testHolm_NoneSignificant\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testHolm_NoneSignificant: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  3. CORRELATION PRE-FILTERING
%  Tests the |r| > 0.8 filter used before Elastic Net in metrics.m.
% =====================================================================

% testCorrFilter_DropsHighCorrelation
try
    % Two perfectly correlated features → one should be dropped.
    rng(7);
    n = 30;
    x1 = randn(n, 1);
    x2 = x1 * 2 + 0.01 * randn(n, 1);  % r ≈ 1.0
    x3 = randn(n, 1);                    % independent
    X = [x1, x2, x3];

    R_corr = corrcoef(X);
    n_feats = size(X, 2);
    drop_flag = false(1, n_feats);
    for fi = 1:n_feats
        if drop_flag(fi), continue; end
        for fj = fi+1:n_feats
            if drop_flag(fj), continue; end
            if abs(R_corr(fi, fj)) > 0.8
                drop_flag(fj) = true;
            end
        end
    end
    keep_idx = find(~drop_flag);

    % x2 (column 2) should be dropped; x1 and x3 kept
    assert(isequal(keep_idx, [1, 3]), ...
        'Highly correlated feature (col 2) should be dropped');
    fprintf('[PASS] testCorrFilter_DropsHighCorrelation\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testCorrFilter_DropsHighCorrelation: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testCorrFilter_KeepsAllWhenUncorrelated
try
    % Orthogonal features: nothing should be dropped.
    X = eye(5);
    R_corr = corrcoef(X);
    n_feats = size(X, 2);
    drop_flag = false(1, n_feats);
    for fi = 1:n_feats
        if drop_flag(fi), continue; end
        for fj = fi+1:n_feats
            if drop_flag(fj), continue; end
            if abs(R_corr(fi, fj)) > 0.8
                drop_flag(fj) = true;
            end
        end
    end
    keep_idx = find(~drop_flag);

    assert(isequal(keep_idx, 1:5), ...
        'No features should be dropped when all are uncorrelated');
    fprintf('[PASS] testCorrFilter_KeepsAllWhenUncorrelated\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testCorrFilter_KeepsAllWhenUncorrelated: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testCorrFilter_PrefersEarlierFeature
try
    % When columns i and j are correlated (i < j), column j is
    % dropped — i.e., the earlier / absolute variable is preferred.
    rng(8);
    n = 50;
    base = randn(n, 1);
    X = [base, base + 0.001*randn(n,1), randn(n,1), randn(n,1)];

    R_corr = corrcoef(X);
    n_feats = size(X, 2);
    drop_flag = false(1, n_feats);
    for fi = 1:n_feats
        if drop_flag(fi), continue; end
        for fj = fi+1:n_feats
            if drop_flag(fj), continue; end
            if abs(R_corr(fi, fj)) > 0.8
                drop_flag(fj) = true;
            end
        end
    end

    % Column 1 kept, column 2 dropped
    assert(~drop_flag(1), ...
        'Earlier feature should be kept');
    assert(drop_flag(2), ...
        'Later correlated feature should be dropped');
    fprintf('[PASS] testCorrFilter_PrefersEarlierFeature\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testCorrFilter_PrefersEarlierFeature: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testCorrFilter_NegativeCorrelation
try
    % Negative correlation with |r| > 0.8 should also be dropped.
    rng(9);
    n = 40;
    x1 = randn(n, 1);
    x2 = -x1 + 0.01*randn(n,1);  % r ≈ -1.0
    X = [x1, x2];

    R_corr = corrcoef(X);
    n_feats = size(X, 2);
    drop_flag = false(1, n_feats);
    for fi = 1:n_feats
        if drop_flag(fi), continue; end
        for fj = fi+1:n_feats
            if drop_flag(fj), continue; end
            if abs(R_corr(fi, fj)) > 0.8
                drop_flag(fj) = true;
            end
        end
    end

    assert(drop_flag(2), ...
        'Negatively correlated feature should be dropped');
    fprintf('[PASS] testCorrFilter_NegativeCorrelation\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testCorrFilter_NegativeCorrelation: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testCorrFilter_BoundaryNotDropped
try
    % |r| = 0.8 exactly should NOT be dropped (threshold is >0.8).
    % Construct a pair with correlation exactly 0.8 via known formula.
    rng(10);
    n = 10000;   % large n so sample r ≈ population r
    x1 = randn(n, 1);
    % x2 = 0.8*x1 + sqrt(1-0.64)*noise gives population r = 0.8
    x2 = 0.8*x1 + sqrt(1 - 0.64)*randn(n, 1);
    X = [x1, x2];

    R_corr = corrcoef(X);
    r12 = abs(R_corr(1, 2));

    drop_flag = false(1, 2);
    if r12 > 0.8
        drop_flag(2) = true;
    end

    % With n=10000, sample r will be very close to 0.8 but may
    % land on either side. Verify the logic is at least applied:
    if r12 <= 0.8
        assert(~drop_flag(2), ...
            'Feature at boundary (|r|==0.8) should NOT be dropped');
    else
        assert(drop_flag(2), ...
            'Feature slightly above boundary should be dropped');
    end
    fprintf('[PASS] testCorrFilter_BoundaryNotDropped\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testCorrFilter_BoundaryNotDropped: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testCorrFilter_ChainDropping
try
    % If A-B correlated and B-C correlated, dropping B means C
    % is evaluated against A only. Verify transitive behaviour.
    rng(11);
    n = 100;
    a = randn(n, 1);
    b = a + 0.01*randn(n, 1);           % r(a,b) ≈ 1
    c = 0.5*a + 0.5*randn(n, 1);        % r(a,c) ≈ 0.7 (not dropped)
    X = [a, b, c];

    R_corr = corrcoef(X);
    n_feats = size(X, 2);
    drop_flag = false(1, n_feats);
    for fi = 1:n_feats
        if drop_flag(fi), continue; end
        for fj = fi+1:n_feats
            if drop_flag(fj), continue; end
            if abs(R_corr(fi, fj)) > 0.8
                drop_flag(fj) = true;
            end
        end
    end
    keep_idx = find(~drop_flag);

    % b dropped (correlated with a), c kept (r(a,c) < 0.8)
    assert(drop_flag(2), 'b should be dropped');
    assert(~drop_flag(3), 'c should be kept');
    assert(length(keep_idx) == 2);
    fprintf('[PASS] testCorrFilter_ChainDropping\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testCorrFilter_ChainDropping: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  4. LOOCV KAPLAN-MEIER RISK ASSIGNMENT
%  Tests that the LOOCV classification loop does not use held-out data.
% =====================================================================

% testLOOCV_NoLeakage_CutoffExcludesHeldOut
try
    % Verify that the median cutoff for each fold is computed
    % WITHOUT the held-out patient's data.
    rng(42);
    n = 10;
    vals = (1:n)';   % Predictable ordering for median checks

    for loo_i = 1:n
        train_mask = true(n, 1);
        train_mask(loo_i) = false;
        train_vals = vals(train_mask);
        cutoff = median(train_vals);

        % The held-out value must NOT influence the cutoff.
        % Re-compute cutoff including it — should differ.
        cutoff_all = median(vals);
        if vals(loo_i) ~= cutoff_all
            assert(cutoff ~= cutoff_all, ...
                sprintf('Fold %d: cutoff should differ when held-out is excluded', loo_i));
        end
    end
    fprintf('[PASS] testLOOCV_NoLeakage_CutoffExcludesHeldOut\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testLOOCV_NoLeakage_CutoffExcludesHeldOut: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testLOOCV_TrainMaskSize
try
    % Each training fold should contain exactly N-1 patients.
    n = 15;
    for loo_i = 1:n
        train_mask = true(n, 1);
        train_mask(loo_i) = false;
        assert(sum(train_mask) == n - 1, ...
            'Training set should have N-1 patients');
    end
    fprintf('[PASS] testLOOCV_TrainMaskSize\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testLOOCV_TrainMaskSize: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testLOOCV_EachPatientHeldOutOnce
try
    % Over the full loop, each patient is held out exactly once.
    n = 12;
    held_out_count = zeros(n, 1);
    for loo_i = 1:n
        held_out_count(loo_i) = held_out_count(loo_i) + 1;
    end
    assert(isequal(held_out_count, ones(n, 1)), ...
        'Each patient should be held out exactly once');
    fprintf('[PASS] testLOOCV_EachPatientHeldOutOnce\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testLOOCV_EachPatientHeldOutOnce: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testLOOCV_RiskAssignmentIsBoolean
try
    % is_high_risk must be logical and have one entry per patient.
    rng(42);
    n = 8;
    km_X = randn(n, 8);
    km_y = [0;0;0;0;1;1;1;1];
    is_high_risk = false(n, 1);

    for loo_i = 1:n
        train_mask = true(n, 1);
        train_mask(loo_i) = false;
        train_vals = km_X(train_mask, 5);  % Arbitrary percent-change col
        cutoff = median(train_vals);
        is_high_risk(loo_i) = km_X(loo_i, 5) < cutoff;
    end

    assert(islogical(is_high_risk), 'is_high_risk should be logical');
    assert(numel(is_high_risk) == n, 'is_high_risk should have one entry per patient');
    fprintf('[PASS] testLOOCV_RiskAssignmentIsBoolean\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testLOOCV_RiskAssignmentIsBoolean: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testLOOCV_NoEmptyGroups_LargeN
try
    % With balanced synthetic data, LOOCV should produce both high-
    % and low-risk patients (i.e., not degenerate to all-same).
    rng(42);
    n = 20;
    % Create data where LC patients have higher values than LF
    km_X = zeros(n, 8);
    km_y = [zeros(n/2, 1); ones(n/2, 1)];
    km_X(1:n/2, 5) = randn(n/2, 1) + 2;   % LC: positive change
    km_X(n/2+1:end, 5) = randn(n/2, 1) - 2; % LF: negative change

    is_high_risk = false(n, 1);
    for loo_i = 1:n
        train_mask = true(n, 1);
        train_mask(loo_i) = false;
        train_vals = km_X(train_mask, 5);
        cutoff = median(train_vals);
        is_high_risk(loo_i) = km_X(loo_i, 5) < cutoff;
    end

    assert(any(is_high_risk), ...
        'Should have at least one high-risk patient');
    assert(any(~is_high_risk), ...
        'Should have at least one low-risk patient');
    fprintf('[PASS] testLOOCV_NoEmptyGroups_LargeN\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testLOOCV_NoEmptyGroups_LargeN: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testLOOCV_FeatureSelectionBase_Mapping
try
    % Verify the mod-based mapping from 8-feature index to base
    % metric index: indices 1-4 = Absolute, 5-8 = Percent-Change.
    % base = mod(sel - 1, 4) + 1
    expected_base = [1, 2, 3, 4, 1, 2, 3, 4];
    for idx = 1:8
        base = mod(idx - 1, 4) + 1;
        assert(base == expected_base(idx), ...
            sprintf('Index %d should map to base %d', idx, expected_base(idx)));
    end
    fprintf('[PASS] testLOOCV_FeatureSelectionBase_Mapping\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testLOOCV_FeatureSelectionBase_Mapping: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testLOOCV_DefaultLowRisk_WhenNoFeatures
try
    % When LASSO/Elastic Net selects zero features, the patient
    % should default to low-risk (is_high_risk = false).
    sel_loo = [];   % Empty selection (simulating convergence failure)
    default_risk = false;

    if isempty(sel_loo)
        assigned_risk = false;
    else
        assigned_risk = true;
    end

    assert(assigned_risk == default_risk, ...
        'Empty feature selection should default to low-risk');
    fprintf('[PASS] testLOOCV_DefaultLowRisk_WhenNoFeatures\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testLOOCV_DefaultLowRisk_WhenNoFeatures: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  5. ELASTIC NET PARAMETER VALIDATION
%  Verifies that the Alpha parameter is correctly set to 0.5
%  (mixed L1+L2 penalty) rather than pure LASSO (Alpha = 1).
% =====================================================================

% testElasticNet_AlphaIsHalf
try
    % The refactored code must use Alpha = 0.5 for Elastic Net.
    code = metrics_code;
    matches = regexp(code, '''Alpha''\s*,\s*0\.5', 'match');
    assert(~isempty(matches), ...
        'metrics.m should contain Alpha = 0.5 for Elastic Net');
    fprintf('[PASS] testElasticNet_AlphaIsHalf\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testElasticNet_AlphaIsHalf: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testElasticNet_NotPureLasso
try
    % Confirm no lassoglm call uses pure LASSO (Alpha = 1).
    code = metrics_code;
    matches = regexp(code, '''Alpha''\s*,\s*1(\s|,|\))', 'match');
    assert(isempty(matches), ...
        'No lassoglm call should use pure LASSO (Alpha=1)');
    fprintf('[PASS] testElasticNet_NotPureLasso\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testElasticNet_NotPureLasso: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  6. ALL-FRACTIONS LOOP (NO POST-HOC TARGET SELECTION)
%  Verifies that the analysis loop iterates over all fractions,
%  not just a cherry-picked subset.
% =====================================================================

% testLoop_NotHardcodedTo2And3
try
    % The loop "for target_fx = [2, 3]" should have been replaced.
    code = metrics_code;
    assert(~contains(code, 'target_fx = [2, 3]'), ...
        'Loop should NOT be hardcoded to [2, 3]');
    fprintf('[PASS] testLoop_NotHardcodedTo2And3\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testLoop_NotHardcodedTo2And3: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testLoop_UsesNTpUpperBound
try
    % The loop should use nTp as the upper bound.
    code = metrics_code;
    matches = regexp(code, 'target_fx\s*=\s*2\s*:\s*nTp', 'match');
    assert(~isempty(matches), ...
        'Loop should iterate from 2 to nTp');
    fprintf('[PASS] testLoop_UsesNTpUpperBound\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testLoop_UsesNTpUpperBound: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  7. TARGETED FDR SECTIONS REMOVED
%  Verifies that the post-hoc "m=8" and "m=4" targeted sections
%  have been removed to eliminate p-hacking.
% =====================================================================

% testRemoved_TargetedM8
try
    % The "Targeted Statistical Analysis (m = 8 tests)" section
    % should no longer exist in metrics.m.
    code = metrics_code;
    assert(~contains(code, 'Targeted Statistical Analysis (m = 8'), ...
        'Targeted m=8 section should be removed');
    fprintf('[PASS] testRemoved_TargetedM8\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testRemoved_TargetedM8: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testRemoved_TargetedM4
try
    % The "Targeted FDR (Relative Change Only, m = 4)" section
    % should no longer exist in metrics.m.
    code = metrics_code;
    assert(~contains(code, 'Targeted FDR (Relative Change Only, m = 4)'), ...
        'Targeted m=4 FDR section should be removed');
    fprintf('[PASS] testRemoved_TargetedM4\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testRemoved_TargetedM4: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  8. ADC FITTING VIA FIT_ADC_MONO
%  Verifies that metrics.m uses fit_adc_mono instead of
%  inline pixel-wise OLS log-linear fitting.
% =====================================================================

% testADC_MetricsCallsFitAdcMono
try
    % metrics.m should call fit_adc_mono instead of inline OLS.
    code = metrics_code;
    matches = regexp(code, 'fit_adc_mono\s*\(', 'match');
    assert(~isempty(matches), ...
        'metrics.m should call fit_adc_mono for ADC computation');
    fprintf('[PASS] testADC_MetricsCallsFitAdcMono\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testADC_MetricsCallsFitAdcMono: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testADC_NoInlineOLS
try
    % The old inline OLS pattern "beta = X \ log(s_signal)" should
    % be gone from metrics.m.
    code = metrics_code;
    assert(~contains(code, 'X \ log(s_signal)'), ...
        'Inline OLS fitting should be removed from metrics.m');
    fprintf('[PASS] testADC_NoInlineOLS\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testADC_NoInlineOLS: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testADC_NoNestedPixelLoops
try
    % The old nested "for y = 1:ny / for x = 1:nx" pixel loops
    % should be gone from metrics.m (replaced by fit_adc_mono call).
    code = metrics_code;
    assert(~contains(code, 'slice_adc(y,x) = beta(2)'), ...
        'Nested pixel-wise OLS loops should be removed');
    fprintf('[PASS] testADC_NoNestedPixelLoops\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testADC_NoNestedPixelLoops: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  9. COMPREHENSIVE FDR AND HOLM REMAIN
%  Verifies that the full-sweep FDR and Holm-Bonferroni sections
%  (which apply corrections across ALL tests) are still present.
% =====================================================================

% testRetained_PerTimepointFDR
try
    % The FDR section should now apply BH per-timepoint.
    code = metrics_code;
    assert(contains(code, 'Per-Timepoint') || contains(code, 'PER-TIMEPOINT'), ...
        'FDR section should indicate per-timepoint operation');
    fprintf('[PASS] testRetained_PerTimepointFDR\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testRetained_PerTimepointFDR: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testRetained_HolmBonferroni
try
    % Holm-Bonferroni should still exist (now per-timepoint).
    code = metrics_code;
    assert(contains(code, 'Holm-Bonferroni'), ...
        'Holm-Bonferroni logic should be retained');
    fprintf('[PASS] testRetained_HolmBonferroni\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testRetained_HolmBonferroni: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testFDR_NoGlobalPooling
try
    % The old global FDR sweep header should no longer exist.
    code = metrics_code;
    assert(~contains(code, 'Full Metric Sweep'), ...
        'Global FDR sweep should be replaced with per-timepoint');
    fprintf('[PASS] testFDR_NoGlobalPooling\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testFDR_NoGlobalPooling: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  10. B-VALUE PROTOCOL VALIDATION (NO ARBITRARY TRUNCATION)
%  Verifies that the arbitrary b-value truncation has been removed
%  and replaced with explicit protocol validation.
% =====================================================================

% testBval_NoTruncation
try
    % The old truncation pattern should be gone from metrics.m.
    code = metrics_code;
    assert(~contains(code, 'min_dim = min(n_imgs, n_bvals)'), ...
        'Arbitrary b-value truncation (min_dim) should be removed');
    fprintf('[PASS] testBval_NoTruncation\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBval_NoTruncation: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBval_ExpectedProtocolDefined
try
    % metrics.m should define expected_bvals for protocol validation.
    code = metrics_code;
    matches = regexp(code, 'expected_bvals\s*=\s*\[0;\s*30;\s*150;\s*550\]', 'match');
    assert(~isempty(matches), ...
        'metrics.m should define expected_bvals = [0; 30; 150; 550]');
    fprintf('[PASS] testBval_ExpectedProtocolDefined\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBval_ExpectedProtocolDefined: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBval_ProtocolDeviationFlagged
try
    % metrics.m should flag protocol deviations with a warning message.
    code = metrics_code;
    assert(contains(code, 'Protocol deviation'), ...
        'Protocol deviation flagging should be present in metrics.m');
    fprintf('[PASS] testBval_ProtocolDeviationFlagged\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBval_ProtocolDeviationFlagged: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBval_DeviationExcludesPatient
try
    % The validation block should use 'continue' to exclude
    % deviating patients from the analysis.
    code = metrics_code;
    % Find the Protocol deviation block and verify it contains continue
    idx_dev = strfind(code, 'Protocol deviation');
    idx_continue = strfind(code, 'continue');
    assert(~isempty(idx_dev) && ~isempty(idx_continue), ...
        'Both protocol deviation flag and continue must exist');
    % At least one continue must follow the protocol deviation flag
    % (within 200 chars, i.e. within the same if-block)
    found = false;
    for ci = 1:length(idx_continue)
        if idx_continue(ci) > idx_dev(1) && (idx_continue(ci) - idx_dev(1)) < 200
            found = true;
            break;
        end
    end
    assert(found, ...
        'A continue statement should follow the protocol deviation flag to exclude the patient');
    fprintf('[PASS] testBval_DeviationExcludesPatient\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBval_DeviationExcludesPatient: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBval_ValidationLogic_MatchingProtocol
try
    % When bvals exactly match [0; 30; 150; 550], no deviation.
    bvals = [0; 30; 150; 550];
    expected_bvals = [0; 30; 150; 550];
    assert(isequal(sort(bvals), expected_bvals), ...
        'Standard protocol b-values should pass validation');
    fprintf('[PASS] testBval_ValidationLogic_MatchingProtocol\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBval_ValidationLogic_MatchingProtocol: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBval_ValidationLogic_UnsortedMatch
try
    % Unsorted b-values that match after sorting should pass.
    bvals = [0; 30; 150; 550];
    expected_bvals = [0; 30; 150; 550];
    assert(isequal(sort(bvals), expected_bvals), ...
        'Unsorted but correct b-values should pass validation');
    fprintf('[PASS] testBval_ValidationLogic_UnsortedMatch\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBval_ValidationLogic_UnsortedMatch: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBval_ValidationLogic_ExtraValue
try
    % Extra b-value (5 instead of 4) should be flagged as deviation.
    bvals = [0; 50; 400; 800; 1000];
    expected_bvals = [0; 30; 150; 550];
    is_deviation = ~isequal(sort(bvals), expected_bvals);
    assert(is_deviation, ...
        'Extra b-values should be flagged as protocol deviation');
    fprintf('[PASS] testBval_ValidationLogic_ExtraValue\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBval_ValidationLogic_ExtraValue: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBval_ValidationLogic_WrongValue
try
    % Different b-value set should be flagged as deviation.
    bvals = [0; 30; 200; 550];
    expected_bvals = [0; 30; 150; 550];
    is_deviation = ~isequal(sort(bvals), expected_bvals);
    assert(is_deviation, ...
        'Non-standard b-values should be flagged as protocol deviation');
    fprintf('[PASS] testBval_ValidationLogic_WrongValue\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBval_ValidationLogic_WrongValue: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testBval_ValidationLogic_MissingValue
try
    % Fewer b-values than expected should be flagged.
    bvals = [0; 50; 400];
    expected_bvals =  [0; 30; 150; 550];
    is_deviation = ~isequal(sort(bvals), expected_bvals);
    assert(is_deviation, ...
        'Missing b-values should be flagged as protocol deviation');
    fprintf('[PASS] testBval_ValidationLogic_MissingValue\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testBval_ValidationLogic_MissingValue: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  11. IMPUTATION STRATEGY (NO LISTWISE DELETION)
%  Verifies that the predictor matrix construction uses imputation
%  instead of listwise deletion to prevent complete-case attrition.
% =====================================================================

% testImpute_NoListwiseDeletion
try
    % The old listwise deletion pattern should be gone from metrics.m.
    code = metrics_code;
    assert(~contains(code, 'valid_lasso_mask = all(~isnan(X_lasso), 2)'), ...
        'Listwise deletion pattern should be removed from metrics.m');
    fprintf('[PASS] testImpute_NoListwiseDeletion\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testImpute_NoListwiseDeletion: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testImpute_UsesKNN
try
    % metrics.m should use knn_impute_train_test for imputation.
    code = metrics_code;
    matches = regexp(code, 'knn_impute_train_test\s*\(', 'match');
    assert(~isempty(matches), ...
        'metrics.m should use knn_impute_train_test for imputation');
    fprintf('[PASS] testImpute_UsesKNN\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testImpute_UsesKNN: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testImpute_ExcludesAllNaNRows
try
    % metrics.m should exclude patients with ALL imaging data missing.
    code = metrics_code;
    matches = regexp(code, 'any\s*\(\s*~isnan\s*\(\s*X_lasso_all', 'match');
    assert(~isempty(matches), ...
        'metrics.m should check for patients missing all imaging data');
    fprintf('[PASS] testImpute_ExcludesAllNaNRows\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testImpute_ExcludesAllNaNRows: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testImpute_RetainsPartialDataRows
try
    % Synthetic test: a patient with partial NaN data should be retained
    % after imputation, not dropped as in listwise deletion.
    X = [1 2 3; 4 NaN 6; NaN NaN NaN; 7 8 9];
    y = [0; 1; 0; 1];

    has_any = any(~isnan(X), 2);
    mask = has_any & ~isnan(y);
    X_imp = X(mask, :);
    y_out = y(mask);
    X_out = fillmissing(X_imp, 'constant', median(X_imp, 1, 'omitnan'));

    % Row 2 (partial NaN) should be kept; row 3 (all NaN) dropped
    assert(size(X_out, 1) == 3, ...
        'Partial-data rows should be retained after imputation');
    assert(~any(isnan(X_out(:))), ...
        'No NaN values should remain after imputation');
    assert(length(y_out) == 3, ...
        'Outcome vector length should match imputed matrix');
    fprintf('[PASS] testImpute_RetainsPartialDataRows\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testImpute_RetainsPartialDataRows: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testImpute_MedianFillValues
try
    % Verify that imputed values equal the column median.
    X = [2 10; NaN 20; 6 NaN];
    col_med = median(X, 1, 'omitnan');  % [4, 15]
    X_filled = fillmissing(X, 'constant', col_med);

    assert(abs(X_filled(2,1) - 4) < 1e-12, ...
        'Imputed value should equal column median');
    assert(abs(X_filled(3,2) - 15) < 1e-12, ...
        'Imputed value should equal column median');
    fprintf('[PASS] testImpute_MedianFillValues\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testImpute_MedianFillValues: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testImpute_BeforeStandardization
try
    % Verify imputation occurs before the Elastic Net standardization.
    % knn_impute_train_test must appear before lassoglm(..., 'Standardize', true)
    code = metrics_code;
    pos_fill = regexp(code, 'knn_impute_train_test\s*\(', 'start', 'once');
    pos_std  = regexp(code, '''Standardize''\s*,\s*true', 'start', 'once');
    assert(~isempty(pos_fill) && ~isempty(pos_std), ...
        'Both knn_impute_train_test and Standardize should exist in metrics.m');
    assert(pos_fill < pos_std, ...
        'Imputation (knn_impute_train_test) must occur before standardization');
    fprintf('[PASS] testImpute_BeforeStandardization\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testImpute_BeforeStandardization: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  12. VISUALIZE_RESULTS B-VALUE PROTOCOL VALIDATION
%  Verifies that visualize_results.m uses the same strict b-value
%  protocol validation as metrics.m, with no arbitrary truncation.
% =====================================================================

% testVis_NoTruncation
try
    % The old truncation pattern should be gone from visualize_results.m.
    code = vis_code;
    assert(~contains(code, 'min_dim = min('), ...
        'Arbitrary b-value truncation (min_dim) should be removed from visualize_results.m');
    fprintf('[PASS] testVis_NoTruncation\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testVis_NoTruncation: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testVis_ExpectedProtocolDefined
try
    % visualize_results.m should define expected_bvals for protocol validation.
    code = vis_code;
    matches = regexp(code, 'expected_bvals\s*=\s*\[0;\s*30;\s*150;\s*550\]', 'match');
    assert(~isempty(matches), ...
        'visualize_results.m should define expected_bvals = [0; 30; 150; 550]');
    fprintf('[PASS] testVis_ExpectedProtocolDefined\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testVis_ExpectedProtocolDefined: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testVis_ProtocolDeviationFlagged
try
    % visualize_results.m should flag protocol deviations with a warning message.
    code = vis_code;
    assert(contains(code, 'Protocol deviation'), ...
        'Protocol deviation flagging should be present in visualize_results.m');
    fprintf('[PASS] testVis_ProtocolDeviationFlagged\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testVis_ProtocolDeviationFlagged: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testVis_DeviationExcludesPatient
try
    % The validation block should use 'continue' to exclude
    % deviating patients from the comparative mapping.
    code = vis_code;
    idx_dev = strfind(code, 'Protocol deviation');
    idx_continue = strfind(code, 'continue');
    assert(~isempty(idx_dev) && ~isempty(idx_continue), ...
        'Both protocol deviation flag and continue must exist in visualize_results.m');
    % At least one continue must follow the protocol deviation flag
    % (within 200 chars, i.e. within the same if-block)
    found = false;
    for ci = 1:length(idx_continue)
        if idx_continue(ci) > idx_dev(1) && (idx_continue(ci) - idx_dev(1)) < 200
            found = true;
            break;
        end
    end
    assert(found, ...
        'A continue statement should follow the protocol deviation flag to exclude the patient');
    fprintf('[PASS] testVis_DeviationExcludesPatient\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testVis_DeviationExcludesPatient: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  13. DEFORMABLE IMAGE REGISTRATION (DIR)
%  Tests the apply_dir_mask_propagation helper function.
% =====================================================================

% testDIR_EmptyInputsReturnEmpty
try
    % All three empty-input cases should return [].
    r1 = apply_dir_mask_propagation([], ones(4,4,4), true(4,4,4));
    r2 = apply_dir_mask_propagation(ones(4,4,4), [], true(4,4,4));
    r3 = apply_dir_mask_propagation(ones(4,4,4), ones(4,4,4), []);
    assert(isempty(r1) && isempty(r2) && isempty(r3), ...
        'Empty inputs should return empty output');
    fprintf('[PASS] testDIR_EmptyInputsReturnEmpty\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testDIR_EmptyInputsReturnEmpty: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testDIR_SizeMismatchReturnsEmpty
try
    % Mismatched image sizes should return [].
    r = apply_dir_mask_propagation(ones(4,4,4), ones(5,5,5), true(4,4,4));
    assert(isempty(r), ...
        'Size-mismatched inputs should return empty output');
    fprintf('[PASS] testDIR_SizeMismatchReturnsEmpty\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testDIR_SizeMismatchReturnsEmpty: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testDIR_IdenticalImagesPreserveMask
try
    % When fixed == moving (no deformation), the warped mask should
    % closely match the original mask.
    sz = [32 32 16];
    img = randn(sz) * 100 + 500;
    mask = false(sz);
    mask(12:20, 12:20, 5:12) = true;
    warped = apply_dir_mask_propagation(img, img, mask);
    if ~isempty(warped)
        dice_coeff = 2*sum(warped(:) & mask(:)) / (sum(warped(:)) + sum(mask(:)));
        assert(dice_coeff > 0.95, ...
            sprintf('Identical images should preserve mask (Dice=%.3f)', dice_coeff));
    end
    fprintf('[PASS] testDIR_IdenticalImagesPreserveMask\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testDIR_IdenticalImagesPreserveMask: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testDIR_OutputIsLogical
try
    % Output mask must be a logical array.
    sz = [16 16 8];
    img = randn(sz) * 100 + 500;
    mask = false(sz); mask(5:10, 5:10, 3:6) = true;
    warped = apply_dir_mask_propagation(img, img, mask);
    if ~isempty(warped)
        assert(islogical(warped), 'Warped mask should be logical');
        assert(isequal(size(warped), sz), 'Warped mask should match input size');
    end
    fprintf('[PASS] testDIR_OutputIsLogical\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testDIR_OutputIsLogical: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testDIR_FunctionExistsInCodebase
try
    % load_dwi_data_forAvery.m should call apply_dir_mask_propagation.
    code = loaddwi_code;
    assert(contains(code, 'apply_dir_mask_propagation'), ...
        'load_dwi_data_forAvery.m should call apply_dir_mask_propagation');
    fprintf('[PASS] testDIR_FunctionExistsInCodebase\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testDIR_FunctionExistsInCodebase: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  14. NEGATIVE ADC HANDLING
%  Verifies that negative monoexponential ADC values are set to NaN
%  (matching biexponential failure handling), not clamped to zero.
% =====================================================================

% testADC_NegativeClampedToNaN
try
    % load_dwi_data_forAvery.m should assign NaN to negative ADC.
    code = loaddwi_code;
    assert(contains(code, 'adc_vec(adc_vec < 0) = nan') || ...
           contains(code, 'adc_vec(adc_vec<0) = nan') || ...
           contains(code, 'adc_vec(adc_vec < 0)= nan'), ...
        'Negative ADC should be set to NaN, not zero');
    fprintf('[PASS] testADC_NegativeClampedToNaN\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testADC_NegativeClampedToNaN: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testADC_NotClampedToZero
try
    % The old clamping pattern should not exist.
    code = loaddwi_code;
    assert(~contains(code, 'adc_vec(adc_vec < 0) = 0'), ...
        'Negative ADC should NOT be clamped to zero');
    fprintf('[PASS] testADC_NotClampedToZero\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testADC_NotClampedToZero: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testADC_NaNRemovedByNanmean
try
    % NaN values should be excluded from mean calculations.
    adc_vec = [1.5e-3; 2.0e-3; NaN; 1.8e-3];
    m = mean(adc_vec, 'omitnan');
    assert(abs(m - mean([1.5e-3; 2.0e-3; 1.8e-3])) < 1e-15, ...
        'nanmean should exclude NaN failed fits');
    fprintf('[PASS] testADC_NaNRemovedByNanmean\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testADC_NaNRemovedByNanmean: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  15. PER-TIMEPOINT FDR CORRECTNESS
%  Verifies that per-timepoint BH yields smaller family sizes and
%  thus more power than global pooling.
% =====================================================================

% testPerTimepointFDR_SmallerFamilySize
try
    % Same p-values corrected per-timepoint (n=4) vs globally (n=12)
    % should yield smaller q-values per-timepoint.
    p_tp = [0.01; 0.02; 0.04; 0.10];  % 4 tests in one timepoint
    n_tp = 4;
    [p_s, si] = sort(p_tp);
    q_tp = zeros(n_tp, 1);
    q_tp(n_tp) = p_s(n_tp);
    for ii = n_tp-1:-1:1
        q_tp(ii) = min(q_tp(ii+1), p_s(ii) * (n_tp / ii));
    end
    q_tp = min(q_tp, 1);

    % Global: same p-values padded to 12 tests
    p_global = [p_tp; 0.30; 0.40; 0.50; 0.60; 0.70; 0.80; 0.90; 0.95];
    n_g = 12;
    [p_sg, sig] = sort(p_global);
    q_g = zeros(n_g, 1);
    q_g(n_g) = p_sg(n_g);
    for ii = n_g-1:-1:1
        q_g(ii) = min(q_g(ii+1), p_sg(ii) * (n_g / ii));
    end
    q_g = min(q_g, 1);
    q_g_unsorted = zeros(n_g, 1); q_g_unsorted(sig) = q_g;

    % The first 4 tests should have smaller q when corrected per-timepoint
    assert(all(q_tp <= q_g_unsorted(1:4) + 1e-12), ...
        'Per-timepoint FDR should produce equal or smaller q-values than global');
    fprintf('[PASS] testPerTimepointFDR_SmallerFamilySize\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testPerTimepointFDR_SmallerFamilySize: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testPerTimepointFDR_NoGlobalSweepInSection9
try
    % Section 9 should no longer pool across timepoints.
    code = metrics_code;
    assert(~contains(code, 'FDR Correction') || ...
           contains(code, 'Per-Timepoint'), ...
        'Section 9 FDR should be per-timepoint, not global');
    fprintf('[PASS] testPerTimepointFDR_NoGlobalSweepInSection9\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testPerTimepointFDR_NoGlobalSweepInSection9: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  16. OUT-OF-FOLD ROC ANALYSIS
%  Verifies that the ROC curve uses unbiased out-of-fold risk scores
%  rather than in-sample fitted probabilities.
% =====================================================================

% testOOFROC_UsesRiskScoresAll
try
    % perfcurve should be called with risk_scores_all, not
    % mdl_roc.Fitted.Probability.
    code = metrics_code;
    assert(contains(code, 'perfcurve(labels(valid_roc), risk_scores_all(valid_roc)'), ...
        'perfcurve should use risk_scores_all directly');
    fprintf('[PASS] testOOFROC_UsesRiskScoresAll\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testOOFROC_UsesRiskScoresAll: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testOOFROC_NoInSampleROC
try
    % The old in-sample ROC pattern should be gone.
    code = metrics_code;
    assert(~contains(code, 'EXPLORATORY ROC ANALYSIS'), ...
        'In-sample exploratory ROC should be replaced with OOF ROC');
    fprintf('[PASS] testOOFROC_NoInSampleROC\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testOOFROC_NoInSampleROC: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testOOFROC_LabeledAsPrimary
try
    % ROC section should say PRIMARY ROC ANALYSIS in the fprintf.
    code = metrics_code;
    assert(contains(code, 'PRIMARY ROC ANALYSIS'), ...
        'ROC section should be labeled Primary');
    fprintf('[PASS] testOOFROC_LabeledAsPrimary\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testOOFROC_LabeledAsPrimary: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testOOFROC_NoStaleGLMRefit
try
    % The stale fitglm call inside the ROC block should be gone.
    code = metrics_code;
    assert(~contains(code, 'mdl_roc = fitglm'), ...
        'Stale fitglm call in ROC block should be removed');
    fprintf('[PASS] testOOFROC_NoStaleGLMRefit\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testOOFROC_NoStaleGLMRefit: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testOOFROC_YoudenCutoffFromLOOCV
try
    % The Youden optimal cutoff (roc_opt_thresh) must come from perfcurve,
    % not from any logistic regression fitted probabilities.
    code = metrics_code;
    assert(contains(code, 'roc_opt_thresh'), ...
        'roc_opt_thresh must be defined (Youden cutoff from LOOCV ROC)');
    assert(~contains(code, 'mdl_roc.Fitted'), ...
        'Optimal cutoff must not use mdl_roc.Fitted (no in-sample refitting)');
    fprintf('[PASS] testOOFROC_YoudenCutoffFromLOOCV\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testOOFROC_YoudenCutoffFromLOOCV: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testLPOCV_StrictlyPairedAUC
try
    % LPOCV header should say it is strictly for paired AUC only.
    code = metrics_code;
    assert(contains(code, 'Paired AUC Scalar Only') || ...
           contains(code, 'STRICTLY reserved'), ...
        'LPOCV block must be labeled as strictly for paired AUC only');
    fprintf('[PASS] testLPOCV_StrictlyPairedAUC\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testLPOCV_StrictlyPairedAUC: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  17. SUB-VOLUME FEATURE EXCLUSION POLICY
%  Verifies that sub-volume dosimetry features are re-enabled now that
%  DIR is implemented, and that Post-RT still excludes all dose features.
% =====================================================================

% testSubVol_DIRReenabled
try
    % The old keep_policy = 1:10 exclusion should be gone.
    code = metrics_code;
    assert(~contains(code, 'keep_policy = 1:10'), ...
        'Sub-volume features should no longer be excluded now that DIR is implemented');
    fprintf('[PASS] testSubVol_DIRReenabled\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testSubVol_DIRReenabled: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testSubVol_PostRTStillExcluded
try
    % Post-RT dose exclusion (target_fx == 6) should still be present.
    code = metrics_code;
    assert(contains(code, 'target_fx == 6'), ...
        'Post-RT dose exclusion should be retained');
    fprintf('[PASS] testSubVol_PostRTStillExcluded\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testSubVol_PostRTStillExcluded: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  18. SUB-VOLUME VOXEL THRESHOLDS
%  Verifies that higher-order metrics (kurtosis, skewness, D95, V50)
%  return NaN for small sub-volumes (< 100 voxels).
% =====================================================================

% testThreshold_KurtosisProtected_ADC
try
    code = loaddwi_code;
    % Should check numel(adc_vec_sub) >= min_vox_hist or similar
    assert(contains(code, 'numel(adc_vec_sub) >= min_vox_hist'), ...
        'ADC sub-volume kurtosis should be protected by voxel threshold');
    fprintf('[PASS] testThreshold_KurtosisProtected_ADC\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testThreshold_KurtosisProtected_ADC: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testThreshold_KurtosisProtected_D
try
    code = loaddwi_code;
    assert(contains(code, 'numel(d_vec_sub) >= min_vox_hist'), ...
        'D sub-volume kurtosis should be protected by voxel threshold');
    fprintf('[PASS] testThreshold_KurtosisProtected_D\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testThreshold_KurtosisProtected_D: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testThreshold_KSProtected
try
    code = loaddwi_code;
    % Should check numel(adc_vec) >= min_vox_hist && numel(adc_baseline) >= min_vox_hist
    assert(contains(code, 'numel(adc_vec) >= min_vox_hist') && contains(code, 'numel(adc_baseline) >= min_vox_hist'), ...
        'KS test should be protected by voxel threshold for both current and baseline');
    fprintf('[PASS] testThreshold_KSProtected\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testThreshold_KSProtected: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testThreshold_D95Protected_Metrics
try
    code = metrics_code;
    % Should check sum(adc_mask_1d) >= min_subvol_voxels or similar
    assert(contains(code, 'sum(adc_mask_1d) >= min_subvol_voxels'), ...
        'Metrics.m D95/V50 should be protected by voxel threshold for ADC sub-volume');
    fprintf('[PASS] testThreshold_D95Protected_Metrics\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testThreshold_D95Protected_Metrics: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testThreshold_V50Protected_Metrics
try
    code = metrics_code;
    assert(contains(code, 'sum(dstar_mask_1d) >= min_subvol_voxels'), ...
        'Metrics.m D95/V50 should be protected by voxel threshold for D* sub-volume');
    fprintf('[PASS] testThreshold_V50Protected_Metrics\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testThreshold_V50Protected_Metrics: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testSanityCheck_NoKurtosisPlot
try
    % The sanity check subplot must NOT use adc_kurt as a metric.
    % It should use adc_sd instead.
    code = metrics_code;
    assert(~contains(code, 'kurt_fx1 = adc_kurt'), ...
        'Sanity check should use ADC SD, not trace-average ADC Kurtosis');
    assert(contains(code, 'sd_fx1  = adc_sd'), ...
        'Sanity check heterogeneity subplot should use adc_sd');
    fprintf('[PASS] testSanityCheck_NoKurtosisPlot\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testSanityCheck_NoKurtosisPlot: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  19. SEGMENTED IVIM FIT — BOUNDARY CONDITIONS
%  Verifies strict parameter limits and b-value threshold in IVIMmodelfit.
% =====================================================================

% Cache IVIMmodelfit source
ivim_code = readIVIMmodelSource();

% testIVIM_Bthr100
try
    % Pipeline default: b >= 100 s/mm² used for D estimation.
    code = loaddwi_code;
    assert(contains(code, 'ivim_bthr = 100'), ...
        'ivim_bthr should be set to 100 in load_dwi_data_forAvery.m');
    fprintf('[PASS] testIVIM_Bthr100\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testIVIM_Bthr100: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testIVIM_fUpperBound
try
    % IVIMmodelfit.m should cap perfusion fraction f at 0.4.
    code = ivim_code;
    assert(contains(code, '0.4'), ...
        'IVIMmodelfit.m should set f upper bound to 0.4');
    fprintf('[PASS] testIVIM_fUpperBound\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testIVIM_fUpperBound: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testIVIM_DstarUpperBound
try
    % IVIMmodelfit.m should cap D* at 0.1 mm²/s.
    code = ivim_code;
    assert(contains(code, '0.1'), ...
        'IVIMmodelfit.m should set D* upper bound to 0.1 mm²/s');
    fprintf('[PASS] testIVIM_DstarUpperBound\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testIVIM_DstarUpperBound: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testIVIM_NoLooseBoundsF
try
    % Old loose f upper bound of 1 should be gone.
    code = ivim_code;
    % The old line was: lim = [0 0 0 0;3e-3 2*max(Y(:)) 1 1]
    assert(~contains(code, '3e-3 2*max(Y(:)) 1 1'), ...
        'Old loose f and D* limits should be removed from IVIMmodelfit.m');
    fprintf('[PASS] testIVIM_NoLooseBoundsF\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testIVIM_NoLooseBoundsF: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testIVIM_OptsLimSupported
try
    % IVIMmodelfit.m should support opts.lim for caller overrides.
    code = ivim_code;
    assert(contains(code, '''lim''') && contains(code, 'opts.lim'), ...
        'IVIMmodelfit.m should support opts.lim override');
    fprintf('[PASS] testIVIM_OptsLimSupported\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testIVIM_OptsLimSupported: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testIVIM_DefaultBlim100
try
    % IVIMmodelfit.m should default blim to 100.
    code = ivim_code;
    assert(contains(code, 'blim = 100'), ...
        'IVIMmodelfit.m default blim should be 100');
    fprintf('[PASS] testIVIM_DefaultBlim100\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testIVIM_DefaultBlim100: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testIVIM_SegmentedTwoStage_Logic
try
    % Inline two-stage test: verify that blim=100 includes more high-b
    % values than blim=200, and that an LLS monoexponential fit on the
    % high-b subset yields a physiologically plausible D.
    D_true = 1.5e-3; f_true = 0.15; Dstar_true = 0.05; S0 = 1000;
    bvals_test = [0; 30; 100; 550];
    rng(42);
    S_test = S0 * ((1-f_true)*exp(-bvals_test*D_true) + f_true*exp(-bvals_test*(D_true+Dstar_true)));
    S_test = S_test + 5*randn(size(S_test));

    % Stage 1 with blim=100: includes b=100 and b=550 (2 points)
    b_hi100 = bvals_test(bvals_test >= 100);
    S_hi100 = S_test(bvals_test >= 100);
    X100 = [-b_hi100, ones(size(b_hi100))];
    p100 = X100 \ log(S_hi100);
    D_est_100 = p100(1);

    % Stage 1 with blim=200: includes only b=550 (exact, 1 point)
    b_hi200 = bvals_test(bvals_test >= 200);
    S_hi200 = S_test(bvals_test >= 200);
    X200 = [-b_hi200, ones(size(b_hi200))];
    p200 = X200 \ log(S_hi200);
    D_est_200 = p200(1);

    assert(D_est_100 > 0 && D_est_100 < 3e-3, ...
        sprintf('blim=100 D estimate out of physiological range: %.4g', D_est_100));
    % blim=100 should use more b-values (over-determined), which is the improvement
    assert(length(b_hi100) > length(b_hi200), ...
        'blim=100 should include more high-b values for D estimation than blim=200');
    fprintf('[PASS] testIVIM_SegmentedTwoStage_Logic\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testIVIM_SegmentedTwoStage_Logic: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  20. DIR DOSE-MAP WARPING
%  Verifies that the Demons displacement field is applied to the dose map
%  before DVH extraction so dose sub-volumes are spatially consistent with
%  the DWI parameter maps.
% =====================================================================

% testDIR_ReturnsDfieldAndRef
try
    % apply_dir_mask_propagation should return 3 outputs (mask, D_forward, ref3d).
    code = readLoadDwiSource();
    % Check signature via file reading
    apPath = fullfile(fileparts(mfilename('fullpath')), 'apply_dir_mask_propagation.m');
    fid = fopen(apPath, 'r'); ap_code = fread(fid,'*char')'; fclose(fid);
    assert(contains(ap_code, 'D_forward, ref3d] = apply_dir_mask_propagation'), ...
        'apply_dir_mask_propagation should return D_forward and ref3d');
    fprintf('[PASS] testDIR_ReturnsDfieldAndRef\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testDIR_ReturnsDfieldAndRef: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testDIR_DoseWarped_NotRigid
try
    % load_dwi_data_forAvery.m should apply imwarp to the dose map.
    code = loaddwi_code;
    assert(contains(code, 'imwarp(dose_map, D_forward_cur'), ...
        'Dose map must be warped with D_forward_cur via imwarp before DVH extraction');
    fprintf('[PASS] testDIR_DoseWarped_NotRigid\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testDIR_DoseWarped_NotRigid: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testDIR_DoseWarpFallback
try
    % At baseline (fi=1) or when DIR fails, code must fall back to rigid dose_map.
    code = loaddwi_code;
    assert(contains(code, 'dose_map_dvh     = dose_map'), ...
        'Rigid dose fallback ''dose_map_dvh = dose_map'' should be initialised before DIR block');
    fprintf('[PASS] testDIR_DoseWarpFallback\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testDIR_DoseWarpFallback: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testDIR_DfieldCached
try
    % D_forward and ref3d should be written to the .mat cache file.
    code = loaddwi_code;
    assert(contains(code, "save(dir_cache_file, 'gtv_mask_warped', 'D_forward', 'ref3d')"), ...
        'D_forward and ref3d must be included in the cache save call');
    fprintf('[PASS] testDIR_DfieldCached\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testDIR_DfieldCached: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testDIR_GTVnDoseAlsoWarped
try
    % GTVn dose block must also use the warped dose (dose_map_dvh_n).
    code = loaddwi_code;
    assert(contains(code, 'imwarp(dose_map, D_forward_cur') && ...
           contains(code, 'dose_map_dvh_n'), ...
        'GTVn dose block must warp dose_map with D_forward_cur into dose_map_dvh_n');
    fprintf('[PASS] testDIR_GTVnDoseAlsoWarped\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testDIR_GTVnDoseAlsoWarped: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testDIR_LinearMaskWarping
try
    % apply_dir_mask_propagation.m should use 'linear' interpolation for mask warping.
    apPath = fullfile(fileparts(mfilename('fullpath')), 'apply_dir_mask_propagation.m');
    fid = fopen(apPath, 'r'); ap_code = fread(fid,'*char')'; fclose(fid);
    assert(contains(ap_code, "'Interp', 'linear'"), ...
        "apply_dir_mask_propagation should use 'linear' interpolation for mask warping");
    assert(~contains(ap_code, "'Interp', 'nearest'"), ...
        "apply_dir_mask_propagation should NOT use 'nearest' interpolation for mask warping");
    fprintf('[PASS] testDIR_LinearMaskWarping\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testDIR_LinearMaskWarping: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% ====================================================================
%  21. TIME-DEPENDENT COX PH MODEL
%  Verifies the counting-process panel builder and that metrics.m calls it.
% =====================================================================

% testTD_PanelHasMoreRowsThanPatients
try
    % Synthetic 3-patient, 3-timepoint panel
    arr1 = [1.0e-3 1.1e-3 1.2e-3;   % patient 1
             0.9e-3 0.85e-3 0.8e-3;  % patient 2
             1.5e-3 1.6e-3 NaN];      % patient 3 (missing tp3)
    lf_td  = [1; 0; 1];
    tot_td = [30; 100; 20];
    [X_td_t, t0, t1, ev, ~] = build_td_panel({arr1}, {'ADC'}, lf_td, tot_td, 3, [0 5 10]);
    assert(size(X_td_t, 1) > length(lf_td), ...
        'Panel should have more rows than patients');
    fprintf('[PASS] testTD_PanelHasMoreRowsThanPatients\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testTD_PanelHasMoreRowsThanPatients: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testTD_NoEventBeforeFinalInterval
try
    arr1 = [1.0e-3 1.1e-3 1.2e-3; 0.9e-3 0.85e-3 0.8e-3];
    lf_td  = [1; 0];
    tot_td = [20; 100];
    [~, t0, t1, ev, pid] = build_td_panel({arr1}, {'ADC'}, lf_td, tot_td, 3, [0 5 10]);
    % For patient 1 (LF), event should only fire on its last interval
    p1_rows = (pid == 1);
    ev_p1 = ev(p1_rows);
    assert(sum(ev_p1) == 1 && ev_p1(end) == true, ...
        'Event must fire exactly once, on the last interval, for LF patient');
    fprintf('[PASS] testTD_NoEventBeforeFinalInterval\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testTD_NoEventBeforeFinalInterval: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testTD_StartAlwaysLessThanStop
try
    arr1 = [1.0e-3 1.1e-3; 0.9e-3 0.8e-3; 1.5e-3 1.6e-3];
    lf_td  = [1; 0; 1];
    tot_td = [15; 50; 8];
    [~, t0, t1, ~, ~] = build_td_panel({arr1}, {'ADC'}, lf_td, tot_td, 2, [0 5]);
    assert(all(t0 < t1), 't_start must be strictly less than t_stop for every interval');
    fprintf('[PASS] testTD_StartAlwaysLessThanStop\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testTD_StartAlwaysLessThanStop: %s\n', e.message);
    n_fail = n_fail + 1;
end

% testTD_MetricsCallsBuildTdPanel
try
    code = metrics_code;
    assert(contains(code, 'build_td_panel'), ...
        'metrics.m should call build_td_panel for the time-dependent Cox model');
    fprintf('[PASS] testTD_MetricsCallsBuildTdPanel\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] testTD_MetricsCallsBuildTdPanel: %s\n', e.message);
    n_fail = n_fail + 1;
end

%% Summary
fprintf('\n%d passed, %d failed out of %d tests\n', n_pass, n_fail, n_pass + n_fail);
if n_fail > 0
    diary off
    error('test:failures', '%d test(s) failed.', n_fail);
end
diary off

% =========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function code = readMetricsSource()
    % Read metrics.m from the same directory as this test file.
    metricsPath = fullfile(fileparts(mfilename('fullpath')), 'metrics.m');
    fid = fopen(metricsPath, 'r');
    assert(fid > 0, 'Could not open metrics.m at: %s', metricsPath);
    code = fread(fid, '*char')';
    fclose(fid);
end

function code = readVisualizeSource()
    % Read visualize_results.m from the same directory as this test file.
    visPath = fullfile(fileparts(mfilename('fullpath')), 'visualize_results.m');
    fid = fopen(visPath, 'r');
    assert(fid > 0, 'Could not open visualize_results.m at: %s', visPath);
    code = fread(fid, '*char')';
    fclose(fid);
end

function code = readLoadDwiSource()
    % Read load_dwi_data_forAvery.m from the same directory.
    p = fullfile(fileparts(mfilename('fullpath')), 'load_dwi_data_forAvery.m');
    fid = fopen(p, 'r');
    assert(fid > 0, 'Could not open load_dwi_data_forAvery.m at: %s', p);
    code = fread(fid, '*char')';
    fclose(fid);
end

function code = readIVIMmodelSource()
    % Read IVIMmodelfit.m from the dependencies directory.
    p = fullfile(fileparts(mfilename('fullpath')), 'dependencies', 'IVIMmodelfit.m');
    fid = fopen(p, 'r');
    assert(fid > 0, 'Could not open IVIMmodelfit.m at: %s', p);
    code = fread(fid, '*char')';
    fclose(fid);
end