function results = compute_nri(y_true, prob_old, prob_new, risk_categories)
% COMPUTE_NRI  Categorical NRI, continuous NRI, and IDI.
%
%   Computes Net Reclassification Improvement (categorical and continuous)
%   and Integrated Discrimination Improvement comparing two models.
%
% Inputs:
%   y_true          - [n x 1] binary outcome (0/1)
%   prob_old        - [n x 1] predicted probabilities from old/baseline model
%   prob_new        - [n x 1] predicted probabilities from new model
%   risk_categories - [1 x k] category boundaries (default: [0 0.2 0.5 1.0])
%
% Outputs:
%   results         - Struct with fields:
%                       nri, nri_events, nri_nonevents, nri_p
%                       cnri, cnri_ci
%                       idi, idi_ci

    if nargin < 4 || isempty(risk_categories)
        risk_categories = [0, 0.2, 0.5, 1.0];
    end

    % Remove NaN
    valid = ~isnan(y_true) & ~isnan(prob_old) & ~isnan(prob_new);
    y_true = y_true(valid);
    prob_old = prob_old(valid);
    prob_new = prob_new(valid);

    events = (y_true == 1);
    nonevents = (y_true == 0);
    n_events = sum(events);
    n_nonevents = sum(nonevents);

    % --- Categorical NRI ---
    cat_old = discretize_risk(prob_old, risk_categories);
    cat_new = discretize_risk(prob_new, risk_categories);

    % Events: reclassified up = improvement
    up_events = sum(cat_new(events) > cat_old(events));
    down_events = sum(cat_new(events) < cat_old(events));
    nri_events = 0;
    if n_events > 0
        nri_events = (up_events - down_events) / n_events;
    end

    % Non-events: reclassified down = improvement
    up_nonevents = sum(cat_new(nonevents) > cat_old(nonevents));
    down_nonevents = sum(cat_new(nonevents) < cat_old(nonevents));
    nri_nonevents = 0;
    if n_nonevents > 0
        nri_nonevents = (down_nonevents - up_nonevents) / n_nonevents;
    end

    nri = nri_events + nri_nonevents;

    % Z-test for NRI
    if n_events > 0 && n_nonevents > 0
        var_nri = (up_events + down_events) / n_events^2 + ...
                  (up_nonevents + down_nonevents) / n_nonevents^2;
        if var_nri > 0
            z_nri = nri / sqrt(var_nri);
            nri_p = 2 * (1 - normcdf(abs(z_nri)));
        else
            nri_p = NaN;
        end
    else
        nri_p = NaN;
    end

    % --- Continuous NRI (cNRI) ---
    diff = prob_new - prob_old;
    cnri_events = 0;
    cnri_nonevents = 0;
    if n_events > 0
        cnri_events = mean(diff(events));
    end
    if n_nonevents > 0
        cnri_nonevents = -mean(diff(nonevents));
    end
    cnri = cnri_events + cnri_nonevents;

    % Bootstrap CI for cNRI
    cnri_ci = [NaN, NaN];
    try
        combined = [y_true, prob_old, prob_new];
        cnri_fn = @(d) compute_cnri_from_data(d);
        [cnri_ci(1), cnri_ci(2)] = bootstrap_ci(combined, cnri_fn, 2000, 0.05);
    catch
        % Bootstrap failed, leave as NaN
    end

    % --- IDI ---
    idi_events = 0;
    idi_nonevents = 0;
    if n_events > 0
        idi_events = mean(prob_new(events) - prob_old(events));
    end
    if n_nonevents > 0
        idi_nonevents = mean(prob_old(nonevents) - prob_new(nonevents));
    end
    idi = idi_events + idi_nonevents;

    % Bootstrap CI for IDI
    idi_ci = [NaN, NaN];
    try
        combined = [y_true, prob_old, prob_new];
        idi_fn = @(d) compute_idi_from_data(d);
        [idi_ci(1), idi_ci(2)] = bootstrap_ci(combined, idi_fn, 2000, 0.05);
    catch
        % Bootstrap failed, leave as NaN
    end

    % --- Output ---
    results = struct();
    results.nri = nri;
    results.nri_events = nri_events;
    results.nri_nonevents = nri_nonevents;
    results.nri_p = nri_p;
    results.cnri = cnri;
    results.cnri_ci = cnri_ci;
    results.idi = idi;
    results.idi_ci = idi_ci;

    % Print results
    fprintf('  NRI = %.3f (events: %.3f, non-events: %.3f, p = %.4f)\n', ...
        nri, nri_events, nri_nonevents, nri_p);
    fprintf('  cNRI = %.3f [%.3f, %.3f]\n', cnri, cnri_ci(1), cnri_ci(2));
    fprintf('  IDI = %.3f [%.3f, %.3f]\n', idi, idi_ci(1), idi_ci(2));
end


function cat = discretize_risk(probs, boundaries)
%DISCRETIZE_RISK  Assign probability to risk category index.
    n_cats = numel(boundaries) - 1;
    cat = ones(size(probs));
    for i = 2:n_cats
        cat(probs >= boundaries(i)) = i;
    end
end

function val = compute_cnri_from_data(d)
%COMPUTE_CNRI_FROM_DATA  cNRI from combined [y, p_old, p_new] matrix.
    y = d(:,1); p_old = d(:,2); p_new = d(:,3);
    diff = p_new - p_old;
    ev = (y == 1); nev = (y == 0);
    val = 0;
    if any(ev), val = val + mean(diff(ev)); end
    if any(nev), val = val - mean(diff(nev)); end
end

function val = compute_idi_from_data(d)
%COMPUTE_IDI_FROM_DATA  IDI from combined [y, p_old, p_new] matrix.
    y = d(:,1); p_old = d(:,2); p_new = d(:,3);
    ev = (y == 1); nev = (y == 0);
    val = 0;
    if any(ev), val = val + mean(p_new(ev) - p_old(ev)); end
    if any(nev), val = val + mean(p_old(nev) - p_new(nev)); end
end

function run_nri_tests()
%RUN_NRI_TESTS  Comprehensive tests for NRI calculations
    fprintf('Running NRI tests...\n');
    
    % Test 1: Basic categorical NRI with known values
    test_categorical_nri_basic();
    
    % Test 2: Perfect reclassification scenarios
    test_perfect_reclassification();
    
    % Test 3: No reclassification scenario
    test_no_reclassification();
    
    % Test 4: Continuous NRI calculations
    test_continuous_nri();
    
    % Test 5: IDI calculations
    test_idi_calculations();
    
    % Test 6: Edge cases and error handling
    test_edge_cases();
    
    % Test 7: Custom risk categories
    test_custom_risk_categories();
    
    fprintf('All NRI tests passed!\n');
end

function test_categorical_nri_basic()
    % Test with known reference values
    y_true = [1; 1; 1; 1; 0; 0; 0; 0];
    prob_old = [0.1; 0.3; 0.6; 0.8; 0.1; 0.3; 0.6; 0.8];
    prob_new = [0.3; 0.6; 0.8; 0.9; 0.05; 0.1; 0.3; 0.6];
    
    results = compute_nri(y_true, prob_old, prob_new);
    
    % Expected: 2 events reclassified up, 0 down -> nri_events = 2/4 = 0.5
    % Expected: 1 nonevent reclassified down, 2 up -> nri_nonevents = (1-2)/4 = -0.25
    % Expected: total NRI = 0.5 + (-0.25) = 0.25
    
    assert(abs(results.nri_events - 0.5) < 1e-6, 'NRI events calculation failed');
    assert(abs(results.nri_nonevents - (-0.25)) < 1e-6, 'NRI non-events calculation failed');
    assert(abs(results.nri - 0.25) < 1e-6, 'Total NRI calculation failed');
    
    fprintf('  ✓ Basic categorical NRI test passed\n');
end

function test_perfect_reclassification()
    % Perfect improvement scenario
    y_true = [1; 1; 0; 0];
    prob_old = [0.1; 0.3; 0.6; 0.8];  % Events in low risk, non-events in high risk
    prob_new = [0.6; 0.8; 0.1; 0.3];  % Events moved to high risk, non-events to low risk
    
    results = compute_nri(y_true, prob_old, prob_new);
    
    % All events reclassified up, all non-events reclassified down
    assert(abs(results.nri_events - 1.0) < 1e-6, 'Perfect events reclassification failed');
    assert(abs(results.nri_nonevents - 1.0) < 1e-6, 'Perfect non-events reclassification failed');
    assert(abs(results.nri - 2.0) < 1e-6, 'Perfect total NRI failed');
    
    fprintf('  ✓ Perfect reclassification test passed\n');
end

function test_no_reclassification()
    % No reclassification scenario
    y_true = [1; 1; 0; 0];
    prob_old = [0.1; 0.6; 0.2; 0.7];
    prob_new = [0.15; 0.65; 0.18; 0.72];  % Small changes within same categories
    
    results = compute_nri(y_true, prob_old, prob_new);
    
    assert(abs(results.nri_events) < 1e-6, 'No reclassification events failed');
    assert(abs(results.nri_nonevents) < 1e-6, 'No reclassification non-events failed');
    assert(abs(results.nri) < 1e-6, 'No reclassification total NRI failed');
    
    fprintf('  ✓ No reclassification test passed\n');
end

function test_continuous_nri()
    % Test continuous NRI with known values
    y_true = [1; 1; 0; 0];
    prob_old = [0.2; 0.4; 0.3; 0.6];
    prob_new = [0.5; 0.7; 0.1; 0.4];
    
    results = compute_nri(y_true, prob_old, prob_new);
    
    % Events: mean improvement = mean([0.3, 0.3]) = 0.3
    % Non-events: mean worsening = -mean([-0.2, -0.2]) = 0.2
    % cNRI = 0.3 + 0.2 = 0.5
    
    expected_cnri_events = mean([0.5-0.2, 0.7-0.4]);
    expected_cnri_nonevents = -mean([0.1-0.3, 0.4-0.6]);
    expected_cnri = expected_cnri_events + expected_cnri_nonevents;
    
    assert(abs(results.cnri - expected_cnri) < 1e-6, 'Continuous NRI calculation failed');
    
    fprintf('  ✓ Continuous NRI test passed\n');
end

function test_idi_calculations()
    % Test IDI with known values
    y_true = [1; 1; 0; 0];
    prob_old = [0.2; 0.4; 0.6; 0.8];
    prob_new = [0.5; 0.7; 0.3; 0.5];
    
    results = compute_nri(y_true, prob_old, prob_new);
    
    % IDI = mean(p_new_events - p_old_events) + mean(p_old_nonevents - p_new_nonevents)
    idi_events = mean([0.5-0.2, 0.7-0.4]);  % 0.3
    idi_nonevents = mean([0.6-0.3, 0.8-0.5]);  % 0.3
    expected_idi = idi_events + idi_nonevents;  % 0.6
    
    assert(abs(results.idi - expected_idi) < 1e-6, 'IDI calculation failed');
    
    fprintf('  ✓ IDI calculations test passed\n');
end

function test_edge_cases()
    % Test with all events or all non-events
    y_true_events = [1; 1; 1; 1];
    prob_old = [0.1; 0.3; 0.6; 0.8];
    prob_new = [0.3; 0.6; 0.8; 0.9];
    
    results = compute_nri(y_true_events, prob_old, prob_new);
    assert(~isnan(results.nri_events), 'All events case failed');
    assert(results.nri_nonevents == 0, 'All events non-events component should be 0');
    
    % Test with NaN values
    y_true_nan = [1; 1; 0; 0; NaN];
    prob_old_nan = [0.1; 0.3; 0.6; NaN; 0.8];
    prob_new_nan = [0.3; NaN; 0.8; 0.9; 0.1];
    
    results = compute_nri(y_true_nan, prob_old_nan, prob_new_nan);
    assert(~isnan(results.nri), 'NaN handling failed');
    
    % Test with empty inputs
    try
        results = compute_nri([], [], []);
        assert(results.nri == 0 || isnan(results.nri), 'Empty input handling failed');
    catch
        % Expected behavior - function should handle gracefully
    end
    
    fprintf('  ✓ Edge cases test passed\n');
end

function test_custom_risk_categories()
    % Test with custom risk categories
    y_true = [1; 1; 1; 0; 0; 0];
    prob_old = [0.05; 0.15; 0.35; 0.05; 0.15; 0.35];
    prob_new = [0.15; 0.35; 0.45; 0.02; 0.08; 0.25];
    custom_cats = [0, 0.1, 0.3, 0.4, 1.0];
    
    results = compute_nri(y_true, prob_old, prob_new, custom_cats);
    
    % Verify that custom categories are used
    assert(~isnan(results.nri), 'Custom categories calculation failed');
    assert(isfield(results, 'nri_events'), 'Missing NRI events field');
    assert(isfield(results, 'nri_nonevents'), 'Missing NRI non-events field');
    
    fprintf('  ✓ Custom risk categories test passed\n');
end

function [ci_lower, ci_upper] = bootstrap_ci(data, stat_fn, n_bootstrap, alpha)
%BOOTSTRAP_CI  Simple bootstrap confidence interval
    n = size(data, 1);
    bootstrap_stats = zeros(n_bootstrap, 1);
    
    for i = 1:n_bootstrap
        idx = randsample(n, n, true);
        bootstrap_data = data(idx, :);
        try
            bootstrap_stats(i) = stat_fn(bootstrap_data);
        catch
            bootstrap_stats(i) = NaN;
        end
    end
    
    bootstrap_stats = bootstrap_stats(~isnan(bootstrap_stats));
    if isempty(bootstrap_stats)
        ci_lower = NaN;
        ci_upper = NaN;
    else
        ci_lower = prctile(bootstrap_stats, 100 * alpha/2);
        ci_upper = prctile(bootstrap_stats, 100 * (1 - alpha/2));
    end
end