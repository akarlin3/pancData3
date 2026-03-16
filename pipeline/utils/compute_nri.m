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
