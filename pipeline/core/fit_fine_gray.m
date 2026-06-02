function results = fit_fine_gray(survival_data, cox_results, config)
% FIT_FINE_GRAY  Fine-Gray subdistribution hazard model for competing risks.
%   Estimates the subdistribution hazard ratio (sHR) which accounts for
%   competing events by keeping competing-event subjects in the risk set
%   with decreasing weights w(t) = G(t)/G(t_comp), where G is the KM
%   estimate of the censoring distribution.
%
%   References
%   ----------
%   Fine JP, Gray RJ. A proportional hazards model for the subdistribution
%   of a competing risk. JASA. 1999;94(446):496-509.

results = struct('success', false, 'message', '');

if ~survival_data.has_sufficient_data
    results.message = 'Insufficient data';
    fprintf('  Insufficient data for Fine-Gray model.\n');
    return;
end

% Count event types
n_competing = sum(survival_data.event == 2);
n_primary = sum(survival_data.event == 1);

if n_primary < 3
    results.message = sprintf('Insufficient primary events (%d)', n_primary);
    fprintf('  Insufficient primary events (%d) for Fine-Gray model.\n', n_primary);
    return;
end

if n_competing == 0
    results.message = 'No competing events: Fine-Gray equivalent to CSH — skipped';
    fprintf('  No competing events detected (0 event=2). Fine-Gray reduces to CSH — skipping.\n');
    return;
end

fprintf('  Fitting Fine-Gray model: %d primary events, %d competing events.\n', ...
    n_primary, n_competing);

% --- Stage 1: Compute per-patient terminal outcomes ---
[unique_pats, ~, pat_group] = unique(survival_data.pat_id);
n_pats_fg = length(unique_pats);
pat_terminal_time = zeros(n_pats_fg, 1);
pat_terminal_event = zeros(n_pats_fg, 1);
for pi = 1:n_pats_fg
    pat_mask = (pat_group == pi);
    pat_rows = find(pat_mask);
    [~, max_idx] = max(survival_data.t_stop(pat_rows));
    terminal_row = pat_rows(max_idx);
    pat_terminal_time(pi) = survival_data.t_stop(terminal_row);
    pat_terminal_event(pi) = survival_data.event(terminal_row);
end

% --- Stage 2: Compute censoring KM G(t) ---
% For censoring KM: "event" = administrative censoring (event==0),
% all actual events (1 or 2) are "censored" in this reversed view.
cens_indicator = double(pat_terminal_event == 0);
[sorted_times, sort_idx] = sort(pat_terminal_time);
sorted_cens = cens_indicator(sort_idx);

uniq_times_fg = unique(sorted_times);
G_values = ones(length(uniq_times_fg), 1);
G_current = 1.0;
for k = 1:length(uniq_times_fg)
    tk = uniq_times_fg(k);
    n_at_risk = sum(sorted_times >= tk);
    n_censored = sum(sorted_times == tk & sorted_cens == 1);
    if n_at_risk > 0
        G_current = G_current * (1 - n_censored / n_at_risk);
    end
    G_values(k) = max(G_current, 0.01);  % floor to avoid div-by-zero
end

% --- Stage 3: Compute subdistribution weights ---
fg_weights = ones(size(survival_data.event));

for pi = 1:n_pats_fg
    if pat_terminal_event(pi) ~= 2
        continue;  % Only competing-event patients need modified weights
    end

    t_comp = pat_terminal_time(pi);
    G_t_comp = interp_G(t_comp, uniq_times_fg, G_values);

    pat_rows = find(pat_group == pi);
    for ri = 1:length(pat_rows)
        row_idx = pat_rows(ri);
        t_row = survival_data.t_stop(row_idx);
        if t_row > t_comp
            G_t = interp_G(t_row, uniq_times_fg, G_values);
            fg_weights(row_idx) = G_t / G_t_comp;
        end
    end
end

% --- Stage 4: Fit weighted Cox model ---
% Reuse CSH scaling if available
if cox_results.success && isfield(cox_results, 'X_scaled')
    X_fg = cox_results.X_scaled;
else
    X_fg = scale_td_panel(survival_data.X, survival_data.feat_names, ...
        survival_data.pat_id, survival_data.t_start, unique(survival_data.pat_id), 'baseline');
end

% Primary event (1) stays as event; competing (2) and censored (0)
% are both "censored" for the subdistribution hazard. Competing-event
% subjects remain in the risk set via the subdistribution weights.
event_fg = survival_data.event;
event_fg(event_fg == 2) = 0;

[X_fg_clean, keep_cols_fg] = remove_constant_columns(X_fg);
if size(X_fg_clean, 2) == 0
    results.message = 'All covariate columns are constant';
    fprintf('  All covariate columns constant — cannot fit Fine-Gray model.\n');
    return;
end

% Apply weights via Frequency (discretized, same pattern as IPCW)
fg_scale = 100;
fg_freq = max(1, round(fg_weights * fg_scale));

ws = warning('off', 'stats:coxphfit:FitWarning');
cleanupWarn = onCleanup(@() warning(ws));

try
    T_fg = [survival_data.t_start, survival_data.t_stop];
    is_censored_fg = (event_fg == 0);

    [b_fg, logl_fg, ~, stats_fg] = coxphfit(X_fg_clean, T_fg, ...
        'Censoring', is_censored_fg, 'Ties', 'efron', 'Frequency', fg_freq);

    % Sandwich SE for weighted model
    sandwich_se_fg = compute_sandwich_se_cox(b_fg, X_fg_clean, T_fg, ...
        is_censored_fg, fg_weights, survival_data.pat_id, 'efron');

    if ~isempty(sandwich_se_fg) && all(isfinite(sandwich_se_fg)) && all(sandwich_se_fg > 0)
        stats_fg.se = sandwich_se_fg;
        fprintf('  Using robust sandwich SE for Fine-Gray model.\n');
    end

    stats_fg.p = 2 * (1 - normcdf(abs(b_fg ./ stats_fg.se)));

    % Map to full feature space
    b_full_fg = zeros(survival_data.n_feat, 1);
    b_full_fg(keep_cols_fg) = b_fg;
    se_full_fg = nan(survival_data.n_feat, 1);
    p_full_fg = nan(survival_data.n_feat, 1);
    se_full_fg(keep_cols_fg) = stats_fg.se;
    p_full_fg(keep_cols_fg) = stats_fg.p;

    results.success = true;
    results.coefficients = b_full_fg;
    results.se = se_full_fg;
    results.p_values = p_full_fg;
    results.hazard_ratios = exp(b_full_fg);  % subdistribution HR (sHR)
    results.loglik = logl_fg;
    results.n_competing = n_competing;
    results.n_primary = n_primary;
    results.fg_weights = fg_weights;
catch ME_fit
    results.message = sprintf('Fine-Gray fit failed: %s', ME_fit.message);
    fprintf('  Fine-Gray model did not converge: %s\n', ME_fit.message);
    return;
end

% --- Stage 5: Print comparison table and CIF plot ---
print_fine_gray_comparison(cox_results, results, survival_data.feat_names);

if ~isempty(config.output_folder)
    try
        plot_cumulative_incidence(survival_data, config);
    catch ME_cif
        fprintf('  ⚠️  CIF plot failed: %s\n', ME_cif.message);
    end
end

fprintf('  Fine-Gray model completed: %d competing events, %d primary events.\n', ...
    n_competing, n_primary);
end
