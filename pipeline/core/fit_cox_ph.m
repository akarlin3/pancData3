function cox_results = fit_cox_ph(survival_data, config)
% FIT_COX_PH Fit Cox Proportional Hazards model with time-dependent covariates
% Uses robust sandwich (Huber-White) variance estimator for IPCW-weighted models.

if ~survival_data.has_sufficient_data
    cox_results = struct('success', false, 'message', 'Insufficient data');
    return;
end

% Scale covariates
X_scaled = scale_td_panel(survival_data.X, survival_data.feat_names, ...
    survival_data.pat_id, survival_data.t_start, unique(survival_data.pat_id), 'baseline');

% Prepare for cause-specific hazard
event_csh = survival_data.event;
event_csh(event_csh == 2) = 0;  % Competing risks → censored

% Compute IPCW weights
ipcw_weights = compute_ipcw_weights(survival_data.event, survival_data.t_start, ...
    survival_data.t_stop, X_scaled, survival_data.pat_id);

% Check for tied event times
event_times = survival_data.t_stop(event_csh == 1);
[unique_times, ~, time_idx] = unique(event_times);
has_tied_times = any(arrayfun(@(i) sum(time_idx == i) > 1, 1:length(unique_times)));
if has_tied_times
    ties_method = 'efron';
else
    ties_method = 'breslow';
end

if has_tied_times
    fprintf('  ⚠️  Detected tied event times. Using Efron approximation.\n');
end

% Fit Cox model with robust sandwich variance estimation
try
    T_matrix = [survival_data.t_start, survival_data.t_stop];
    is_censored = (event_csh == 0);

    [X_clean, keep_cols] = remove_constant_columns(X_scaled);
    if size(X_clean, 2) == 0
        error('Cox:NoVariableColumns', 'All covariate columns are constant.');
    end

    % --- Step 1: Fit unweighted Cox model to get initial estimates ---
    % We first fit the unweighted model, then compute the sandwich variance
    % using the IPCW weights and score residuals.

    % Suppress only the specific coxphfit warnings, with guaranteed restoration
    ws = warning('off', 'stats:coxphfit:FitWarning');
    cleanupWarn = onCleanup(@() warning(ws));

    % Apply IPCW via frequency weights for point estimation
    ipcw_scale = 100;
    ipcw_freq = max(1, round(ipcw_weights * ipcw_scale));

    [b_short, logl, H_out, stats_short] = coxphfit(X_clean, T_matrix, ...
        'Censoring', is_censored, 'Ties', ties_method, 'Frequency', ipcw_freq);

    % --- Step 2: Compute robust sandwich (Huber-White) standard errors ---
    % The sandwich variance is: V_robust = A^{-1} B A^{-1}
    % where A = negative Hessian (observed information), B = meat matrix
    % B = sum over clusters (patients) of (score_i * w_i)' * (score_i * w_i)

    fprintf('  Computing robust sandwich variance estimator for IPCW weights.\n');

    sandwich_se = compute_sandwich_se_cox(b_short, X_clean, T_matrix, ...
        is_censored, ipcw_weights, survival_data.pat_id, ties_method);

    if ~isempty(sandwich_se) && all(isfinite(sandwich_se)) && all(sandwich_se > 0)
        % Use sandwich SE
        stats_short.se = sandwich_se;
        fprintf('  Using robust sandwich SE estimator.\n');
    else
        % Fallback: use clustered bootstrap SE stratified by patient ID
        fprintf('  Sandwich SE computation failed; falling back to clustered bootstrap.\n');
        bootstrap_se = compute_bootstrap_se_cox(X_clean, T_matrix, is_censored, ...
            ipcw_weights, survival_data.pat_id, ties_method, config.n_bootstrap);

        if ~isempty(bootstrap_se) && all(isfinite(bootstrap_se)) && all(bootstrap_se > 0)
            stats_short.se = bootstrap_se;
            fprintf('  Using clustered bootstrap SE (%d replicates).\n', config.n_bootstrap);
        else
            % Last resort: keep model-based SE with a warning
            fprintf('  ⚠️  Both sandwich and bootstrap SE failed; using model-based SE (interpret with caution).\n');
        end
    end

    % Recompute p-values with corrected SE
    stats_short.p = 2 * (1 - normcdf(abs(b_short ./ stats_short.se)));

    % Map back to full feature space
    b_full = zeros(survival_data.n_feat, 1);
    b_full(keep_cols) = b_short;
    se_full = nan(survival_data.n_feat, 1);
    p_full = nan(survival_data.n_feat, 1);
    se_full(keep_cols) = stats_short.se;
    p_full(keep_cols) = stats_short.p;

    cox_results = struct();
    cox_results.success = true;
    cox_results.coefficients = b_full;
    cox_results.se = se_full;
    cox_results.p_values = p_full;
    cox_results.hazard_ratios = exp(b_full);
    cox_results.loglik = logl;
    cox_results.keep_cols = keep_cols;
    cox_results.X_scaled = X_scaled;
    cox_results.event_csh = event_csh;
    cox_results.ipcw_weights = ipcw_weights;

    % Print results
    print_cox_results(cox_results, survival_data.feat_names);

catch ME
    if contains(ME.identifier, 'FitWarning') || contains(ME.identifier, 'IterationLimit')
        fprintf('  ⚠️  Cox model convergence issue: %s\n', ME.message);
        cox_results = create_failed_cox_results(survival_data.n_feat);
    else
        rethrow(ME);
    end
end
end
