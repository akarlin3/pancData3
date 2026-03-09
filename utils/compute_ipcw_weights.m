function ipcw_weights = compute_ipcw_weights(event_td, t_start_td, t_stop_td, X_td_global, pat_id_td)
% COMPUTE_IPCW_WEIGHTS — Inverse probability of censoring weights for Cox PH.
%
%   Corrects for informative administrative censoring using a Cox-based
%   censoring model.  Competing events are excluded from the censoring model
%   (violating IPCW independence would bias hazard ratios).
%
%   Parameters
%   ----------
%   event_td : numeric vector
%       Event indicator: 0=censored, 1=primary event, 2=competing event.
%   t_start_td : numeric vector
%       Interval start times (counting process format).
%   t_stop_td : numeric vector
%       Interval stop times.
%   X_td_global : numeric matrix
%       Covariate matrix (rows x features).
%   pat_id_td : numeric or cell vector
%       Patient IDs for each row (to identify terminal intervals).
%
%   Returns
%   -------
%   ipcw_weights : numeric vector
%       Stabilised IPCW weights (mean-normalised, truncated at 0.05).
%       Returns ones() if no administrative censoring or estimation fails.
%
%   References
%   ----------
%   Robins JM, Finkelstein DM. Correcting for noncompliance and dependent
%   censoring in an AIDS clinical trial with IPCW log-rank tests.
%   Biometrics. 2000;56(3):779-788.
%
%   See also: metrics_survival, remove_constant_columns

    td_n_feat = size(X_td_global, 2);
    ipcw_weights = ones(size(event_td));  % default: unweighted

    % Identify terminal intervals (last row per patient).
    is_terminal_td = false(size(event_td));
    [~, last_idx_td] = unique(pat_id_td, 'last');
    is_terminal_td(last_idx_td) = true;
    is_admin_cens_event = is_terminal_td & (event_td == 0);
    has_admin_cens = any(is_admin_cens_event);

    if ~has_admin_cens
        return;
    end

    try
        % Restrict IPCW model to rows that are NOT competing events.
        not_competing = (event_td ~= 2);
        is_admin_cens_subset = double(is_admin_cens_event(not_competing));
        T_cens = [t_start_td(not_competing), t_stop_td(not_competing)];
        X_cens_subset = X_td_global(not_competing, :);
        [X_cens_clean, keep_ipcw] = remove_constant_columns(X_cens_subset);
        if size(X_cens_clean, 2) == 0
            error('IPCW:NoVariableColumns', 'All covariate columns are constant.');
        end
        w_ipcw = warning('off', 'all');
        [b_cens_short, ~, ~, stats_cens] = coxphfit(X_cens_clean, T_cens, ...
            'Censoring', (is_admin_cens_subset == 0), 'Ties', 'breslow');
        warning(w_ipcw);
        b_cens = zeros(td_n_feat, 1);
        b_cens(keep_ipcw) = b_cens_short;

        % Compute cumulative censoring survival function G(t|X) via Cox.
        lp_cens = X_td_global * b_cens;
        lp_cens_sub = X_cens_subset * b_cens;
        t_start_sub = t_start_td(not_competing);
        t_stop_sub  = t_stop_td(not_competing);
        uniq_times = sort(unique(t_stop_sub));

        % Step 1: Cumulative baseline hazard H0(t) at each event time.
        h0_increments = zeros(length(uniq_times), 1);
        for ui = 1:length(uniq_times)
            t_u = uniq_times(ui);
            at_risk_sub  = (t_start_sub < t_u) & (t_stop_sub >= t_u);
            events_u_sub = at_risk_sub & (t_stop_sub == t_u) & (is_admin_cens_subset == 1);
            if any(events_u_sub) && any(at_risk_sub)
                h0_increments(ui) = sum(events_u_sub) / sum(exp(lp_cens_sub(at_risk_sub)));
            end
        end
        H0_cumulative = cumsum(h0_increments);

        % Step 2: G(t_start | X) using cumulative baseline hazard.
        G_hat = ones(size(t_stop_td));
        for ri = 1:length(t_stop_td)
            ts = t_start_td(ri);
            idx = find(uniq_times <= ts, 1, 'last');
            if ~isempty(idx)
                G_hat(ri) = exp(-H0_cumulative(idx) * exp(lp_cens(ri)));
            end
        end

        % Stabilised weights: truncated and mean-normalised
        G_hat = max(G_hat, 0.05);
        ipcw_weights = 1 ./ G_hat;
        ipcw_weights = ipcw_weights / mean(ipcw_weights);
        fprintf('  IPCW weights applied (admin censoring model, %d competing events excluded). Range: [%.2f, %.2f]\n', ...
            sum(event_td == 2), min(ipcw_weights), max(ipcw_weights));
    catch ME_ipcw
        fprintf('  ⚠️  IPCW weight estimation failed (%s). Proceeding unweighted.\n', ME_ipcw.message);
        ipcw_weights = ones(size(event_td));
    end
end
