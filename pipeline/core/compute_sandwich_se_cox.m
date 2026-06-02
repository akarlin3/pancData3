function sandwich_se = compute_sandwich_se_cox(beta, X, T_matrix, is_censored, ipcw_weights, pat_id, ties_method)
% COMPUTE_SANDWICH_SE_COX Compute robust sandwich (Huber-White) SE for IPCW-weighted Cox model.
%
% The sandwich estimator accounts for the IPCW weighting by computing:
%   V_sandwich = A^{-1} * B * A^{-1}
% where A is the observed information (negative Hessian) and B is the "meat"
% matrix formed from weighted score residuals clustered by patient.

try
    n = size(X, 1);
    p = size(X, 2);
    t_stop = T_matrix(:, 2);
    event = double(~is_censored);  % 1 = event, 0 = censored

    % Compute linear predictor and risk scores
    eta = X * beta;
    risk = exp(eta);
    w = ipcw_weights(:);  % per-observation IPCW weights

    % Get unique ordered event times
    event_times = sort(unique(t_stop(event == 1)));
    n_events = length(event_times);

    if n_events == 0
        sandwich_se = [];
        return;
    end

    % --- Compute score residuals for each observation ---
    % Score residual for observation i:
    %   U_i = delta_i * w_i * (x_i - xbar(t_i))
    %         - w_i * sum_{j: t_j <= t_i, delta_j=1} [ w_j * risk_i * (x_i - xbar(t_j)) / S0(t_j) ] * I(i in risk set at t_j)
    % where S0(t) = sum_{k in R(t)} w_k * risk_k
    %       xbar(t) = sum_{k in R(t)} w_k * risk_k * x_k / S0(t)

    score_resid = zeros(n, p);

    % Precompute risk set quantities at each event time
    S0 = zeros(n_events, 1);
    S1 = zeros(n_events, p);

    for j = 1:n_events
        tj = event_times(j);
        in_risk = (T_matrix(:, 1) < tj) & (t_stop >= tj);
        if ~any(in_risk), continue; end

        w_risk = w(in_risk) .* risk(in_risk);
        S0(j) = sum(w_risk);
        S1(j, :) = sum(bsxfun(@times, X(in_risk, :), w_risk), 1);
    end

    % Avoid division by zero
    S0(S0 == 0) = eps;
    xbar = bsxfun(@rdivide, S1, S0);  % n_events x p

    % For each observation, accumulate score residual contributions
    for i = 1:n
        ti = t_stop(i);

        % Event contribution (if this observation had an event)
        if event(i) == 1
            % Find which event time this corresponds to
            eidx = find(event_times == ti, 1);
            if ~isempty(eidx)
                score_resid(i, :) = score_resid(i, :) + w(i) * (X(i, :) - xbar(eidx, :));
            end
        end

        % At-risk contribution: for each event time t_j where i is in the risk set
        for j = 1:n_events
            tj = event_times(j);
            if T_matrix(i, 1) < tj && t_stop(i) >= tj
                % Number of events at this time (for ties)
                n_events_at_tj = sum(event == 1 & t_stop == tj);

                % Contribution from being in risk set
                contrib = w(i) * risk(i) * n_events_at_tj * (X(i, :) - xbar(j, :)) / S0(j);
                score_resid(i, :) = score_resid(i, :) - w(i) * contrib;
            end
        end
    end

    % --- Compute the Hessian (observed information matrix A) ---
    % A = sum over event times of [ S2(t)/S0(t) - (S1(t)/S0(t))' * (S1(t)/S0(t)) ]
    % weighted by number of events at each time
    A = zeros(p, p);
    for j = 1:n_events
        tj = event_times(j);
        in_risk = (T_matrix(:, 1) < tj) & (t_stop >= tj);
        if ~any(in_risk), continue; end

        w_risk = w(in_risk) .* risk(in_risk);
        X_risk = X(in_risk, :);

        % S2 matrix
        S2 = (bsxfun(@times, X_risk, w_risk))' * X_risk;

        n_events_at_tj = sum(event == 1 & t_stop == tj);

        A = A + n_events_at_tj * (S2 / S0(j) - xbar(j, :)' * xbar(j, :));
    end

    % --- Compute meat matrix B, clustered by patient ---
    unique_pats = unique(pat_id);
    B = zeros(p, p);
    for k = 1:length(unique_pats)
        idx_k = (pat_id == unique_pats(k));
        U_k = sum(score_resid(idx_k, :), 1);  % 1 x p: sum score residuals within cluster
        B = B + U_k' * U_k;
    end

    % --- Sandwich covariance: A^{-1} * B * A^{-1} ---
    A_inv = pinv(A);  % Use pseudoinverse for numerical stability
    V_sandwich = A_inv * B * A_inv;

    % Extract SE from diagonal
    sandwich_var = diag(V_sandwich);

    if any(sandwich_var < 0)
        fprintf('  ⚠️  Negative sandwich variance detected; attempting nearest PSD correction.\n');
        % Force positive semi-definiteness
        [V_eig, D_eig] = eig(V_sandwich);
        D_eig = diag(max(diag(D_eig), 0));
        V_sandwich = V_eig * D_eig * V_eig';
        sandwich_var = diag(V_sandwich);
    end

    sandwich_se = sqrt(sandwich_var);

    if any(~isfinite(sandwich_se))
        sandwich_se = [];
    end

catch ME
    fprintf('  Sandwich SE computation error: %s\n', ME.message);
    sandwich_se = [];
end
end
