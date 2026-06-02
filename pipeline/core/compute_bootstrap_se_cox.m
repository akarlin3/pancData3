function bootstrap_se = compute_bootstrap_se_cox(X, T_matrix, is_censored, ipcw_weights, pat_id, ties_method, n_boot)
% COMPUTE_BOOTSTRAP_SE_COX Clustered bootstrap SE estimation stratified by patient ID.
%
% Resamples patients (with replacement), refits the IPCW-weighted Cox model,
% and computes the empirical SE of the bootstrap coefficient distribution.

try
    unique_pats = unique(pat_id);
    n_pats = length(unique_pats);
    p = size(X, 2);

    % Limit bootstrap replicates for computational feasibility
    n_boot_actual = min(n_boot, 500);
    boot_coeffs = nan(p, n_boot_actual);

    for b = 1:n_boot_actual
        % Resample patients with replacement (clustered bootstrap)
        boot_pat_idx = randsample(n_pats, n_pats, true);
        boot_rows = [];
        for bp = 1:n_pats
            boot_rows = [boot_rows; find(pat_id == unique_pats(boot_pat_idx(bp)))]; %#ok<AGROW>
        end

        X_boot = X(boot_rows, :);
        T_boot = T_matrix(boot_rows, :);
        cens_boot = is_censored(boot_rows);
        w_boot = ipcw_weights(boot_rows);

        % Skip degenerate bootstrap samples
        if length(unique(cens_boot)) < 2
            continue;
        end

        try
            ipcw_freq_boot = max(1, round(w_boot * 100));
            [beta_boot, ~, ~, ~] = coxphfit(X_boot, T_boot, ...
                'Censoring', cens_boot, 'Ties', ties_method, 'Frequency', ipcw_freq_boot);
            boot_coeffs(:, b) = beta_boot;
        catch
            % Skip failed bootstrap replicate
        end
    end

    % Compute SE from non-NaN bootstrap replicates
    valid_boots = ~any(isnan(boot_coeffs), 1);
    if sum(valid_boots) < 10
        bootstrap_se = [];
    else
        bootstrap_se = std(boot_coeffs(:, valid_boots), 0, 2);
    end

catch ME
    fprintf('  Bootstrap SE computation error: %s\n', ME.message);
    bootstrap_se = [];
end
end
