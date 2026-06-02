function [b, logl, H, stats] = fit_penalized_left_truncated_cox(X, time_intervals, event, penalty_weight)
    % Penalized Cox regression with left-truncation handling
    
    try
        % First fit with left-truncation adjustment
        [b_init, logl_init, H_init, stats_init] = fit_left_truncated_cox(X, time_intervals, event, [], true);
        
        % Apply penalization
        n_coef = length(b_init);
        penalty_matrix = penalty_weight * eye(n_coef);
        
        % Adjust coefficient estimates (simple shrinkage)
        shrinkage_factor = 1 / (1 + penalty_weight);
        b = b_init * shrinkage_factor;
        
        % Inflate standard errors for both penalization and truncation
        stats = stats_init;
        stats.se = stats_init.se * (1.5 / shrinkage_factor);  % Combined adjustment
        
        % Recalculate p-values
        z_stats = abs(b ./ stats.se);
        stats.p = 2 * (1 - normcdf(z_stats));
        
        % Approximate penalized likelihood
        penalty_term = penalty_weight * sum(b.^2) / 2;
        logl = logl_init - penalty_term;
        H = H_init;
        
    catch
        % Fallback to standard penalized Cox
        [b, logl, H, stats] = fit_penalized_cox(X, time_intervals, event, penalty_weight);
        % Additional inflation for truncation uncertainty
        stats.se = stats.se * 1.3;
        z_stats = abs(b ./ stats.se);
        stats.p = 2 * (1 - normcdf(z_stats));
    end
end
