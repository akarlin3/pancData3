function [b, logl, H, stats] = fit_left_truncated_cox(X, time_intervals, event, strat_group, is_interaction)
    % Fit Cox model with proper left-truncation handling
    
    t_start = time_intervals(:, 1);
    t_stop = time_intervals(:, 2);
    
    try
        % Adjust likelihood for left-truncation using counting process approach
        % This is an approximation - full implementation would require custom likelihood
        
        % Weight observations by inverse of entry probability
        entry_weights = calculate_left_truncation_weights(t_start, t_stop);
        
        if ~isempty(strat_group)
            % Stratified model via per-stratum fitting
            strata_lvls = unique(strat_group);
            n_coef = size(X, 2);
            b = zeros(n_coef, 1); se_acc = zeros(n_coef, 1); tw = zeros(n_coef, 1);
            logl = 0; H = [];
            for ssi = 1:numel(strata_lvls)
                sm = strat_group == strata_lvls(ssi);
                if sum(event(sm)) < 2, continue; end
                [b_s, ll_s, H_s, st_s] = coxphfit(X(sm, :), [t_start(sm), t_stop(sm)], ...
                    'Censoring', event(sm) == 0, 'Ties', 'breslow', ...
                    'Frequency', entry_weights(sm));
                ws = 1 ./ (st_s.se.^2 + eps);
                b = b + b_s .* ws; tw = tw + ws; logl = logl + ll_s;
                if isempty(H), H = H_s; end
            end
            ok = tw > 0; b(ok) = b(ok) ./ tw(ok);
            se_acc(ok) = 1 ./ sqrt(tw(ok));
            z_acc = abs(b) ./ (se_acc + eps);
            stats = struct('se', se_acc, 'p', 2 * (1 - normcdf(z_acc)));
        else
            % Standard model with weights
            [b, logl, H, stats] = coxphfit(X, [t_start, t_stop], ...
                'Censoring', event == 0, 'Ties', 'breslow', ...
                'Frequency', entry_weights);
        end
        
        % Adjust standard errors for left-truncation
        if is_interaction && isfield(stats, 'se')
            % Inflate SEs to account for truncation uncertainty
            truncation_inflation = 1 + 0.2 * (std(t_start) / (mean(t_stop) + eps));
            stats.se = stats.se * truncation_inflation;
            
            % Recalculate p-values with adjusted SEs
            z_stats = abs(b ./ stats.se);
            stats.p = 2 * (1 - normcdf(z_stats));
        end
        
    catch ME
        % Fallback to standard Cox if left-truncation adjustment fails
        warning('Left-truncation adjustment failed: %s. Using standard Cox.', ME.message);
        % Fallback: ignore stratification, fit pooled model
        [b, logl, H, stats] = coxphfit(X, [t_start, t_stop], ...
            'Censoring', event == 0, 'Ties', 'breslow');
    end
end
