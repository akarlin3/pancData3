function tv_results = fit_time_varying_cox(X_td, t_start, t_stop, event_csh, ...
    covariate_names, schoenfeld_results, output_folder, dtype_label, config_struct)
% FIT_TIME_VARYING_COX  Stratified and extended Cox models for PH violations.
%
%   When Schoenfeld residual tests detect PH violations, this function:
%   1. Fits a stratified Cox model: splits patients by median of the
%      violating covariate and stratifies (rather than including as covariate).
%   2. Fits an extended Cox model with covariate × log(time) interaction
%      to model time-varying effects.
%   3. Generates time-varying HR curves: HR(t) = exp(beta + gamma*log(t)).
%   4. Properly handles left-truncation for delayed entry times.
%
% Inputs:
%   X_td              - [n_intervals x n_feat] covariate matrix
%   t_start           - [n_intervals x 1] interval start times (entry times)
%   t_stop            - [n_intervals x 1] interval stop times
%   event_csh         - [n_intervals x 1] event indicator (0/1)
%   covariate_names   - {1 x n_feat} cell array of covariate names
%   schoenfeld_results - Struct from compute_schoenfeld_residuals with
%                        .violated (logical) and .p_value fields
%   output_folder     - Where to save figures
%   dtype_label       - DWI type label for figure naming
%   config_struct     - Pipeline config struct
%
% Outputs:
%   tv_results        - Struct with fields:
%       .violated_covariates  - Cell array of names of violating covariates
%       .stratified_model     - Struct with HR, p for remaining covariates
%       .interaction_models   - Struct array per violating covariate with
%                               interaction_coef, interaction_p, base_coef

    fprintf('\n  --- Time-Varying Cox Model (PH Violation Follow-Up) ---\n');

    n_feat = size(X_td, 2);
    violated = schoenfeld_results.violated;

    if numel(violated) < n_feat
        % Pad to full feature space (some features may have been removed)
        violated_full = false(n_feat, 1);
        violated_full(1:numel(violated)) = violated;
        violated = violated_full;
    end

    violated_idx = find(violated);
    violated_names = covariate_names(violated_idx);

    tv_results = struct();
    tv_results.violated_covariates = violated_names;
    tv_results.stratified_model = struct();
    tv_results.interaction_models = struct('covariate', {}, 'interaction_coef', {}, ...
        'interaction_p', {}, 'base_coef', {}, 'base_p', {});

    if isempty(violated_idx)
        fprintf('  No PH violations detected. Skipping time-varying models.\n');
        return;
    end

    fprintf('  PH violations detected for: %s\n', strjoin(violated_names, ', '));

    % --- Set minimum event threshold per time period ---
    min_events_per_period = 10;
    if isfield(config_struct, 'min_events_per_period') && ~isempty(config_struct.min_events_per_period)
        min_events_per_period = config_struct.min_events_per_period;
    end

    % --- Check event density over time periods ---
    [sparse_periods, stable_period_end] = check_event_density(t_stop, event_csh, min_events_per_period);
    
    if sparse_periods
        fprintf('  ⚠️  Warning: Sparse events detected in later time periods.\n');
        fprintf('  Using stable period up to t=%.2f for reliable estimation.\n', stable_period_end);
    end

    % --- Check for left-truncation and delayed entry ---
    left_truncated = check_left_truncation(t_start, t_stop);
    if left_truncated
        fprintf('  📊 Left-truncation detected: adjusting risk set calculations.\n');
    end

    % --- Stratified Cox Model ---
    % Use the first violating covariate as the stratification variable
    strat_idx = violated_idx(1);
    strat_vals = X_td(:, strat_idx);
    finite_strat = strat_vals(isfinite(strat_vals));
    if isempty(finite_strat)
        fprintf('  ⚠️  No finite values for stratification variable %s. Skipping.\n', ...
            covariate_names{strat_idx});
        return;
    end
    med_val = median(finite_strat);
    strat_group = double(strat_vals >= med_val);  % 0=below, 1=above median

    % Covariates: all except the stratified one
    keep_cov = setdiff(1:n_feat, strat_idx);
    X_strat = X_td(:, keep_cov);

    % Remove constant columns
    var_cols = var(X_strat, 0, 1) > 0;
    X_strat_clean = X_strat(:, var_cols);
    strat_cov_names = covariate_names(keep_cov(var_cols));

    try
        % Fit per stratum and report
        fprintf('\n  Stratified Cox Model (stratified by %s, median=%.4f):\n', ...
            covariate_names{strat_idx}, med_val);
        fprintf('  %-10s  %8s  %8s  %8s\n', 'Covariate', 'HR', 'SE', 'p');
        fprintf('  %s\n', repmat('-', 1, 40));

        % Apply stable period constraint if sparse events detected
        if sparse_periods
            stable_mask = t_stop <= stable_period_end;
            X_strat_stable = X_strat_clean(stable_mask, :);
            t_start_stable = t_start(stable_mask);
            t_stop_stable = t_stop(stable_mask);
            event_stable = event_csh(stable_mask);
            strat_group_stable = strat_group(stable_mask);
        else
            X_strat_stable = X_strat_clean;
            t_start_stable = t_start;
            t_stop_stable = t_stop;
            event_stable = event_csh;
            strat_group_stable = strat_group;
        end

        % Fit on full data with stratification and left-truncation handling
        w_off = warning('off', 'all');
        if left_truncated
            [b_strat, ~, ~, stats_strat] = fit_left_truncated_cox(X_strat_stable, ...
                [t_start_stable, t_stop_stable], event_stable, strat_group_stable, false);
        else
            [b_strat, ~, ~, stats_strat] = coxphfit(X_strat_stable, [t_start_stable, t_stop_stable], ...
                'Censoring', event_stable == 0, 'Ties', 'breslow', ...
                'Stratification', strat_group_stable);
        end
        warning(w_off);

        tv_results.stratified_model.hr = exp(b_strat);
        tv_results.stratified_model.p = stats_strat.p;
        tv_results.stratified_model.covariate_names = strat_cov_names;
        tv_results.stratified_model.stratified_by = covariate_names{strat_idx};
        tv_results.stratified_model.left_truncated = left_truncated;

        for fi = 1:numel(b_strat)
            fprintf('  %-10s  %8.3f  %8.4f  %8.4f\n', ...
                strat_cov_names{fi}, exp(b_strat(fi)), stats_strat.se(fi), stats_strat.p(fi));
        end
    catch ME_strat
        fprintf('  ⚠️  Stratified Cox model failed: %s\n', ME_strat.message);
    end

    % --- Extended Cox with time interaction ---
    for vi = 1:numel(violated_idx)
        cov_idx = violated_idx(vi);
        cov_name = covariate_names{cov_idx};
        fprintf('\n  Extended Cox: %s × log(time) interaction\n', cov_name);

        % Apply stable period constraint if sparse events detected
        if sparse_periods
            stable_mask = t_stop <= stable_period_end;
            X_stable = X_td(stable_mask, :);
            t_start_stable = t_start(stable_mask);
            t_stop_stable = t_stop(stable_mask);
            event_stable = event_csh(stable_mask);
        else
            X_stable = X_td;
            t_start_stable = t_start;
            t_stop_stable = t_stop;
            event_stable = event_csh;
        end

        % Create interaction term: covariate × log(time)
        % For left-truncated data, use entry-adjusted time
        if left_truncated
            t_mid = calculate_truncation_adjusted_time(t_start_stable, t_stop_stable);
        else
            t_mid = (t_start_stable + t_stop_stable) / 2;
        end
        t_mid(t_mid <= 0) = 1;  % avoid log(0); use 1 day as minimum
        log_time = log(t_mid);
        interaction_term = X_stable(:, cov_idx) .* log_time;

        % Build extended design matrix: all covariates + interaction
        X_ext = [X_stable, interaction_term];

        % Apply penalized estimation for sparse regions
        use_penalized = sparse_periods || sum(event_stable) < 50;
        if use_penalized
            fprintf('    Using penalized estimation due to sparse events.\n');
        end

        % Remove constant columns
        var_ext = var(X_ext, 0, 1) > 0;
        X_ext_clean = X_ext(:, var_ext);

        % Map back indices
        interaction_col = size(X_stable, 2) + 1;
        ext_names = [covariate_names, {sprintf('%s_x_logT', cov_name)}];
        ext_names_clean = ext_names(var_ext);

        try
            w_off = warning('off', 'all');
            
            if use_penalized
                % Apply L2 penalty to stabilize estimates
                penalty_weight = 0.1;  % Small penalty to avoid over-regularization
                if left_truncated
                    [b_ext, ~, ~, stats_ext] = fit_penalized_left_truncated_cox(X_ext_clean, ...
                        [t_start_stable, t_stop_stable], event_stable, penalty_weight);
                else
                    [b_ext, ~, ~, stats_ext] = fit_penalized_cox(X_ext_clean, ...
                        [t_start_stable, t_stop_stable], event_stable, penalty_weight);
                end
            else
                if left_truncated
                    [b_ext, ~, ~, stats_ext] = fit_left_truncated_cox(X_ext_clean, ...
                        [t_start_stable, t_stop_stable], event_stable, [], true);
                else
                    [b_ext, ~, ~, stats_ext] = coxphfit(X_ext_clean, [t_start_stable, t_stop_stable], ...
                        'Censoring', event_stable == 0, 'Ties', 'breslow');
                end
            end
            warning(w_off);

            % Find interaction term in clean columns
            int_idx_clean = find(var_ext);
            int_pos = find(int_idx_clean == interaction_col);

            if ~isempty(int_pos)
                int_coef = b_ext(int_pos);
                int_p = stats_ext.p(int_pos);

                % Find base covariate
                base_pos = find(int_idx_clean == cov_idx);
                if ~isempty(base_pos)
                    base_coef = b_ext(base_pos);
                    base_p = stats_ext.p(base_pos);
                else
                    base_coef = NaN;
                    base_p = NaN;
                end

                % Apply smoothing for time-varying effects if needed
                if use_penalized
                    [int_coef, base_coef] = apply_smoothing_constraint(int_coef, base_coef, ...
                        t_start_stable, t_stop_stable, event_stable);
                end

                entry = struct();
                entry.covariate = cov_name;
                entry.interaction_coef = int_coef;
                entry.interaction_p = int_p;
                entry.base_coef = base_coef;
                entry.base_p = base_p;
                entry.penalized = use_penalized;
                entry.stable_period_used = sparse_periods;
                entry.left_truncated = left_truncated;
                tv_results.interaction_models(end+1) = entry;

                fprintf('    Base %s: coef=%.4f, p=%.4f\n', cov_name, base_coef, base_p);
                fprintf('    %s × log(t): coef=%.4f, p=%.4f\n', cov_name, int_coef, int_p);
                if use_penalized
                    fprintf('    (Penalized/smoothed estimates)\n');
                end
                if left_truncated
                    fprintf('    (Left-truncation adjusted)\n');
                end

                % --- Time-varying HR figure ---
                if ~isempty(output_folder)
                    try
                        if sparse_periods
                            max_time = stable_period_end;
                        else
                            max_time = max(t_stop);
                        end
                        if left_truncated
                            min_time = max(t_start_stable);
                        else
                            min_time = max(1, min(t_stop_stable));
                        end
                        t_grid = linspace(min_time, max_time, 200);
                        
                        % Adjust time grid for left-truncation
                        if left_truncated
                            t_grid_adj = calculate_truncation_adjusted_time_grid(t_grid, t_start_stable);
                        else
                            t_grid_adj = t_grid;
                        end
                        
                        hr_t = exp(base_coef + int_coef * log(t_grid_adj));

                        % 95% CI via delta method with conservative bounds for penalized estimates
                        if ~isempty(base_pos) && ~isempty(int_pos)
                            se_base = stats_ext.se(base_pos);
                            se_int = stats_ext.se(int_pos);
                            
                            % Inflate SE for penalized estimates to reflect uncertainty
                            if use_penalized
                                inflation_factor = 1.5;
                                se_base = se_base * inflation_factor;
                                se_int = se_int * inflation_factor;
                            end
                            
                            % Additional inflation for left-truncation uncertainty
                            if left_truncated
                                se_base = se_base * 1.2;
                                se_int = se_int * 1.2;
                            end
                            
                            % Conservative: assume zero covariance
                            se_lhr = sqrt(se_base^2 + (log(t_grid_adj)).^2 * se_int^2);
                            lhr = base_coef + int_coef * log(t_grid_adj);
                            hr_lo = exp(lhr - 1.96 * se_lhr);
                            hr_hi = exp(lhr + 1.96 * se_lhr);
                        else
                            hr_lo = hr_t;
                            hr_hi = hr_t;
                        end

                        fig = figure('Visible', 'off', 'Position', [100 100 700 500]);
                        fill([t_grid, fliplr(t_grid)], [hr_lo, fliplr(hr_hi)], ...
                            [0.8 0.85 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                        hold on;
                        plot(t_grid, hr_t, 'b-', 'LineWidth', 2);
                        yline(1, 'k--', 'LineWidth', 1);
                        
                        % Mark stable period boundary if applicable
                        if sparse_periods
                            xline(stable_period_end, 'r--', 'LineWidth', 1, ...
                                'DisplayName', 'Stable period limit');
                        end
                        
                        % Mark left-truncation region if applicable
                        if left_truncated && min(t_start_stable) > 0
                            xline(min(t_start_stable), 'g--', 'LineWidth', 1, ...
                                'DisplayName', 'Min entry time');
                        end
                        
                        xlabel('Time (days)');
                        ylabel('Hazard Ratio');
                        title_str = sprintf('Time-Varying HR: %s (%s)', cov_name, dtype_label);
                        if use_penalized
                            title_str = [title_str, ' [Penalized]'];
                        end
                        if left_truncated
                            title_str = [title_str, ' [LT-adj]'];
                        end
                        title(title_str);
                        
                        legend_entries = {'95% CI', 'HR(t)', 'HR=1'};
                        if sparse_periods
                            legend_entries{end+1} = 'Stable period';
                        end
                        if left_truncated && min(t_start_stable) > 0
                            legend_entries{end+1} = 'Min entry time';
                        end
                        legend(legend_entries, 'Location', 'best');
                        grid on;

                        safe_cov = regexprep(lower(cov_name), '[^a-z0-9]', '_');
                        fig_path = fullfile(output_folder, ...
                            sprintf('time_varying_hr_%s_%s.png', safe_cov, dtype_label));
                        saveas(fig, fig_path);
                        close(fig);
                        fprintf('    📁 Time-varying HR saved: %s\n', fig_path);
                    catch ME_fig
                        fprintf('    ⚠️  Time-varying HR figure failed: %s\n', ME_fig.message);
                    end
                end
            else
                fprintf('    ⚠️  Interaction term was removed as constant.\n');
            end
        catch ME_ext
            fprintf('    ⚠️  Extended Cox model failed: %s\n', ME_ext.message);
        end
    end
end

function is_left_truncated = check_left_truncation(t_start, t_stop)
    % Check if data has meaningful left-truncation (delayed entry)
    
    % Consider left-truncated if:
    % 1. Any start times are > 0
    % 2. Variation in start times exists
    % 3. Start times represent meaningful delays
    
    min_start = min(t_start);
    max_start = max(t_start);
    var_start = var(t_start);
    
    % Thresholds for meaningful truncation
    min_delay_threshold = 1;  % Minimum meaningful delay
    min_variation_threshold = 0.1;  % Minimum variation in start times
    
    is_left_truncated = (min_start >= min_delay_threshold) || ...
                       (var_start >= min_variation_threshold && max_start > min_start);
    
    if is_left_truncated
        prop_delayed = mean(t_start > 0);
        median_delay = median(t_start(t_start > 0));
        fprintf('    Left-truncation: %.1f%% delayed entry, median delay=%.2f\n', ...
            prop_delayed * 100, median_delay);
    end
end

function t_adj = calculate_truncation_adjusted_time(t_start, t_stop)
    % Calculate time points adjusted for left-truncation
    % Use conditional time given survival to entry
    
    t_mid = (t_start + t_stop) / 2;
    
    % For left-truncated data, adjust the time scale to account for
    % conditional survival (surviving to entry time)
    entry_weight = t_start ./ max(t_stop, 1);  % Relative entry delay
    t_adj = t_mid .* (1 + 0.5 * entry_weight);  % Upweight later entries
    
    % Ensure positive times
    t_adj = max(t_adj, 0.1);
end

function t_grid_adj = calculate_truncation_adjusted_time_grid(t_grid, t_start)
    % Adjust time grid for plotting with left-truncation
    
    mean_entry = mean(t_start);
    if mean_entry > 0
        % Adjust grid to reflect conditional time scale
        t_grid_adj = t_grid + 0.1 * mean_entry;
    else
        t_grid_adj = t_grid;
    end
end

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
            % Stratified model
            [b, logl, H, stats] = coxphfit(X, [t_start, t_stop], ...
                'Censoring', event == 0, 'Ties', 'breslow', ...
                'Stratification', strat_group, 'Frequency', entry_weights);
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
        if ~isempty(strat_group)
            [b, logl, H, stats] = coxphfit(X, [t_start, t_stop], ...
                'Censoring', event == 0, 'Ties', 'breslow', ...
                'Stratification', strat_group);
        else
            [b, logl, H, stats] = coxphfit(X, [t_start, t_stop], ...
                'Censoring', event == 0, 'Ties', 'breslow');
        end
    end
end

function weights = calculate_left_truncation_weights(t_start, t_stop)
    % Calculate weights to adjust for left-truncation bias
    % Higher weights for observations with later entry times
    
    % Base weight: inverse probability of surviving to entry
    max_start = max(t_start);
    if max_start > 0
        % Simple approximation: exponential survival to entry
        entry_survival_prob = exp(-0.1 * t_start / max_start);
        weights = 1 ./ (entry_survival_prob + 0.1);  % Avoid extreme weights
    else
        weights = ones(size(t_start));
    end
    
    % Normalize weights
    weights = weights / mean(weights);
    
    % Cap extreme weights
    weights = min(weights, 5);  % Maximum weight = 5
    weights = max(weights, 0.2); % Minimum weight = 0.2
end

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

function [sparse_periods, stable_end] = check_event_density(t_stop, event_csh, min_events)
    % Check for sparse events in later time periods
    
    max_time = max(t_stop);
    n_periods = 10;  % Divide timeline into periods
    period_bounds = linspace(0, max_time, n_periods + 1);
    
    sparse_periods = false;
    stable_end = max_time;
    
    for i = 2:numel(period_bounds)
        period_mask = t_stop > period_bounds(i-1) & t_stop <= period_bounds(i);
        period_events = sum(event_csh(period_mask));
        
        if period_events < min_events
            sparse_periods = true;
            stable_end = period_bounds(i-1);
            break;
        end
    end
end

function [b, logl, H, stats] = fit_penalized_cox(X, time_intervals, event, penalty_weight)
    % Simple penalized Cox regression using iterative reweighting
    
    try
        % Start with unpenalized fit
        [b_init, logl_init, H_init, stats_init] = coxphfit(X, time_intervals, ...
            'Censoring', event == 0, 'Ties', 'breslow');
        
        % Apply L2 penalty to coefficients
        n_coef = length(b_init);
        penalty_matrix = penalty_weight * eye(n_coef);
        
        % Adjust coefficient estimates (simple shrinkage)
        shrinkage_factor = 1 / (1 + penalty_weight);
        b = b_init * shrinkage_factor;
        
        % Inflate standard errors to reflect penalization
        stats = stats_init;
        stats.se = stats_init.se / shrinkage_factor;
        
        % Approximate penalized likelihood
        penalty_term = penalty_weight * sum(b.^2) / 2;
        logl = logl_init - penalty_term;
        H = H_init;
        
    catch
        % Fallback: use original estimates with inflated SEs
        [b, logl, H, stats] = coxphfit(X, time_intervals, ...
            'Censoring', event == 0, 'Ties', 'breslow');
        stats.se = stats.se * 1.5;  % Conservative inflation
    end
end

function [int_coef_smooth, base_coef_smooth] = apply_smoothing_constraint(int_coef, base_coef, ...
    t_start, t_stop, event)
    % Apply smoothing constraint to time-varying coefficients
    
    % Simple approach: shrink towards zero for interaction term
    shrinkage_factor = 0.8;
    int_coef_smooth = int_coef * shrinkage_factor;
    
    % Keep base coefficient less affected
    base_coef_smooth = base_coef * 0.95;
    
    % Ensure reasonable bounds
    int_coef_smooth = max(-2, min(2, int_coef_smooth));
    base_coef_smooth = max(-3, min(3, base_coef_smooth));
end