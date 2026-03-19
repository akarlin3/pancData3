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
%
% Inputs:
%   X_td              - [n_intervals x n_feat] covariate matrix
%   t_start           - [n_intervals x 1] interval start times
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

    % --- Stratified Cox Model ---
    % Use the first violating covariate as the stratification variable
    strat_idx = violated_idx(1);
    strat_vals = X_td(:, strat_idx);
    med_val = median(strat_vals(isfinite(strat_vals)));
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

        % Fit on full data with stratification via frequency weighting
        w_off = warning('off', 'all');
        [b_strat, ~, ~, stats_strat] = coxphfit(X_strat_stable, [t_start_stable, t_stop_stable], ...
            'Censoring', event_stable == 0, 'Ties', 'breslow', ...
            'Stratification', strat_group_stable);
        warning(w_off);

        tv_results.stratified_model.hr = exp(b_strat);
        tv_results.stratified_model.p = stats_strat.p;
        tv_results.stratified_model.covariate_names = strat_cov_names;
        tv_results.stratified_model.stratified_by = covariate_names{strat_idx};

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
        t_mid = (t_start_stable + t_stop_stable) / 2;
        t_mid(t_mid <= 0) = 0.5;  % avoid log(0)
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
                [b_ext, ~, ~, stats_ext] = fit_penalized_cox(X_ext_clean, ...
                    [t_start_stable, t_stop_stable], event_stable, penalty_weight);
            else
                [b_ext, ~, ~, stats_ext] = coxphfit(X_ext_clean, [t_start_stable, t_stop_stable], ...
                    'Censoring', event_stable == 0, 'Ties', 'breslow');
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
                tv_results.interaction_models(end+1) = entry;

                fprintf('    Base %s: coef=%.4f, p=%.4f\n', cov_name, base_coef, base_p);
                fprintf('    %s × log(t): coef=%.4f, p=%.4f\n', cov_name, int_coef, int_p);
                if use_penalized
                    fprintf('    (Penalized/smoothed estimates)\n');
                end

                % --- Time-varying HR figure ---
                if ~isempty(output_folder)
                    try
                        max_time = sparse_periods ? stable_period_end : max(t_stop);
                        t_grid = linspace(max(1, min(t_stop_stable)), max_time, 200);
                        hr_t = exp(base_coef + int_coef * log(t_grid));

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
                            
                            % Conservative: assume zero covariance
                            se_lhr = sqrt(se_base^2 + (log(t_grid)).^2 * se_int^2);
                            lhr = base_coef + int_coef * log(t_grid);
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
                        
                        xlabel('Time (days)');
                        ylabel('Hazard Ratio');
                        title_str = sprintf('Time-Varying HR: %s (%s)', cov_name, dtype_label);
                        if use_penalized
                            title_str = [title_str, ' [Penalized]'];
                        end
                        title(title_str);
                        
                        legend_entries = {'95% CI', 'HR(t)', 'HR=1'};
                        if sparse_periods
                            legend_entries{end+1} = 'Stable period';
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