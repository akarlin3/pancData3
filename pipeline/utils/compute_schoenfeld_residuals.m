function results = compute_schoenfeld_residuals(X_td, t_start, t_stop, event_csh, beta, feat_names, output_folder, dtype_label)
% COMPUTE_SCHOENFELD_RESIDUALS  Test the proportional hazards assumption.
%
%   Computes scaled Schoenfeld residuals for each covariate in a fitted
%   Cox PH model, regresses them against time, and reports the p-value
%   for the correlation (null hypothesis: PH holds).
%
% Inputs:
%   X_td         - Covariate matrix (n_intervals x p)
%   t_start      - Interval start times
%   t_stop       - Interval stop times
%   event_csh    - CSH event indicator (1=event, 0=censored)
%   beta         - Fitted Cox coefficients (p x 1)
%   feat_names   - Cell array of covariate names
%   output_folder - Directory for saving diagnostic figure (optional)
%   dtype_label  - DWI type label for file naming (optional)
%
% Outputs:
%   results      - Struct with fields:
%                  .rho       - Correlation coefficients (p x 1)
%                  .chi2      - Chi-squared statistics (p x 1)
%                  .p_value   - P-values for PH test (p x 1)
%                  .violated  - Logical vector (p x 1), true if PH violated
%                  .residuals - Schoenfeld residuals matrix (n_events x p)
%                  .event_times - Times at which events occurred

    if nargin < 7, output_folder = ''; end
    if nargin < 8, dtype_label = ''; end

    p = size(X_td, 2);
    event_idx = find(event_csh == 1);
    n_events = length(event_idx);

    results = struct();
    results.rho = nan(p, 1);
    results.chi2 = nan(p, 1);
    results.p_value = nan(p, 1);
    results.violated = false(p, 1);
    results.residuals = [];
    results.event_times = [];

    if n_events < 3 || any(~isfinite(beta))
        fprintf('  Schoenfeld test skipped: insufficient events (%d) or non-finite coefficients.\n', n_events);
        return;
    end

    % Compute Schoenfeld residuals
    % For each event time, the Schoenfeld residual for covariate j is:
    %   x_ij - E[x_j | R_i] = x_ij - sum_{k in R_i} x_kj * w_k / sum_{k in R_i} w_k
    % where w_k = exp(X_k * beta) and R_i is the risk set at event time t_i.

    schoenfeld = nan(n_events, p);
    event_times = nan(n_events, 1);

    lp = X_td * beta;  % linear predictor

    for e = 1:n_events
        idx = event_idx(e);
        t_event = t_stop(idx);
        event_times(e) = t_event;

        % Risk set: intervals that span this event time
        at_risk = (t_start < t_event) & (t_stop >= t_event);
        if sum(at_risk) < 2
            continue;
        end

        % Weights
        w = exp(lp(at_risk));
        w_sum = sum(w);
        if w_sum == 0, continue; end

        % Expected covariate value under the model
        X_risk = X_td(at_risk, :);
        E_x = (w' * X_risk) / w_sum;

        % Schoenfeld residual = observed - expected
        schoenfeld(e, :) = X_td(idx, :) - E_x;
    end

    results.residuals = schoenfeld;
    results.event_times = event_times;

    % Scale residuals by variance (Grambsch-Therneau scaling)
    % and test correlation with time
    valid_events = ~any(isnan(schoenfeld), 2) & isfinite(event_times);
    if sum(valid_events) < 3
        fprintf('  Schoenfeld test: too few valid events for correlation test.\n');
        return;
    end

    t_valid = event_times(valid_events);
    S_valid = schoenfeld(valid_events, :);

    fprintf('\n  --- Schoenfeld Residuals: PH Assumption Test ---\n');
    fprintf('  %-10s  %8s  %8s  %8s  %12s\n', 'Covariate', 'rho', 'chi2', 'p-value', 'PH_violated');
    fprintf('  %s\n', repmat('-', 1, 56));

    any_violated = false;
    for j = 1:p
        resid_j = S_valid(:, j);
        if all(resid_j == 0) || std(resid_j) == 0
            results.rho(j) = 0;
            results.chi2(j) = 0;
            results.p_value(j) = 1;
            results.violated(j) = false;
        else
            % Spearman rank correlation with time
            rho_j = corr(t_valid, resid_j, 'type', 'Spearman');
            results.rho(j) = rho_j;
            % Chi-squared test: chi2 = n * rho^2
            n_v = sum(valid_events);
            results.chi2(j) = n_v * rho_j^2;
            results.p_value(j) = 1 - chi2cdf(results.chi2(j), 1);
            results.violated(j) = results.p_value(j) < 0.05;
        end

        violated_str = '';
        if results.violated(j)
            violated_str = '***';
            any_violated = true;
        end
        if j <= length(feat_names)
            name = feat_names{j};
        else
            name = sprintf('X%d', j);
        end
        fprintf('  %-10s  %8.4f  %8.4f  %8.4f  %12s\n', ...
            name, results.rho(j), results.chi2(j), results.p_value(j), violated_str);
    end

    if any_violated
        violated_names = feat_names(results.violated);
        fprintf('  ⚠️  PH assumption violated for: %s\n', strjoin(violated_names, ', '));
    else
        fprintf('  ✅  PH assumption holds for all covariates (alpha=0.05).\n');
    end

    % Generate diagnostic figure
    if ~isempty(output_folder) && sum(valid_events) >= 3
        try
            fig = figure('Visible', 'off', 'Position', [100 100 1200 300*ceil(p/2)]);
            for j = 1:p
                subplot(ceil(p/2), 2, j);
                resid_j = S_valid(:, j);
                scatter(t_valid, resid_j, 20, 'filled', 'MarkerFaceAlpha', 0.6);
                hold on;

                % LOWESS smoother
                if length(t_valid) >= 5
                    [t_sorted, sort_idx] = sort(t_valid);
                    r_sorted = resid_j(sort_idx);
                    smoothed = lowess_smooth(t_sorted, r_sorted, 0.5);
                    plot(t_sorted, smoothed, 'r-', 'LineWidth', 2);
                end

                yline(0, '--k', 'LineWidth', 0.5);
                if j <= length(feat_names)
                    name = feat_names{j};
                else
                    name = sprintf('X%d', j);
                end
                title(sprintf('%s (p=%.4f)', name, results.p_value(j)));
                xlabel('Time (days)');
                ylabel('Schoenfeld Residual');
                hold off;
            end
            sgtitle('Schoenfeld Residuals vs Time (PH Diagnostics)');

            fig_path = fullfile(output_folder, sprintf('schoenfeld_residuals_%s.png', dtype_label));
            saveas(fig, fig_path);
            close(fig);
            fprintf('  📁 Schoenfeld diagnostic figure saved: %s\n', fig_path);
        catch ME
            fprintf('  ⚠️  Could not generate Schoenfeld figure: %s\n', ME.message);
        end
    end
end


function ys = lowess_smooth(x, y, span)
% Simple LOWESS smoother implementation
    n = length(x);
    ys = zeros(n, 1);
    h = max(3, ceil(span * n));  % At least 3 points for stable regression

    for i = 1:n
        dists = abs(x - x(i));
        [~, idx] = sort(dists);
        idx = idx(1:min(h, n));

        max_dist = max(dists(idx));
        if max_dist < eps
            % All neighbours are at the same x-location; regression is
            % degenerate, so return the weighted mean directly.
            ys(i) = mean(y(idx));
        else
            w = (1 - (dists(idx) / (max_dist + eps)).^3).^3;
            w = max(w, 0);

            if sum(w) == 0
                ys(i) = y(i);
            else
                % Weighted linear regression
                X_local = [ones(length(idx), 1), x(idx)];
                W = diag(w);
                A = X_local' * W * X_local;
                if rcond(A) < eps
                    % Ill-conditioned: fall back to weighted mean
                    ys(i) = sum(w .* y(idx)) / sum(w);
                else
                    b_local = A \ (X_local' * W * y(idx));
                    ys(i) = [1, x(i)] * b_local;
                end
            end
        end
    end
end
