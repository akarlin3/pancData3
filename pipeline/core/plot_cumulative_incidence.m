function plot_cumulative_incidence(survival_data, config)
% PLOT_CUMULATIVE_INCIDENCE  Non-parametric CIF plot via Aalen-Johansen.

% Get per-patient terminal outcomes
[unique_pats, ~, pat_group] = unique(survival_data.pat_id);
n_pats = length(unique_pats);
pat_time = zeros(n_pats, 1);
pat_event = zeros(n_pats, 1);
for pi = 1:n_pats
    rows = find(pat_group == pi);
    [~, max_idx] = max(survival_data.t_stop(rows));
    pat_time(pi) = survival_data.t_stop(rows(max_idx));
    pat_event(pi) = survival_data.event(rows(max_idx));
end

[sorted_t, si] = sort(pat_time);
sorted_ev = pat_event(si);
uniq_t = unique(sorted_t);

% Aalen-Johansen estimator
S_prev = 1.0;
CIF_primary = zeros(length(uniq_t), 1);
CIF_competing = zeros(length(uniq_t), 1);
cum_primary = 0;
cum_competing = 0;

for k = 1:length(uniq_t)
    tk = uniq_t(k);
    n_risk = sum(sorted_t >= tk);
    d_primary = sum(sorted_t == tk & sorted_ev == 1);
    d_competing = sum(sorted_t == tk & sorted_ev == 2);

    if n_risk > 0
        h_primary = d_primary / n_risk;
        h_competing = d_competing / n_risk;
        cum_primary = cum_primary + S_prev * h_primary;
        cum_competing = cum_competing + S_prev * h_competing;
        S_prev = S_prev * (1 - (d_primary + d_competing) / n_risk);
    end

    CIF_primary(k) = cum_primary;
    CIF_competing(k) = cum_competing;
end

fig = figure('Visible', 'off', 'Position', [100 100 700 500]);
stairs(uniq_t, CIF_primary, 'b-', 'LineWidth', 2);
hold on;
stairs(uniq_t, CIF_competing, 'r--', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Cumulative Incidence');
title(sprintf('Cumulative Incidence Functions (%s)', config.dtype_label));
legend('Local Failure', 'Competing Events', 'Location', 'northwest');
grid on;
ylim([0, 1]);

fig_path = fullfile(config.output_folder, sprintf('cif_plot_%s.png', config.dtype_label));
saveas(fig, fig_path);
close(fig);
fprintf('  📁 CIF plot saved: %s\n', fig_path);
end
