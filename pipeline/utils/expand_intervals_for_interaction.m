function [X_out, t0_out, t1_out, ev_out] = expand_intervals_for_interaction( ...
    X_in, t_start, t_stop, event)
    % Expand counting-process intervals by splitting at event-time
    % quantiles.  This ensures the time-dependent interaction term
    % X × log(t_mid) uses shared interval boundaries rather than each
    % subject's own event/censoring time, preventing endogeneity bias.

    n_events = sum(event);
    if n_events < 4
        X_out = X_in; t0_out = t_start; t1_out = t_stop; ev_out = event;
        return;
    end

    event_ts = sort(t_stop(event == 1));
    n_cuts = min(8, max(2, floor(n_events / 5)));
    cut_q = linspace(0, 1, n_cuts + 2);
    cuts = unique(quantile(event_ts, cut_q(2:end-1)));

    if isempty(cuts)
        X_out = X_in; t0_out = t_start; t1_out = t_stop; ev_out = event;
        return;
    end

    n_obs = numel(event);
    X_cells  = cell(n_obs, 1);
    t0_cells = cell(n_obs, 1);
    t1_cells = cell(n_obs, 1);
    ev_cells = cell(n_obs, 1);

    for ii = 1:n_obs
        s = t_start(ii);  e = t_stop(ii);
        interior = cuts(cuts > s & cuts < e);
        bd = unique([s; interior(:); e]);
        k = numel(bd) - 1;
        X_cells{ii}  = repmat(X_in(ii, :), k, 1);
        t0_cells{ii} = bd(1:end-1);
        t1_cells{ii} = bd(2:end);
        ev_sub = zeros(k, 1);
        ev_sub(end) = event(ii);   % event only in last sub-interval
        ev_cells{ii} = ev_sub;
    end

    X_out  = vertcat(X_cells{:});
    t0_out = vertcat(t0_cells{:});
    t1_out = vertcat(t1_cells{:});
    ev_out = vertcat(ev_cells{:});
end
