function G_val = interp_G(t, uniq_times, G_values)
% INTERP_G  Step-function lookup for censoring survival G(t).
    idx = find(uniq_times <= t, 1, 'last');
    if isempty(idx)
        G_val = 1.0;
    else
        G_val = G_values(idx);
    end
end
