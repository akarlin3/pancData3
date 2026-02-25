% test_landmark_cindex_mock.m
% Generate some mock data for the Landmark validation logic
rng(42);
clc;
n_vp = 40;
valid_pts = true(n_vp, 1);
td_n_feat = 4;
td_scan_days = [0, 5, 10, 15, 20, 90]; 
nTp = length(td_scan_days);

ADC_abs = rand(n_vp, nTp) * 2;
D_abs   = rand(n_vp, nTp) * 1.5;
f_abs   = rand(n_vp, nTp);
Dstar_abs = rand(n_vp, nTp) * 0.5;

td_feat_arrays = {ADC_abs, D_abs, f_abs, Dstar_abs};
td_feat_names  = {'ADC', 'D', 'f', 'D*'};
td_tot_time = randi([10, 200], n_vp, 1);
td_lf = randi([0, 1], n_vp, 1);
m_id_list = arrayfun(@(x) sprintf('Pt%d', x), 1:n_vp, 'UniformOutput', false);

% minimal provenance
dl_provenance = struct();
dtype = 1;

[X_td, t_start_td, t_stop_td, event_td, pat_id_td, frac_td] = ...
    build_td_panel(td_feat_arrays, td_feat_names, td_lf, td_tot_time, nTp, td_scan_days);

%% ---------- Time-Dependent IPCW Concordance ----------
fprintf('\n--- TIME-DEPENDENT IPCW CONCORDANCE (Out-of-Fold Validation) ---\n');

td_n_feat = numel(td_feat_arrays);

surv_time_all = td_tot_time; 
surv_event_all = td_lf;
valid_idx_pts = find(valid_pts);
n_vp = length(valid_idx_pts);

oof_risk_history = cell(n_vp, 1);

for p_idx = 1:n_vp
    pt_id = m_id_list{valid_idx_pts(p_idx)};
    
    is_leaky = false;
    if dtype == 2 && isfield(dl_provenance, 'dncnn_train_ids') && any(strcmp(dl_provenance.dncnn_train_ids, pt_id))
        is_leaky = true;
    elseif dtype == 3 && isfield(dl_provenance, 'ivimnet_train_ids') && any(strcmp(dl_provenance.ivimnet_train_ids, pt_id))
        is_leaky = true;
    end
    if is_leaky
        continue;
    end
    
    train_mask = (pat_id_td ~= p_idx);
    
    % Scale specifically for this LOOCV fold
    train_pat_ids = unique(pat_id_td(train_mask));
    X_td_scaled = scale_td_panel(X_td, td_feat_names, pat_id_td, frac_td, train_pat_ids);

    X_train = X_td_scaled(train_mask, :);
    T_start_train = t_start_td(train_mask);
    T_stop_train = t_stop_td(train_mask);
    E_train = event_td(train_mask);
    
    test_mask = (pat_id_td == p_idx);
    X_test = X_td_scaled(test_mask, :);
    T_start_test = t_start_td(test_mask);
    T_stop_test = t_stop_td(test_mask);
    
    if sum(test_mask) == 0
        continue;
    end
    
    w_state = warning('off', 'all');
    try
        [b_loo, ~] = coxphfit(X_train, [T_start_train, T_stop_train], 'Censoring', ~E_train, 'Ties', 'breslow', 'Options', statset('MaxIter', 100));
    catch
        try
            unique_train_pats = unique(pat_id_td(train_mask));
            X_last = nan(length(unique_train_pats), td_n_feat);
            E_last = nan(length(unique_train_pats), 1);
            T_train_pat = pat_id_td(train_mask);
            for ui = 1:length(unique_train_pats)
                rows = find(T_train_pat == unique_train_pats(ui));
                X_last(ui, :) = X_train(rows(end), :);
                E_last(ui) = E_train(rows(end));
            end
            mdl_firth = fitglm(X_last, E_last, 'Distribution', 'binomial', 'LikelihoodPenalty', 'jeffreys-prior');
            b_loo = mdl_firth.Coefficients.Estimate(2:end);
        catch
            b_loo = nan(td_n_feat, 1);
        end
    end
    warning(w_state);
    
    oof_risk_history{p_idx} = struct('risk', X_test * b_loo, 't_start', T_start_test, 't_stop', T_stop_test);
    
    train_pat_mask = true(n_vp, 1);
    train_pat_mask(p_idx) = false;
    train_times = surv_time_all(train_pat_mask);
    train_cens = surv_event_all(train_pat_mask);
    
    % Covariate-adjusted conditional censoring model
    Z_train = ADC_abs(train_pat_mask, 1);
    Z_test  = ADC_abs(p_idx, 1);
    cens_model_events = (train_cens == 0); % 1 if censored, 0 if failed
    
    w_state = warning('off', 'all');
    try
        [b_cens, ~, H_cens, ~] = coxphfit(Z_train, train_times, 'Censoring', ~cens_model_events, 'Ties', 'breslow', 'Options', statset('MaxIter', 100));
        if isempty(H_cens)
            error('Empty hazard');
        end
    catch
        b_cens = 0;
        [F_cens_km, x_cens_km] = ecdf(train_times, 'Censoring', train_cens);
        H_cens = [x_cens_km, -log(1 - F_cens_km)];
    end
    warning(w_state);
    
    oof_risk_history{p_idx}.cens_b = b_cens;
    oof_risk_history{p_idx}.cens_H = H_cens;
    oof_risk_history{p_idx}.Z_test = Z_test;
end

concordant = 0; discordant = 0; weights_sum = 0;

for i = 1:n_vp
    if surv_event_all(i) == 1 && ~isempty(oof_risk_history{i})
        T_i = surv_time_all(i);
        
        risk_i = NaN;
        idx_i = find(T_i > oof_risk_history{i}.t_start & T_i <= oof_risk_history{i}.t_stop, 1, 'last');
        if isempty(idx_i) && T_i == oof_risk_history{i}.t_start(1)
            idx_i = 1;
        end
        if ~isempty(idx_i)
            risk_i = oof_risk_history{i}.risk(idx_i);
        end
        
        if isnan(risk_i)
            continue;
        end
        
        H_c = oof_risk_history{i}.cens_H;
        b_c = oof_risk_history{i}.cens_b;
        Z_test_i = oof_risk_history{i}.Z_test;
        
        idx_G = find(H_c(:,1) <= T_i, 1, 'last');
        if isempty(idx_G)
            G_Ti = 1;
        else
            H_0_t = H_c(idx_G, 2);
            G_Ti = exp(-H_0_t * exp(Z_test_i * b_c));
        end
        
        if G_Ti <= 0 || isnan(G_Ti)
            G_Ti = 1e-5;
        end
        W_i = 1 / (G_Ti^2);
        
        for j = 1:n_vp
            if i ~= j && surv_time_all(j) > T_i && ~isempty(oof_risk_history{j})
                risk_j = NaN;
                idx_j = find(T_i > oof_risk_history{j}.t_start & T_i <= oof_risk_history{j}.t_stop, 1, 'last');
                if isempty(idx_j) && T_i == oof_risk_history{j}.t_start(1)
                    idx_j = 1;
                end
                if ~isempty(idx_j)
                    risk_j = oof_risk_history{j}.risk(idx_j);
                end
                
                if ~isnan(risk_j)
                    if risk_i > risk_j
                        concordant = concordant + W_i;
                    elseif risk_i < risk_j
                        discordant = discordant + W_i;
                    else
                        concordant = concordant + 0.5 * W_i;
                        discordant = discordant + 0.5 * W_i;
                    end
                    weights_sum = weights_sum + W_i;
                end
            end
        end
    end
end

if weights_sum > 0
    IPCW_C = concordant / weights_sum;
    fprintf('  Continuous Time-Dependent IPCW C-index: %.4f\n\n', IPCW_C);
else
    fprintf('  Cannot compute IPCW C-index (weights=0 or no comparable pairs).\n\n');
end
