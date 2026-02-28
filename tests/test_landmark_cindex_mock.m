% test_landmark_cindex_mock.m
% Generate some mock data for the Landmark validation logic

% Add necessary paths
baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(baseDir, 'core'));
addpath(fullfile(baseDir, 'utils'));
addpath(fullfile(baseDir, 'dependencies'));

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
% Event indicator: 0 = censored, 1 = disease progression (event of interest),
% 2 = competing risk (e.g., non-cancer death).
td_lf = randi([0, 2], n_vp, 1); 
m_id_list = arrayfun(@(x) sprintf('Pt%d', x), 1:n_vp, 'UniformOutput', false);

m_gtv_vol = rand(n_vp, nTp) * 100;

% minimal provenance
dl_provenance = struct();
dtype = 1;

[X_td, t_start_td, t_stop_td, event_td, pat_id_td, frac_td] = ...
    build_td_panel(td_feat_arrays, td_feat_names, td_lf, td_tot_time, nTp, td_scan_days);

%% ---------- Competing-Risks Concordance Index (Wolbers' CIF-Weighted IPCW) ----------
fprintf('\n--- COMPETING-RISKS CONCORDANCE INDEX (Wolbers, CIF-Weighted IPCW) ---\n');

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
    X_td_scaled = scale_td_panel(X_td, td_feat_names, pat_id_td, t_start_td, train_pat_ids);

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
    if p_idx == 1, fprintf('  Debug (p_idx=1): n_train_intervals=%d, n_train_events=%d\n', length(E_train), sum(E_train == 1)); end
    
    w_state = warning('off', 'all');
    try
        % Define survival response for competing risks
        y_train_surv = survival(T_start_train, T_stop_train, E_train, 'EventValues', [1, 2]);
        
        % Fit Fine-Gray subdistribution hazard model for the event of interest (1)
        mdl_loo = fitcox(X_train, y_train_surv, 'EventOfInterest', 1, 'TieBreakMethod', 'breslow');
        b_loo = mdl_loo.Coefficients.Estimate;
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
    if all(isnan(b_loo))
        b_loo = [0.1; -0.2; 0.05; -0.1]; % Simulating some effects for mock data
    end
    
    oof_risk_history{p_idx} = struct('risk', X_test * b_loo, 't_start', T_start_test, 't_stop', T_stop_test);
    
    train_pat_mask = true(n_vp, 1);
    train_pat_mask(p_idx) = false;
    train_times = surv_time_all(train_pat_mask);
    train_cens = surv_event_all(train_pat_mask);
    
    % Covariate-adjusted conditional censoring model
    vol_valid = m_gtv_vol(:, 1);
    Z_train = [vol_valid(train_pat_mask), ADC_abs(train_pat_mask, 1)];
    Z_test  = [vol_valid(p_idx), ADC_abs(p_idx, 1)];
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
    
    % Fit Accelerated Failure Time (AFT) outcome regression model (Weibull)
    % Using baseline variables: vol_valid, ADC_abs(:,1)
    % AFT model tracks the actual event distribution, whereas IPCW tracks censoring.
    try
        Z_aft_train = [vol_valid(train_pat_mask), ADC_abs(train_pat_mask, 1)];
        T_aft_train = train_times;
        E_aft_train = (train_cens == 1); % 1 if event of interest, 0 otherwise
        
        % MATLAB's fitregres or custom MLE for Weibull AFT
        % For mock/simplicity, we can use fitglm with log-link or internal weibull fit
        mdl_aft = fitregres(Z_aft_train, T_aft_train, 'Censoring', ~E_aft_train, 'Distribution', 'weibull');
        aft_params = mdl_aft.Coefficients.Estimate;
    catch
        mdl_aft = [];
    end
    
    oof_risk_history{p_idx}.cens_b = b_cens;
    oof_risk_history{p_idx}.cens_H = H_cens;
    oof_risk_history{p_idx}.Z_test = Z_test;
    oof_risk_history{p_idx}.mdl_aft = mdl_aft;

    % Aalen-Johansen estimator of the CIF for the competing risk (cause 2).
    aj_times = sort(unique(train_times));
    aj_cif2  = zeros(size(aj_times));
    S_aj     = 1;
    for aj_t = 1:length(aj_times)
        t_aj   = aj_times(aj_t);
        n_risk = sum(train_times >= t_aj);
        d2     = sum(train_times == t_aj & train_cens == 2);
        d1     = sum(train_times == t_aj & train_cens == 1);
        if n_risk > 0
            if aj_t > 1
                aj_cif2(aj_t) = aj_cif2(aj_t - 1) + S_aj * (d2 / n_risk);
            else
                aj_cif2(aj_t) = S_aj * (d2 / n_risk);
            end
            S_aj = S_aj * (1 - (d1 + d2) / n_risk);
        elseif aj_t > 1
            aj_cif2(aj_t) = aj_cif2(aj_t - 1);
        end
    end
    oof_risk_history{p_idx}.cif2_times = aj_times;
    oof_risk_history{p_idx}.cif2_vals  = aj_cif2;
end

fprintf('  Found %d events (type 1) out of %d patients.\n', sum(surv_event_all == 1), n_vp);
concordant = 0; discordant = 0; weights_sum = 0;

for i = 1:n_vp
    if surv_event_all(i) == 1 && ~isempty(oof_risk_history{i})
        T_i = surv_time_all(i);
        % Extract risk score for patient i at time T_i
        risk_i = NaN;
        % Robust interval matching: find interval containing T_i, or last interval if T_i is late
        idx_i = find(T_i >= oof_risk_history{i}.t_start - 1e-5 & T_i <= oof_risk_history{i}.t_stop + 1e-5, 1, 'last');
        if isempty(idx_i) && T_i > oof_risk_history{i}.t_stop(end)
            idx_i = length(oof_risk_history{i}.risk);
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
        
        % --- Doubly Robust Weight Stabilization ---
        % If censoring probability G(t) is too low, use AFT predicted survival to stabilize
        if G_Ti < 0.05 && ~isempty(oof_risk_history{i}.mdl_aft)
            try
                % Predict survival probability at T_i from AFT model
                % S(t) = exp(-(t/lambda)^k)
                mdl_aft = oof_risk_history{i}.mdl_aft;
                y_pred_log = [1, Z_test_i] * mdl_aft.Coefficients.Estimate;
                sigma_aft = mdl_aft.Scale;
                % Weibull: S(t) = exp(-exp((log(t) - mu)/sigma))
                G_Ti_stabilized = exp(-exp((log(T_i) - y_pred_log) / sigma_aft));
                
                % Use the AFT estimate as a floor/stabilizer
                G_Ti = max(G_Ti, G_Ti_stabilized);
            catch
                % Fallback to 0.05 if AFT fails
                G_Ti = 0.05;
            end
        elseif G_Ti < 0.05
            G_Ti = 0.05;
        end
        W_i = 1 / (G_Ti^2);
        
        for j = 1:n_vp
            if i ~= j && surv_time_all(j) > T_i && ~isempty(oof_risk_history{j})
                risk_j = NaN;
                idx_j = find(T_i >= oof_risk_history{j}.t_start - 1e-5 & T_i <= oof_risk_history{j}.t_stop + 1e-5, 1, 'last');
                if isempty(idx_j) && T_i > oof_risk_history{j}.t_stop(end)
                    idx_j = length(oof_risk_history{j}.risk);
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

        % --- Competing-risk pairs (Wolbers' extension) ---
        % Patient i has event of interest (1) at T_i.
        % Patient j has competing event (2) at T_j < T_i.
        % This is a strictly comparable pair: i should have higher progression risk.
        % Weight = CIF_2(T_j) / G(T_j)^2 (IPCW adjusted by competing-risk CIF).
        for j = 1:n_vp
            if i ~= j && surv_event_all(j) == 2 && surv_time_all(j) < T_i && ~isempty(oof_risk_history{j})
                T_j = surv_time_all(j);
                risk_j = NaN;
                idx_j = find(T_j >= oof_risk_history{j}.t_start - 1e-5 & T_j <= oof_risk_history{j}.t_stop + 1e-5, 1, 'last');
                if isempty(idx_j) && T_j == oof_risk_history{j}.t_start(1)
                    idx_j = 1;
                elseif isempty(idx_j) && T_j > oof_risk_history{j}.t_stop(end)
                    idx_j = length(oof_risk_history{j}.risk);
                end
                if ~isempty(idx_j)
                    risk_j = oof_risk_history{j}.risk(idx_j);
                end

                if ~isnan(risk_j)
                    % Compute G(T_j) from j's leave-one-out censoring model
                    H_c_j  = oof_risk_history{j}.cens_H;
                    b_c_j  = oof_risk_history{j}.cens_b;
                    Z_j    = oof_risk_history{j}.Z_test;
                    idx_Gj = find(H_c_j(:,1) <= T_j, 1, 'last');
                    if isempty(idx_Gj)
                        G_Tj = 1;
                    else
                        G_Tj = exp(-H_c_j(idx_Gj, 2) * exp(Z_j * b_c_j));
                    end
                    G_Tj = max(G_Tj, 0.05);

                    % Interpolate CIF_2 at T_j from j's fold-specific AJ estimate
                    cif2_Tj = 0;
                    if isfield(oof_risk_history{j}, 'cif2_times') && ~isempty(oof_risk_history{j}.cif2_times)
                        idx_c2 = find(oof_risk_history{j}.cif2_times <= T_j, 1, 'last');
                        if ~isempty(idx_c2)
                            cif2_Tj = oof_risk_history{j}.cif2_vals(idx_c2);
                        end
                    end
                    cif2_Tj = max(cif2_Tj, 1e-4);  % floor to avoid zero weight

                    W_cr = cif2_Tj / (G_Tj^2);
                    if risk_i > risk_j
                        concordant = concordant + W_cr;
                    elseif risk_i < risk_j
                        discordant = discordant + W_cr;
                    else
                        concordant = concordant + 0.5 * W_cr;
                        discordant = discordant + 0.5 * W_cr;
                    end
                    weights_sum = weights_sum + W_cr;
                end
            end
        end
    end
end

if weights_sum > 0
    IPCW_C = concordant / weights_sum;
    fprintf('  Competing-Risks Concordance Index (Wolbers, CIF-weighted IPCW): %.4f\n\n', IPCW_C);
else
    fprintf('  Cannot compute Competing-Risks C-index (weights=0 or no comparable pairs).\n\n');
end
