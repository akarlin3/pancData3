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

td_tot_time = randi([10, 200], n_vp, 1);
td_lf = randi([0, 1], n_vp, 1);
m_id_list = arrayfun(@(x) sprintf('Pt%d', x), 1:n_vp, 'UniformOutput', false);

% minimal provenance
dl_provenance = struct();
dtype = 1;

fprintf('\n--- MOCK LANDMARK ANALYSIS (Out-of-Fold Validation) ---\n');

surv_time_all = td_tot_time; 
surv_event_all = td_lf;
valid_idx_pts = find(valid_pts);
landmarks = [30, 90, 180];
n_lm = length(landmarks);

for lm_idx = 1:n_lm
    L = landmarks(lm_idx);
    fprintf('  Evaluating Landmark L = %d days...\n', L);
    at_risk = (surv_time_all > L);
    if sum(at_risk) < 10
        continue;
    end
    
    X_L = nan(n_vp, td_n_feat);
    for i = 1:n_vp
        if at_risk(i)
            for f_idx = 1:td_n_feat
                feat_mat = td_feat_arrays{f_idx};
                for t = length(td_scan_days):-1:1
                    if td_scan_days(t) <= L && ~isnan(feat_mat(i, t))
                        X_L(i, f_idx) = feat_mat(i, t);
                        break;
                    end
                end
            end
        end
    end
    
    complete_data = at_risk & ~any(isnan(X_L), 2);
    surv_time_L = surv_time_all - L; 
    surv_event_L = surv_event_all;
    
    oof_risk = nan(n_vp, 1);
    patients_to_test = find(complete_data);
    
    for p_idx = 1:length(patients_to_test)
        loo_i = patients_to_test(p_idx);
        
        train_mask = complete_data;
        train_mask(loo_i) = false;
        
        X_train = X_L(train_mask, :);
        T_train = surv_time_L(train_mask);
        E_train = surv_event_L(train_mask);
        X_test = X_L(loo_i, :);
        
        w_state = warning('off', 'all');
        try
            b_loo = coxphfit(X_train, T_train, 'Censoring', ~E_train, 'Options', statset('MaxIter', 100));
            oof_risk(loo_i) = X_test * b_loo;
        catch
            try
                mdl_firth = fitglm(X_train, E_train, 'Distribution', 'binomial', 'LikelihoodPenalty', 'jeffreys-prior');
                b_loo = mdl_firth.Coefficients.Estimate(2:end);
                oof_risk(loo_i) = X_test * b_loo;
            catch
                oof_risk(loo_i) = NaN;
            end
        end
        warning(w_state);
    end
    
    valid_oof = ~isnan(oof_risk) & complete_data;
    if sum(valid_oof) < 5 || sum(surv_event_L(valid_oof)) < 2
        continue;
    end
    
    T_valid = surv_time_L(valid_oof);
    E_valid = surv_event_L(valid_oof);
    R_valid = oof_risk(valid_oof);
    
    [F_cens, x_cens] = ecdf(surv_time_L(complete_data), 'Censoring', surv_event_L(complete_data));
    G_cens = 1 - F_cens;
    
    n_v = length(T_valid);
    concordant = 0; discordant = 0; weights = 0;
    
    for i = 1:n_v
        if E_valid(i) == 1
            idx_G = find(x_cens <= T_valid(i), 1, 'last');
            if isempty(idx_G) || G_cens(idx_G) == 0
                G_Ti = 1;
            else
                G_Ti = G_cens(idx_G);
            end
            
            W_i = 1 / (G_Ti^2);
            
            for j = 1:n_v
                if T_valid(j) > T_valid(i)
                    if R_valid(i) > R_valid(j) 
                        concordant = concordant + W_i;
                    elseif R_valid(i) < R_valid(j)
                        discordant = discordant + W_i;
                    else
                        concordant = concordant + 0.5 * W_i;
                        discordant = discordant + 0.5 * W_i;
                    end
                    weights = weights + W_i;
                end
            end
        end
    end
    
    if weights > 0
        Uno_C = concordant / weights;
        fprintf('    Uno''s C-statistic (L=%d): %.4f\n\n', L, Uno_C);
    else
        fprintf('    Cannot compute Uno''s C (weights=0).\n\n');
    end
end
disp('DONE MOCK TEST');
