function [survival_data, config] = prepare_survival_data_dispatch(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ...
    m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label, ...
    m_gtv_vol, actual_scan_days, config_struct_in, output_folder)
% PREPARE_SURVIVAL_DATA_DISPATCH Orchestrate data preparation by delegating
% to independently testable utility functions:
%   build_survival_features  — feature array assembly
%   prepare_outcome_data     — censoring/competing risk logic
%   select_landmark_day      — landmark day selection

% --- 1. Scan days ---
if ~isempty(actual_scan_days)
    td_scan_days = actual_scan_days;
    fprintf('  Using provided scan days [%s].\n', num2str(td_scan_days));
else
    td_scan_days = [0, 5, 10, 15, 20, 90];
    fprintf('  ⚠️  CAUTION: Using default scan days [%s].\n', num2str(td_scan_days));
    fprintf('      Pass actual DICOM-derived scan days to avoid immortal time bias.\n');
end

% --- 2. Build feature arrays (delegated to build_survival_features) ---
[td_feat_arrays, td_feat_names] = build_survival_features(valid_pts, ...
    ADC_abs, D_abs, f_abs, Dstar_abs, m_gtv_vol, nTp);

% --- 3. Prepare outcome data with censoring logic (delegated to prepare_outcome_data) ---
[td_lf, td_tot_time] = prepare_outcome_data(valid_pts, m_lf, m_total_time, m_total_follow_up_time);

% --- 4. Parse half-life config ---
if isfield(config_struct_in, 'td_halflife_days') && ~isempty(config_struct_in.td_halflife_days)
    td_halflife_days = config_struct_in.td_halflife_days;
    fprintf('  Using configured td_halflife_days = %.1f (from config).\n', td_halflife_days);
else
    td_halflife_days = 18;
    fprintf('  Using default td_halflife_days = %.1f.\n', td_halflife_days);
    fprintf('      Set config.td_halflife_days to customize the decay half-life for time-dependent feature weighting.\n');
end

% --- 5. Build time-dependent panel ---
[X_td_def, t_start_td_def, t_stop_td_def, event_td_def, pat_id_td_def, frac_td_def] = ...
    build_td_panel(td_feat_arrays, td_feat_names, td_lf, td_tot_time, nTp, td_scan_days, td_halflife_days);

% --- 6. Check initial data sufficiency ---
td_n_feat = numel(td_feat_arrays);
td_ok = (sum(event_td_def == 1) >= 3) && (size(X_td_def, 1) > td_n_feat + 1);

% --- 7. Select landmark day (delegated to select_landmark_day) ---
landmark_day = select_landmark_day(td_scan_days, config_struct_in);

% --- 8. Apply landmark filtering ---
lm_keep = (t_start_td_def >= landmark_day);

if any(lm_keep)
    X_td = X_td_def(lm_keep, :);
    t_start_td = t_start_td_def(lm_keep);
    t_stop_td = t_stop_td_def(lm_keep);
    event_td = event_td_def(lm_keep);
    pat_id_td = pat_id_td_def(lm_keep);
    frac_td = frac_td_def(lm_keep);

    n_events_post_lm = sum(event_td == 1);
    td_ok = td_ok && (n_events_post_lm >= 3) && (size(X_td, 1) > td_n_feat + 1);

    fprintf('  Landmark at day %d: %d events, %d patients\n', ...
        landmark_day, n_events_post_lm, numel(unique(pat_id_td)));
else
    td_ok = false;
end

% --- 9. Package results ---
survival_data = struct();
survival_data.has_sufficient_data = td_ok;
if td_ok
    survival_data.X = X_td;
    survival_data.t_start = t_start_td;
    survival_data.t_stop = t_stop_td;
    survival_data.event = event_td;
    survival_data.pat_id = pat_id_td;
    survival_data.frac = frac_td;
    survival_data.feat_names = td_feat_names;
    survival_data.n_feat = td_n_feat;
    survival_data.landmark_day = landmark_day;
    survival_data.scan_days = td_scan_days;
    survival_data.feat_arrays = td_feat_arrays;
    survival_data.lf_status = td_lf;
    survival_data.total_time = td_tot_time;
end

config = struct();
config.half_life_grid = [3, 6, 12, 18, 24];
config.n_bootstrap = 1000;
config.output_folder = output_folder;
config.dtype_label = dtype_label;
config.td_halflife_days = td_halflife_days;

% Propagate user config flags for optional analyses
if isfield(config_struct_in, 'compute_fine_gray')
    config.compute_fine_gray = config_struct_in.compute_fine_gray;
else
    config.compute_fine_gray = true;
end
if isfield(config_struct_in, 'fit_time_varying_cox')
    config.fit_time_varying_cox = config_struct_in.fit_time_varying_cox;
else
    config.fit_time_varying_cox = true;
end
end
