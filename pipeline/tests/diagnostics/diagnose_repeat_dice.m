function diagnose_repeat_dice(output_folder)
% DIAGNOSE_REPEAT_DICE  Investigate why dice_rpt_* arrays are all NaN.
%
% Usage:
%   diagnose_repeat_dice();                              % auto-detect most recent saved_files_*
%   diagnose_repeat_dice('saved_files_20260418_194711'); % explicit output folder
%
% Checks, for each patient in summary_metrics_Standard.mat:
%   1. Whether gtv_locations{j, 1, rpi} points to real files for Fx1 repeats.
%   2. Whether each Fx1 repeat .mat file contains a variable named 'Stvol3d'.
%      (safe_load_mask returns [] if the variable is named something else.)
%   3. How many patients have rp_count >= 2 (required for spatial repeatability).
%
% Usage (from MATLAB):
%     cd C:\Projects\pancData3
%     addpath('pipeline/utils');
%     diagnose_repeat_dice();
%
% Output: a printed table showing per-patient Fx1 repeat status plus a summary.

    if nargin < 1, output_folder = ''; end

    % Make sure utils (normalize_path_preserving_roots, safe_load_mask,
    % compute_spatial_repeatability) are on the path.
    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', '..', 'utils'));

    cfg = parse_config('config.json');
    dataloc = cfg.dataloc;

    fprintf('\n=== Repeat Dice Diagnostic ===\n');
    fprintf('dataloc: %s\n', dataloc);

    % Auto-detect the most recent saved_files_* folder when none was provided.
    % compute_summary_metrics saves the checkpoint under dataloc, but the
    % canonical output is under <project>/saved_files_<timestamp>/.
    if isempty(output_folder)
        % Search both current working directory and dataloc's parent.
        search_roots = {pwd};
        if ~isempty(dataloc); search_roots{end+1} = dataloc; end
        candidates = [];
        for r = 1:length(search_roots)
            d = dir(fullfile(search_roots{r}, 'saved_files_*'));
            d = d([d.isdir]);
            for i = 1:length(d)
                d(i).fullpath = fullfile(search_roots{r}, d(i).name);
            end
            if ~isempty(d); candidates = [candidates; d]; end %#ok<AGROW>
        end
        if ~isempty(candidates)
            [~, idx] = max([candidates.datenum]);
            output_folder = candidates(idx).fullpath;
            fprintf('output_folder (auto): %s\n', output_folder);
        end
    else
        fprintf('output_folder: %s\n', output_folder);
    end

    % Locate summary_metrics and dwi_vectors, trying several candidate names.
    sm_candidates = { ...
        fullfile(dataloc, 'summary_metrics_Standard.mat'), ...
        fullfile(dataloc, 'summary_metrics.mat') };
    dv_candidates = { ...
        fullfile(dataloc, 'dwi_vectors_Standard.mat'), ...
        fullfile(dataloc, 'dwi_vectors.mat') };
    if ~isempty(output_folder)
        sm_candidates = [ ...
            {fullfile(output_folder, 'summary_metrics_Standard.mat'), ...
             fullfile(output_folder, 'Standard', 'summary_metrics_Standard.mat')}, ...
            sm_candidates];
        dv_candidates = [ ...
            {fullfile(output_folder, 'dwi_vectors_Standard.mat'), ...
             fullfile(output_folder, 'Standard', 'dwi_vectors_Standard.mat')}, ...
            dv_candidates];
    end

    sm_path = '';
    for i = 1:length(sm_candidates)
        if exist(sm_candidates{i}, 'file'); sm_path = sm_candidates{i}; break; end
    end
    dv_path = '';
    for i = 1:length(dv_candidates)
        if exist(dv_candidates{i}, 'file'); dv_path = dv_candidates{i}; break; end
    end

    if isempty(sm_path)
        fprintf('\n❌ summary_metrics not found at any of:\n');
        for i = 1:length(sm_candidates); fprintf('     %s\n', sm_candidates{i}); end

        fprintf('\nFiles in dataloc matching summary_metrics*.mat:\n');
        d1 = dir(fullfile(dataloc, 'summary_metrics*.mat'));
        if isempty(d1); fprintf('  (none)\n'); else
            for i = 1:length(d1); fprintf('  %s  (%s)\n', d1(i).name, d1(i).date); end
        end

        if ~isempty(output_folder)
            fprintf('\nFiles in output_folder matching summary_metrics*.mat (recursive):\n');
            d2 = dir(fullfile(output_folder, '**', 'summary_metrics*.mat'));
            if isempty(d2); fprintf('  (none)\n'); else
                for i = 1:length(d2)
                    fprintf('  %s  (%s)\n', fullfile(d2(i).folder, d2(i).name), d2(i).date);
                end
            end
        end

        error(['Please pass the correct output folder: ' ...
               'diagnose_repeat_dice(''C:\\Projects\\pancData3\\saved_files_YYYYMMDD_HHMMSS'')']);
    end
    if isempty(dv_path)
        fprintf('\n❌ dwi_vectors not found at any of:\n');
        for i = 1:length(dv_candidates); fprintf('     %s\n', dv_candidates{i}); end
        error('dwi_vectors missing.');
    end

    fprintf('summary_metrics: %s\n', sm_path);
    fprintf('dwi_vectors:     %s\n\n', dv_path);

    sm = load(sm_path, 'summary_metrics');
    sm = sm.summary_metrics;
    dv = load(dv_path);
    dwi_vectors_gtvp = dv.data_vectors_gtvp;

    if isfield(sm, 'gtv_locations')
        gtv_locs = sm.gtv_locations;
    else
        error('summary_metrics has no gtv_locations field');
    end

    nPat = length(sm.id_list);
    nRpt = size(dwi_vectors_gtvp, 3);
    fprintf('Cohort: %d patients, %d repeat slots, gtv_locations size = [%s]\n\n', ...
        nPat, nRpt, num2str(size(gtv_locs)));

    % Tally
    has_path        = zeros(nPat, 1);   % #Fx1 repeats with non-empty raw path
    file_exists_raw = zeros(nPat, 1);   % #Fx1 repeats where raw path exists
    file_exists_norm= zeros(nPat, 1);   % #Fx1 repeats where normalized path exists
    has_stvol3d     = zeros(nPat, 1);   % #Fx1 repeats whose .mat has 'Stvol3d'
    has_path_effective = zeros(nPat, 1);% post-shared-Fx1 fallback path count
    voxel_match     = zeros(nPat, 1);   % #Fx1 repeats where mask n_gtv==numel(adc_vec)
    non_stvol_var_list = containers.Map('KeyType','char','ValueType','double');
    rp_count_vec = zeros(nPat, 1);

    for j = 1:nPat
        % Mirror the rp_count logic in compute_summary_metrics (Standard DWI).
        rp_count = 0;
        for rpi = 1:nRpt
            adc_vec = dwi_vectors_gtvp(j, 1, rpi).adc_vector;
            if ~isempty(adc_vec); rp_count = rp_count + 1; end
        end
        rp_count_vec(j) = rp_count;

        % --- Shared Fx1 fallback (mirrors compute_spatial_repeatability) ---
        shared_fx1 = '';
        for ri = 1:nRpt
            c = gtv_locs{j, 1, ri};
            if ~isempty(c); shared_fx1 = c; break; end
        end

        for rpi = 1:nRpt
            p_raw = gtv_locs{j, 1, rpi};
            if ~isempty(p_raw); has_path(j) = has_path(j) + 1; end

            % Effective path = raw or shared-Fx1 fallback.
            p_eff = p_raw;
            if isempty(p_eff); p_eff = shared_fx1; end
            if isempty(p_eff); continue; end
            has_path_effective(j) = has_path_effective(j) + 1;

            % Raw-path exist check (pre-normalization).
            if exist(p_eff, 'file') == 2
                file_exists_raw(j) = file_exists_raw(j) + 1;
            end

            % Normalized-path exist check (what compute_spatial_repeatability
            % actually uses).
            p_norm = normalize_path_preserving_roots(p_eff);
            if exist(p_norm, 'file') == 2
                file_exists_norm(j) = file_exists_norm(j) + 1;
                try
                    vars = whos('-file', p_norm);
                    names = {vars.name};
                    if any(strcmp(names, 'Stvol3d'))
                        has_stvol3d(j) = has_stvol3d(j) + 1;
                        % Voxel-count match check: load mask and compare to
                        % adc_vector length for this repeat.
                        try
                            m = safe_load_mask(p_norm, 'Stvol3d');
                            adc_vec = dwi_vectors_gtvp(j, 1, rpi).adc_vector;
                            n_gtv = sum(m(:) == 1);
                            if ~isempty(adc_vec) && numel(adc_vec) == n_gtv
                                voxel_match(j) = voxel_match(j) + 1;
                            end
                        catch
                        end
                    else
                        for nn = 1:numel(names)
                            key = names{nn};
                            if isKey(non_stvol_var_list, key)
                                non_stvol_var_list(key) = non_stvol_var_list(key) + 1;
                            else
                                non_stvol_var_list(key) = 1;
                            end
                        end
                    end
                catch
                end
            end
        end
    end

    % Per-patient Fx1 diagnostic (first 15 patients).
    fprintf(['%-20s %-8s %-8s %-11s %-12s %-12s %-8s %-10s\n'], ...
        'Patient','rp_cnt','raw_p','eff_path','exist_raw','exist_norm','Stvol3d','vox_match');
    fprintf('%s\n', repmat('-', 1, 100));
    for j = 1:min(nPat, 15)
        fprintf('%-20s %-8d %-8d %-11d %-12d %-12d %-8d %-10d\n', ...
            sm.id_list{j}, rp_count_vec(j), has_path(j), has_path_effective(j), ...
            file_exists_raw(j), file_exists_norm(j), has_stvol3d(j), voxel_match(j));
    end
    if nPat > 15
        fprintf('... (%d more patients)\n', nPat - 15);
    end

    % Cohort summary.
    fprintf('\n=== Summary ===\n');
    fprintf('Patients with rp_count >= 2 (required for Dice):      %d / %d\n', ...
        sum(rp_count_vec >= 2), nPat);
    fprintf('Patients with >=2 raw Fx1 GTV paths:                  %d / %d\n', ...
        sum(has_path >= 2), nPat);
    fprintf('Patients with >=2 paths post-shared-Fx1 fallback:     %d / %d\n', ...
        sum(has_path_effective >= 2), nPat);
    fprintf('Patients with >=2 files that exist (raw path):        %d / %d\n', ...
        sum(file_exists_raw >= 2), nPat);
    fprintf('Patients with >=2 files that exist (normalized path): %d / %d\n', ...
        sum(file_exists_norm >= 2), nPat);
    fprintf('Patients with >=2 Fx1 GTV files containing Stvol3d:   %d / %d\n', ...
        sum(has_stvol3d >= 2), nPat);
    fprintf('Patients with >=2 voxel-count matches:                %d / %d\n', ...
        sum(voxel_match >= 2), nPat);

    % Actually call compute_spatial_repeatability for first 3 rp_count>=2
    % patients and show the result (this will also write a per-patient
    % trace to <tempdir>/repeat_dice_trace.txt).
    fprintf('\n=== Live compute_spatial_repeatability for first 3 eligible patients ===\n');
    trace_file = fullfile(tempdir, 'repeat_dice_trace.txt');
    if exist(trace_file, 'file') == 2
        try; delete(trace_file); catch; end
    end
    try; clear compute_spatial_repeatability; catch; end
    eligible = find(rp_count_vec >= 2);
    nTry = min(3, numel(eligible));
    adc_thresh  = cfg.adc_thresh;
    d_thresh    = cfg.d_thresh;
    f_thresh    = cfg.f_thresh;
    dstar_thresh= cfg.dstar_thresh;
    morph_se    = strel('disk', 1);
    morph_min_cc= 3;
    for kk = 1:nTry
        j = eligible(kk);
        last_m = ''; last_mask = [];
        try
            [dice_a, ~, ~, dice_d, ~, ~, dice_f, ~, ~, dice_ds, ~, ~] = ...
                compute_spatial_repeatability(dwi_vectors_gtvp, j, 1, gtv_locs, ...
                    adc_thresh, d_thresh, f_thresh, dstar_thresh, ...
                    morph_se, morph_min_cc, last_m, last_mask);
            fprintf('  j=%d (%s): dice_adc=%g dice_d=%g dice_f=%g dice_dstar=%g\n', ...
                j, sm.id_list{j}, dice_a, dice_d, dice_f, dice_ds);
        catch err
            fprintf('  j=%d (%s): ERROR %s\n', j, sm.id_list{j}, err.message);
        end
    end
    if exist(trace_file, 'file') == 2
        fprintf('\nPer-patient trace written to: %s\n', trace_file);
    end

    if non_stvol_var_list.Count > 0
        fprintf('\n=== Non-Stvol3d variable names found in Fx1 GTV files ===\n');
        keys_all = non_stvol_var_list.keys;
        for i = 1:length(keys_all)
            fprintf('  %-40s  (%d files)\n', keys_all{i}, non_stvol_var_list(keys_all{i}));
        end
        fprintf(['\nNOTE: safe_load_mask defaults to varname="Stvol3d".  If most GTV files ' ...
                 'use a different variable name, compute_spatial_repeatability will receive ' ...
                 'empty masks and return NaN.\n']);
    else
        fprintf('\nAll Fx1 GTV files that exist contain Stvol3d.\n');
    end

    % Sanity check: print actual dice_rpt_adc column for dwi_type=1.
    if isfield(sm, 'dice_rpt_adc')
        col = sm.dice_rpt_adc(:, 1);
        fprintf('\nFinite dice_rpt_adc (Standard) entries: %d / %d\n', ...
            sum(~isnan(col)), numel(col));
        if any(~isnan(col))
            finite_vals = col(~isnan(col));
            fprintf('Mean = %.3f, range = [%.3f, %.3f]\n', ...
                mean(finite_vals), min(finite_vals), max(finite_vals));
        end
    end
end
