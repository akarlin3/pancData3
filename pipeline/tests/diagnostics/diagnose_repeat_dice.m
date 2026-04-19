function diagnose_repeat_dice()
% DIAGNOSE_REPEAT_DICE  Investigate why dice_rpt_* arrays are all NaN.
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

    cfg = parse_config();
    dataloc = cfg.dataloc;
    sm_path = fullfile(dataloc, 'summary_metrics_Standard.mat');
    dv_path = fullfile(dataloc, 'dwi_vectors_Standard.mat');

    fprintf('\n=== Repeat Dice Diagnostic ===\n');
    fprintf('summary_metrics: %s\n', sm_path);
    fprintf('dwi_vectors:     %s\n\n', dv_path);

    if ~exist(sm_path, 'file'); error('summary_metrics not found at %s', sm_path); end
    if ~exist(dv_path, 'file'); error('dwi_vectors not found at %s', dv_path); end

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
    has_path     = zeros(nPat, 1);   % #Fx1 repeats with a non-empty path
    file_exists  = zeros(nPat, 1);   % #Fx1 repeats where the file exists on disk
    has_stvol3d  = zeros(nPat, 1);   % #Fx1 repeats whose .mat has 'Stvol3d'
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

        for rpi = 1:nRpt
            p = gtv_locs{j, 1, rpi};
            if ~isempty(p)
                has_path(j) = has_path(j) + 1;
                if exist(p, 'file')
                    file_exists(j) = file_exists(j) + 1;
                    try
                        vars = whos('-file', p);
                        names = {vars.name};
                        if any(strcmp(names, 'Stvol3d'))
                            has_stvol3d(j) = has_stvol3d(j) + 1;
                        else
                            % Collect other variable names in this mask file.
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
                        % leave counters alone
                    end
                end
            end
        end
    end

    % Per-patient Fx1 diagnostic (first 15 patients).
    fprintf('%-20s %-10s %-10s %-12s %-12s\n', 'Patient', 'rp_count', 'paths', 'files exist', 'Stvol3d');
    fprintf('%s\n', repmat('-', 1, 68));
    for j = 1:min(nPat, 15)
        fprintf('%-20s %-10d %-10d %-12d %-12d\n', ...
            sm.id_list{j}, rp_count_vec(j), has_path(j), file_exists(j), has_stvol3d(j));
    end
    if nPat > 15
        fprintf('... (%d more patients)\n', nPat - 15);
    end

    % Cohort summary.
    fprintf('\n=== Summary ===\n');
    fprintf('Patients with rp_count >= 2 (required for Dice): %d / %d\n', ...
        sum(rp_count_vec >= 2), nPat);
    fprintf('Patients with >=2 Fx1 GTV paths:                  %d / %d\n', ...
        sum(has_path >= 2), nPat);
    fprintf('Patients with >=2 Fx1 GTV files that exist:       %d / %d\n', ...
        sum(file_exists >= 2), nPat);
    fprintf('Patients with >=2 Fx1 GTV files containing Stvol3d: %d / %d\n', ...
        sum(has_stvol3d >= 2), nPat);

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
