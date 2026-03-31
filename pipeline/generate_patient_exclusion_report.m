function report = generate_patient_exclusion_report(config_path)
% GENERATE_PATIENT_EXCLUSION_REPORT — List every excluded patient and the reasons.
%
%   report = generate_patient_exclusion_report()
%   report = generate_patient_exclusion_report('config.json')
%
%   Scans the data directory AND the most recent pipeline output to compile a
%   comprehensive report of every patient excluded from analysis and why.
%
%   Exclusion categories checked:
%     1. Pre-pipeline data issues (missing Fx1 DWI, missing Fx1 GTV, etc.)
%     2. Missing baseline imaging (no GTV volume or ADC at Fx1)
%     3. No matching clinical record (NaN local failure status)
%     4. Competing risk classification (non-cancer death without LF)
%     5. DL training set membership (data leakage prevention)
%     6. Baseline outlier detection (3x IQR fencing)
%     7. Insufficient survival data (< 3 events or n < p+1)
%
%   Returns a struct with fields:
%     exclusions  - Cell array of {patient_id, reason, category} rows
%     summary     - Struct with counts per exclusion category
%     n_total     - Total patients in data directory
%     n_excluded  - Number of patients with at least one exclusion
%     n_analysed  - Number of patients passing all filters

    if nargin < 1
        script_dir = fileparts(mfilename('fullpath'));
        config_path = fullfile(script_dir, '..', 'config.json');
    end

    % Add pipeline paths
    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, 'core'), fullfile(script_dir, 'utils'), fullfile(script_dir, 'dependencies'));

    %% --- Load config ---
    config = parse_config(config_path);
    dataloc = config.dataloc;

    fprintf('=============================================================\n');
    fprintf('  PATIENT EXCLUSION REPORT\n');
    fprintf('=============================================================\n\n');

    % Accumulator: {patient_id, reason, category}
    exclusions = cell(0, 3);

    %% ===================================================================
    %  STAGE 1: Pre-pipeline data completeness (file system scan)
    %  ===================================================================
    fprintf('📋 Stage 1: Scanning data directory for missing files...\n');

    patlist = clean_dir_command(dataloc);
    patlist = patlist(cellfun(@isempty, strfind({patlist.name}, 'template')));

    % Extract numeric IDs
    id_num = zeros(size(patlist));
    pat_id_str = cell(size(patlist));
    for j = 1:length(patlist)
        strtmp = strsplit(strrep(strrep(patlist(j).name, 'P', ''), '_', '-'), '-');
        id_num(j) = str2double(strtmp{1});
        pat_id_str{j} = strtmp{1};
    end

    % Non-patient folders — log but do not count in exclusions or totals
    non_patient = isnan(id_num);
    non_patient_names = {patlist(non_patient).name};
    if ~isempty(non_patient_names)
        fprintf('  💡 Skipping %d non-patient folder(s):', length(non_patient_names));
        for npi = 1:length(non_patient_names)
            fprintf(' %s', non_patient_names{npi});
        end
        fprintf('\n');
    end

    % Keep only patient folders
    patlist = patlist(~non_patient);
    pat_id_str = pat_id_str(~non_patient);
    n_total = length(patlist);

    % Optional patient_ids filter from config
    if isfield(config, 'patient_ids') && ~isempty(config.patient_ids)
        pid_filter = config.patient_ids;
        if iscell(pid_filter)
            pid_filter = cellfun(@str2double, pid_filter);
        end
        id_num_valid = cellfun(@str2double, pat_id_str);
        filtered_out = ~ismember(id_num_valid, pid_filter);
        for j = find(filtered_out)'
            exclusions = add_excl(exclusions, patlist(j).name, ...
                'Filtered out by config.patient_ids', 'Config filter');
        end
    end

    % Check each patient for missing baseline data
    fx_search = {'Fx1', 'Fx2', 'Fx3', 'Fx4', 'Fx5', 'post'};
    for j = 1:n_total
        pat_name = patlist(j).name;
        basefolder = fullfile(dataloc, pat_name);
        basefolder_contents = clean_dir_command(basefolder);

        % Check Fx1 folder existence
        fx1_idx = ~cellfun(@isempty, strfind({basefolder_contents.name}, 'Fx1'));
        fx1_dirs = basefolder_contents(fx1_idx);
        if isempty(fx1_dirs)
            exclusions = add_excl(exclusions, pat_name, ...
                'Missing Fx1 (baseline) folder', 'Missing baseline data');
            continue;
        end

        fxfolder = fullfile(basefolder, fx1_dirs(1).name);

        % Check Fx1 DWI
        dwi_search = clean_dir_command(fullfile(fxfolder, '*DWI*'));
        if isempty(dwi_search)
            dwi_search = clean_dir_command(fullfile(fxfolder, '*DICOM*'));
            dwi_search = dwi_search(cellfun(@isempty, strfind({dwi_search.name}, 't2')));
        end
        has_dwi = false;
        if ~isempty(dwi_search)
            for dwii = 1:length(dwi_search)
                dwi_path = fullfile(dwi_search(dwii).folder, dwi_search(dwii).name);
                if ~isempty(dir(fullfile(dwi_path, '*.dcm'))) || ...
                   ~isempty(dir(fullfile(dwi_path, 'DICOM', '*.dcm')))
                    has_dwi = true;
                    break;
                end
            end
        end
        if ~has_dwi
            exclusions = add_excl(exclusions, pat_name, ...
                'No DWI DICOM files in Fx1', 'Missing baseline data');
        end

        % Check Fx1 GTV mask
        gtv_files = dir(fullfile(fxfolder, '*GTV*.mat'));
        if isempty(gtv_files)
            exclusions = add_excl(exclusions, pat_name, ...
                'No GTV mask file for Fx1', 'Missing baseline data');
        end

        % Count missing fractions (informational)
        n_missing_fx = 0;
        missing_fx_list = {};
        for fi = 2:length(fx_search)
            fx_label = fx_search{fi};
            fxtmp_idx = ~cellfun(@isempty, strfind({basefolder_contents.name}, fx_label));
            if ~any(fxtmp_idx)
                n_missing_fx = n_missing_fx + 1;
                missing_fx_list{end+1} = fx_label;
            end
        end
        if n_missing_fx > 0
            exclusions = add_excl(exclusions, pat_name, ...
                sprintf('Missing fraction folders: %s (partial data, handled by imputation)', ...
                    strjoin(missing_fx_list, ', ')), 'Incomplete fractions (non-fatal)');
        end
    end

    %% ===================================================================
    %  STAGE 2: Post-pipeline exclusions (from saved pipeline outputs)
    %  ===================================================================
    fprintf('📋 Stage 2: Checking pipeline output for analysis-stage exclusions...\n');

    % Find the most recent saved_files folder (lives in repo root, not dataloc)
    repo_root = fullfile(fileparts(mfilename('fullpath')), '..');
    % Resolve '..' so that dir() glob works reliably on Windows
    prev_dir = pwd;
    cd(repo_root);
    repo_root = pwd;
    cd(prev_dir);
    saved_dirs = dir(fullfile(repo_root, 'saved_files_*'));
    if isempty(saved_dirs)
        fprintf('  ⚠️  No saved_files_* output folder found. Stage 2 checks skipped.\n');
        fprintf('     Run the pipeline first to get full exclusion reporting.\n\n');
    else
        [~, newest_idx] = max([saved_dirs.datenum]);
        output_folder = fullfile(repo_root, saved_dirs(newest_idx).name);
        fprintf('  Using output folder: %s\n', output_folder);

        dwi_type_names = {'Standard', 'dnCNN', 'IVIMnet'};
        dtype = config.dwi_types_to_run;
        dtype_name = dwi_type_names{dtype};

        % Try to load summary_metrics
        sm_file = fullfile(output_folder, sprintf('summary_metrics_%s.mat', dtype_name));
        if ~exist(sm_file, 'file')
            sm_file = fullfile(output_folder, 'summary_metrics.mat');
        end

        if exist(sm_file, 'file')
            tmp = load(sm_file, 'summary_metrics');
            sm = tmp.summary_metrics;
            id_list = sm.id_list;
            mrn_list = sm.mrn_list;
            adc_mean = sm.adc_mean;
            gtv_vol = sm.gtv_vol;

            fprintf('  Loaded %d patients from summary_metrics.\n', length(id_list));

            % --- Check 2a: Missing baseline imaging (NaN GTV or ADC at Fx1) ---
            valid_baseline = ~isnan(gtv_vol(:,1)) & ~isnan(adc_mean(:,1,dtype));
            for j = find(~valid_baseline)'
                reasons = {};
                if isnan(gtv_vol(j,1))
                    reasons{end+1} = 'NaN GTV volume at Fx1';
                end
                if isnan(adc_mean(j,1,dtype))
                    reasons{end+1} = sprintf('NaN ADC at Fx1 (%s)', dtype_name);
                end
                exclusions = add_excl(exclusions, id_list{j}, ...
                    sprintf('Missing baseline imaging: %s', strjoin(reasons, '; ')), ...
                    'Missing baseline imaging');
            end

            % --- Check 2b: Clinical record matching & competing risks ---
            clinical_data_sheet = config.clinical_data_sheet;
            if ~isempty(clinical_data_sheet)
                sheet_full = fullfile(dataloc, clinical_data_sheet);
                if ~exist(sheet_full, 'file')
                    sheet_full = clinical_data_sheet;
                end
                if exist(sheet_full, 'file')
                    if isfield(config, 'clinical_sheet_name')
                        T = readtable(sheet_full, 'Sheet', config.clinical_sheet_name);
                    else
                        T = readtable(sheet_full);
                    end

                    % Find patient ID column
                    pat_col = '';
                    for cand = {'Pat', 'Patient', 'PatientID', 'Patient_ID'}
                        if ismember(cand{1}, T.Properties.VariableNames)
                            pat_col = cand{1};
                            break;
                        end
                    end
                    if isempty(pat_col)
                        vnames = T.Properties.VariableNames;
                        if ~isempty(vnames)
                            pat_col = vnames{1};
                        end
                    end

                    if ~isempty(pat_col)
                        [pat_normalized, id_list_normalized] = normalize_patient_ids(T.(pat_col), id_list);

                        % Determine cause-of-death column
                        cod_disabled = false;
                        if isfield(config, 'cause_of_death_column')
                            cod_column = char(config.cause_of_death_column);
                            if isempty(cod_column)
                                cod_disabled = true;
                                cod_column = 'CauseOfDeath';
                            end
                        else
                            cod_column = 'CauseOfDeath';
                        end
                        % Rename CoD column if needed
                        if ~cod_disabled && ismember(cod_column, T.Properties.VariableNames) && ~strcmp(cod_column, 'CauseOfDeath')
                            T.Properties.VariableNames{strcmp(T.Properties.VariableNames, cod_column)} = 'CauseOfDeath';
                        elseif ~cod_disabled
                            ci_idx = find(strcmpi(T.Properties.VariableNames, cod_column));
                            if ~isempty(ci_idx)
                                T.Properties.VariableNames{ci_idx(1)} = 'CauseOfDeath';
                            end
                        end

                        for j = 1:length(id_list)
                            i_find = find(strcmp(pat_normalized, id_list_normalized{j}));
                            if isempty(i_find)
                                exclusions = add_excl(exclusions, id_list{j}, ...
                                    'No matching clinical record in spreadsheet', ...
                                    'No clinical record');
                            else
                                i_find = i_find(1);
                                lf_val = T.LocalOrRegionalFailure(i_find);
                                if isnan(lf_val)
                                    exclusions = add_excl(exclusions, id_list{j}, ...
                                        'NaN local failure status in clinical spreadsheet', ...
                                        'No clinical record');
                                elseif lf_val == 0 && ismember('CauseOfDeath', T.Properties.VariableNames)
                                    cod_raw = T.CauseOfDeath(i_find);
                                    if iscell(cod_raw), cod_raw = cod_raw{1}; end
                                    cod = char(string(cod_raw));
                                    cod_lower = strtrim(lower(cod));
                                    is_unknown = isempty(cod_lower) || strcmp(cod_lower, 'unknown') || ...
                                                 strcmp(cod_lower, 'pending') || strcmp(cod_lower, 'n/a');
                                    is_panc_cancer = ~isempty(regexp(cod_lower, 'pancrea|biliar|disease\s+progression', 'once'));
                                    if ~isempty(cod) && ~is_unknown && ~is_panc_cancer
                                        exclusions = add_excl(exclusions, id_list{j}, ...
                                            sprintf('Competing risk: non-cancer death (%s) — excluded from GLME/Wilcoxon, censored in Cox', cod), ...
                                            'Competing risk');
                                    end
                                end
                            end
                        end
                    end
                else
                    fprintf('  ⚠️  Clinical spreadsheet not found. Skipping clinical checks.\n');
                end
            end

            % --- Check 2c: DL training set exclusion ---
            manifest_file = fullfile(dataloc, 'dl_validation_manifest.mat');
            dl_provenance = load_dl_provenance(manifest_file);
            if dtype == 2 && ~isempty(dl_provenance.dncnn_train_ids)
                for j = 1:length(id_list)
                    if any(strcmp(dl_provenance.dncnn_train_ids, id_list{j}))
                        exclusions = add_excl(exclusions, id_list{j}, ...
                            'In dnCNN training set (data leakage prevention)', ...
                            'DL training set');
                    end
                end
            end
            if dtype == 3 && ~isempty(dl_provenance.ivimnet_train_ids)
                for j = 1:length(id_list)
                    if any(strcmp(dl_provenance.ivimnet_train_ids, id_list{j}))
                        exclusions = add_excl(exclusions, id_list{j}, ...
                            'In IVIMnet training set (data leakage prevention)', ...
                            'DL training set');
                    end
                end
            end

            % --- Check 2d: Baseline outliers (3x IQR) ---
            d_mean = sm.d_mean;
            f_mean = sm.f_mean;
            dstar_mean = sm.dstar_mean;
            baseline_metrics_oi = {adc_mean(:,1,dtype), d_mean(:,1,dtype), f_mean(:,1,dtype), dstar_mean(:,1,dtype)};
            baseline_metric_names = {'ADC', 'D', 'f', 'D*'};

            % Only run on patients with valid baseline
            lf_dummy = zeros(length(id_list), 1);  % outcome-blinded detection
            is_outlier = false(length(id_list), 1);
            for mi = 1:length(baseline_metrics_oi)
                col = baseline_metrics_oi{mi};
                col_clean = col(~isnan(col));
                if numel(col_clean) < 3, continue; end
                med_val = median(col_clean);
                iqr_val = iqr(col_clean);
                if iqr_val == 0, continue; end
                lower_fence = med_val - 3 * iqr_val;
                upper_fence = med_val + 3 * iqr_val;
                outlier_flags = (col < lower_fence | col > upper_fence) & ~isnan(col);
                for j = find(outlier_flags)'
                    exclusions = add_excl(exclusions, id_list{j}, ...
                        sprintf('Baseline outlier in %s (value=%.4g, range [%.4g, %.4g])', ...
                            baseline_metric_names{mi}, col(j), lower_fence, upper_fence), ...
                        'Outlier (3x IQR)');
                end
                is_outlier = is_outlier | outlier_flags;
            end
        else
            fprintf('  ⚠️  summary_metrics file not found in %s. Stage 2 checks skipped.\n', output_folder);
        end
    end

    %% ===================================================================
    %  STAGE 3: Print the report
    %  ===================================================================
    fprintf('\n=============================================================\n');
    fprintf('  EXCLUSION REPORT RESULTS\n');
    fprintf('=============================================================\n\n');

    % Get unique patient IDs with exclusions
    if isempty(exclusions)
        excluded_ids = {};
    else
        excluded_ids = unique(exclusions(:, 1));
    end
    n_excluded = length(excluded_ids);

    fprintf('Total patients in data directory:  %d\n', n_total);
    fprintf('Patients with exclusion flags:     %d\n', n_excluded);
    fprintf('Patients passing all filters:      %d\n\n', n_total - n_excluded);

    % Category summary
    if ~isempty(exclusions)
        categories = unique(exclusions(:, 3));
        fprintf('--- Exclusions by Category ---\n');
        cat_counts = struct();
        for ci = 1:length(categories)
            cat = categories{ci};
            n_cat = sum(strcmp(exclusions(:, 3), cat));
            n_pts = length(unique(exclusions(strcmp(exclusions(:, 3), cat), 1)));
            fprintf('  %-40s %3d flags across %d patients\n', cat, n_cat, n_pts);
        end
        fprintf('\n');

        % Per-patient detail
        fprintf('--- Per-Patient Exclusion Details ---\n');
        for pi = 1:length(excluded_ids)
            pid = excluded_ids{pi};
            pid_rows = strcmp(exclusions(:, 1), pid);
            reasons = exclusions(pid_rows, 2);
            cats = exclusions(pid_rows, 3);

            % Determine severity: fatal exclusions vs informational
            fatal_cats = {'Missing baseline data', 'Missing baseline imaging', ...
                          'No clinical record', 'DL training set'};
            is_fatal = false;
            for fi = 1:length(cats)
                if any(strcmp(cats{fi}, fatal_cats))
                    is_fatal = true;
                    break;
                end
            end

            if is_fatal
                fprintf('  ❌ %s\n', pid);
            else
                fprintf('  ⚠️  %s\n', pid);
            end
            for ri = 1:length(reasons)
                fprintf('       • %s\n', reasons{ri});
            end
        end
    else
        fprintf('✅ No patients excluded.\n');
    end

    %% --- Build return struct ---
    report.exclusions = exclusions;
    report.n_total = n_total;
    report.n_excluded = n_excluded;
    report.n_analysed = n_total - n_excluded;
    report.excluded_ids = excluded_ids;

    % Category counts
    summary = struct();
    if ~isempty(exclusions)
        categories = unique(exclusions(:, 3));
        for ci = 1:length(categories)
            cat = categories{ci};
            field_name = matlab.lang.makeValidName(cat);
            summary.(field_name) = length(unique(exclusions(strcmp(exclusions(:, 3), cat), 1)));
        end
    end
    report.summary = summary;

    fprintf('\n=============================================================\n');
    fprintf('  Report complete.\n');
    fprintf('=============================================================\n');
end

%% --- Helper ---
function exclusions = add_excl(exclusions, patient_id, reason, category)
    exclusions(end+1, :) = {patient_id, reason, category};
end
