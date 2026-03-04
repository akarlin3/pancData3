function report = patient_data_check(config_path)
% PATIENT_DATA_CHECK Pre-pipeline check for patient data issues.
%
%   report = patient_data_check('config.json')
%
%   Scans the data directory specified in config.json and reports:
%     - Missing or malformed patient folders
%     - Missing fraction folders (Fx1–Fx5, post)
%     - Missing DWI DICOM files
%     - Missing GTV mask files
%     - Missing RT dose folders
%     - dcm2niix availability
%     - Clinical spreadsheet accessibility
%
%   Returns a struct with fields:
%     n_patients       - Number of patient folders found
%     patients         - Table of per-patient issue counts
%     issues           - Cell array of {severity, patient, message} rows
%     n_errors         - Count of critical issues (missing DWI data)
%     n_warnings       - Count of non-critical issues (missing optional data)

    if nargin < 1
        config_path = 'config.json';
    end

    addpath('core', 'utils', 'dependencies');

    %% --- Load and validate configuration ---
    fprintf('📋 Patient Data Check\n');
    fprintf('=====================\n\n');

    config = parse_config(config_path);
    issues = cell(0, 3);  % {severity, patient, message}

    %% --- Check infrastructure ---
    fprintf('⚙️  Checking infrastructure...\n');

    % Data directory
    if ~isfield(config, 'dataloc') || isempty(config.dataloc)
        issues = add_issue(issues, 'ERROR', 'config', 'dataloc not set in config');
    elseif ~isfolder(config.dataloc)
        issues = add_issue(issues, 'ERROR', 'config', ...
            sprintf('Data directory not found: %s', config.dataloc));
    end

    % dcm2niix
    if isfield(config, 'dcm2nii_call') && ~isempty(config.dcm2nii_call)
        [exit_code, ~] = system(sprintf('%s --version', escape_shell_arg(config.dcm2nii_call)));
        if exit_code ~= 0
            issues = add_issue(issues, 'WARNING', 'config', ...
                sprintf('dcm2niix not runnable at: %s', config.dcm2nii_call));
        end
    else
        issues = add_issue(issues, 'WARNING', 'config', 'dcm2nii_call not set in config');
    end

    % Clinical spreadsheet
    if isfield(config, 'clinical_data_sheet') && ~isempty(config.clinical_data_sheet)
        sheet_path = config.clinical_data_sheet;
        if isfield(config, 'dataloc') && ~isempty(config.dataloc)
            sheet_full = fullfile(config.dataloc, sheet_path);
            if ~isfile(sheet_full) && ~isfile(sheet_path)
                issues = add_issue(issues, 'WARNING', 'config', ...
                    sprintf('Clinical spreadsheet not found: %s', sheet_path));
            end
        end
    else
        issues = add_issue(issues, 'WARNING', 'config', 'clinical_data_sheet not set in config');
    end

    %% --- Scan patient directories ---
    if ~isfield(config, 'dataloc') || ~isfolder(config.dataloc)
        fprintf('❌ Cannot scan patients: data directory inaccessible.\n');
        report = build_report(issues, 0);
        return;
    end

    fprintf('⚙️  Scanning patient directories in %s ...\n', config.dataloc);

    patlist = clean_dir_command(config.dataloc);
    patlist = patlist(cellfun(@isempty, strfind({patlist.name}, 'template')));

    % Filter to valid patient folders (numeric ID prefix)
    id_num = zeros(size(patlist));
    for j = 1:length(patlist)
        strtmp = strsplit(strrep(strrep(patlist(j).name, 'P', ''), '_', '-'), '-');
        id_num(j) = str2double(strtmp{1});
    end
    valid_idx = ~isnan(id_num);

    n_skipped = sum(~valid_idx);
    if n_skipped > 0
        issues = add_issue(issues, 'INFO', 'global', ...
            sprintf('%d non-patient folders skipped', n_skipped));
    end

    patlist = patlist(valid_idx);
    [~, id_sort] = sort(id_num(valid_idx), 'ascend');
    patlist = patlist(id_sort);

    n_patients = length(patlist);
    fprintf('   Found %d patient folders.\n\n', n_patients);

    % Optional: filter to requested patient_ids
    if isfield(config, 'patient_ids') && ~isempty(config.patient_ids)
        keep = false(n_patients, 1);
        for j = 1:n_patients
            strtmp = strsplit(strrep(strrep(patlist(j).name, 'P', ''), '_', '-'), '-');
            pid = str2double(strtmp{1});
            if ismember(pid, config.patient_ids)
                keep(j) = true;
            end
        end
        patlist = patlist(keep);
        n_patients = length(patlist);
        fprintf('   Filtered to %d patients per config.patient_ids.\n\n', n_patients);
    end

    fx_search = {'Fx1', 'Fx2', 'Fx3', 'Fx4', 'Fx5', 'post'};

    % Per-patient counters
    pat_names = cell(n_patients, 1);
    n_fx_found = zeros(n_patients, 1);
    n_dwi_found = zeros(n_patients, 1);
    n_gtv_found = zeros(n_patients, 1);
    n_dose_found = zeros(n_patients, 1);
    n_issues_per_pat = zeros(n_patients, 1);

    for j = 1:n_patients
        pat_name = patlist(j).name;
        pat_names{j} = pat_name;
        basefolder = fullfile(config.dataloc, pat_name);
        basefolder_contents = clean_dir_command(basefolder);

        for fi = 1:length(fx_search)
            fx_label = fx_search{fi};
            fxtmp_idx = ~cellfun(@isempty, strfind({basefolder_contents.name}, fx_label));
            fxtmp = basefolder_contents(fxtmp_idx);

            if isempty(fxtmp)
                % Fx1 is critical; others are warnings
                if fi == 1
                    issues = add_issue(issues, 'ERROR', pat_name, ...
                        sprintf('Missing %s folder (baseline)', fx_label));
                    n_issues_per_pat(j) = n_issues_per_pat(j) + 1;
                else
                    issues = add_issue(issues, 'INFO', pat_name, ...
                        sprintf('No %s folder', fx_label));
                end
                continue;
            end

            n_fx_found(j) = n_fx_found(j) + 1;
            fxfolder = fullfile(basefolder, fxtmp(1).name);

            % --- Check DWI data ---
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
                        n_dwi_found(j) = n_dwi_found(j) + 1;
                    end
                end
            end

            if ~has_dwi
                sev = 'ERROR';
                if fi > 1; sev = 'WARNING'; end
                issues = add_issue(issues, sev, pat_name, ...
                    sprintf('No DWI DICOM files in %s', fx_label));
                n_issues_per_pat(j) = n_issues_per_pat(j) + 1;
            end

            % --- Check GTV masks ---
            gtv_files = dir(fullfile(fxfolder, '*GTV*.mat'));
            if isempty(gtv_files)
                sev = 'WARNING';
                if fi == 1; sev = 'ERROR'; end
                issues = add_issue(issues, sev, pat_name, ...
                    sprintf('No GTV mask found for %s', fx_label));
                n_issues_per_pat(j) = n_issues_per_pat(j) + 1;
            else
                n_gtv_found(j) = n_gtv_found(j) + 1;
            end

            % --- Check RT dose (fractions only, not post) ---
            if fi <= 5
                dose_search = clean_dir_command(fullfile(fxfolder, '*rtdose*'));
                has_dose = false;
                if ~isempty(dose_search)
                    for di = 1:length(dose_search)
                        if ~isempty(dir(fullfile(dose_search(di).folder, dose_search(di).name, '*.dcm')))
                            has_dose = true;
                            break;
                        end
                    end
                end
                if has_dose
                    n_dose_found(j) = n_dose_found(j) + 1;
                else
                    issues = add_issue(issues, 'INFO', pat_name, ...
                        sprintf('No RT dose folder for %s', fx_label));
                end
            end
        end
    end

    %% --- Print summary ---
    fprintf('\n📊 Summary\n');
    fprintf('==========\n');
    fprintf('Patients scanned:  %d\n', n_patients);

    n_errors = sum(strcmp(issues(:,1), 'ERROR'));
    n_warnings = sum(strcmp(issues(:,1), 'WARNING'));
    n_info = sum(strcmp(issues(:,1), 'INFO'));

    if n_errors > 0
        fprintf('❌ Errors:   %d\n', n_errors);
    end
    if n_warnings > 0
        fprintf('⚠️  Warnings: %d\n', n_warnings);
    end
    fprintf('💡 Info:     %d\n\n', n_info);

    % Print errors and warnings
    if n_errors + n_warnings > 0
        fprintf('--- Issues ---\n');
        for i = 1:size(issues, 1)
            sev = issues{i,1};
            if strcmp(sev, 'ERROR')
                fprintf('❌ [%s] %s\n', issues{i,2}, issues{i,3});
            elseif strcmp(sev, 'WARNING')
                fprintf('⚠️  [%s] %s\n', issues{i,2}, issues{i,3});
            end
        end
        fprintf('\n');
    end

    % Per-patient overview (patients with issues only)
    problem_pats = find(n_issues_per_pat > 0);
    if ~isempty(problem_pats)
        fprintf('--- Patients with issues ---\n');
        fprintf('%-20s %s/%s/%s/%s\n', 'Patient', 'Fractions', 'DWI', 'GTV', 'Dose');
        for idx = problem_pats'
            fprintf('%-20s %d/6      %d    %d    %d\n', ...
                pat_names{idx}, n_fx_found(idx), n_dwi_found(idx), ...
                n_gtv_found(idx), n_dose_found(idx));
        end
        fprintf('\n');
    end

    if n_errors == 0
        fprintf('✅ No critical issues found. Pipeline should be able to run.\n');
    else
        fprintf('❌ %d critical issue(s) found. Fix these before running the pipeline.\n', n_errors);
    end

    %% --- Build return struct ---
    report = build_report(issues, n_patients);
    if n_patients > 0
        report.patients = table(pat_names, n_fx_found, n_dwi_found, n_gtv_found, ...
            n_dose_found, n_issues_per_pat, ...
            'VariableNames', {'Patient', 'Fractions', 'DWI', 'GTV', 'Dose', 'Issues'});
    end
end

%% --- Helper functions ---

function issues = add_issue(issues, severity, patient, message)
    issues(end+1, :) = {severity, patient, message};
end

function report = build_report(issues, n_patients)
    report.n_patients = n_patients;
    report.issues = issues;
    report.n_errors = sum(strcmp(issues(:,1), 'ERROR'));
    report.n_warnings = sum(strcmp(issues(:,1), 'WARNING'));
end
