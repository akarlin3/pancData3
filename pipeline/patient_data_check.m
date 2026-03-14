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
%
% [ANALYTICAL RATIONALE — WHY THIS PRE-FLIGHT CHECK EXISTS]:
% The DWI analysis pipeline processes a large patient cohort and takes hours
% to complete. Discovering missing data mid-run (e.g., a patient missing their
% baseline Fx1 DWI folder, or a GTV mask file that was never exported from the
% treatment planning system) wastes compute time and produces incomplete results
% that may silently bias downstream statistics (e.g., survivorship bias if
% patients with missing data are systematically different from those with
% complete data).
%
% This function performs a comprehensive audit of the data directory structure
% BEFORE committing to a pipeline run, categorizing issues by severity:
%   ERROR   — Data gaps that will cause the pipeline to fail or produce invalid
%             results (e.g., missing baseline DWI, missing baseline GTV mask).
%             These MUST be resolved before running.
%   WARNING — Data gaps that degrade but don't prevent analysis (e.g., missing
%             late-fraction DWI, missing RT dose for dosimetry sub-analysis).
%   INFO    — Expected or benign observations (e.g., non-patient folders in the
%             data directory, patients with fewer than 5 fractions).
%
% [DATA DIRECTORY STRUCTURE EXPECTED]:
% The pipeline expects a specific directory hierarchy reflecting the clinical
% workflow of pancreatic cancer radiotherapy:
%   dataloc/
%     PatientID_Name/           (one folder per patient, numeric ID prefix)
%       Fx1_date/               (baseline fraction — REQUIRED)
%         *DWI*/                (DWI DICOM series folder)
%           *.dcm               (DICOM files)
%         *GTV*.mat             (GTV mask from treatment planning, MATLAB format)
%         *rtdose*/             (RT dose DICOM folder)
%       Fx2_date/ ... Fx5_date/ (treatment fractions 2-5)
%       post_date/              (post-treatment follow-up scan)

    if nargin < 1
        script_dir = fileparts(mfilename('fullpath'));
        config_path = fullfile(script_dir, '..', 'config.json');
    end

    % Add pipeline modules to MATLAB path so that parse_config, escape_shell_arg,
    % and clean_dir_command are accessible from any working directory.
    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, 'core'), fullfile(script_dir, 'utils'), fullfile(script_dir, 'dependencies'));

    %% --- Load and validate configuration ---
    fprintf('📋 Patient Data Check\n');
    fprintf('=====================\n\n');

    config = parse_config(config_path);
    % Issue accumulator: each row is {severity, patient_name, message}.
    % Severity levels: 'ERROR' (pipeline will fail), 'WARNING' (degraded
    % analysis), 'INFO' (informational/benign observation).
    issues = cell(0, 3);

    %% --- Check infrastructure ---
    % Verify that the three critical infrastructure components are accessible
    % before scanning patient data. These are prerequisites for any pipeline run:
    %   1. Data directory (dataloc) — where all patient imaging data lives
    %   2. dcm2niix — external DICOM-to-NIFTI converter needed by load_dwi_data
    %   3. Clinical spreadsheet — contains outcome data (survival, local failure)
    %      needed by metrics_baseline and all downstream statistical analyses
    fprintf('⚙️  Checking infrastructure...\n');

    % Data directory — the root folder containing all patient subdirectories
    if ~isfield(config, 'dataloc') || isempty(config.dataloc)
        issues = add_issue(issues, 'ERROR', 'config', 'dataloc not set in config');
    elseif ~isfolder(config.dataloc)
        issues = add_issue(issues, 'ERROR', 'config', ...
            sprintf('Data directory not found: %s', config.dataloc));
    end

    % dcm2niix — the DICOM-to-NIFTI converter (from MRIcroGL).
    % This is essential for the load step: raw DICOM DWI files from the MRI
    % scanner must be converted to NIFTI format before voxel-wise model fitting.
    % dcm2niix handles the complex DICOM header parsing (b-value extraction,
    % slice ordering, orientation matrices) that would be error-prone to
    % implement manually. The escape_shell_arg call prevents command injection
    % if the dcm2niix path contains special characters.
    % Classified as WARNING (not ERROR) because dcm2niix is only needed for
    % the load step — if the user has already converted data and is rerunning
    % only metrics, dcm2niix is not required.
    if isfield(config, 'dcm2nii_call') && ~isempty(config.dcm2nii_call)
        [exit_code, ~] = system(sprintf('%s --version', escape_shell_arg(config.dcm2nii_call)));
        if exit_code ~= 0
            issues = add_issue(issues, 'WARNING', 'config', ...
                sprintf('dcm2niix not runnable at: %s', config.dcm2nii_call));
        end
    else
        issues = add_issue(issues, 'WARNING', 'config', 'dcm2nii_call not set in config');
    end

    % Clinical spreadsheet — contains patient outcome data (local failure
    % status, survival time, follow-up duration, treatment details) and
    % clinical covariates (GTV volume, dose coverage metrics) that are
    % joined with imaging-derived metrics in metrics_baseline.
    % Without this spreadsheet, the pipeline can still compute diffusion
    % parameters but cannot perform any outcome correlation, statistical
    % comparison, predictive modeling, or survival analysis — hence WARNING.
    % The spreadsheet is searched both relative to dataloc and as an absolute
    % path, since clinical data may be stored separately from imaging data
    % for PHI management reasons.
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

    % List all subdirectories in the data root, excluding template folders
    % (which contain reference data or contour templates, not patient data).
    % clean_dir_command is a wrapper around dir() that excludes '.' and '..'
    % entries and handles cross-platform path separators.
    patlist = clean_dir_command(config.dataloc);
    patlist = patlist(cellfun(@isempty, strfind({patlist.name}, 'template')));

    % Filter to valid patient folders by extracting the numeric patient ID
    % from the folder name. Patient folders follow the naming convention:
    %   "123-LastName" or "P123_LastName" or "123_LastName-FirstName"
    % The parsing strips the optional 'P' prefix and splits on '-' or '_'
    % to isolate the leading numeric ID. Folders that do not start with a
    % numeric ID (e.g., "config_backup", "scripts") are excluded — they are
    % infrastructure folders, not patient data.
    id_num = zeros(size(patlist));
    for j = 1:length(patlist)
        strtmp = strsplit(strrep(strrep(patlist(j).name, 'P', ''), '_', '-'), '-');
        id_num(j) = str2double(strtmp{1});
    end
    % Folders whose names do not start with a numeric ID produce NaN from
    % str2double and are classified as non-patient (infrastructure) folders.
    valid_idx = ~isnan(id_num);

    n_skipped = sum(~valid_idx);
    if n_skipped > 0
        issues = add_issue(issues, 'INFO', 'global', ...
            sprintf('%d non-patient folders skipped', n_skipped));
    end

    % Sort patients by numeric ID for deterministic, human-readable output.
    patlist = patlist(valid_idx);
    [~, id_sort] = sort(id_num(valid_idx), 'ascend');
    patlist = patlist(id_sort);

    n_patients = length(patlist);
    fprintf('   Found %d patient folders.\n\n', n_patients);

    % Optional: filter to requested patient_ids from config.json.
    % This supports targeted data checks when the researcher wants to verify
    % data completeness for a specific subset of patients (e.g., newly added
    % patients, or patients who failed in a previous pipeline run) without
    % scanning the entire cohort.
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

    % [FRACTION STRUCTURE — PANCREATIC CANCER RT WORKFLOW]:
    % Pancreatic cancer patients undergoing stereotactic body radiation therapy
    % (SBRT) or hypofractionated RT typically receive 5 treatment fractions
    % (Fx1-Fx5) over 1-2 weeks, with DWI MRI acquired at each fraction to
    % monitor intra-treatment response. A post-treatment scan captures the
    % early post-RT state (typically 4-6 weeks after completion).
    %
    % Fx1 (baseline) is the most critical — it establishes pre-treatment
    % diffusion characteristics against which all changes are measured.
    % Missing Fx1 data is classified as ERROR because without baseline,
    % percent change metrics are undefined and the patient cannot contribute
    % to longitudinal or predictive analyses.
    % Missing Fx2-Fx5 or post is classified as INFO (not ERROR) because
    % many patients have incomplete fraction data due to clinical scheduling
    % constraints, patient withdrawal, or scanner availability — this is
    % expected and handled via KNN imputation in the metrics pipeline.
    fx_search = {'Fx1', 'Fx2', 'Fx3', 'Fx4', 'Fx5', 'post'};

    % Per-patient counters track data completeness across the four
    % data types needed for full analysis: fraction folders (temporal
    % coverage), DWI DICOM files (diffusion imaging data), GTV masks
    % (tumor delineation), and RT dose files (dosimetry analysis).
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

        % Check each expected fraction folder within this patient's directory.
        % Uses substring matching (strfind) rather than exact match because
        % fraction folders typically include the date suffix (e.g., "Fx1_20230315").
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
            % Search for DWI DICOM series folders. The primary search looks
            % for folders containing "DWI" in the name (standard naming from
            % most MRI scanner export protocols). The fallback searches for
            % "DICOM" folders while excluding T2-weighted series — some
            % institutions export all MRI series into a generic "DICOM"
            % folder structure rather than protocol-named subfolders.
            % T2 series are excluded because they are anatomical images,
            % not diffusion-weighted, and would produce meaningless IVIM fits.
            dwi_search = clean_dir_command(fullfile(fxfolder, '*DWI*'));
            if isempty(dwi_search)
                dwi_search = clean_dir_command(fullfile(fxfolder, '*DICOM*'));
                dwi_search = dwi_search(cellfun(@isempty, strfind({dwi_search.name}, 't2')));
            end

            % Verify that at least one DWI series folder contains actual DICOM
            % files (*.dcm). An empty DWI folder (created by a partial export
            % or interrupted transfer) would cause dcm2niix to fail silently,
            % producing no NIFTI output and ultimately NaN-filled parameter maps.
            % The nested DICOM/ subfolder check handles the common case where
            % the scanner exports into DWI_series_name/DICOM/*.dcm.
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

            % Missing DWI at Fx1 (baseline) is ERROR — the patient cannot
            % contribute to any analysis without baseline diffusion data.
            % Missing DWI at later fractions is WARNING — the patient can
            % still contribute baseline metrics and partial longitudinal data.
            if ~has_dwi
                sev = 'ERROR';
                if fi > 1; sev = 'WARNING'; end
                issues = add_issue(issues, sev, pat_name, ...
                    sprintf('No DWI DICOM files in %s', fx_label));
                n_issues_per_pat(j) = n_issues_per_pat(j) + 1;
            end

            % --- Check GTV masks ---
            % GTV (Gross Tumor Volume) masks define the 3D tumor boundary
            % within which voxel-wise diffusion parameters are extracted.
            % These are typically exported from the radiation treatment
            % planning system (e.g., Eclipse, RayStation) as MATLAB .mat
            % files containing a binary 3D volume registered to the DWI
            % image space.
            %
            % Missing GTV at Fx1 is ERROR: without knowing WHERE the tumor
            % is at baseline, voxel extraction is impossible and no diffusion
            % metrics can be computed for this patient.
            % Missing GTV at later fractions is WARNING: the pipeline can
            % propagate the baseline GTV mask to subsequent fractions via
            % deformable image registration (apply_dir_mask_propagation.m),
            % though this introduces registration uncertainty.
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
            % RT dose DICOM files contain the 3D dose distribution delivered
            % to the patient at each fraction. These are needed by
            % metrics_dosimetry to compute dose-volume metrics (D95, V50)
            % within diffusion-defined tumor sub-volumes.
            %
            % Only checked for Fx1-Fx5 (fi <= 5), not post-treatment,
            % because RT dose is only delivered during treatment fractions.
            % The post-treatment scan has no associated dose delivery.
            %
            % Missing dose is classified as INFO (not ERROR or WARNING)
            % because dosimetry is an optional sub-analysis — the core
            % pipeline (diffusion fitting, longitudinal analysis, survival)
            % can operate without dose data. Dose coverage features are
            % simply absent from the predictive feature set.
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
    % The summary report is designed for rapid triage: the researcher can
    % immediately see whether the cohort is pipeline-ready (zero errors) or
    % requires data collection/export before proceeding. The per-patient
    % table for problem patients shows the data completeness pattern at a
    % glance (e.g., a patient with 3/6 fractions and 0 GTV masks clearly
    % needs mask export from the treatment planning system).
    fprintf('\n📊 Summary\n');
    fprintf('==========\n');
    fprintf('Patients scanned:  %d\n', n_patients);

    % Tally issues by severity level for the summary header.
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
    % The return struct provides machine-readable output for programmatic use
    % (e.g., automated pre-pipeline checks in execute_all_workflows, or
    % generating data completeness reports for IRB submissions).
    % The patients table provides a compact per-patient overview with counts
    % out of the maximum possible (6 fractions, 6 DWI series, 6 GTV masks,
    % 5 dose files). Low counts relative to the maximum indicate incomplete
    % data that may require follow-up with the clinical data management team.
    report = build_report(issues, n_patients);
    if n_patients > 0
        report.patients = table(pat_names, n_fx_found, n_dwi_found, n_gtv_found, ...
            n_dose_found, n_issues_per_pat, ...
            'VariableNames', {'Patient', 'Fractions', 'DWI', 'GTV', 'Dose', 'Issues'});
    end
end

%% --- Helper functions ---

function issues = add_issue(issues, severity, patient, message)
% Append a single issue row to the issues cell array accumulator.
    issues(end+1, :) = {severity, patient, message};
end

function report = build_report(issues, n_patients)
% Construct the machine-readable return struct from the accumulated issues.
    report.n_patients = n_patients;
    report.issues = issues;
    report.n_errors = sum(strcmp(issues(:,1), 'ERROR'));
    report.n_warnings = sum(strcmp(issues(:,1), 'WARNING'));
end
