function [id_list, mrn_list, fx_dates, dwi_locations, rtdose_locations, gtv_locations, gtvn_locations] = discover_patient_files(dataloc)
% DISCOVER_PATIENT_FILES — Locates DWI, GTV, and RT dose files on the network
% Part of the load_dwi_data.m refactoring.
%
% Inputs:
%   dataloc           - Base directory path where patient data is stored
%
% Outputs:
%   id_list           - Cell array of patient folder names
%   mrn_list          - Cell array of patient medical record numbers
%   fx_dates          - Cell matrix of study dates for each fraction
%   dwi_locations     - Cell matrix of absolute paths to DWI DICOM folders
%   rtdose_locations  - Cell matrix of absolute paths to RT Dose DICOM folders
%   gtv_locations     - Cell matrix of absolute paths to primary GTV files
%   gtvn_locations    - Cell matrix of absolute paths to nodal GTV files
%
% ANALYTICAL RATIONALE — AUTOMATED FILE DISCOVERY
%   Clinical MRI research data is stored in a hierarchical folder structure:
%     dataloc/P##-ID/Fx#/DWI#/   (DWI DICOM files)
%     dataloc/P##-ID/Fx#/        (GTV mask .mat files, rtdose DICOM folder)
%   This function automates the mapping from folder structure to the
%   patient x fraction x repeat indexing scheme used throughout the pipeline.
%   Automated discovery is preferred over manual lookup tables because:
%     1. New patients can be added by simply creating properly named folders
%     2. Missing data (e.g., patient without Fx3 scan) is naturally handled
%        as empty cells, which propagate as NaN in downstream computations
%     3. Repeat scans at Fx1 (for repeatability/reproducibility analysis)
%        are automatically detected by enumerating DWI subfolders
%
%   The function also extracts MRN (Medical Record Number) and StudyDate
%   from DICOM headers, linking the imaging data to clinical outcome data
%   in the spreadsheet. This avoids manual data entry errors that could
%   corrupt patient-level correlations.

% List all patient folders, excluding any 'template' directories
% (template folders contain example data structures, not patient data)
% clean_dir_command wraps dir() and removes '.' and '..' entries.
patlist = clean_dir_command(dataloc);
% Exclude template directories (contain example folder structures, not real data)
patlist = patlist(cellfun(@isempty, strfind({patlist.name},'template')));

% Sort patients by numeric ID extracted from folder names.
% Folder naming convention is "P##-descriptor" (e.g., "P12-ABC", "P03_XYZ").
% The numeric prefix determines processing order, which also determines
% the row index in all output arrays. Consistent ordering is essential
% because downstream analysis uses row indices to cross-reference between
% imaging data (dwi_locations, gtv_locations) and clinical data (lf, immuno).
id_num = zeros(size(patlist));
for j=1:length(patlist)
    % Strip 'P' prefix and normalize separators, then parse leading number
    strtmp = strsplit(strrep(strrep(patlist(j).name,'P',''),'_','-'),'-');
    id_num(j) = str2double(strtmp{1});
end

% Filter out non-patient folders (e.g., "README", "scripts") whose names
% do not start with a numeric ID (str2double returns NaN for these).
valid_idx = ~isnan(id_num);
patlist = patlist(valid_idx);
id_num = id_num(valid_idx);

% Sort ascending by numeric ID so row index in output arrays is deterministic
[~,id_sort] = sort(id_num,'ascend');
patlist = patlist(id_sort);

%% --- Locate DWI, GTV mask, and RT dose paths for every patient x fraction ---

% Fractions to search for: 5 on-treatment fractions + 1 post-treatment scan.
% Pancreatic SBRT typically delivers 5 fractions over 1-2 weeks. Each
% fraction has an associated DWI scan for longitudinal biomarker tracking.
% The "post" scan is acquired 4-6 weeks after RT completion to assess
% early treatment response. This temporal sampling captures the evolution
% of diffusion parameters during and shortly after radiotherapy.
fx_search = {'Fx1','Fx2','Fx3','Fx4','Fx5','post'};

% Pre-allocate identifier and date storage
id_list = cell(length(patlist),1);       % patient folder names (used as analysis IDs)
mrn_list = cell(length(patlist),1);      % medical record numbers (links imaging to clinical records)
fx_dates = cell(length(patlist),6);      % DICOM StudyDate per fraction (for temporal spacing analysis)

% Cell arrays to hold file paths: (patient, fraction, repeat-scan index).
% The 3rd dimension (up to 6 repeats) accommodates Fx1 repeatability
% acquisitions where the patient is scanned multiple times in the same
% session. These repeats enable within-session coefficient of variation
% (wCV) calculation, which quantifies measurement reproducibility and
% helps distinguish true treatment-induced changes from measurement noise.
dwi_locations = cell(length(patlist),6,6);
rtdose_locations = cell(length(patlist),6);   % RT dose only for Fx1-Fx5 (no dose for post-RT)
gtv_locations = cell(length(patlist),6,6);    % primary tumor (GTVp) contour masks
gtvn_locations = cell(length(patlist),6,6);   % nodal disease (GTVn) masks — present for a subset of cases

% --- Main discovery loop: iterate over patients × fractions ---
n_pat_discover = length(patlist);
for j=1:n_pat_discover
    text_progress_bar(j, n_pat_discover, 'Discovering patient files');
    basefolder = fullfile(dataloc, patlist(j).name);
    basefolder_contents = clean_dir_command(basefolder);  % list all subdirectories (Fx1, Fx2, ...)
    id_list{j} = patlist(j).name;
    have_mrn = 0;    % flag: MRN already extracted for this patient (only need one per patient)

    for fi=1:length(fx_search)
        have_fx_date = 0;   % flag: study date already extracted for this fraction

        % Anchor the fraction-folder name at the start and require a
        % non-alphanumeric separator (space, hyphen, underscore, …) or end
        % of string after it. This prevents 'Fx1' matching 'Fx10'/'Fx11'
        % while still matching real-world variants like
        % 'Fx1 - repeatability' that some sites use for the repeat-scan
        % session. An earlier version used '\b' here, but some MATLAB
        % runtimes did not treat a following space as a word boundary for
        % this pattern in practice, causing every 'Fx1 - repeatability'
        % folder to be silently skipped and cascading into all-NaN
        % dice_rpt_* across the cohort. Case-insensitive via regexpi.
        fxtmp_idx = ~cellfun(@isempty, regexpi({basefolder_contents.name}, ['^' fx_search{fi} '([^A-Za-z0-9]|$)'], 'once'));
        fxtmp = basefolder_contents(fxtmp_idx);

        if isempty(fxtmp)
            fprintf('%s, no %s folder\n',patlist(j).name,fx_search{fi});
        else
            fxfolder = fullfile(basefolder, fxtmp(1).name);
            % ok, now find the DWI data for that fraction (there may be 1
            % or more for repeatability analysis)
            % Search for DWI DICOM folders. The naming convention varies across
            % institutions and scanner vendors: some use "DWI_1", "DWI_2" for
            % repeat scans, others use a generic "DICOM" folder. The fallback
            % to *DICOM* handles sites that did not use DWI-specific naming.
            % T2-weighted folders are excluded from the DICOM fallback because
            % they contain structural (not diffusion) data and would cause
            % model fitting to fail silently with meaningless parameter values.
            dwi_search = clean_dir_command(fullfile(fxfolder, '*DWI*'));
            % Fallback: if no DWI-specific folder, try generic DICOM folders
            % but exclude T2-weighted dirs (structural, not diffusion data)
            if isempty(dwi_search), dwi_search = clean_dir_command(fullfile(fxfolder, '*DICOM*')); dwi_search = dwi_search(cellfun(@isempty, strfind({dwi_search.name}, 't2'))); end

            % Verify each candidate folder actually contains .dcm files.
            % Empty folders or folders with only metadata (e.g., DICOMDIR)
            % would cause dcm2niix to fail or produce empty NIfTI volumes.
            if ~isempty(dwi_search)
                contains_dicom = zeros(size(dwi_search));
                for dwii=1:length(dwi_search)
                    if ~isempty(dir(fullfile(dwi_search(dwii).folder, dwi_search(dwii).name, '*.dcm')))
                        contains_dicom(dwii) = 1;
                    else
                        % Some institutions nest DICOMs in a subdirectory.
                        % This two-level search handles both flat and nested layouts.
                        if ~isempty(dir(fullfile(dwi_search(dwii).folder, dwi_search(dwii).name, 'DICOM', '*.dcm')))
                            contains_dicom(dwii) = 1;
                            dwi_search(dwii).name = fullfile(dwi_search(dwii).name, 'DICOM');
                        end
                    end
                end
                dwi_search = dwi_search(contains_dicom==1);  % keep only valid DICOM dirs
            end

            if isempty(dwi_search)
                fprintf('%s, no DWI data found for %s\n',patlist(j).name,fx_search{fi});
            else
                % --- Store DWI path & extract DICOM metadata for each repeat scan ---
                for dwii=1:length(dwi_search)
                    dwi_locations{j,fi,dwii} = fullfile(dwi_search(dwii).folder, dwi_search(dwii).name);

                    if have_fx_date==0
                        % Extract MRN and StudyDate from a DICOM header.
                        % Using the 5th file (or last available) rather than
                        % the 1st avoids incomplete headers that some scanners
                        % write for the initial localizer/scout images.
                        % The MRN links imaging data to the clinical outcomes
                        % spreadsheet; StudyDate enables temporal analysis
                        % (e.g., days between fractions, treatment duration).
                        dicom_files = dir(fullfile(dwi_locations{j,fi,dwii}, '*.dcm'));
                        if isempty(dicom_files)
                            fprintf('  ⚠️  No .dcm files found in %s — skipping MRN/date extraction.\n', dwi_locations{j,fi,dwii});
                            continue;
                        end
                        dcm_idx = min(5, length(dicom_files));
                        try
                            pat_data = dicominfo(fullfile(dicom_files(dcm_idx).folder, dicom_files(dcm_idx).name));
                            if have_mrn==0
                                mrn_list{j} = pat_data.PatientID;
                                id_list{j} = patlist(j).name;
                                have_mrn=1;
                            end
                            fx_dates{j,fi} = pat_data.StudyDate;
                            have_fx_date=1;
                        catch ME_dcm
                            fprintf('  ⚠️  Failed to read DICOM header from %s: %s\n', ...
                                dicom_files(dcm_idx).name, ME_dcm.message);
                            continue;
                        end
                    end

                    % --- Locate GTV mask .mat files for this DWI repeat ---
                    % find_gtv_files searches for GTVp (primary tumor) and
                    % GTVn (nodal disease) mask files with flexible naming
                    % patterns (e.g., "gtv_1.mat", "GTV_p.mat", etc.)
                    [gtvp_path, gtvn_path] = find_gtv_files(fxfolder, dwii, patlist(j).name);

                    if ~isempty(gtvp_path)
                        gtv_locations{j,fi,dwii} = gtvp_path;
                    end
                    if ~isempty(gtvn_path)
                        gtvn_locations{j,fi,dwii} = gtvn_path;
                    end
                end

            end

            % --- Locate RT dose DICOM folder for this fraction ---
            % RT dose files contain the spatial distribution of radiation
            % dose delivered to the patient at each fraction. These are
            % essential for dose-response analysis: correlating local dose
            % (e.g., D95, V50Gy) with diffusion parameter changes to
            % determine if higher radiation doses produce greater changes
            % in tumor cellularity (as reflected by ADC/D changes).
            dose_search = clean_dir_command(fullfile(fxfolder, '*rtdose*'));
            if ~isempty(dose_search)
                % Verify the rtdose folder contains DICOM files
                contains_dicom = zeros(size(dose_search));
                for di=1:length(dose_search)
                    if ~isempty(dir(fullfile(dose_search(di).folder, dose_search(di).name, '*.dcm')))
                        contains_dicom(di) = 1;
                    end
                end
                dose_search = dose_search(contains_dicom==1);
            end
            if ~isempty(dose_search)
                fprintf('%s/%s: found rtdose\n',patlist(j).name,fx_search{fi});
                rtdose_locations{j,fi} = fullfile(fxfolder, dose_search(1).name);
            end
        end
    end
end
end
