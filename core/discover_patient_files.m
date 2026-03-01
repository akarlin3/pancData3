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

% List all patient folders, excluding any 'template' directories
patlist = clean_dir_command(dataloc);
patlist = patlist(~contains({patlist.name},'template'));

% sort by pat id — extract the numeric patient ID from folder names
% (e.g., "P12-ABC" → 12) and sort in ascending order
id_num = zeros(size(patlist));
for j=1:length(patlist)
    strtmp = strsplit(strrep(strrep(patlist(j).name,'P',''),'_','-'),'-');
    id_num(j) = str2double(strtmp{1});
end

% Filter out non-patient folders (which evaluate to NaN)
valid_idx = ~isnan(id_num);
patlist = patlist(valid_idx);
id_num = id_num(valid_idx);

[~,id_sort] = sort(id_num,'ascend');
patlist = patlist(id_sort);

%% --- Locate DWI, GTV mask, and RT dose paths for every patient × fraction ---

% Fractions to search for: 5 on-treatment fractions + 1 post-treatment scan
fx_search = {'Fx1','Fx2','Fx3','Fx4','Fx5','post'};

% Pre-allocate identifier and date storage
id_list = cell(length(patlist),1);       % patient folder names
mrn_list = cell(length(patlist),1);      % medical record numbers (from DICOM headers)
fx_dates = cell(length(patlist),6);      % study dates per fraction (from DICOM headers)

% Cell arrays to hold file paths: (patient, fraction, repeat-scan index)
% Up to 6 fractions and 6 repeat scans (Fx1 repeatability acquisitions)
dwi_locations = cell(length(patlist),6,6);
rtdose_locations = cell(length(patlist),6);   % RT dose only for Fx1–Fx5
gtv_locations = cell(length(patlist),6,6);
gtvn_locations = cell(length(patlist),6,6); % a few cases have a nodal gtv as well

% --- Main discovery loop: iterate over patients × fractions ---
for j=1:length(patlist)
    basefolder = fullfile(dataloc, patlist(j).name);
    have_mrn = 0;    % flag: MRN already extracted for this patient
    for fi=1:length(fx_search)
        have_fx_date = 0;   % flag: study date already extracted for this fraction
        fxtmp = clean_dir_command(fullfile(basefolder, ['*' fx_search{fi} '*']));
        if isempty(fxtmp)
            fprintf('%s, no %s folder\n',patlist(j).name,fx_search{fi});
        else
            fxfolder = fullfile(basefolder, fxtmp(1).name);
            % ok, now find the DWI data for that fraction (there may be 1
            % or more for repeatability analysis)
            dwi_search = clean_dir_command(fullfile(fxfolder, '*DWI*'));
            % Fallback: some sites stored DICOMs under a generic 'DICOM' folder
            if isempty(dwi_search), dwi_search = clean_dir_command(fullfile(fxfolder, '*DICOM*')); dwi_search = dwi_search(~contains({dwi_search.name}, 't2')); end

            % Verify each candidate folder actually contains .dcm files
            if ~isempty(dwi_search)
                contains_dicom = zeros(size(dwi_search));
                for dwii=1:length(dwi_search)
                    if ~isempty(dir(fullfile(dwi_search(dwii).folder, dwi_search(dwii).name, '*.dcm')))
                        contains_dicom(dwii) = 1;
                    else
                        % Some folders nest DICOMs in a subdirectory called 'DICOM'
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
                        % get the MRN and study date from an arbitrary DICOM header
                        dicom_files = dir(fullfile(dwi_locations{j,fi,dwii}, '*.dcm'));
                        pat_data = dicominfo(fullfile(dicom_files(5).folder, dicom_files(5).name));
                        if have_mrn==0
                            mrn_list{j} = pat_data.PatientID;
                            id_list{j} = patlist(j).name;
                            have_mrn=1;
                        end
                        fx_dates{j,fi} = pat_data.StudyDate;
                        have_fx_date=1;
                    end

                    % --- Locate GTV mask .mat files for this DWI repeat ---
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
            % find rtdose
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
