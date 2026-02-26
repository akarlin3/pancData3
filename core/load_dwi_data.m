function [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = load_dwi_data(config_struct)
%% load_dwi_data_forAvery.m
% =========================================================================
% PURPOSE
%   Main data loading, processing, and analysis pipeline for a pancreatic
%   cancer DWI (Diffusion-Weighted Imaging) research study. This script
%   discovers patient imaging data on a network share, converts DICOM to
%   NIfTI, fits diffusion models, extracts voxel-level biomarkers within
%   GTV (Gross Tumor Volume) contours, and computes longitudinal summary
%   metrics for downstream statistical analysis.
%
% WORKFLOW (5 sections — run sequentially or start from Section 4)
%   Set skip_to_reload = true (in the USER OPTIONS section below) to skip
%   Sections 1-3 and start directly from Section 4.
%   Section 1 — File Discovery:
%       Scans the network share for patient folders, locating DWI DICOM
%       directories, GTV masks (.mat), nodal GTV masks, and RT dose
%       folders for each patient × fraction combination.
%   Section 2 — DICOM-to-NIfTI Conversion & Model Fitting:
%       Converts DWI DICOMs to NIfTI via dcm2niix, saves GTV masks as
%       NIfTI, resamples RT dose onto DWI geometry, loads volumes, fits
%       ADC (monoexponential) and IVIM (segmented biexponential) models,
%       extracts voxel-level biomarker vectors within GTV masks, computes
%       DVH (Dose-Volume Histogram) parameters, and handles denoised
%       (DnCNN) and IVIMnet pipeline variants.
%   Section 3 — Save:
%       Backs up any previous save file and writes the workspace to a
%       .mat file.
%   Section 4 — Reload (entry point for analysis):
%       Reloads the saved .mat file so downstream scripts can skip
%       Sections 1-3.
%   Section 5 — Longitudinal Summary Metrics:
%       Computes mean, kurtosis, skewness, SD of ADC and IVIM parameters
%       (D, f, D*) per patient × timepoint × DWI-type. Also computes
%       sub-volume metrics, histogram distributions, KS test statistics
%       vs baseline, corruption flags, and repeatability metrics from
%       Fx1 repeat scans.
%
% KEY OUTPUTS
%   data_vectors_gtvp  — struct array (patient × fraction × repeat) with
%                         voxel-level ADC, D, f, D*, dose vectors for the
%                         primary GTV (GTVp).
%   data_vectors_gtvn  — same structure for the nodal GTV (GTVn), when
%                         available.
%   Summary arrays     — adc_mean, d_mean, f_mean, dstar_mean, etc.
%                         (patient × timepoint × dwi_type).
%
% DEPENDENCIES
%   External executables:
%     - dcm2niix (MRIcroGL)      — DICOM-to-NIfTI conversion
%   MATLAB toolboxes / functions:
%     - clean_dir_command         — wrapper around dir() that removes '.' entries
%     - load_untouch_nii          — NIfTI reader (NIfTI toolbox, Jimmy Shen)
%     - IVIMmodelfit              — segmented biexponential IVIM fitting
%     - sample_rtdose_on_image    — resamples RT dose grid onto DWI geometry
%     - dvh                       — computes dose-volume histogram parameters
%     - Image Processing Toolbox  — mat2gray (used for DnCNN normalisation)
%     - Statistics Toolbox        — kurtosis, skewness, kstest2
%
% NOTES
%   - The network path \\pensmph6\... must be mapped / accessible.
%   - Voxel volumes assume isotropic in-plane resolution of ~1.98 mm with
%     5 mm slice thickness (0.19792 × 0.19792 × 0.5 cm³).
%   - Three DWI processing pipelines are compared throughout:
%       dwi_type 1 = standard (raw DWI)
%       dwi_type 2 = DnCNN-denoised DWI + conventional IVIM fit
%       dwi_type 3 = IVIMnet (deep-learning IVIM fit on raw DWI)
% =========================================================================

%% USER OPTIONS
% Set skip_to_reload = true to skip Sections 1-3 (file discovery, conversion,
% and saving) and jump directly to Section 4 (reload saved data).
% This is useful when the .mat save file already exists and you only want
% to re-run the downstream analysis (Section 5).
skip_to_reload = config_struct.skip_to_reload;

% b-value threshold (s/mm²) used for segmented IVIM fitting.
% b-values >= ivim_bthr are used to estimate D (true diffusion) via a 
% monoexponential fit. b-values < ivim_bthr are used to estimate f and D*
% (pseudo-diffusion), with D held fixed from Stage 1.
% Set to 100 (includes b=0, 30 in low set; b=100, 150, 550 in high set).
ivim_bthr = config_struct.ivim_bthr;

%% ========================================================================
fprintf('\n--- SECTION 1: File Discovery ---\n');
%  SECTION 1 — FILE DISCOVERY

if ~skip_to_reload

% find the data (here: W:\ is the mapped location for
% \\pensmph6\mpcsresearch1\aliottae\

% Root path to the pancreas DWI data on the network share
dataloc = config_struct.dataloc;

% List all patient folders, excluding any 'template' directories
patlist= clean_dir_command(dataloc);
patlist = patlist(~contains({patlist.name},'template'));

% sort by pat id — extract the numeric patient ID from folder names
% (e.g., "P12-ABC" → 12) and sort in ascending order
id_num = zeros(size(patlist));
for j=1:length(patlist)
    strtmp = strsplit(strrep(strrep(patlist(j).name,'P',''),'_','-'),'-');
    id_num(j) = str2double(strtmp{1});
end
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
    basefolder = [dataloc patlist(j).name '/'];
    have_mrn = 0;    % flag: MRN already extracted for this patient
    for fi=1:length(fx_search)
        have_fx_date = 0;   % flag: study date already extracted for this fraction
        fxtmp = clean_dir_command([basefolder '*' fx_search{fi} '*']);
        if isempty(fxtmp)
            fprintf('%s, no %s folder\n',patlist(j).name,fx_search{fi});
        else
            fxfolder = [basefolder fxtmp(1).name '/'];
            % ok, now find the DWI data for that fraction (there may be 1
            % or more for repeatability analysis)
            dwi_search = clean_dir_command([fxfolder '*DWI*']);
            % Fallback: some sites stored DICOMs under a generic 'DICOM' folder
            if isempty(dwi_search), dwi_search = clean_dir_command([fxfolder '*DICOM*']); dwi_search = dwi_search(~contains({dwi_search.name}, 't2')); end

            % Verify each candidate folder actually contains .dcm files
            if ~isempty(dwi_search)
                contains_dicom = zeros(size(dwi_search));
                for dwii=1:length(dwi_search)
                    if ~isempty(dir([dwi_search(dwii).folder '/' dwi_search(dwii).name '/*.dcm']))
                        contains_dicom(dwii) = 1;
                    else
                        % Some folders nest DICOMs in a subdirectory called 'DICOM'
                        if ~isempty(dir([dwi_search(dwii).folder '/' dwi_search(dwii).name '/DICOM/*.dcm']))
                            contains_dicom(dwii) = 1;
                            dwi_search(dwii).name = [dwi_search(dwii).name '/DICOM/'];
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
                    dwi_locations{j,fi,dwii} = [dwi_search(dwii).folder '/' dwi_search(dwii).name '/'];

                    if have_fx_date==0
                        % get the MRN and study date from an arbitrary DICOM header
                        dicom_files = dir([dwi_locations{j,fi,dwii} '*.dcm']);
                        pat_data = dicominfo([dicom_files(5).folder '/' dicom_files(5).name]);
                        if have_mrn==0
                            mrn_list{j} = pat_data.PatientID;
                            id_list{j} = patlist(j).name;
                            have_mrn=1;
                        end
                        fx_dates{j,fi} = pat_data.StudyDate;
                        have_fx_date=1;
                    end

                    % --- Locate GTV mask .mat files for this DWI repeat ---
                    % use special logic for the 3 cases with nodal gtvs
                    % ('two' in name indicates patient has both GTVp and GTVn)
                    if ~contains(patlist(j).name,'two')
                        % find associated GTV (need to avoid using the date in
                        % the filenames to query)
                        gtv_search = dir([fxfolder '*GTV*' int2str(dwii) '*.mat']);
                        single_gtv_search = dir([fxfolder '*GTV*.mat']);
                        % this can be simpler if only one GTV exists
                        % If only one GTV mask exists, use it directly
                        if isscalar(single_gtv_search)
                            gtv_locations{j,fi,dwii} = [fxfolder '/' single_gtv_search.name];
                        else
                            % Multiple GTV masks exist — match by repeat index (dwii)
                            % using the second underscore-delimited token in the filename
                            gtv_names = {gtv_search.name};
                            gtv_search_result = zeros(size(gtv_search));
                            for gi=1:length(gtv_names)
                                gtmp = strsplit(gtv_names{gi},'_');
                                gtmp = gtmp{2};
                                regex_pattern = regexptranslate('wildcard', ['*GTV*' int2str(dwii)]);
                                isfound = regexp(gtmp, regex_pattern);
                                if isfound==1
                                    gtv_search_result(gi) = 1;
                                end
                            end
                            % a second format was also used for some cases...
                            % (repeat index as trailing numeric token after last '_')
                            if sum(gtv_search_result)==0
                                for gi=1:length(gtv_names)
                                    gtmp = strsplit(gtv_names{gi},'_');
                                    gtmp = strrep(gtmp{end},'.mat','');
                                    if str2double(gtmp)==dwii
                                        gtv_search_result(gi) = 1;
                                    end
                                end
                            end

                            if sum(gtv_search_result)==0
                                fprintf('%s/%s: No GTV%d found\n',patlist(j).name,fx_search{fi},dwii);
                            end
                            if sum(gtv_search_result)>1
                                fprintf('%s/%s: Redundant GTV%ds found\n',patlist(j).name,fx_search{fi},dwii);
                            end
                            if sum(gtv_search_result)==1
                                gtv_locations{j,fi,dwii} = [fxfolder '/' gtv_search(gtv_search_result==1).name];
                            end
                        end
                    else
                        % --- Special logic for patients with both GTVp and GTVn ---
                        % Search multiple naming conventions for the primary
                        % pancreatic GTV (GTVp): GTV_MR, GTVp, GTV_panc
                        % special logic for GTVn cases.
                        gtv_search1 = dir([fxfolder '*GTV_MR' int2str(dwii) '*.mat']);
                        gtv_search2 = dir([fxfolder '*GTVp' int2str(dwii) '*.mat']);
                        gtv_search3 = dir([fxfolder '*GTV_panc*' int2str(dwii) '*.mat']);
                        single_gtv_search1 = dir([fxfolder '*GTV_MR*.mat']);
                        single_gtv_search2 = dir([fxfolder '*GTVp*.mat']);
                        single_gtv_search3 = dir([fxfolder '*GTV_panc*.mat']);

                        % Combine results from all naming conventions
                        gtv_search = cat(1,gtv_search1,gtv_search2,gtv_search3);
                        single_gtv_search = cat(1,single_gtv_search1,single_gtv_search2,single_gtv_search3);

                        % Search for nodal GTV masks: GTV_LN, GTVn, GTV_node
                        gtvn_search1 = dir([fxfolder '*GTV*LN' int2str(dwii) '*.mat']);
                        gtvn_search2 = dir([fxfolder '*GTVn' int2str(dwii) '*.mat']);
                        gtvn_search3 = dir([fxfolder '*GTV_node*' int2str(dwii) '*.mat']);
                        single_gtvn_search1 = dir([fxfolder '*GTV*LN*.mat']);
                        single_gtvn_search2 = dir([fxfolder '*GTVn*.mat']);
                        single_gtvn_search3 = dir([fxfolder '*GTV_node*.mat']);

                        gtvn_search = cat(1,gtvn_search1,gtvn_search2,gtvn_search3);
                        single_gtvn_search = cat(1,single_gtvn_search1,single_gtvn_search2,single_gtvn_search3);

                        % Assign GTVp path (prefer unique match; fall back to repeat-indexed)
                        if isscalar(single_gtv_search)
                            gtv_locations{j,fi,dwii} = [fxfolder '/' single_gtv_search.name];
                        else
                            if isscalar(gtv_search)
                                gtv_name = gtv_search.name;
                                gtv_locations{j,fi,dwii} = [fxfolder '/' gtv_name];
                            end
                        end

                        % Assign GTVn path (same logic as GTVp)
                        if isscalar(single_gtvn_search)
                            gtvn_locations{j,fi,dwii} = [fxfolder '/' single_gtvn_search.name];
                        else
                            if isscalar(gtvn_search)
                                gtvn_name = gtvn_search.name;
                                gtvn_locations{j,fi,dwii} = [fxfolder '/' gtvn_name];
                            end
                        end
                    end

                end

            end

            % --- Locate RT dose DICOM folder for this fraction ---
            % find rtdose
            dose_search = clean_dir_command([fxfolder '*rtdose*']);
            if ~isempty(dose_search)
                % Verify the rtdose folder contains DICOM files
                contains_dicom = zeros(size(dose_search));
                for di=1:length(dose_search)
                    if ~isempty(dir([dose_search(di).folder '/' dose_search(di).name '/*.dcm']))
                        contains_dicom(di) = 1;
                    end
                end
                dose_search = dose_search(contains_dicom==1);
            end
            if ~isempty(dose_search)
                fprintf('%s/%s: found rtdose\n',patlist(j).name,fx_search{fi});
                rtdose_locations{j,fi} = [fxfolder '/' dose_search(1).name '/'];
            end
        end
    end
end

%% ========================================================================
fprintf('\n--- SECTION 2: DICOM-to-NIFTI Conversion, Model Fitting, & Data Extraction ---\n');
%  SECTION 2 — DICOM-TO-NIFTI CONVERSION, MODEL FITTING & DATA EXTRACTION

% I use two executables (non matlab functions_ which need to be at the
% following two locations. You can adjust these two directories to wherever
% they are installed on your computer.
dcm2nii_call = config_struct.dcm2nii_call;   % dcm2niix path
if isfield(config_struct, 'conv_call'), conv_call = config_struct.conv_call; end  % RT struct converter (unused below)

dataloc = config_struct.dataloc;

% load clinical data (local failure and immunotherapy status per patient)
clinical_data_sheet = fullfile(dataloc, config_struct.clinical_data_sheet);
T = readtable(clinical_data_sheet,'Sheet','Clin List');

% load dwi data locations (optionally reload from a previous run)
data_file = [dataloc 'adc_vectors.mat'];
%load(data_file,'dwi_locations','rtdose_locations','gtv_locations','gtvn_locations','mrn_list','id_list','fx_dates');

fx_search = {'Fx1','Fx2','Fx3','Fx4','Fx5','post'};

% Initialise output struct arrays for GTVp (primary) and GTVn (nodal)
data_vectors_gtvp = struct;
data_vectors_gtvn = struct;

% Pre-allocate summary metric arrays (patient × fraction × repeat)
adc_mean = nan(size(dwi_locations));
adc_kurtosis = nan(size(dwi_locations));

% IVIM true diffusion coefficient D — mean values per pipeline variant
d_mean = nan(size(dwi_locations));         % standard IVIM fit
d_mean_dncnn = nan(size(dwi_locations));   % DnCNN-denoised IVIM fit
d_mean_ivimnet = nan(size(dwi_locations)); % IVIMnet deep-learning fit

d_kurtosis = nan(size(dwi_locations));
d_skewness = nan(size(dwi_locations));

% DVH parameters within GTVp and GTVn (patient × fraction)
% Sized to 6 columns to accommodate on-treatment fractions + Post-RT scan
dmean_gtvp = nan(length(mrn_list), 6);  % mean dose in GTVp (Gy)
dmean_gtvn = nan(length(mrn_list), 6);  % mean dose in GTVn (Gy)

d95_gtvp = nan(length(mrn_list), 6);    % D95% — dose received by 95% of GTV (Gy)
d95_gtvn = nan(length(mrn_list), 6);

v50gy_gtvp = nan(length(mrn_list), 6);  % V50Gy — fraction of GTV receiving ≥50 Gy
v50gy_gtvn = nan(length(mrn_list), 6);

% Clinical outcome flags (per patient)
lf = zeros(size(mrn_list));      % local failure (1 = yes)
immuno = zeros(size(mrn_list));  % received immunotherapy (1 = yes)

% Track problematic DWI acquisitions for manual review
bad_dwi_locations = cell(1);
bad_dwi_count = 0;

h = waitbar(0,'loading dwi data');

% --- Main processing loop: iterate over patients ---
for j = 1:length(mrn_list)
    mrn = mrn_list{j};
    fprintf('\n******* MRN: %s\n',mrn);

    % Per-patient DIR reference: populated at Fx1, reused at Fx2+
    b0_fx1_ref        = [];   % b=0 volume at baseline fraction
    gtv_mask_fx1_ref  = [];   % GTVp mask at baseline fraction
    gtvn_mask_fx1_ref = [];   % GTVn mask at baseline fraction (when present)

    % read clinical data (local failure, immunotherapy) from spreadsheet
    i_pat = find(contains(strrep(T.Pat,'_','-'),strrep(id_list{j},'_','-')));
    immuno(j) = T.Immuno(i_pat(1));
    lf(j) = T.LF(i_pat(1));

    % --- Loop over fractions (fi) and repeat acquisitions (rpi) ---
    for fi=1:size(dwi_locations,2)
        for rpi = 1:size(dwi_locations,3)

            basefolder = [dataloc id_list{j} '/'];
            fxtmp = clean_dir_command([basefolder '*' fx_search{fi} '*']);

            if isempty(fxtmp)
                fprintf('%s, no %s folder\n',id_list{j},fx_search{fi});
            else
                fxfolder = [basefolder fxtmp(1).name '/'];
            end

            % Retrieve previously discovered file paths for this combination
            dicomloc = dwi_locations{j,fi,rpi};
            struct_file = gtv_locations{j,fi,rpi};
            struct_file_gtvn = gtvn_locations{j,fi,rpi};
            if fi<=size(rtdose_locations,2)
                dicomdoseloc = rtdose_locations{j,fi};  % RT dose only for Fx1–Fx5
            else
                dicomdoseloc = [];  % no dose for post-treatment scan
            end

            outloc = [basefolder '/nii/'];   % output directory for NIfTI files

            % get the dwi file list (DICOM)
            dicom_files = dir([dicomloc '*.dcm']);
            % keep track of dwi issues
            bad_dwi_found = 0;

            % Build standardised naming IDs for this scan
            if fi<=size(rtdose_locations,2)
                fx_id = ['fx' int2str(fi)];
            else
                fx_id = 'post';
            end

            scanID = [fx_id '_dwi' int2str(rpi)]; %'fx1_dwi1';
            gtvname = [fx_id '_gtv' int2str(rpi)]; %'fx1_gtv1';
            gtvn_name = [fx_id '_gtvn' int2str(rpi)]; %'fx1_gtvn1';
            dosename = [fx_id '_dose_on_dwi' int2str(rpi)]; %'fx1_dose_on_dwi1';

            if ~isfolder(outloc), mkdir(outloc); end

            % --- Convert DWI DICOMs to NIfTI using dcm2niix ---
            % save dwi as nii.gz
            if ~isempty(dicomloc)
                if ~exist([outloc scanID '.nii.gz'],'file')
                    % keep track of number of files generated by this
                    % command (make sure exactly 3 are generated:
                    % .nii.gz, .bval, .bvec, plus the dir entry = 4 delta)
                    nfiles_before = length(dir(outloc));
                    nii_cmd = ['"' dcm2nii_call '" -z y -f ' scanID ' -o "' outloc '" "' dicomloc '"'];
                    system(nii_cmd);
                    nfiles_after = length(dir(outloc));
                    if (nfiles_after - nfiles_before) ~=4
                        fprintf('!!!! incorrect number of DWI files generated for %s... need to fix!\n',fx_id);
                        bad_dwi_count = bad_dwi_count+1;
                        bad_dwi_locations{bad_dwi_count} = dicomloc;
                        bad_dwi_found = 1;
                    end
                end
            end

            % --- Save GTVp mask as NIfTI for consistency with DWI volumes ---
            % save the mask as nii.gz (just for consistency)
            % Stvol3d is the 3-D binary mask variable stored in the .mat file
            if ~isempty(struct_file)
                if ~exist([outloc gtvname '.nii.gz'],'file')
                    load(struct_file);
                    gtv_mask = Stvol3d;
                    niftiwrite(rot90(double(gtv_mask),-1),[outloc gtvname],'Compressed',true);
                end
            end

            % --- Save GTVn (nodal) mask as NIfTI, if present ---
            if ~isempty(struct_file_gtvn)
                if ~exist([outloc gtvn_name '.nii.gz'],'file')
                    load(struct_file_gtvn);
                    gtvn_mask = Stvol3d;
                    niftiwrite(rot90(double(gtvn_mask),-1),[outloc gtvn_name],'Compressed',true);
                end
            end

            % --- Resample RT dose onto DWI geometry and save as NIfTI ---
            % sample RTdose on the DWI geometry, save as nii.gz
            % but first, see if we actually have a dose map to load....
            if ~isempty(dicomdoseloc) && ~isempty(dicomloc)
                if ~exist([outloc dosename '.nii.gz'],'file')
                    % extract only b=0 images — these define the DWI spatial
                    % geometry used as the resampling target for the dose grid
                    % extract only b0 images for the dwi sampling
                    isb0 = zeros(length(dicom_files),1);
                    b0list = cell(1);
                    b0count = 0;
                    for bi = 1:length(isb0)
                        data_tmp = dicominfo([dicom_files(bi).folder '/' dicom_files(bi).name], 'UseDictionaryVR', true);
                        if data_tmp.DiffusionBValue == 0
                            isb0(bi) = 1;
                            b0count = b0count+1;
                            b0list{b0count,1} = [dicom_files(bi).folder '/' dicom_files(bi).name];
                        end
                    end

                    rtdose_dicom = dir([dicomdoseloc '*.dcm']);
                    rtdosefile = [rtdose_dicom.folder '/' rtdose_dicom.name];

                    % Resample RT dose onto b=0 image grid and write NIfTI
                    dose_sampled = sample_rtdose_on_image(b0list,rtdosefile);
                    niftiwrite(rot90(dose_sampled,-1),[outloc dosename],'Compressed',true);
                end
            end

            % --- Load NIfTI DWI volume and extract b-values ---
            % load the data
            havedwi = 0;
            if exist([outloc scanID '.nii.gz'],'file')
                dwi_dat = load_untouch_nii([outloc scanID '.nii.gz']);
                dwi = rot90(dwi_dat.img);            % orient to standard view
                dwi_dims = dwi_dat.hdr.dime.pixdim(2:4);  % voxel dimensions in mm
                dwi_vox_vol = prod(dwi_dims*0.1);    % voxel volume in cc (mm→cm)
                fprintf('...Loaded %s. ',[outloc scanID '.nii.gz']);
                havedwi = 1;
                % Read b-values from the sidecar .bval file produced by dcm2niix
                % now extract bvalues
                bval_file = [outloc scanID '.bval'];
                if exist(bval_file,'file')
                    fid = fopen(bval_file);
                    tline = fgetl(fid);
                    bval_data = tline;
                    fclose(fid);

                    bvalues = strsplit(bval_data,' ');
                    bvalues = cellfun(@str2num,bvalues');

                    % Sort b-values ascending and reorder the 4-D DWI volume
                    [b_sort,i_sort] = sort(bvalues,'ascend');

                    bvalues = bvalues(i_sort);
                    dwi = double(dwi(:,:,:,i_sort));
                    fprintf('loaded bvalues\n');
                else
                    fprintf('bvalue file not found!\n');
                    havedwi = 0;
                    if bad_dwi_found==0
                        bad_dwi_count = bad_dwi_count+1;
                        bad_dwi_locations{bad_dwi_count} = dicomloc;
                        bad_dwi_found = 1;
                    end
                end
                % Expect exactly 4 b-values (e.g., 0, 50, 400, 800 s/mm²)
                if size(dwi,4)~=4
                    fprintf('DWI does not have expected dimensions found: %s skipping\n',mat2str(size(dwi)))
                    havedwi = 0;
                    if bad_dwi_found==0
                        bad_dwi_count = bad_dwi_count+1;
                        bad_dwi_locations{bad_dwi_count} = dicomloc;
                        bad_dwi_found = 1;
                    end
                end
            end

            % --- Load DnCNN-denoised DWI (deep learning denoising) ---
            havedenoised = 0;
            if havedwi==1
                dncnnid = [scanID '_dncnn.nii.gz'];
                dncnn_file = [basefolder '/dncnn/' dncnnid];
                if exist(dncnn_file,'file')
                    dwi_dat_dncnn = load_untouch_nii(dncnn_file);
                    dwi_dncnn = rot90(dwi_dat_dncnn.img);
                    % Normalise denoised signal to [0,1] for IVIM fitting
                    dwi_dncnn = double(mat2gray(dwi_dncnn(:,:,:,i_sort)));
                    havedenoised=1;
                end
            end

            % --- Load IVIMnet deep-learning fit results (pre-computed) ---
            haveivimnet = 0;
            if havedwi==1
                ivimid = [scanID '_ivimnet.mat'];
                ivimnet_file = [basefolder '/ivimnet/' ivimid];
                if exist(ivimnet_file,'file')
                    % Loads D_ivimnet, f_ivimnet, Dstar_ivimnet, S0_ivimnet
                    load(ivimnet_file);
                    haveivimnet=1;
                end
            end

            % --- Load GTVp mask and validate spatial dimensions ---
            havegtvp = 0;
            if exist([outloc gtvname '.nii.gz'],'file')
                gtv_dat = load_untouch_nii([outloc gtvname '.nii.gz']);
                gtv_mask = rot90(gtv_dat.img);
                fprintf('...Loaded %s\n',[outloc gtvname '.nii.gz']);
                havegtvp = 1;

                % Ensure GTV mask dimensions match the DWI spatial dims
                gtv_size = size(gtv_mask);
                dwi_size = size(dwi);
                if sum(gtv_size ~= dwi_size(1:3))>0
                    havegtvp = 0;
                    fprintf('size mismatch. excluding gtvp\n');
                end
            end

            % --- Load GTVn (nodal) mask and validate dimensions ---
            havegtvn = 0;
            if exist([outloc gtvn_name '.nii.gz'],'file')
                gtvn_dat = load_untouch_nii([outloc gtvn_name '.nii.gz']);
                gtvn_mask = rot90(gtvn_dat.img);
                fprintf('... *NODAL* Loaded %s\n',[outloc gtvn_name '.nii.gz']);
                havegtvn = 1;

                gtv_size = size(gtvn_mask);
                dwi_size = size(dwi);
                if sum(gtv_size ~= dwi_size(1:3))>0
                    havegtvn = 0;
                    fprintf('size mismatch. excluding gtvn\n');
                end
            end

            % --- Load resampled RT dose map ---
            havedose = 0;
            if exist([outloc dosename '.nii.gz'],'file')
                dose_dat = load_untouch_nii([outloc dosename '.nii.gz']);
                dose_map = rot90(dose_dat.img);
                fprintf('...Loaded %s\n',[outloc dosename '.nii.gz']);
                havedose = 1;
            end

            % --- Fit ADC (monoexponential) and IVIM (biexponential) models ---
            % fit adc and IVIM models
            if havedwi && (havegtvn || havegtvp)
                % only fit IVIM model within mask(s) to save time
                % (avoids fitting thousands of background voxels)
                mask_ivim = false(size(dwi,1),size(dwi,2),size(dwi,3));
                if havegtvp, mask_ivim = logical(mask_ivim + logical(gtv_mask)); end
                if havegtvn, mask_ivim = logical(mask_ivim + logical(gtvn_mask)); end

                % Segmented IVIM fit: b < bthr separates perfusion from diffusion.
                % b >= bthr used to estimate D (true tissue diffusion),
                % b < bthr used to estimate f and D* (pseudo-diffusion).
                opts = [];
                opts.bthr = ivim_bthr; % set in USER OPTIONS (default 200 s/mm2; clinical standard for abdominal IVIM)

                sz3 = [size(dwi,1), size(dwi,2), size(dwi,3)];
                has_sufficient_bvalues = sum(bvalues >= opts.bthr) >= 2;
                if has_sufficient_bvalues
                    ivim_fit = IVIMmodelfit(dwi,bvalues,"seg",mask_ivim,opts);
                    % 4th dimension of ivim_fit is: adc, s0, f, dstar
                    d_map = ivim_fit(:,:,:,1);      % D — true diffusion coefficient
                    f_map = ivim_fit(:,:,:,3);      % f — perfusion fraction
                    dstar_map = ivim_fit(:,:,:,4);  % D* — pseudo-diffusion coefficient

                    % Replace zero-fit voxels with NaN (failed fits)
                    zero_mask = d_map==0;
                    d_map(zero_mask) = nan;
                    f_map(zero_mask) = nan;
                    dstar_map(zero_mask) = nan;
                else
                    fprintf('Insufficient b-values >= %d for IVIM fit; skipping IVIM (maps set to NaN)\n', opts.bthr);
                    d_map = nan(sz3);
                    f_map = nan(sz3);
                    dstar_map = nan(sz3);
                end

                % Monoexponential ADC fit — log-linear OLS on active voxels only.
                % Pre-filtering to voxels with all-positive signal prevents the
                % log-linear operation from freezing on large zero-signal backgrounds.
                adc_sz  = [size(dwi,1), size(dwi,2), size(dwi,3)];
                S_2d    = reshape(double(dwi), [prod(adc_sz), length(bvalues)]);
                adc_idx = find(all(S_2d > 0, 2));
                adc_vec = zeros(prod(adc_sz), 1);
                if ~isempty(adc_idx)
                    S_a = S_2d(adc_idx, :);
                    adc_vec(adc_idx) = (-bvalues(2:end) \ ...
                        permute(log(S_a(:,2:end) ./ S_a(:,1)), [2 1]))';
                    adc_vec(adc_vec < 0) = nan;  % clamp noise-driven negative ADC estimates
                end
                adc_map = reshape(adc_vec, adc_sz);

                % Repeat IVIM + ADC fitting on DnCNN-denoised data
                if havedenoised==1
                    if has_sufficient_bvalues
                        ivim_fit_dncnn = IVIMmodelfit(dwi_dncnn,bvalues,"seg",mask_ivim,opts);
                        % 4th dimension of ivim_fit is: adc, s0, f, dstar
                        d_map_dncnn = ivim_fit_dncnn(:,:,:,1);
                        f_map_dncnn = ivim_fit_dncnn(:,:,:,3);
                        dstar_map_dncnn = ivim_fit_dncnn(:,:,:,4);

                        zero_mask_dncnn = d_map_dncnn==0;
                        d_map_dncnn(zero_mask_dncnn) = nan;
                        f_map_dncnn(zero_mask_dncnn) = nan;
                        dstar_map_dncnn(zero_mask_dncnn) = nan;
                    else
                        d_map_dncnn = nan(sz3);
                        f_map_dncnn = nan(sz3);
                        dstar_map_dncnn = nan(sz3);
                    end

                    % Monoexponential ADC fit on DnCNN-denoised data (pre-filtered)
                    adc_sz_d   = [size(dwi_dncnn,1), size(dwi_dncnn,2), size(dwi_dncnn,3)];
                    S_2d_d     = reshape(double(dwi_dncnn), [prod(adc_sz_d), length(bvalues)]);
                    adc_idx_d  = find(all(S_2d_d > 0, 2));
                    adc_vec_d  = zeros(prod(adc_sz_d), 1);
                    if ~isempty(adc_idx_d)
                        S_a_d = S_2d_d(adc_idx_d, :);
                        adc_vec_d(adc_idx_d) = (-bvalues(2:end) \ ...
                            permute(log(S_a_d(:,2:end) ./ S_a_d(:,1)), [2 1]))';
                        adc_vec_d(adc_vec_d < 0) = nan;  % drop noise-driven negative ADC estimates (failed fit)
                    end
                    adc_map_dncnn = reshape(adc_vec_d, adc_sz_d);
                end
            end

            % --- Deformable Image Registration (DIR): propagate GTV mask and ---
            % --- compute displacement field for inter-fraction alignment       ---
            % Performed BEFORE biomarker extraction so that the DIR inverse field
            % is available to warp native-space DnCNN parameter maps into the
            % baseline geometry prior to radiomic feature extraction.
            % The Demons displacement field is cached to disk to avoid redundant
            % computation on re-runs.
            gtv_mask_for_dvh = gtv_mask;  % default: raw mask (updated below for fi>1)
            dose_map_dvh     = dose_map;  % default: rigid dose (updated below for fi>1)
            D_forward_cur    = [];        % Demons field for this fraction
            ref3d_cur        = [];        % imref3d object for this fraction

            if havedwi && havegtvp
                b0_current = dwi(:,:,:,1); % b=0 is always the first volume
                if fi == 1
                    % Cache the Fx1 references for downstream DIR calls
                    b0_fx1_ref        = b0_current;
                    gtv_mask_fx1_ref  = gtv_mask;
                    if havegtvn
                        gtvn_mask_fx1_ref = gtvn_mask;
                    end
                elseif fi > 1 && ~isempty(b0_fx1_ref) && ~isempty(gtv_mask_fx1_ref)
                    dir_cache_file = fullfile([dataloc id_list{j} '/nii/'], ...
                        sprintf('dir_field_rpi%d_fx%d.mat', rpi, fi));
                    if exist(dir_cache_file, 'file')
                        tmp_dir = load(dir_cache_file, 'gtv_mask_warped', 'D_forward', 'ref3d');
                        gtv_mask_for_dvh = tmp_dir.gtv_mask_warped;
                        if isfield(tmp_dir, 'D_forward'),  D_forward_cur = tmp_dir.D_forward;  end
                        if isfield(tmp_dir, 'ref3d'),      ref3d_cur     = tmp_dir.ref3d;      end
                        fprintf('  [DIR] Loaded cached warped mask + D_forward for Fx%d rpi%d\n', fi, rpi);
                    else
                        fprintf('  [DIR] Running imregdemons for Fx%d rpi%d...\n', fi, rpi);
                        [gtv_mask_warped, D_forward_cur, ref3d_cur] = ...
                            apply_dir_mask_propagation(b0_fx1_ref, b0_current, gtv_mask_fx1_ref);
                        if ~isempty(gtv_mask_warped)
                            gtv_mask_for_dvh = gtv_mask_warped;
                            % Cache mask + displacement field to disk
                            D_forward = D_forward_cur;  ref3d = ref3d_cur; %#ok<NASGU>
                            save(dir_cache_file, 'gtv_mask_warped', 'D_forward', 'ref3d');
                            fprintf('  [DIR] Done. Warped mask + D_forward saved.\n');
                        else
                            fprintf('  [DIR] Registration failed for Fx%d rpi%d. Falling back to rigid dose/mask.\n', fi, rpi);
                        end
                    end

                    % --- DO NOT warp the dose map ---
                    % The dose map remains rigidly aligned. Sub-volume dose metrics (V50, D95)
                    % are computed by sampling this static beam grid against the deformed daily
                    % anatomy (gtv_mask_for_dvh).
                end
            end

            % --- Warp native-space DnCNN parameter maps to baseline geometry ---
            % The DnCNN network is applied to the raw, unwarped fractional MRI to
            % preserve the physical Rician noise distribution. After fitting IVIM
            % on the native-space denoised signal, the resulting parameter maps are
            % warped into baseline geometry using the DIR inverse field (-D_forward)
            % so that radiomic features are extracted in a consistent anatomical
            % reference frame. -D_forward is the first-order approximation of the
            % inverse displacement field (current-fraction to baseline) for the
            % symmetric Demons registration used in apply_dir_mask_propagation.
            if havedenoised && fi > 1 && ~isempty(D_forward_cur) && ~isempty(b0_fx1_ref)
                ref3d_bl = imref3d(size(b0_fx1_ref));
                d_map_dncnn     = imwarp(d_map_dncnn,     -D_forward_cur, 'Interp', 'linear', ...
                    'OutputView', ref3d_bl, 'FillValues', nan);
                f_map_dncnn     = imwarp(f_map_dncnn,     -D_forward_cur, 'Interp', 'linear', ...
                    'OutputView', ref3d_bl, 'FillValues', nan);
                dstar_map_dncnn = imwarp(dstar_map_dncnn, -D_forward_cur, 'Interp', 'linear', ...
                    'OutputView', ref3d_bl, 'FillValues', nan);
                adc_map_dncnn   = imwarp(adc_map_dncnn,   -D_forward_cur, 'Interp', 'linear', ...
                    'OutputView', ref3d_bl, 'FillValues', nan);
                fprintf('  [DnCNN] Warped native-space parameter maps to baseline geometry.\n');
            end

            % --- Extract voxel-level biomarkers within GTVp ---
            if havegtvp
                % Compute summary statistics within the primary GTV
                adc_mean(j,fi,rpi) = nanmean(adc_map(gtv_mask==1));
                % NOTE: Histogram kurtosis of trace-average ADC — NOT valid DKI.
                adc_kurtosis(j,fi,rpi) = kurtosis(adc_map(gtv_mask==1));

                d_mean(j,fi,rpi) = nanmean(d_map(gtv_mask==1));
                % NOTE: Histogram kurtosis of trace-average map — NOT valid DKI.
                d_kurtosis(j,fi,rpi) = kurtosis(d_map(gtv_mask==1));

                % Store voxel-level vectors and metadata in the output struct
                data_vectors_gtvp(j,fi,rpi).adc_vector = adc_map(gtv_mask==1);
                data_vectors_gtvp(j,fi,rpi).d_vector = d_map(gtv_mask==1);
                data_vectors_gtvp(j,fi,rpi).f_vector = f_map(gtv_mask==1);
                data_vectors_gtvp(j,fi,rpi).dstar_vector = dstar_map(gtv_mask==1);    
                data_vectors_gtvp(j,fi,rpi).ID = id_list{j};
                data_vectors_gtvp(j,fi,rpi).MRN = mrn_list{j};
                data_vectors_gtvp(j,fi,rpi).LF = lf(j);
                data_vectors_gtvp(j,fi,rpi).Immuno = immuno(j);
                data_vectors_gtvp(j,fi,rpi).Fraction = fi;
                data_vectors_gtvp(j,fi,rpi).Repeatability_index = rpi;
                data_vectors_gtvp(j,fi,rpi).vox_vol = dwi_vox_vol;

                % Store DnCNN-denoised biomarker vectors (pipeline variant 2).
                % For fi > 1: the parameter maps have been warped to baseline
                % geometry; extract within the baseline GTV mask for spatial
                % consistency. For fi == 1: maps are already in native/baseline
                % space; use the native mask directly.
                if havedenoised
                    dncnn_mask_p = gtv_mask;
                    if fi > 1 && ~isempty(gtv_mask_fx1_ref)
                        dncnn_mask_p = gtv_mask_fx1_ref;
                    end
                    d_mean_dncnn(j,fi,rpi) = nanmean(d_map_dncnn(dncnn_mask_p==1));
                    data_vectors_gtvp(j,fi,rpi).d_vector_dncnn    = d_map_dncnn(dncnn_mask_p==1);
                    data_vectors_gtvp(j,fi,rpi).f_vector_dncnn    = f_map_dncnn(dncnn_mask_p==1);
                    data_vectors_gtvp(j,fi,rpi).dstar_vector_dncnn = dstar_map_dncnn(dncnn_mask_p==1);
                    data_vectors_gtvp(j,fi,rpi).adc_vector_dncnn  = adc_map_dncnn(dncnn_mask_p==1);
                end

                % Store IVIMnet deep-learning fit vectors (pipeline variant 3)
                if haveivimnet
%                     mdic = {"D_ivimnet": D_out, "f_ivimnet":f_out, "Dstar_ivimnet":Dstar_out, "S0_ivimnet":S0_out}
                    D_ivimnet = rot90(D_ivimnet);
                    f_ivimnet = rot90(f_ivimnet);
                    Dstar_ivimnet = rot90(Dstar_ivimnet);
                    S0_ivimnet = rot90(S0_ivimnet);
                    d_mean_ivimnet(j,fi,rpi) = nanmean(D_ivimnet(gtv_mask==1));

                    data_vectors_gtvp(j,fi,rpi).d_vector_ivimnet = D_ivimnet(gtv_mask==1);
                    data_vectors_gtvp(j,fi,rpi).f_vector_ivimnet = f_ivimnet(gtv_mask==1);
                    data_vectors_gtvp(j,fi,rpi).dstar_vector_ivimnet = Dstar_ivimnet(gtv_mask==1);
                end
            end

            if havedose && havegtvp
                dmean_gtvp(j,fi) = nanmean(dose_map_dvh(gtv_mask_for_dvh==1));
                dwi_dims = dwi_dat.hdr.dime.pixdim(2:4);
                % DVH uses the strictly rigid dose against the DIR-warped daily GTV tissue mask
                % so that dose correctly reflects the true static beam delivered to deformed anatomy.
                [dvhparams, dvh_values] = dvh(dose_map_dvh, gtv_mask_for_dvh, dwi_dims, 2000, 'Dperc',95,'Vperc',50,'Normalize',true);
                d95_gtvp(j,fi) = dvhparams.("D95% (Gy)");
                v50gy_gtvp(j,fi) = dvhparams.("V50Gy (%)");

                data_vectors_gtvp(j,fi,rpi).dose_vector = dose_map_dvh(gtv_mask_for_dvh==1);
                data_vectors_gtvp(j,fi,rpi).dvh = dvh_values;
                data_vectors_gtvp(j,fi,rpi).d95 = dvhparams.("D95% (Gy)");
                data_vectors_gtvp(j,fi,rpi).v50gy = dvhparams.("V50Gy (%)");
            end

            % --- Extract voxel-level biomarkers within GTVn (nodal GTV) ---
            % now collect gtvn info
            if havegtvn
                data_vectors_gtvn(j,fi,rpi).adc_vector = adc_map(gtvn_mask==1);
                data_vectors_gtvn(j,fi,rpi).d_vector = d_map(gtvn_mask==1);
                data_vectors_gtvn(j,fi,rpi).f_vector = f_map(gtvn_mask==1);
                data_vectors_gtvn(j,fi,rpi).dstar_vector = dstar_map(gtvn_mask==1);
                data_vectors_gtvn(j,fi,rpi).ID = id_list{j};
                data_vectors_gtvn(j,fi,rpi).MRN = mrn_list{j};
                data_vectors_gtvn(j,fi,rpi).LF = lf(j);
                data_vectors_gtvn(j,fi,rpi).Immuno = immuno(j);
                data_vectors_gtvn(j,fi,rpi).Fraction = fi;
                data_vectors_gtvn(j,fi,rpi).Repeatability_index = rpi;
                data_vectors_gtvn(j,fi,rpi).vox_vol = dwi_vox_vol;
            end

            % Store DnCNN-denoised vectors for GTVn.
            % For fi > 1 the parameter maps have been warped to baseline geometry;
            % extract within the baseline GTVn mask for spatial consistency.
            if havedenoised && havegtvn
                dncnn_mask_n = gtvn_mask;
                if fi > 1 && ~isempty(gtvn_mask_fx1_ref)
                    dncnn_mask_n = gtvn_mask_fx1_ref;
                end
                data_vectors_gtvn(j,fi,rpi).d_vector_dncnn    = d_map_dncnn(dncnn_mask_n==1);
                data_vectors_gtvn(j,fi,rpi).f_vector_dncnn    = f_map_dncnn(dncnn_mask_n==1);
                data_vectors_gtvn(j,fi,rpi).dstar_vector_dncnn = dstar_map_dncnn(dncnn_mask_n==1);
                data_vectors_gtvn(j,fi,rpi).adc_vector_dncnn  = adc_map_dncnn(dncnn_mask_n==1);
            end

            % Store IVIMnet vectors for GTVn
            if haveivimnet && havegtvn
                %                     mdic = {"D_ivimnet": D_out, "f_ivimnet":f_out, "Dstar_ivimnet":Dstar_out, "S0_ivimnet":S0_out}
                D_ivimnet = rot90(D_ivimnet);
                f_ivimnet = rot90(f_ivimnet);
                Dstar_ivimnet = rot90(Dstar_ivimnet);
                S0_ivimnet = rot90(S0_ivimnet);

                data_vectors_gtvn(j,fi,rpi).d_vector_ivimnet = D_ivimnet(gtvn_mask==1);
                data_vectors_gtvn(j,fi,rpi).f_vector_ivimnet = f_ivimnet(gtvn_mask==1);
                data_vectors_gtvn(j,fi,rpi).dstar_vector_ivimnet = Dstar_ivimnet(gtvn_mask==1);
            end

            % Compute DVH parameters within GTVn
            % Sample rigidly aligned dose map against the GTVn mask.
            if havedose && havegtvn
                dose_map_dvh_n = dose_map;  % rigidly aligned dose
                dmean_gtvn(j,fi) = nanmean(dose_map_dvh_n(gtvn_mask==1));
                dwi_dims = dwi_dat.hdr.dime.pixdim(2:4);
                [dvhparams, dvh_values] = dvh(dose_map_dvh_n, gtvn_mask, dwi_dims, 2000, 'Dperc',95,'Vperc',50,'Normalize',true);
                d95_gtvn(j,fi) = dvhparams.("D95% (Gy)");
                v50gy_gtvn(j,fi) = dvhparams.("V50Gy (%)");

                data_vectors_gtvn(j,fi,rpi).dose_vector = dose_map_dvh_n(gtvn_mask==1);
                data_vectors_gtvn(j,fi,rpi).dvh = dvh_values;
                data_vectors_gtvn(j,fi,rpi).d95 = dvhparams.("D95% (Gy)");
                data_vectors_gtvn(j,fi,rpi).v50gy = dvhparams.("V50Gy (%)");
            end
        end
    end
    h = waitbar(j/length(mrn_list));  % update progress bar
end

delete(h);  % close waitbar

%% ========================================================================
fprintf('\n--- SECTION 3: Save Results ---\n');
%  SECTION 3 — SAVE RESULTS

datasave = fullfile(dataloc, 'dwi_vectors.mat');
% Create a date-stamped backup before overwriting
if exist(datasave,'file')
    dt = datetime('now');
    dateString = char(dt, 'yyyy_MMM_dd');
    newfilename = cat(2,dataloc, 'dwi_vectors_', dateString, '.mat');
    copyfile(datasave,newfilename);
    fprintf('backed up existing save to %s\n',newfilename);
end
save(datasave,'data_vectors_gtvn','data_vectors_gtvp','lf','immuno','mrn_list','id_list','fx_dates','dwi_locations','rtdose_locations','gtv_locations','gtvn_locations','dmean_gtvp','dmean_gtvn','d95_gtvp','d95_gtvn','v50gy_gtvp','v50gy_gtvn','bad_dwi_locations','bad_dwi_count');
fprintf('saved %s\n',datasave);

end % if ~skip_to_reload

%% ========================================================================
fprintf('\n--- SECTION 4: Reload Saved Data ---\n');
%  SECTION 4 — RELOAD SAVED DATA

% Set data path from configuration
dataloc = config_struct.dataloc;

datasave = fullfile(dataloc, 'dwi_vectors.mat');
load(datasave);

%% ========================================================================
fprintf('\n--- SECTION 5: Longitudinal Summary Metrics ---\n');
%  SECTION 5 — LONGITUDINAL SUMMARY METRICS

% ADC threshold for identifying "restricted diffusion" sub-volume
% (1.15×10⁻³ mm²/s, per Muraoka et al. 2013)
adc_thresh = config_struct.adc_thresh; % https://pubmed.ncbi.nlm.nih.gov/23545001/

% Secondary ADC threshold for identifying "high ADC" sub-volume
high_adc_thresh = config_struct.high_adc_thresh;

% Minimum voxel threshold for higher-order histogram metrics
% (kurtosis, skewness, KS test). Returns NaN for smaller volumes to
% prevent unstable estimates.
min_vox_hist = config_struct.min_vox_hist;

nTp = 6;   % number of timepoints (Fx1–Fx5 + post)
nRpt = 6;  % max number of repeat scans at Fx1

% ADC: last entry: 1 - normal, 2: dncnn + ivim fit
% --- Pre-allocate sub-volume metric arrays (patient × timepoint × pipeline) ---
% "sub" = voxels with ADC below adc_thresh (restricted diffusion sub-volume)
adc_sub_vol_pc = nan(length(id_list),nTp,2);  % sub-volume as fraction of GTV
adc_sub_vol = nan(length(id_list),nTp,2);     % sub-volume in cc
adc_sub_mean = nan(length(id_list),nTp,2);
adc_sub_kurt = nan(length(id_list),nTp,2);
adc_sub_skew = nan(length(id_list),nTp,2);

high_adc_sub_vol = nan(length(id_list),nTp,2); % volume of voxels above high_adc_thresh

ivim_sub_vol = nan(length(id_list),nTp,2); % voxels with D<0.001 AND f<0.1

gtv_vol = nan(length(id_list),nTp);  % total GTV volume in cc

% --- Whole-GTV summary statistics for ADC ---
adc_mean = nan(length(id_list),nTp,2);
adc_kurt = nan(length(id_list),nTp,2);
adc_skew = nan(length(id_list),nTp,2);
adc_sd = nan(length(id_list),nTp,2);

% IVIM: last entry: 1 - normal, 2: dncnn + ivim fit, 3: ivimnet fit
% --- Whole-GTV summary statistics for IVIM parameters (D, f, D*) ---
d_mean = nan(length(id_list),nTp,3);
d_kurt = nan(length(id_list),nTp,3);
d_skew = nan(length(id_list),nTp,3);
d_sd = nan(length(id_list),nTp,3);

% D sub-volume statistics (restricted sub-volume only)
d_sub_mean = nan(length(id_list),nTp,3);
d_sub_kurt = nan(length(id_list),nTp,3);
d_sub_skew = nan(length(id_list),nTp,3);

% Perfusion fraction (f) summary statistics
f_mean = nan(length(id_list),nTp,3);
f_kurt = nan(length(id_list),nTp,3);
f_skew = nan(length(id_list),nTp,3);

% Pseudo-diffusion coefficient (D*) summary statistics
dstar_mean = nan(length(id_list),nTp,3);
dstar_kurt = nan(length(id_list),nTp,3);
dstar_skew = nan(length(id_list),nTp,3);

% --- Histogram and KS-test arrays for longitudinal distribution comparison ---
bin_edges = 0:0.5e-4:3e-3; % to compute histograms and histogram distances
adc_histograms = nan(length(id_list),nTp,length(bin_edges)-1,3);  % normalised ADC histograms
d_histograms = nan(length(id_list),nTp,length(bin_edges)-1,3);    % normalised D histograms
ks_stats_adc = nan(length(id_list),nTp,3);  % KS statistic: ADC at timepoint k vs baseline
ks_pvals_adc = nan(length(id_list),nTp,3);
ks_stats_d = nan(length(id_list),nTp,3);    % KS statistic: D at timepoint k vs baseline
ks_pvals_d = nan(length(id_list),nTp,3);


% --- Motion corruption flag ---
% ADC: last entry: 1 - normal, 2: dncnn + ivim fit
% Fraction of GTV voxels exceeding adc_max — high values suggest motion artefact
fx_corrupted = nan(length(id_list),nTp,2); % count voxels with adc>3e-3 (indicative of motion corruption)
adc_max = config_struct.adc_max;  % corruption threshold (mm²/s)

% --- Repeatability arrays (Fx1 repeat scans only, for wCV calculation) ---
adc_mean_rpt = nan(length(id_list),nRpt,2);
adc_sub_rpt = nan(length(id_list),nRpt,2);

fx_corrupted_rpt = nan(length(id_list),nRpt,2);

d_mean_rpt = nan(length(id_list),nRpt,3);
f_mean_rpt = nan(length(id_list),nRpt,3);
dstar_mean_rpt = nan(length(id_list),nRpt,3);

% --- Pooled voxel vectors across all patients (for population-level analysis) ---
adc_vec_all = [];
d_vec_all = [];
f_vec_all = [];
fx_vec_all = [];     % fraction index for each pooled voxel

% Pooled vectors stratified by local failure (lf) vs complete response (cr)
adc_vec_lf = [];
d_vec_lf= [];
f_vec_lf = [];
fx_vec_lf = [];

adc_vec_cr = [];
d_vec_cr = [];
f_vec_cr = [];
fx_vec_cr = [];

n_rpt = nan(length(id_list),1);  % number of valid repeat scans per patient

% --- Main analysis loop: patient × timepoint × DWI pipeline ---
for j=1:length(id_list)
    for k=1:nTp
        vox_vol = data_vectors_gtvp(j,k,1).vox_vol;
        for dwi_type=1:3 % 1- normal, 2-dncnn, 3-ivimnet

            % Select the appropriate voxel vectors depending on pipeline
            switch dwi_type
                case 1  % Standard (raw DWI)
                    adc_vec = data_vectors_gtvp(j,k,1).adc_vector;
                    d_vec = data_vectors_gtvp(j,k,1).d_vector;
                    f_vec = data_vectors_gtvp(j,k,1).f_vector;
                    dstar_vec = data_vectors_gtvp(j,k,1).dstar_vector;

                    adc_baseline = data_vectors_gtvp(j,1,1).adc_vector;  % Fx1 as baseline
                    d_baseline = data_vectors_gtvp(j,1,1).d_vector;
                case 2  % DnCNN-denoised + conventional IVIM fit
                    adc_vec = data_vectors_gtvp(j,k,1).adc_vector_dncnn;
                    d_vec = data_vectors_gtvp(j,k,1).d_vector_dncnn;
                    f_vec = data_vectors_gtvp(j,k,1).f_vector_dncnn;
                    dstar_vec = data_vectors_gtvp(j,k,1).dstar_vector_dncnn;

                    adc_baseline = data_vectors_gtvp(j,1,1).adc_vector_dncnn;
                    d_baseline = data_vectors_gtvp(j,1,1).d_vector_dncnn;
                case 3  % IVIMnet deep-learning IVIM fit (ADC uses standard pipeline)
                    adc_vec = data_vectors_gtvp(j,k,1).adc_vector;
                    d_vec = data_vectors_gtvp(j,k,1).d_vector_ivimnet;
                    f_vec = data_vectors_gtvp(j,k,1).f_vector_ivimnet;
                    dstar_vec = data_vectors_gtvp(j,k,1).dstar_vector_ivimnet;

                    adc_baseline = data_vectors_gtvp(j,1,1).adc_vector;
                    d_baseline = data_vectors_gtvp(j,1,1).d_vector_ivimnet;
            end

            % --- Compute ADC summary metrics for this patient/timepoint ---
            if ~isempty(adc_vec)
%                 adc_vec = data_vectors_gtvp(j,k,1).adc_vector;
                gtv_vol(j,k) = numel(adc_vec)*vox_vol;   % GTV volume (cc)
                adc_mean(j,k,dwi_type) = nanmean(adc_vec);
                if numel(adc_vec) >= min_vox_hist
                    % NOTE: Histogram kurtosis of trace-average ADC — NOT valid DKI.
                    % Retained for archival completeness; do not use in primary feature pool.
                    adc_kurt(j,k,dwi_type) = kurtosis(adc_vec);
                    adc_skew(j,k,dwi_type) = skewness(adc_vec);
                end
                adc_sd(j,k,dwi_type) = nanstd(adc_vec);
                
                % Sub-volume: voxels with restricted diffusion (ADC < threshold)
                adc_vec_sub = adc_vec(adc_vec<adc_thresh);
                adc_vec_high_sub = adc_vec(adc_vec>high_adc_thresh);
                % adc_vec_sub(adc_vec_sub>adc_thresh) = nan;

                adc_sub_vol(j,k,dwi_type) = numel(adc_vec_sub)*vox_vol;
                adc_sub_vol_pc(j,k,dwi_type) = adc_sub_vol(j,k,dwi_type)/gtv_vol(j,k);
                adc_sub_mean(j,k,dwi_type) = nanmean(adc_vec_sub);
                if numel(adc_vec_sub) >= min_vox_hist
                    % NOTE: Histogram kurtosis of trace-average map — NOT valid DKI.
                    adc_sub_kurt(j,k,dwi_type) = kurtosis(adc_vec_sub);
                    adc_sub_skew(j,k,dwi_type) = skewness(adc_vec_sub);
                end

                % Normalised histogram (replace zeros with eps for log-safety)
                % store histograms
                [c1, ~] = histcounts(adc_vec, bin_edges);
                p1 = c1 / numel(adc_vec); p1(p1==0)=eps;
                adc_histograms(j,k,:,dwi_type) = p1;
                % Two-sample KS test: current timepoint vs baseline (Fx1)
                % compare difference from baseline (if available)
                if ~isempty(adc_baseline) && numel(adc_vec) >= min_vox_hist && numel(adc_baseline) >= min_vox_hist
                    [~,p,ks2stat] = kstest2(adc_vec,adc_baseline);
                    ks_stats_adc(j,k,dwi_type) = ks2stat;
                    ks_pvals_adc(j,k,dwi_type) = p;
                end

                high_adc_sub_vol(j,k,dwi_type) = numel(adc_vec_high_sub)*vox_vol;

                % Corruption metric: fraction of GTV voxels above adc_max
                fx_corrupted(j,k,dwi_type) = numel(adc_vec(adc_vec>adc_max))/numel(adc_vec);

                % Pool standard-pipeline voxels across all patients/timepoints
                if dwi_type==1
                    adc_vec_all = cat(1,adc_vec_all,adc_vec);
                end
            end

            % --- Compute IVIM summary metrics (D, f, D*) ---
            if ~isempty(d_vec)
                %                 d_vec = data_vectors_gtvp(j,k,1).d_vector;
                %                 f_vec = data_vectors_gtvp(j,k,1).f_vector;
                f_vec(f_vec==0) = nan;  % zero perfusion fraction = failed fit
                %                 dstar_vec = data_vectors_gtvp(j,k,1).dstar_vector;

                % Pool voxels for population-level analysis (standard pipeline only)
                if dwi_type==1
                    d_vec_all = cat(1,d_vec_all,d_vec);
                    f_vec_all = cat(1,f_vec_all,f_vec);
                    fx_vec_all = cat(1,fx_vec_all,ones(size(f_vec))*k);

                    % Stratify pooled voxels by clinical outcome
                    if lf(j)==1  % local failure
                        adc_vec_lf = cat(1,adc_vec_lf,adc_vec);
                        d_vec_lf = cat(1,d_vec_lf,d_vec);
                        f_vec_lf = cat(1,f_vec_lf,f_vec);
                        fx_vec_lf = cat(1,fx_vec_lf,ones(size(f_vec))*k);
                    else         % complete response / no local failure
                        adc_vec_cr = cat(1,adc_vec_cr,adc_vec);
                        d_vec_cr = cat(1,d_vec_cr,d_vec);
                        f_vec_cr = cat(1,f_vec_cr,f_vec);
                        fx_vec_cr = cat(1,fx_vec_cr,ones(size(f_vec))*k);
                    end
                end

                % Joint D–f sub-volume: voxels with low D AND low f
                % come up with a 2D metric for identifying subvolumes
                ivim_vec_sub = d_vec(d_vec<0.001 & f_vec<0.1);
                ivim_sub_vol(j,k,dwi_type) = numel(ivim_vec_sub)*vox_vol;

                % D sub-volume restricted to ADC-thresholded voxels
                d_vec_sub = d_vec(adc_vec<adc_thresh);

                % Whole-GTV D (true diffusion) statistics
                d_mean(j,k,dwi_type) = nanmean(d_vec);
                if numel(d_vec) >= min_vox_hist
                    % NOTE: Histogram kurtosis of trace-average map — NOT valid DKI.
                    d_kurt(j,k,dwi_type) = kurtosis(d_vec);
                    d_skew(j,k,dwi_type) = skewness(d_vec);
                end
                d_sd(j,k,dwi_type) = nanstd(d_vec);

                % Normalised D histogram
                [c1, ~] = histcounts(d_vec, bin_edges);
                p1 = c1 / numel(d_vec); p1(p1==0)=eps;
                d_histograms(j,k,:,dwi_type) = p1;
                % Two-sample KS test: D distribution vs baseline (Fx1)
                % compare difference from baseline (if available)
                if ~isempty(d_baseline) && numel(d_vec) >= min_vox_hist && numel(d_baseline) >= min_vox_hist
                    [~,p,ks2stat] = kstest2(d_vec,d_baseline);
                    ks_stats_d(j,k,dwi_type) = ks2stat;
                    ks_pvals_d(j,k,dwi_type) = p;
                end

                % D sub-volume statistics (restricted region only)
                d_sub_mean(j,k,dwi_type) = nanmean(d_vec_sub);
                if numel(d_vec_sub) >= min_vox_hist
                    % NOTE: Histogram kurtosis of trace-average map — NOT valid DKI.
                    d_sub_kurt(j,k,dwi_type) = kurtosis(d_vec_sub);
                    d_sub_skew(j,k,dwi_type) = skewness(d_vec_sub);
                end

                % Perfusion fraction (f) statistics
                f_mean(j,k,dwi_type) = nanmean(f_vec);
                if numel(f_vec) >= min_vox_hist
                    % NOTE: Histogram kurtosis of trace-average map — NOT valid DKI.
                    f_kurt(j,k,dwi_type) = kurtosis(f_vec);
                    f_skew(j,k,dwi_type) = skewness(f_vec);
                end

                % Pseudo-diffusion coefficient (D*) statistics
                dstar_mean(j,k,dwi_type) = nanmean(dstar_vec);
                if numel(dstar_vec) >= min_vox_hist
                    % NOTE: Histogram kurtosis of trace-average map — NOT valid DKI.
                    dstar_kurt(j,k,dwi_type) = kurtosis(dstar_vec);
                    dstar_skew(j,k,dwi_type) = skewness(dstar_vec);
                end
            end

            % --- Repeatability analysis: extract metrics from Fx1 repeat scans ---
            % Only at k==1 (Fx1) where repeat acquisitions exist
            if k==1
                rp_count = 0;  % count valid repeats for wCV computation
                for rpi=1:nRpt
                    % Select vectors for the appropriate pipeline variant
                    switch dwi_type
                        case 1  % standard
                            adc_vec = data_vectors_gtvp(j,k,rpi).adc_vector;
                            d_vec = data_vectors_gtvp(j,k,rpi).d_vector;
                            f_vec = data_vectors_gtvp(j,k,rpi).f_vector;
                            dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector;
                        case 2  % DnCNN
                            adc_vec = data_vectors_gtvp(j,k,rpi).adc_vector_dncnn;
                            d_vec = data_vectors_gtvp(j,k,rpi).d_vector_dncnn;
                            f_vec = data_vectors_gtvp(j,k,rpi).f_vector_dncnn;
                            dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector_dncnn;
                        case 3  % IVIMnet (no ADC variant; uses standard ADC)
                            adc_vec = [];
                            d_vec = data_vectors_gtvp(j,k,rpi).d_vector_ivimnet;
                            f_vec = data_vectors_gtvp(j,k,rpi).f_vector_ivimnet;
                            dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector_ivimnet;
                    end

                    % Store per-repeat mean ADC and corruption metric
                    if ~isempty(adc_vec)
                        rp_count = rp_count+1;
                        adc_mean_rpt(j,rpi,dwi_type) = nanmean(adc_vec);
                        fx_corrupted_rpt(j,rpi,dwi_type) = numel(adc_vec(adc_vec>adc_max))/numel(adc_vec);

                        adc_vec_sub = adc_vec(adc_vec<adc_thresh);
                        adc_sub_rpt(j,rpi,dwi_type) = nanmean(adc_vec_sub);
                    end

                    % Store per-repeat mean IVIM parameters
                    if ~isempty(d_vec)
                        d_mean_rpt(j,rpi,dwi_type) = nanmean(d_vec);
                        f_mean_rpt(j,rpi,dwi_type) = nanmean(f_vec);
                        dstar_mean_rpt(j,rpi,dwi_type) = nanmean(dstar_vec);
                    end
                end
                % Record number of valid repeats (standard pipeline only)
                if dwi_type==1
                    n_rpt(j) = rp_count;
                end
            end
        end
    end
end
summary_metrics = struct('adc_mean', adc_mean, 'adc_kurt', adc_kurt, 'adc_skew', adc_skew, 'adc_sd', adc_sd, 'd_mean', d_mean, 'f_mean', f_mean, 'dstar_mean', dstar_mean, ...
    'ivim_sub_vol', ivim_sub_vol, ...
    'adc_sub_vol', adc_sub_vol, 'adc_sub_vol_pc', adc_sub_vol_pc, 'high_adc_sub_vol', high_adc_sub_vol, 'd_kurt', d_kurt, 'd_skew', d_skew, 'd_sd', d_sd, ...
    'f_kurt', f_kurt, 'f_skew', f_skew, 'dstar_kurt', dstar_kurt, 'dstar_skew', dstar_skew, 'd_sub_mean', d_sub_mean, 'd_sub_kurt', d_sub_kurt, ...
    'd_sub_skew', d_sub_skew, 'adc_histograms', adc_histograms, 'd_histograms', d_histograms, 'ks_stats_adc', ks_stats_adc, 'ks_pvals_adc', ks_pvals_adc, ...
    'ks_stats_d', ks_stats_d, 'ks_pvals_d', ks_pvals_d, 'fx_corrupted', fx_corrupted, 'gtv_vol', gtv_vol, ...
    'id_list', {id_list}, 'mrn_list', {mrn_list}, 'd95_gtvp', d95_gtvp, 'v50gy_gtvp', v50gy_gtvp, 'lf', lf, 'immuno', immuno, ...
    'adc_mean_rpt', adc_mean_rpt, 'adc_sub_rpt', adc_sub_rpt, 'd_mean_rpt', d_mean_rpt, 'f_mean_rpt', f_mean_rpt, 'dstar_mean_rpt', dstar_mean_rpt, ...
    'n_rpt', n_rpt, 'dmean_gtvp', dmean_gtvp, 'gtv_locations', {gtv_locations}, 'dwi_locations', {dwi_locations});
end