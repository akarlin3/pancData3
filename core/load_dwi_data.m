function [data_vectors_gtvp, data_vectors_gtvn, summary_metrics] = load_dwi_data(config_struct)
%% load_dwi_data_forAvery.m
% Author: Avery Karlin
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
%     - niftiread                 — Native MATLAB NIfTI reader
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
% [CHECKPOINTING STRATEGY]:
% If skip_to_reload is TRUE, the script bypasses raw data discovery and processing,
% assuming 'dwi_vectors.mat' is present. This allows rapid iteration on
% statistical metrics (Section 5) without re-running expensive DICOM conversions.
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

[id_list, mrn_list, fx_dates, dwi_locations, rtdose_locations, gtv_locations, gtvn_locations] = discover_patient_files(config_struct.dataloc);

if isfield(config_struct, 'patient_ids') && ~isempty(config_struct.patient_ids)
    keep_idx = find(ismember(mrn_list, config_struct.patient_ids));
    if ~isempty(keep_idx)
        id_list = id_list(keep_idx);
        mrn_list = mrn_list(keep_idx);
        fx_dates = fx_dates(keep_idx);
        dwi_locations = dwi_locations(keep_idx, :, :);
        if ~isempty(rtdose_locations)
            rtdose_locations = rtdose_locations(keep_idx, :);
        end
        gtv_locations = gtv_locations(keep_idx, :, :);
        gtvn_locations = gtvn_locations(keep_idx, :, :);
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
try
    T = readtable(clinical_data_sheet,'Sheet','Clin List');
catch ME
    if exist('OCTAVE_VERSION', 'builtin')
        warning('Octave could not read %s. Creating dummy clinical data table.', clinical_data_sheet);
        T = struct();
        T.Pat = id_list;
        T.MRN = mrn_list;
        T.LF = zeros(size(id_list));
        T.IO_start = repmat({''}, size(id_list));
    else
        rethrow(ME);
    end
end

% load dwi data locations (optionally reload from a previous run)
data_file = fullfile(dataloc, 'adc_vectors.mat');

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
bad_dwi_locations_per_patient = cell(length(mrn_list), 1);

% Checkpoint directory setup
checkpoint_dir = fullfile(dataloc, 'processed_patients');
if ~isfolder(checkpoint_dir)
    mkdir(checkpoint_dir);
end

% Scan for existing checkpoints to identify completed patients
patient_completed = false(size(mrn_list));
for j = 1:length(mrn_list)
    mrn = mrn_list{j};
    patient_id = id_list{j};
    checkpoint_file = fullfile(checkpoint_dir, sprintf('patient_%03d_%s.mat', j, patient_id));
    if exist(checkpoint_file, 'file')
        patient_completed(j) = true;
    end
end

% --- Main processing loop: iterate over patients ---
% [PARALLELIZATION STRATEGY]:
% Iterates over patients (j) in parallel using parfor. Since patients are
% statistically independent units, this is embarrassingly parallel.
% Output variables (data_vectors_gtvp/n, summary arrays) are sliced variables,
% meaning each worker writes to a specific index (j,:,:), preventing race conditions.
% 'bad_dwi_locations_per_patient' accumulates errors per worker and is flattened later.

% [PERFORMANCE OPTIMIZATION]:
% Pre-compute normalized strings outside the loop to avoid redundant strrep calls.
% Ensure T.Pat is a cellstr before passing to strrep
if exist('OCTAVE_VERSION', 'builtin')
    % iscategorical is missing or mocked, T.Pat might be char array, need to make sure it's cellstr
    if isfield(T, 'Pat')
        T_Pat_cell_tmp = T.Pat;
        if ischar(T_Pat_cell_tmp)
            if size(T_Pat_cell_tmp, 1) > 1
                % if multiple rows
                T_Pat_cell = {};
                for i_pat_row = 1:size(T_Pat_cell_tmp, 1)
                    T_Pat_cell{i_pat_row} = strtrim(T_Pat_cell_tmp(i_pat_row, :));
                end
            else
                T_Pat_cell = {T_Pat_cell_tmp};
            end
        elseif isnumeric(T_Pat_cell_tmp)
            T_Pat_cell = {};
        elseif iscell(T_Pat_cell_tmp)
            T_Pat_cell = T_Pat_cell_tmp;
        else
            T_Pat_cell = {T_Pat_cell_tmp};
        end
    else
        T_Pat_cell = {};
    end
else
    if iscategorical(T.Pat)
        T_Pat_cell = cellstr(T.Pat);
    else
        T_Pat_cell = T.Pat;
    end
end

if isempty(T_Pat_cell)
    T_Pat_normalized = {};
else
    try
        T_Pat_normalized = strrep(T_Pat_cell, '_', '-');
    catch
        T_Pat_normalized = {};
    end
end

if isempty(id_list)
    id_list_normalized = {};
else
    try
        id_list_normalized = strrep(id_list, '_', '-');
    catch
        id_list_normalized = {};
    end
end

parfor j = 1:length(mrn_list)
    mrn = mrn_list{j};

    if patient_completed(j)
        fprintf('Skipping patient %d/%d (Patient ID %s) - already processed.\n', j, length(mrn_list), patient_id);
        continue;
    end

    fprintf('\n******* MRN: %s\n',mrn);

    max_dwis_j = size(dwi_locations, 2) * size(dwi_locations, 3);
    bad_dwi_list_j = cell(1, max_dwis_j);
    bad_dwi_idx_j = 0;

    % Localized output variables for this patient iteration
    pat_data_vectors_gtvp = struct;
    pat_data_vectors_gtvn = struct;
    pat_dmean_gtvp = nan(1, size(dwi_locations, 2));
    pat_dmean_gtvn = nan(1, size(dwi_locations, 2));
    pat_d95_gtvp = nan(1, size(dwi_locations, 2));
    pat_d95_gtvn = nan(1, size(dwi_locations, 2));
    pat_v50gy_gtvp = nan(1, size(dwi_locations, 2));
    pat_v50gy_gtvn = nan(1, size(dwi_locations, 2));
    pat_adc_mean = nan(1, size(dwi_locations, 2), size(dwi_locations, 3));
    pat_adc_kurtosis = nan(1, size(dwi_locations, 2), size(dwi_locations, 3));
    pat_d_mean = nan(1, size(dwi_locations, 2), size(dwi_locations, 3));
    pat_d_kurtosis = nan(1, size(dwi_locations, 2), size(dwi_locations, 3));
    pat_d_mean_dncnn = nan(1, size(dwi_locations, 2), size(dwi_locations, 3));
    pat_d_mean_ivimnet = nan(1, size(dwi_locations, 2), size(dwi_locations, 3));

    % Initialize potential fields to ensure struct consistency
    for fi=1:size(dwi_locations,2)
        for rpi=1:size(dwi_locations,3)
            % GTVp fields
            pat_data_vectors_gtvp(fi,rpi).adc_vector = [];
            pat_data_vectors_gtvp(fi,rpi).d_vector = [];
            pat_data_vectors_gtvp(fi,rpi).f_vector = [];
            pat_data_vectors_gtvp(fi,rpi).dstar_vector = [];
            pat_data_vectors_gtvp(fi,rpi).dose_vector = [];
            pat_data_vectors_gtvp(fi,rpi).dvh = [];
            pat_data_vectors_gtvp(fi,rpi).d95 = [];
            pat_data_vectors_gtvp(fi,rpi).v50gy = [];
            pat_data_vectors_gtvp(fi,rpi).d_vector_dncnn = [];
            pat_data_vectors_gtvp(fi,rpi).f_vector_dncnn = [];
            pat_data_vectors_gtvp(fi,rpi).dstar_vector_dncnn = [];
            pat_data_vectors_gtvp(fi,rpi).adc_vector_dncnn = [];
            pat_data_vectors_gtvp(fi,rpi).d_vector_ivimnet = [];
            pat_data_vectors_gtvp(fi,rpi).f_vector_ivimnet = [];
            pat_data_vectors_gtvp(fi,rpi).dstar_vector_ivimnet = [];

            % GTVn fields
            pat_data_vectors_gtvn(fi,rpi).adc_vector = [];
            pat_data_vectors_gtvn(fi,rpi).d_vector = [];
            pat_data_vectors_gtvn(fi,rpi).f_vector = [];
            pat_data_vectors_gtvn(fi,rpi).dstar_vector = [];
            pat_data_vectors_gtvn(fi,rpi).dose_vector = [];
            pat_data_vectors_gtvn(fi,rpi).dvh = [];
            pat_data_vectors_gtvn(fi,rpi).d95 = [];
            pat_data_vectors_gtvn(fi,rpi).v50gy = [];
            pat_data_vectors_gtvn(fi,rpi).d_vector_dncnn = [];
            pat_data_vectors_gtvn(fi,rpi).f_vector_dncnn = [];
            pat_data_vectors_gtvn(fi,rpi).dstar_vector_dncnn = [];
            pat_data_vectors_gtvn(fi,rpi).adc_vector_dncnn = [];
            pat_data_vectors_gtvn(fi,rpi).d_vector_ivimnet = [];
            pat_data_vectors_gtvn(fi,rpi).f_vector_ivimnet = [];
            pat_data_vectors_gtvn(fi,rpi).dstar_vector_ivimnet = [];
        end
    end

    % Per-patient DIR reference: populated at Fx1, reused at Fx2+
    b0_fx1_ref        = [];   % b=0 volume at baseline fraction
    gtv_mask_fx1_ref  = [];   % GTVp mask at baseline fraction
    gtvn_mask_fx1_ref = [];   % GTVn mask at baseline fraction (when present)

    % read clinical data (local failure, immunotherapy) from spreadsheet
    i_pat = find(contains(T_Pat_normalized, id_list_normalized{j}));
    pat_immuno = T.Immuno(i_pat(1));
    pat_lf = T.LF(i_pat(1));

    basefolder = fullfile(dataloc, id_list{j});
    basefolder_contents = clean_dir_command(basefolder);

    % --- Loop over fractions (fi) and repeat acquisitions (rpi) ---
    for fi=1:size(dwi_locations,2)
        % Pre-filter fraction folder for this fraction to avoid redundant dir() calls
        fxtmp_idx = contains({basefolder_contents.name}, fx_search{fi});
        fxtmp = basefolder_contents(fxtmp_idx);

        for rpi = 1:size(dwi_locations,3)

            if isempty(fxtmp)
                fprintf('%s, no %s folder\n',id_list{j},fx_search{fi});
            else
                fxfolder = fullfile(basefolder, fxtmp(1).name);
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

            outloc = fullfile(basefolder, 'nii');   % output directory for NIfTI files

            % get the dwi file list (DICOM)
            dicom_files = dir(fullfile(dicomloc, '*.dcm'));
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
                bad_dwi_found_flag = convert_dicom(dicomloc, outloc, scanID, dcm2nii_call, fx_id);
                if bad_dwi_found_flag
                    bad_dwi_idx_j = bad_dwi_idx_j + 1;
                    bad_dwi_list_j{bad_dwi_idx_j} = dicomloc;
                    bad_dwi_found = 1;
                end
            end

            % --- Save GTVp mask as NIfTI for consistency with DWI volumes ---
            % save the mask as nii.gz (just for consistency)
            % Stvol3d is the 3-D binary mask variable stored in the .mat file
            if ~isempty(struct_file)
                if ~exist(fullfile(outloc, [gtvname '.nii.gz']),'file')
                    gtv_mask = safe_load_mask(struct_file, 'Stvol3d');
                    if ~isempty(gtv_mask)
                        niftiwrite(rot90(double(gtv_mask),-1),fullfile(outloc, gtvname),'Compressed',true);
                    else
                         fprintf('Warning: Failed to load GTVp mask from %s\n', struct_file);
                    end
                end
            end

            % --- Save GTVn (nodal) mask as NIfTI, if present ---
            if ~isempty(struct_file_gtvn)
                if ~exist(fullfile(outloc, [gtvn_name '.nii.gz']),'file')
                    gtvn_mask = safe_load_mask(struct_file_gtvn, 'Stvol3d');
                    if ~isempty(gtvn_mask)
                        niftiwrite(rot90(double(gtvn_mask),-1),fullfile(outloc, gtvn_name),'Compressed',true);
                    else
                        fprintf('Warning: Failed to load GTVn mask from %s\n', struct_file_gtvn);
                    end
                end
            end

            % --- Resample RT dose onto DWI geometry and save as NIfTI ---
            % sample RTdose on the DWI geometry, save as nii.gz
            % but first, see if we actually have a dose map to load....
            if ~isempty(dicomdoseloc) && ~isempty(dicomloc)
                if ~exist(fullfile(outloc, [dosename '.nii.gz']),'file')
                    % extract only b=0 images — these define the DWI spatial
                    % geometry used as the resampling target for the dose grid
                    % extract only b0 images for the dwi sampling
                    isb0 = zeros(length(dicom_files),1);
                    b0list = cell(1);
                    b0count = 0;
                    for bi = 1:length(isb0)
                        data_tmp = dicominfo(fullfile(dicom_files(bi).folder, dicom_files(bi).name), 'UseDictionaryVR', true);
                        if data_tmp.DiffusionBValue == 0
                            isb0(bi) = 1;
                            b0count = b0count+1;
                            b0list{b0count,1} = fullfile(dicom_files(bi).folder, dicom_files(bi).name);
                        end
                    end

                    rtdose_dicom = dir(fullfile(dicomdoseloc, '*.dcm'));
                    rtdosefile = fullfile(rtdose_dicom.folder, rtdose_dicom.name);

                    % Resample RT dose onto b=0 image grid and write NIfTI
                    dose_sampled = sample_rtdose_on_image(b0list,rtdosefile);
                    niftiwrite(rot90(dose_sampled,-1),fullfile(outloc, dosename),'Compressed',true);
                end
            end

            % --- Load NIfTI DWI volume and extract b-values ---
            % load the data
            havedwi = 0;
            if exist(fullfile(outloc, [scanID '.nii.gz']),'file')
                dwi_info = niftiinfo(fullfile(outloc, [scanID '.nii.gz']));
                dwi = rot90(niftiread(dwi_info));            % orient to standard view
                dwi_dims = dwi_info.PixelDimensions(1:3);  % voxel dimensions in mm (Native NIfTI)
                dwi_vox_vol = prod(dwi_dims*0.1);    % voxel volume in cc (mm→cm)
                fprintf('...Loaded %s. ',fullfile(outloc, [scanID '.nii.gz']));
                havedwi = 1;
                % Read b-values from the sidecar .bval file produced by dcm2niix
                % now extract bvalues
                bval_file = fullfile(outloc, [scanID '.bval']);
                if exist(bval_file,'file')
                    fid = fopen(bval_file);
                    tline = fgetl(fid);
                    bval_data = tline;
                    fclose(fid);

                    bvalues = sscanf(bval_data, '%f');

                    % Sort b-values ascending and reorder the 4-D DWI volume
                    [b_sort,i_sort] = sort(bvalues,'ascend');

                    bvalues = bvalues(i_sort);
                    dwi = double(dwi(:,:,:,i_sort));
                    fprintf('loaded bvalues\n');
                else
                    fprintf('bvalue file not found!\n');
                    havedwi = 0;
                    if bad_dwi_found==0
                        bad_dwi_idx_j = bad_dwi_idx_j + 1;
                        bad_dwi_list_j{bad_dwi_idx_j} = dicomloc;
                        bad_dwi_found = 1;
                    end
                end
                % Expect exactly 4 b-values (e.g., 0, 50, 400, 800 s/mm²)
                if size(dwi,4)~=4
                    fprintf('DWI does not have expected dimensions found: %s skipping\n',mat2str(size(dwi)))
                    havedwi = 0;
                    if bad_dwi_found==0
                        bad_dwi_idx_j = bad_dwi_idx_j + 1;
                        bad_dwi_list_j{bad_dwi_idx_j} = dicomloc;
                        bad_dwi_found = 1;
                    end
                end
            end

            % --- Load DnCNN-denoised DWI (deep learning denoising) ---
            havedenoised = 0;
            if havedwi==1
                dncnnid = [scanID '_dncnn.nii.gz'];
                dncnn_file = fullfile(basefolder, 'dncnn', dncnnid);
                if exist(dncnn_file,'file')
                    dncnn_info = niftiinfo(dncnn_file);
                    dwi_dncnn = rot90(niftiread(dncnn_info));
                    % Normalise denoised signal to [0,1] for IVIM fitting
                    dwi_dncnn = double(mat2gray(dwi_dncnn(:,:,:,i_sort)));
                    havedenoised=1;
                else
                    % [MODULARIZATION STAGE 2]: CPU Acceleration + Dynamic Fallback
                    % If the cached 'dwi_dncnn' doesn't exist, we compute it on the fly.
                    % We pull the heavy dependency (apply_dncnn_symmetric) into this loop,
                    % but we cast the matrices to single() *first* so MATLAB accelerates it.
                    fprintf('  [DnCNN] Cache missing. Executing deep learning denoising on CPU...\n');

                    % 1. Load the pre-trained neural network (assuming it exists in dependencies)
                    try
                        % NOTE: Replace 'dncnn_model.mat' with the actual model file if known.
                        % For now, we assume a generic 'net' variable is loaded.
                        loaded_model = load(fullfile(config_struct.dataloc, '../dependencies/dncnn_model.mat'), 'net');
                        dncnn_net = loaded_model.net;

                        % 2. Cast the raw 4D DWI volume and 3D GTV mask
                        dwi_cpu = single(dwi);

                        % If we have a primary GTV, use it. Otherwise, use nodal if available.
                        if exist('gtv_mask', 'var') && ~isempty(gtv_mask)
                            mask_cpu = single(gtv_mask);
                        elseif exist('gtvn_mask', 'var') && ~isempty(gtvn_mask)
                            mask_cpu = single(gtvn_mask);
                        else
                            mask_cpu = ones(size(dwi, 1), size(dwi, 2), size(dwi, 3), 'single');
                        end

                        % 3. Pre-allocate the denoised 4D volume
                        dwi_dncnn_cpu = zeros(size(dwi_cpu), 'single');

                        % 4. Apply the black-box dependency slice-by-slice (or volume-by-volume)
                        %    apply_dncnn_symmetric handles 3D, so we loop over the 4th dimension (b-values)
                        for b_idx = 1:size(dwi_cpu, 4)
                            % The dependency is untouched, but it operates on CPU arrays natively!
                            dwi_dncnn_cpu(:,:,:,b_idx) = apply_dncnn_symmetric(dwi_cpu(:,:,:,b_idx), mask_cpu, dncnn_net, 15);
                        end

                        % 5. Format back to double for the rest of the pipeline
                        dwi_dncnn = double(dwi_dncnn_cpu);

                        % Normalise denoised signal to [0,1] for IVIM fitting
                        dwi_dncnn = double(mat2gray(dwi_dncnn(:,:,:,i_sort)));
                        havedenoised = 1;

                        % (Optional: Save the dwi_dncnn back to disk to cache it for the next run)
                        % This would require creating the directories and saving the NIfTI.
                        fprintf('  [DnCNN] Deep learning denoising completed.\n');
                    catch CPU_ME
                        fprintf('  [DnCNN] CPU Computation failed: %s\n', CPU_ME.message);
                        fprintf('  [DnCNN] Proceeding without denoised data.\n');
                        havedenoised = 0;
                    end
                end
            end

            % --- Load IVIMnet deep-learning fit results (pre-computed) ---
            haveivimnet = 0;
            if havedwi==1
                ivimid = [scanID '_ivimnet.mat'];
                ivimnet_file = fullfile(basefolder, 'ivimnet', ivimid);
                if exist(ivimnet_file,'file')
                    % Loads D_ivimnet, f_ivimnet, Dstar_ivimnet, S0_ivimnet
                    tmp = load(ivimnet_file);
                    D_ivimnet = tmp.D_ivimnet;
                    f_ivimnet = tmp.f_ivimnet;
                    Dstar_ivimnet = tmp.Dstar_ivimnet;
                    S0_ivimnet = tmp.S0_ivimnet;
                    haveivimnet=1;
                end
            end

            dwi_size = size(dwi);

            % --- Load GTVp mask and validate spatial dimensions ---
            gtvp_filepath = fullfile(outloc, [gtvname '.nii.gz']);
            [havegtvp, gtv_mask] = load_mask(gtvp_filepath, dwi_size, '', 'gtvp');

            % --- Load GTVn (nodal) mask and validate dimensions ---
            gtvn_filepath = fullfile(outloc, [gtvn_name '.nii.gz']);
            [havegtvn, gtvn_mask] = load_mask(gtvn_filepath, dwi_size, '*NODAL* ', 'gtvn');

            % --- Load resampled RT dose map ---
            havedose = 0;
            if exist(fullfile(outloc, [dosename '.nii.gz']),'file')
                dose_info = niftiinfo(fullfile(outloc, [dosename '.nii.gz']));
                dose_map = rot90(niftiread(dose_info));
                fprintf('...Loaded %s\n',fullfile(outloc, [dosename '.nii.gz']));
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
                % [PHYSICS EXPLANATION]:
                % The Intravoxel Incoherent Motion (IVIM) model describes signal decay as:
                % S/S0 = f * exp(-b * D*) + (1-f) * exp(-b * D)
                % where:
                %   D  = True diffusion coefficient (thermal Brownian motion)
                %   f  = Perfusion fraction (volume fraction of capillaries)
                %   D* = Pseudo-diffusion coefficient (blood flow velocity)
                %
                % Segmented Fitting Approach:
                % 1. Fit D using only high b-values (>= bthr), where perfusion contribution is negligible.
                %    Approximation: S ~ (1-f) * exp(-b * D)  => log(S) linear vs b
                % 2. Fix D, then fit f and D* using all b-values (or low b-values).
                opts = [];
                opts.bthr = ivim_bthr; % set in USER OPTIONS (default 200 s/mm2; clinical standard for abdominal IVIM)

                [d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvalues, mask_ivim, opts);
                % Repeat IVIM + ADC fitting on DnCNN-denoised data
                if havedenoised==1
                    [d_map_dncnn, f_map_dncnn, dstar_map_dncnn, adc_map_dncnn] = fit_models(dwi_dncnn, bvalues, mask_ivim, opts);
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
                    dir_cache_file = fullfile(dataloc, id_list{j}, 'nii', ...
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
                            % Use helper function to save within parfor loop (transparency fix)
                            parsave_dir_cache(dir_cache_file, gtv_mask_warped, D_forward, ref3d);
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
                    'FillValues', nan);
                f_map_dncnn     = imwarp(f_map_dncnn,     -D_forward_cur, 'Interp', 'linear', ...
                    'FillValues', nan);
                dstar_map_dncnn = imwarp(dstar_map_dncnn, -D_forward_cur, 'Interp', 'linear', ...
                    'FillValues', nan);
                adc_map_dncnn   = imwarp(adc_map_dncnn,   -D_forward_cur, 'Interp', 'linear', ...
                    'FillValues', nan);
                fprintf('  [DnCNN] Warped native-space parameter maps to baseline geometry.\n');
            end

            % --- Extract voxel-level biomarkers within GTVp ---
            if havegtvp
                % Compute summary statistics within the primary GTV
                pat_adc_mean(1,fi,rpi) = nanmean(adc_map(gtv_mask==1));
                % NOTE: Histogram kurtosis of trace-average ADC — NOT valid DKI.
                pat_adc_kurtosis(1,fi,rpi) = kurtosis(adc_map(gtv_mask==1));

                pat_d_mean(1,fi,rpi) = nanmean(d_map(gtv_mask==1));
                % NOTE: Histogram kurtosis of trace-average map — NOT valid DKI.
                pat_d_kurtosis(1,fi,rpi) = kurtosis(d_map(gtv_mask==1));

                % Store voxel-level vectors and metadata in the output struct
                pat_data_vectors_gtvp(fi,rpi).adc_vector = adc_map(gtv_mask==1);
                pat_data_vectors_gtvp(fi,rpi).d_vector = d_map(gtv_mask==1);
                pat_data_vectors_gtvp(fi,rpi).f_vector = f_map(gtv_mask==1);
                pat_data_vectors_gtvp(fi,rpi).dstar_vector = dstar_map(gtv_mask==1);
                pat_data_vectors_gtvp(fi,rpi).ID = id_list{j};
                pat_data_vectors_gtvp(fi,rpi).MRN = mrn_list{j};
                pat_data_vectors_gtvp(fi,rpi).LF = pat_lf;
                pat_data_vectors_gtvp(fi,rpi).Immuno = pat_immuno;
                pat_data_vectors_gtvp(fi,rpi).Fraction = fi;
                pat_data_vectors_gtvp(fi,rpi).Repeatability_index = rpi;
                pat_data_vectors_gtvp(fi,rpi).vox_vol = dwi_vox_vol;

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
                    pat_d_mean_dncnn(1,fi,rpi) = nanmean(d_map_dncnn(dncnn_mask_p==1));
                    pat_data_vectors_gtvp(fi,rpi).d_vector_dncnn    = d_map_dncnn(dncnn_mask_p==1);
                    pat_data_vectors_gtvp(fi,rpi).f_vector_dncnn    = f_map_dncnn(dncnn_mask_p==1);
                    pat_data_vectors_gtvp(fi,rpi).dstar_vector_dncnn = dstar_map_dncnn(dncnn_mask_p==1);
                    pat_data_vectors_gtvp(fi,rpi).adc_vector_dncnn  = adc_map_dncnn(dncnn_mask_p==1);
                end

                % Store IVIMnet deep-learning fit vectors (pipeline variant 3)
                if haveivimnet
                    D_ivimnet = rot90(D_ivimnet);
                    f_ivimnet = rot90(f_ivimnet);
                    Dstar_ivimnet = rot90(Dstar_ivimnet);
                    S0_ivimnet = rot90(S0_ivimnet);
                    pat_d_mean_ivimnet(1,fi,rpi) = nanmean(D_ivimnet(gtv_mask==1));

                    pat_data_vectors_gtvp(fi,rpi).d_vector_ivimnet = D_ivimnet(gtv_mask==1);
                    pat_data_vectors_gtvp(fi,rpi).f_vector_ivimnet = f_ivimnet(gtv_mask==1);
                    pat_data_vectors_gtvp(fi,rpi).dstar_vector_ivimnet = Dstar_ivimnet(gtv_mask==1);
                end
            end

            if havedose && havegtvp
                pat_dmean_gtvp(1,fi) = nanmean(dose_map_dvh(gtv_mask_for_dvh==1));
                % DVH uses the strictly rigid dose against the DIR-warped daily GTV tissue mask
                % so that dose correctly reflects the true static beam delivered to deformed anatomy.
                [dvhparams, dvh_values] = dvh(dose_map_dvh, gtv_mask_for_dvh, dwi_dims, 2000, 'Dperc',95,'Vperc',50,'Normalize',true);
                pat_d95_gtvp(1,fi) = dvhparams.("D95% (Gy)");
                pat_v50gy_gtvp(1,fi) = dvhparams.("V50Gy (%)");

                pat_data_vectors_gtvp(fi,rpi).dose_vector = dose_map_dvh(gtv_mask_for_dvh==1);
                pat_data_vectors_gtvp(fi,rpi).dvh = dvh_values;
                pat_data_vectors_gtvp(fi,rpi).d95 = dvhparams.("D95% (Gy)");
                pat_data_vectors_gtvp(fi,rpi).v50gy = dvhparams.("V50Gy (%)");
            end

            % --- Extract voxel-level biomarkers within GTVn (nodal GTV) ---
            % now collect gtvn info
            if havegtvn
                pat_data_vectors_gtvn(fi,rpi).adc_vector = adc_map(gtvn_mask==1);
                pat_data_vectors_gtvn(fi,rpi).d_vector = d_map(gtvn_mask==1);
                pat_data_vectors_gtvn(fi,rpi).f_vector = f_map(gtvn_mask==1);
                pat_data_vectors_gtvn(fi,rpi).dstar_vector = dstar_map(gtvn_mask==1);
                pat_data_vectors_gtvn(fi,rpi).ID = id_list{j};
                pat_data_vectors_gtvn(fi,rpi).MRN = mrn_list{j};
                pat_data_vectors_gtvn(fi,rpi).LF = pat_lf;
                pat_data_vectors_gtvn(fi,rpi).Immuno = pat_immuno;
                pat_data_vectors_gtvn(fi,rpi).Fraction = fi;
                pat_data_vectors_gtvn(fi,rpi).Repeatability_index = rpi;
                pat_data_vectors_gtvn(fi,rpi).vox_vol = dwi_vox_vol;
            end

            % Store DnCNN-denoised vectors for GTVn.
            % For fi > 1 the parameter maps have been warped to baseline geometry;
            % extract within the baseline GTVn mask for spatial consistency.
            if havedenoised && havegtvn
                dncnn_mask_n = gtvn_mask;
                if fi > 1 && ~isempty(gtvn_mask_fx1_ref)
                    dncnn_mask_n = gtvn_mask_fx1_ref;
                end
                pat_data_vectors_gtvn(fi,rpi).d_vector_dncnn    = d_map_dncnn(dncnn_mask_n==1);
                pat_data_vectors_gtvn(fi,rpi).f_vector_dncnn    = f_map_dncnn(dncnn_mask_n==1);
                pat_data_vectors_gtvn(fi,rpi).dstar_vector_dncnn = dstar_map_dncnn(dncnn_mask_n==1);
                pat_data_vectors_gtvn(fi,rpi).adc_vector_dncnn  = adc_map_dncnn(dncnn_mask_n==1);
            end

            % Store IVIMnet vectors for GTVn
            if haveivimnet && havegtvn
                D_ivimnet = rot90(D_ivimnet);
                f_ivimnet = rot90(f_ivimnet);
                Dstar_ivimnet = rot90(Dstar_ivimnet);
                S0_ivimnet = rot90(S0_ivimnet);

                pat_data_vectors_gtvn(fi,rpi).d_vector_ivimnet = D_ivimnet(gtvn_mask==1);
                pat_data_vectors_gtvn(fi,rpi).f_vector_ivimnet = f_ivimnet(gtvn_mask==1);
                pat_data_vectors_gtvn(fi,rpi).dstar_vector_ivimnet = Dstar_ivimnet(gtvn_mask==1);
            end

            % Compute DVH parameters within GTVn
            % Sample rigidly aligned dose map against the GTVn mask.
            if havedose && havegtvn
                dose_map_dvh_n = dose_map;  % rigidly aligned dose
                pat_dmean_gtvn(1,fi) = nanmean(dose_map_dvh_n(gtvn_mask==1));
                [dvhparams, dvh_values] = dvh(dose_map_dvh_n, gtvn_mask, dwi_dims, 2000, 'Dperc',95,'Vperc',50,'Normalize',true);
                pat_d95_gtvn(1,fi) = dvhparams.("D95% (Gy)");
                pat_v50gy_gtvn(1,fi) = dvhparams.("V50Gy (%)");

                pat_data_vectors_gtvn(fi,rpi).dose_vector = dose_map_dvh_n(gtvn_mask==1);
                pat_data_vectors_gtvn(fi,rpi).dvh = dvh_values;
                pat_data_vectors_gtvn(fi,rpi).d95 = dvhparams.("D95% (Gy)");
                pat_data_vectors_gtvn(fi,rpi).v50gy = dvhparams.("V50Gy (%)");
            end
        end
    end
    bad_dwi_list_j = bad_dwi_list_j(1:bad_dwi_idx_j);
    bad_dwi_locations_per_patient{j} = bad_dwi_list_j;

    % Collect output struct for checkpointing
    pat_data_out = struct();
    pat_data_out.data_vectors_gtvp = pat_data_vectors_gtvp;
    pat_data_out.data_vectors_gtvn = pat_data_vectors_gtvn;
    pat_data_out.dmean_gtvp = pat_dmean_gtvp;
    pat_data_out.dmean_gtvn = pat_dmean_gtvn;
    pat_data_out.d95_gtvp = pat_d95_gtvp;
    pat_data_out.d95_gtvn = pat_d95_gtvn;
    pat_data_out.v50gy_gtvp = pat_v50gy_gtvp;
    pat_data_out.v50gy_gtvn = pat_v50gy_gtvn;
    pat_data_out.adc_mean = pat_adc_mean;
    pat_data_out.adc_kurtosis = pat_adc_kurtosis;
    pat_data_out.d_mean = pat_d_mean;
    pat_data_out.d_kurtosis = pat_d_kurtosis;
    pat_data_out.d_mean_dncnn = pat_d_mean_dncnn;
    pat_data_out.d_mean_ivimnet = pat_d_mean_ivimnet;
    pat_data_out.lf = pat_lf;
    pat_data_out.immuno = pat_immuno;
    pat_data_out.bad_dwi_list = bad_dwi_list_j;

    % Save checkpoint
    checkpoint_file = fullfile(checkpoint_dir, sprintf('patient_%03d_%s.mat', j, patient_id));
    parsave_checkpoint(checkpoint_file, pat_data_out);

    fprintf('Finished processing patient %d/%d (MRN: %s)\n', j, length(mrn_list), mrn);
end

% Reconstruct global arrays from checkpoints
for j = 1:length(mrn_list)
    mrn = mrn_list{j};
    checkpoint_file = fullfile(checkpoint_dir, sprintf('patient_%03d_%s.mat', j, mrn));

    if exist(checkpoint_file, 'file')
        % Load checkpoint
        loaded_data = load(checkpoint_file);

        % Assign back to global arrays
        % Struct arrays
        data_vectors_gtvp(j,:,:) = loaded_data.data_vectors_gtvp;
        data_vectors_gtvn(j,:,:) = loaded_data.data_vectors_gtvn;

        % Scalar/Vector arrays (patient x fraction)
        dmean_gtvp(j,:) = loaded_data.dmean_gtvp;
        dmean_gtvn(j,:) = loaded_data.dmean_gtvn;
        d95_gtvp(j,:) = loaded_data.d95_gtvp;
        d95_gtvn(j,:) = loaded_data.d95_gtvn;
        v50gy_gtvp(j,:) = loaded_data.v50gy_gtvp;
        v50gy_gtvn(j,:) = loaded_data.v50gy_gtvn;

        % Summary metrics (patient x fraction x repeat)
        adc_mean(j,:,:) = loaded_data.adc_mean;
        if isfield(loaded_data, 'adc_kurtosis')
            adc_kurtosis(j,:,:) = loaded_data.adc_kurtosis;
        end
        d_mean(j,:,:) = loaded_data.d_mean;
        if isfield(loaded_data, 'd_kurtosis')
            d_kurtosis(j,:,:) = loaded_data.d_kurtosis;
        end
        d_mean_dncnn(j,:,:) = loaded_data.d_mean_dncnn;
        d_mean_ivimnet(j,:,:) = loaded_data.d_mean_ivimnet;

        % Clinical data and tracking
        lf(j) = loaded_data.lf;
        immuno(j) = loaded_data.immuno;
        bad_dwi_locations_per_patient{j} = loaded_data.bad_dwi_list;
    else
        fprintf('Warning: No checkpoint found for patient %d (Patient ID %s) during reconstruction.\n', j, patient_id);
    end
end

% Flatten bad_dwi_locations
bad_dwi_locations = [bad_dwi_locations_per_patient{:}];
bad_dwi_count = length(bad_dwi_locations);

%% ========================================================================
fprintf('\n--- SECTION 3: Save Results ---\n');
%  SECTION 3 — SAVE RESULTS
%  [CHECKPOINT]: Saves the output of the computationally intensive Section 2.
%  This .mat file serves as the input for Section 4, allowing the pipeline
%  to resume from here in future runs.

datasave = fullfile(dataloc, 'dwi_vectors.mat');
% Create a date-stamped backup before overwriting
if exist(datasave,'file')
    dt = datetime('now');
    dateString = char(dt, 'yyyy_MMM_dd');
    newfilename = cat(2,dataloc, 'dwi_vectors_', dateString, '.mat');
    copyfile(datasave,newfilename);
    fprintf('backed up existing save to %s\n',newfilename);
end
if isfield(config_struct, 'dwi_type_name')
    file_prefix = ['_' config_struct.dwi_type_name];
else
    file_prefix = '';
end
datasave = fullfile(dataloc, ['dwi_vectors' file_prefix '.mat']);
save(datasave,'data_vectors_gtvn','data_vectors_gtvp','lf','immuno','mrn_list','id_list','fx_dates','dwi_locations','rtdose_locations','gtv_locations','gtvn_locations','dmean_gtvp','dmean_gtvn','d95_gtvp','d95_gtvn','v50gy_gtvp','v50gy_gtvn','bad_dwi_locations','bad_dwi_count');
fprintf('saved %s\n',datasave);

else
    %% ========================================================================
    fprintf('\n--- SECTION 4: Reload Saved Data ---\n');
    %  SECTION 4 — RELOAD SAVED DATA
    %  [ENTRY POINT]: If skip_to_reload=true, execution begins here.
    %  Loads the pre-processed 'dwi_vectors.mat' containing voxel-level data.

    % Set data path from configuration
    dataloc = config_struct.dataloc;

    if isfield(config_struct, 'dwi_type_name')
        file_prefix = ['_' config_struct.dwi_type_name];
    else
        file_prefix = '';
    end
    datasave = fullfile(dataloc, ['dwi_vectors' file_prefix '.mat']);
    if exist(datasave, 'file')
        load(datasave);
    else
        fallback_datasave = fullfile(dataloc, 'dwi_vectors.mat');
        if exist(fallback_datasave, 'file')
            fprintf('  Specific %s not found. Falling back to %s\n', ['dwi_vectors' file_prefix '.mat'], 'dwi_vectors.mat');
            load(fallback_datasave);
        else
            error('Unable to find file or directory ''%s''.', datasave);
        end
    end

if exist('OCTAVE_VERSION', 'builtin') && ~exist('id_list', 'var')
    warning('id_list not loaded from save file. This may occur during mock tests. Proceeding with dummy data.');
    id_list = {}; mrn_list = {}; lf = []; immuno = {}; gtv_locations = []; dwi_locations = []; dmean_gtvp = []; d95_gtvp = []; v50gy_gtvp = []; data_vectors_gtvp = [];
end
end % if ~skip_to_reload

%% ========================================================================
fprintf('\n--- SECTION 5: Longitudinal Summary Metrics ---\n');
%  SECTION 5 — LONGITUDINAL SUMMARY METRICS

summary_metrics = compute_summary_metrics(config_struct, data_vectors_gtvp, id_list, mrn_list, lf, immuno, gtv_locations, dwi_locations, dmean_gtvp, d95_gtvp, v50gy_gtvp);

end

function parsave_checkpoint(fname, data)
    save(fname, '-struct', 'data');
end

function parsave_dir_cache(fname, gtv_mask_warped, D_forward, ref3d)
    save(fname, 'gtv_mask_warped', 'D_forward', 'ref3d');
end


function filepath = find_gtv_file(folder, patterns, index, pat_name, fx_name)
    % FIND_GTV_FILE Helper function to locate GTV mask files
    %   Handles multiple naming conventions and repeat indices

    if ischar(patterns) || isstring(patterns)
        patterns = {char(patterns)};
    end

    filepath = '';
    gtv_search = [];
    single_gtv_search = [];

    % Search for all possible patterns
    for p = 1:length(patterns)
        pat = patterns{p};
        % Specific repeat index search
        gtv_search = cat(1, gtv_search, dir(fullfile(folder, [pat int2str(index) '*.mat'])));
        % General search (no index)
        single_gtv_search = cat(1, single_gtv_search, dir(fullfile(folder, [pat '*.mat'])));
    end

    % If exactly one general GTV mask exists, use it directly
    if isscalar(single_gtv_search)
        filepath = fullfile(folder, single_gtv_search.name);
        return;
    end

    if isempty(gtv_search)
        fprintf('%s/%s: No GTV%d found\n', pat_name, fx_name, index);
        return;
    end

    % If exactly one specific index mask exists, use it directly
    if isscalar(gtv_search)
        filepath = fullfile(folder, gtv_search.name);
        return;
    end

    % If multiple specific GTV masks exist, try to match by repeat index (index)
    % using the second underscore-delimited token in the filename
    gtv_names = {gtv_search.name};
    gtv_search_result = zeros(size(gtv_search));

    for gi = 1:length(gtv_names)
        gtmp = strsplit(gtv_names{gi}, '_');
        if length(gtmp) >= 2
            gtmp_tok = gtmp{2};
            for p = 1:length(patterns)
                pat = patterns{p};
                regex_pattern = regexptranslate('wildcard', [pat int2str(index)]);
                isfound = regexp(gtmp_tok, regex_pattern);
                if ~isempty(isfound) && isfound(1) == 1
                    gtv_search_result(gi) = 1;
                    break;
                end
            end
        end
    end

    % Fallback: a second format was also used for some cases...
    % (repeat index as trailing numeric token after last '_')
    if sum(gtv_search_result) == 0
        for gi = 1:length(gtv_names)
            gtmp = strsplit(gtv_names{gi}, '_');
            gtmp_end = strrep(gtmp{end}, '.mat', '');
            if str2double(gtmp_end) == index
                gtv_search_result(gi) = 1;
            end
        end
    end

    if sum(gtv_search_result) == 0
        fprintf('%s/%s: No GTV%d found\n', pat_name, fx_name, index);
    elseif sum(gtv_search_result) > 1
        fprintf('%s/%s: Redundant GTV%ds found\n', pat_name, fx_name, index);
    elseif sum(gtv_search_result) == 1
        filepath = fullfile(folder, gtv_search(gtv_search_result == 1).name);
    end
end


function [d_map, f_map, dstar_map] = fit_ivim_model(dwi_data, bvalues, mask, opts)
    % Helper function to fit segmented IVIM model to valid voxels within a mask
    sz3 = [size(dwi_data,1), size(dwi_data,2), size(dwi_data,3)];
    valid_voxels_idx = find(mask);
    n_valid = length(valid_voxels_idx);

    % Preallocate output 1D arrays
    d_vec = nan(n_valid, 1);
    f_vec = nan(n_valid, 1);
    dstar_vec = nan(n_valid, 1);

    has_sufficient_bvalues = sum(bvalues >= opts.bthr) >= 2;
    if has_sufficient_bvalues && n_valid > 0
        % Extract 1D signal decay curves for valid voxels
        dwi_flat = reshape(dwi_data, [prod(sz3), length(bvalues)]);
        dwi_valid = dwi_flat(valid_voxels_idx, :);

        % Pad to even number of elements
        pad_len = mod(2 - mod(n_valid, 2), 2);
        dwi_valid_padded = [dwi_valid; zeros(pad_len, length(bvalues))];
        n_padded = n_valid + pad_len;

        % Reshape for dependency [N, 1, 2, bval]
        dwi_1d_vol = reshape(dwi_valid_padded, [n_padded/2, 1, 2, length(bvalues)]);
        mask_1d_vol = true(n_padded/2, 1, 2);

        % Execute the untouched dependency on the flattened array
        ivim_fit_1d = IVIMmodelfit(dwi_1d_vol, bvalues, "seg", mask_1d_vol, opts);

        % Restructure output back to strictly 1D and snip padding
        ivim_out_flat = reshape(ivim_fit_1d, [n_padded, 4]);
        d_vec = squeeze(ivim_out_flat(1:n_valid, 1));
        f_vec = squeeze(ivim_out_flat(1:n_valid, 3));
        dstar_vec = squeeze(ivim_out_flat(1:n_valid, 4));

        % Replace zero-fit voxels with NaN
        zero_mask = (d_vec == 0);
        d_vec(zero_mask) = nan;
        f_vec(zero_mask) = nan;
        dstar_vec(zero_mask) = nan;
    end

    % Safely reconstruct 1D outputs back into 3D volume geometry
    d_map = nan(sz3);
    f_map = nan(sz3);
    dstar_map = nan(sz3);

    if n_valid > 0
        d_map(valid_voxels_idx) = d_vec;
        f_map(valid_voxels_idx) = f_vec;
        dstar_map(valid_voxels_idx) = dstar_vec;
    end
end % function fit_ivim_model

function [have_mask, mask_data] = load_mask(filepath, dwi_size, message_prefix, mask_name)
    % LOAD_MASK Helper function to load a NIfTI mask and validate dimensions
    have_mask = 0;
    mask_data = [];

    if exist(filepath, 'file')
        info = niftiinfo(filepath);
        mask_data = rot90(niftiread(info));
        fprintf('...%sLoaded %s\n', message_prefix, filepath);
        have_mask = 1;

        mask_size = size(mask_data);
        if sum(mask_size ~= dwi_size(1:3)) > 0
            have_mask = 0;
            mask_data = [];
            fprintf('size mismatch. excluding %s\n', mask_name);
        end
    end
end % function [have_mask, mask_data] = load_mask(...)
