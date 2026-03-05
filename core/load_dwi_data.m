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
if isfield(config_struct, 'clinical_sheet_name')
    sheet_name = config_struct.clinical_sheet_name;
else
    sheet_name = 'Clin List_MR';
end
try
    T = readtable(clinical_data_sheet,'Sheet',sheet_name);
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

% Scan for existing checkpoints to identify completed patients.
% A patient is considered complete ONLY if its checkpoint .mat exists AND
% no .lock sentinel is present (which would indicate an in-progress or
% interrupted write).  This prevents a race where Worker A is still writing
% Patient N's checkpoint while Worker B sees a partial file and skips it.
patient_completed = false(size(mrn_list));
for j = 1:length(mrn_list)
    mrn = mrn_list{j};
    patient_id = id_list{j};
    checkpoint_file = fullfile(checkpoint_dir, sprintf('patient_%03d_%s.mat', j, patient_id));
    lock_file = fullfile(checkpoint_dir, sprintf('patient_%03d_%s.lock', j, patient_id));
    if exist(checkpoint_file, 'file') && ~exist(lock_file, 'file')
        patient_completed(j) = true;
    elseif exist(lock_file, 'file')
        % A .lock file without a completed checkpoint indicates a
        % previously interrupted run.  Delete the stale lock so this
        % patient gets re-processed.
        delete(lock_file);
        if exist(checkpoint_file, 'file')
            delete(checkpoint_file);
            fprintf('Removed stale lock + partial checkpoint for patient %d (%s).\n', j, patient_id);
        end
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
        % Strip leading/trailing single quotes that Excel may embed in text cells
        T_Pat_cell = strrep(T_Pat_cell, '''', '');
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

% --- DEBUG: print spreadsheet vs folder patient IDs for matching diagnosis ---
fprintf('\n--- DEBUG: Clinical spreadsheet Pat column (first 5) ---\n');
for dbg_i = 1:min(5, numel(T_Pat_normalized))
    fprintf('  Spreadsheet[%d]: "%s"\n', dbg_i, T_Pat_normalized{dbg_i});
end
fprintf('--- DEBUG: Folder id_list (first 5) ---\n');
for dbg_i = 1:min(5, numel(id_list_normalized))
    fprintf('  Folder[%d]: "%s"\n', dbg_i, id_list_normalized{dbg_i});
end
fprintf('--- DEBUG: Total spreadsheet entries: %d, Total folder entries: %d ---\n\n', numel(T_Pat_normalized), numel(id_list_normalized));

% --- Progress bar for parallel patient processing ---
dq_progress = [];  % initialise so parfor sees the variable regardless of branch
n_to_process = sum(~patient_completed);
if ~exist('OCTAVE_VERSION', 'builtin') && n_to_process > 0
    dq_progress = parallel.pool.DataQueue;
    afterEach(dq_progress, parfor_progress(n_to_process, 'Processing patients'));
end

parfor j = 1:length(mrn_list)
    mrn = mrn_list{j};
    patient_id = id_list{j};

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
    n_fx = size(dwi_locations,2);
    n_rp = size(dwi_locations,3);
    [pat_data_vectors_gtvp, pat_data_vectors_gtvn] = init_scan_structs(n_fx, n_rp);

    % Per-patient DIR reference: populated at Fx1, reused at Fx2+
    b0_fx1_ref        = [];   % b=0 volume at baseline fraction
    gtv_mask_fx1_ref  = [];   % GTVp mask at baseline fraction
    gtvn_mask_fx1_ref = [];   % GTVn mask at baseline fraction (when present)

    % read clinical data (local failure, immunotherapy) from spreadsheet
    % Use exact matching (strcmp) to prevent substring collisions (e.g. "P12" matching "P123").
    i_pat = find(strcmp(T_Pat_normalized, id_list_normalized{j}));
    if isempty(i_pat)
        warning('load_dwi_data:patientNotFound', 'Patient ID ''%s'' not found in clinical spreadsheet. Skipping.', id_list{j});
        continue;
    end
    if numel(i_pat) > 1
        warning('load_dwi_data:duplicatePatient', 'Patient ID ''%s'' has %d matches in clinical spreadsheet. Using first.', id_list{j}, numel(i_pat));
    end
    pat_immuno = T.Immuno(i_pat(1));
    pat_lf = T.LF(i_pat(1));

    basefolder = fullfile(dataloc, id_list{j});
    basefolder_contents = clean_dir_command(basefolder);

    % --- Loop over fractions (fi) and repeat acquisitions (rpi) ---
    n_rp_dim = size(dwi_locations, 3);
    for fi=1:size(dwi_locations,2)
        % Collect DVH metrics across repeats for correct averaging
        dvh_dmean_gtvp_rp  = nan(1, n_rp_dim);
        dvh_dmean_gtvn_rp  = nan(1, n_rp_dim);
        dvh_d95_gtvp_rp    = nan(1, n_rp_dim);
        dvh_d95_gtvn_rp    = nan(1, n_rp_dim);
        dvh_v50gy_gtvp_rp  = nan(1, n_rp_dim);
        dvh_v50gy_gtvn_rp  = nan(1, n_rp_dim);
        % Pre-filter fraction folder for this fraction to avoid redundant dir() calls
        fxtmp_idx = ~cellfun(@isempty, strfind({basefolder_contents.name}, fx_search{fi}));
        fxtmp = basefolder_contents(fxtmp_idx);

        for rpi = 1:size(dwi_locations,3)

            if isempty(fxtmp)
                fprintf('%s, no %s folder\n',id_list{j},fx_search{fi});
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

            % Build scan context for this fraction × repeat
            scan_ctx = struct();
            scan_ctx.fi = fi;
            scan_ctx.rpi = rpi;
            scan_ctx.dicomloc = dicomloc;
            scan_ctx.struct_file = struct_file;
            scan_ctx.struct_file_gtvn = struct_file_gtvn;
            scan_ctx.dicomdoseloc = dicomdoseloc;
            scan_ctx.basefolder = basefolder;
            scan_ctx.dataloc = dataloc;
            scan_ctx.patient_id = patient_id;
            scan_ctx.id_j = id_list{j};
            scan_ctx.mrn_j = mrn_list{j};
            scan_ctx.pat_lf = pat_lf;
            scan_ctx.pat_immuno = pat_immuno;
            scan_ctx.dcm2nii_call = dcm2nii_call;
            scan_ctx.ivim_bthr = ivim_bthr;
            scan_ctx.n_rtdose_cols = size(rtdose_locations,2);
            scan_ctx.b0_fx1_ref = b0_fx1_ref;
            scan_ctx.gtv_mask_fx1_ref = gtv_mask_fx1_ref;
            scan_ctx.gtvn_mask_fx1_ref = gtvn_mask_fx1_ref;

            % Process this scan and collect results
            [scan_result, b0_ref_upd, gtvp_ref_upd, gtvn_ref_upd] = ...
                process_single_scan(scan_ctx);

            % Update Fx1 reference volumes when returned
            if ~isempty(b0_ref_upd),   b0_fx1_ref = b0_ref_upd;           end
            if ~isempty(gtvp_ref_upd), gtv_mask_fx1_ref = gtvp_ref_upd;   end
            if ~isempty(gtvn_ref_upd), gtvn_mask_fx1_ref = gtvn_ref_upd;  end

            % Warn if baseline reference could not be established at Fx1.
            % Subsequent fractions will silently skip DIR registration,
            % producing spatially misaligned parameter and dose maps.
            if fi == 1 && isempty(b0_fx1_ref)
                warning('load_dwi_data:missingFx1Reference', ...
                    'Baseline b0 reference not available for patient %s. DIR registration will be skipped for all subsequent fractions.', id_list{j});
            end

            % Collect bad DWI locations
            for bi_bad = 1:length(scan_result.bad_dwi_list)
                bad_dwi_idx_j = bad_dwi_idx_j + 1;
                bad_dwi_list_j{bad_dwi_idx_j} = scan_result.bad_dwi_list{bi_bad};
            end

            % Assign scan results back to patient-level arrays
            pat_data_vectors_gtvp(fi,rpi) = scan_result.data_gtvp;
            pat_data_vectors_gtvn(fi,rpi) = scan_result.data_gtvn;
            pat_adc_mean(1,fi,rpi) = scan_result.adc_mean;
            pat_adc_kurtosis(1,fi,rpi) = scan_result.adc_kurtosis;
            pat_d_mean(1,fi,rpi) = scan_result.d_mean;
            pat_d_kurtosis(1,fi,rpi) = scan_result.d_kurtosis;
            pat_d_mean_dncnn(1,fi,rpi) = scan_result.d_mean_dncnn;
            pat_d_mean_ivimnet(1,fi,rpi) = scan_result.d_mean_ivimnet;
            % Collect DVH metrics for this repeat (averaged after loop)
            dvh_dmean_gtvp_rp(rpi)  = scan_result.dmean_gtvp;
            dvh_dmean_gtvn_rp(rpi)  = scan_result.dmean_gtvn;
            dvh_d95_gtvp_rp(rpi)    = scan_result.d95_gtvp;
            dvh_d95_gtvn_rp(rpi)    = scan_result.d95_gtvn;
            dvh_v50gy_gtvp_rp(rpi)  = scan_result.v50gy_gtvp;
            dvh_v50gy_gtvn_rp(rpi)  = scan_result.v50gy_gtvn;
        end
        % Average DVH metrics across repeats (nanmean ignores missing repeats)
        pat_dmean_gtvp(1,fi)  = nanmean(dvh_dmean_gtvp_rp);
        pat_dmean_gtvn(1,fi)  = nanmean(dvh_dmean_gtvn_rp);
        pat_d95_gtvp(1,fi)    = nanmean(dvh_d95_gtvp_rp);
        pat_d95_gtvn(1,fi)    = nanmean(dvh_d95_gtvn_rp);
        pat_v50gy_gtvp(1,fi)  = nanmean(dvh_v50gy_gtvp_rp);
        pat_v50gy_gtvn(1,fi)  = nanmean(dvh_v50gy_gtvn_rp);
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

    % Save checkpoint with lock-file protection to prevent race conditions
    % in parfor: the .lock file signals that a write is in progress.
    checkpoint_file = fullfile(checkpoint_dir, sprintf('patient_%03d_%s.mat', j, patient_id));
    lock_file = fullfile(checkpoint_dir, sprintf('patient_%03d_%s.lock', j, patient_id));
    parsave_checkpoint(checkpoint_file, pat_data_out, lock_file);

    % Signal progress to the client-side progress bar
    if ~exist('OCTAVE_VERSION', 'builtin') && ~isempty(dq_progress)
        send(dq_progress, j);
    end

    fprintf('Finished processing patient %d/%d (MRN: %s)\n', j, length(mrn_list), mrn);
end

% Reconstruct global arrays from checkpoints
n_reconstruct = length(mrn_list);
for j = 1:n_reconstruct
    text_progress_bar(j, n_reconstruct, 'Reconstructing checkpoints');
    mrn = mrn_list{j};
    patient_id = id_list{j};
    checkpoint_file = fullfile(checkpoint_dir, sprintf('patient_%03d_%s.mat', j, patient_id));

    if exist(checkpoint_file, 'file')
        % Load checkpoint with basic corruption detection
        loaded_data = load(checkpoint_file);

        required_fields = {'data_vectors_gtvp', 'data_vectors_gtvn', ...
            'dmean_gtvp', 'dmean_gtvn', 'd95_gtvp', 'd95_gtvn', ...
            'v50gy_gtvp', 'v50gy_gtvn', 'adc_mean', 'd_mean', ...
            'd_mean_dncnn', 'd_mean_ivimnet', 'lf', 'immuno', 'bad_dwi_list'};
        missing_fields = setdiff(required_fields, fieldnames(loaded_data));
        if ~isempty(missing_fields)
            warning('load_dwi_data:corruptCheckpoint', ...
                'Checkpoint for patient %d (%s) is missing fields: %s. Skipping.', ...
                j, patient_id, strjoin(missing_fields, ', '));
            continue;
        end

        % Validate dimensions of checkpoint arrays against expected sizes
        expected_n_fx = size(dwi_locations, 2);
        expected_n_rp = size(dwi_locations, 3);
        cp_adc_size = size(loaded_data.adc_mean);
        if length(cp_adc_size) < 2 || cp_adc_size(2) ~= expected_n_fx
            warning('load_dwi_data:checkpointSizeMismatch', ...
                'Checkpoint for patient %d (%s): adc_mean has %d fractions, expected %d. Skipping.', ...
                j, patient_id, cp_adc_size(2), expected_n_fx);
            continue;
        end
        if length(cp_adc_size) >= 3 && cp_adc_size(3) ~= expected_n_rp
            warning('load_dwi_data:checkpointSizeMismatch', ...
                'Checkpoint for patient %d (%s): adc_mean has %d repeats, expected %d. Skipping.', ...
                j, patient_id, cp_adc_size(3), expected_n_rp);
            continue;
        end

        % Assign back to global arrays
        % Struct arrays
        data_vectors_gtvp = align_and_assign_struct(data_vectors_gtvp, loaded_data.data_vectors_gtvp, j);
        data_vectors_gtvn = align_and_assign_struct(data_vectors_gtvn, loaded_data.data_vectors_gtvn, j);

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
    newfilename = fullfile(dataloc, ['dwi_vectors_', dateString, '.mat']);
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
if ~exist(datasave, 'file') && ~isempty(file_prefix)
    % Fallback to the default (un-typed) file when the variant-specific
    % file has not been created yet (e.g. first run before per-type saves).
    % NOTE: run_dwi_pipeline.m (reload branch, ~line 288) restricts this
    % fallback to Standard (dtype==1) only, preventing cross-type
    % contamination.  This path is only reached during the initial 'load'
    % step or direct calls outside the orchestrator.
    datasave_fallback = fullfile(dataloc, 'dwi_vectors.mat');
    if exist(datasave_fallback, 'file')
        fprintf('💡 %s not found — falling back to %s\n', ...
            ['dwi_vectors' file_prefix '.mat'], 'dwi_vectors.mat');
        datasave = datasave_fallback;
    end
end
if ~exist(datasave, 'file')
    type_label = '';
    if isfield(config_struct, 'dwi_type_name')
        type_label = config_struct.dwi_type_name;
    end
    error('load_dwi_data:fileNotFound', ...
        'Required data file ''%s'' not found. Run the load step for DWI type ''%s'' before reloading.', ...
        datasave, type_label);
end
tmp_data = load(datasave);
data_vectors_gtvn = tmp_data.data_vectors_gtvn; data_vectors_gtvp = tmp_data.data_vectors_gtvp; lf = tmp_data.lf;
immuno = tmp_data.immuno; mrn_list = tmp_data.mrn_list; id_list = tmp_data.id_list; fx_dates = tmp_data.fx_dates;
dwi_locations = tmp_data.dwi_locations; rtdose_locations = tmp_data.rtdose_locations; gtv_locations = tmp_data.gtv_locations;
gtvn_locations = tmp_data.gtvn_locations; dmean_gtvp = tmp_data.dmean_gtvp; dmean_gtvn = tmp_data.dmean_gtvn;
d95_gtvp = tmp_data.d95_gtvp; d95_gtvn = tmp_data.d95_gtvn; v50gy_gtvp = tmp_data.v50gy_gtvp; v50gy_gtvn = tmp_data.v50gy_gtvn;
bad_dwi_locations = tmp_data.bad_dwi_locations; bad_dwi_count = tmp_data.bad_dwi_count;
    
if exist('OCTAVE_VERSION', 'builtin') && ~exist('id_list', 'var')
    warning('id_list not loaded from save file. This may occur during mock tests. Proceeding with dummy data.');
    id_list = {}; mrn_list = {}; lf = []; immuno = {}; gtv_locations = []; dwi_locations = []; dmean_gtvp = []; d95_gtvp = []; v50gy_gtvp = []; data_vectors_gtvp = []; data_vectors_gtvn = [];
end
end % if ~skip_to_reload

%% ========================================================================
fprintf('\n--- SECTION 5: Longitudinal Summary Metrics ---\n');
%  SECTION 5 — LONGITUDINAL SUMMARY METRICS

summary_metrics = compute_summary_metrics(config_struct, data_vectors_gtvp, id_list, mrn_list, lf, immuno, gtv_locations, dwi_locations, dmean_gtvp, d95_gtvp, v50gy_gtvp);

end

function parsave_checkpoint(fname, data, lock_file)
    % Create lock sentinel BEFORE writing to prevent race conditions.
    % The lock is removed only after the .mat write completes successfully.
    if nargin >= 3 && ~isempty(lock_file)
        fid = fopen(lock_file, 'w');
        if fid > 0, fclose(fid); end
    end
    save(fname, '-struct', 'data');
    if nargin >= 3 && ~isempty(lock_file) && exist(lock_file, 'file')
        delete(lock_file);
    end
end

function global_struct = align_and_assign_struct(global_struct, new_struct, index)
    % ALIGN_AND_ASSIGN_STRUCT Helper to assign struct arrays with potentially missing fields

    % Per-patient checkpoint data is stored as nFx × nRp (no patient dim).
    % Reshape to 1 × nFx × nRp so it can be slotted into the global
    % nPatients × nFx × nRp array at (index, :, :).
    if ndims(new_struct) <= 2 && size(new_struct, 1) > 1
        new_struct = reshape(new_struct, [1, size(new_struct)]);
    end

    if isempty(fieldnames(global_struct))
        % Initialise global struct: create a matching-fields template so
        % MATLAB can perform subscripted assignment at any index.
        fields = fieldnames(new_struct);
        empty_vals = repmat({[]}, numel(fields), 1);
        template = cell2struct(empty_vals, fields, 1);
        sz_new = size(new_struct);
        dims = [index, sz_new(2:end)];
        global_struct = repmat(template, dims);
        global_struct(index, :, :) = new_struct;
        return;
    end

    fields_global = fieldnames(global_struct);
    fields_new = fieldnames(new_struct);

    % Add any fields that exist in global_struct but are missing in new_struct
    missing_in_new = setdiff(fields_global, fields_new);
    for i = 1:length(missing_in_new)
        [new_struct.(missing_in_new{i})] = deal([]);
    end

    % Add any fields that exist in new_struct but are missing in global_struct
    missing_in_global = setdiff(fields_new, fields_global);
    for i = 1:length(missing_in_global)
        [global_struct.(missing_in_global{i})] = deal([]);
    end

    % Order the new struct fields to match global_struct
    new_struct = orderfields(new_struct, global_struct);

    % Perform the assignment safely
    global_struct(index, :, :) = new_struct;
end
