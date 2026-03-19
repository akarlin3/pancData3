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
% ANALYTICAL RATIONALE
%   Pancreatic DWI analysis requires fitting physical diffusion models
%   (IVIM biexponential and ADC monoexponential) to multi-b-value MRI
%   signal decay curves at every voxel within the tumor volume. The IVIM
%   model separates true tissue diffusivity (D, reflecting cellularity)
%   from pseudo-diffusion (D*, reflecting capillary blood flow) and the
%   perfusion fraction (f, reflecting microvasculature). These biomarkers
%   are hypothesized to be early predictors of treatment response in
%   pancreatic cancer radiotherapy — changes in cellularity and perfusion
%   during fractionated RT may precede volumetric tumor shrinkage by weeks.
%
%   The pipeline processes data longitudinally across treatment fractions
%   (Fx1 through Fx5 + post-RT) because the temporal trajectory of
%   diffusion parameters is as clinically informative as baseline values.
%   Deformable image registration (DIR) aligns all timepoints to the
%   baseline (Fx1) anatomy so that voxel-level longitudinal comparisons
%   are spatially valid despite inter-fraction patient setup variation
%   and tumor deformation during treatment.
%
%   Three parallel DWI processing pipelines are maintained (Standard,
%   DnCNN-denoised, IVIMnet) to enable methodological comparison: DnCNN
%   reduces Rician noise in DWI images before conventional IVIM fitting,
%   while IVIMnet uses a neural network to estimate IVIM parameters
%   directly, potentially with better noise robustness for the unstable
%   D* parameter.
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
% This two-phase design (compute-once + reload-many) is essential because
% DICOM conversion + model fitting for a typical 30-patient cohort with
% 6 fractions each can take 6-12 hours, while summary metric computation
% takes only minutes.
skip_to_reload = config_struct.skip_to_reload;

% b-value threshold (s/mm²) used for segmented IVIM fitting.
% The IVIM model S(b) = S0*[f*exp(-b*D*) + (1-f)*exp(-b*D)] has two
% exponential components. At low b-values (< ~100 s/mm²), the
% pseudo-diffusion term (D*, typically 10-30x larger than D) dominates
% signal decay, reflecting intravascular blood flow. At high b-values
% (>= ~100 s/mm²), the perfusion signal has decayed away and the
% remaining signal reflects true tissue diffusion (D).
%
% The segmented fitting approach exploits this separation:
%   Stage 1: Fit ln(S) vs b using only high b-values to estimate D
%   Stage 2: With D fixed, fit the full model to estimate f and D*
% This avoids the numerical instability of fitting all 3 parameters
% simultaneously, which is particularly problematic for D* due to its
% high variance in pancreatic tissue.
%
% Default threshold = 100 s/mm². For typical pancreatic DWI protocols
% with b = [0, 30, 100, 150, 550], this places b=0 and b=30 in the
% perfusion-sensitive low set, and b=100, 150, 550 in the diffusion-
% dominated high set.
ivim_bthr = config_struct.ivim_bthr;
use_gpu = config_struct.use_gpu;

% When GPU is requested, validate availability once at startup rather than
% per-scan to avoid redundant checks inside the parfor loop.
if use_gpu
    [gpu_ok, ~] = gpu_available(config_struct.gpu_device);
    if ~gpu_ok
        fprintf('  ⚠️ use_gpu is true but no suitable GPU found. Falling back to CPU.\n');
        use_gpu = false;
    end
end

%% ========================================================================
% Pre-check: if skip_to_reload is requested for a specific DWI type
% (dnCNN, IVIMnet) but the type-specific .mat file does not exist yet
% (e.g., first run from execute_all_workflows), fall through to the full
% load path instead of erroring out later.  This only applies when
% dwi_type_name is set — for untyped/Standard runs, a missing file with
% skip_to_reload=true is still a genuine error (the user explicitly
% requested reload).
if skip_to_reload && isfield(config_struct, 'dwi_type_name') && ~isempty(config_struct.dwi_type_name)
    dataloc_check = config_struct.dataloc;
    file_prefix_check = ['_' config_struct.dwi_type_name];
    datasave_check = fullfile(dataloc_check, ['dwi_vectors' file_prefix_check '.mat']);
    if ~exist(datasave_check, 'file')
        % Before falling through to full load, check if the default
        % (un-typed) dwi_vectors.mat exists as a valid fallback.  The
        % Section 4 reload path already handles this fallback, so we only
        % need to set skip_to_reload = false when neither file exists.
        datasave_fallback_check = fullfile(dataloc_check, 'dwi_vectors.mat');
        if exist(datasave_fallback_check, 'file')
            fprintf('  💡 %s not found — will fallback to dwi_vectors.mat for %s.\n', ...
                ['dwi_vectors' file_prefix_check '.mat'], config_struct.dwi_type_name);
        else
            fprintf('  💡 skip_to_reload requested but %s not found. Running full load for %s.\n', ...
                ['dwi_vectors' file_prefix_check '.mat'], config_struct.dwi_type_name);
            skip_to_reload = false;
        end
    end
end

fprintf('\n--- SECTION 1: File Discovery ---\n');
%  SECTION 1 — FILE DISCOVERY

if ~skip_to_reload

% Discover all patient imaging data on the network share. This traverses
% the directory hierarchy to locate DWI DICOM folders, GTV contour masks,
% RT dose files, and nodal GTV masks for every patient x fraction x repeat
% combination. The structured naming convention (P##/Fx#/DWI#/) allows
% automated discovery without a manual lookup table.
[id_list, mrn_list, fx_dates, dwi_locations, rtdose_locations, gtv_locations, gtvn_locations] = discover_patient_files(config_struct.dataloc);

% Optional patient subset filtering: when config specifies patient_ids,
% restrict processing to only those patients. This is useful for debugging
% individual cases or re-processing specific patients after data corrections.
if isfield(config_struct, 'patient_ids') && ~isempty(config_struct.patient_ids)
    keep_idx = find(ismember(mrn_list, config_struct.patient_ids));
    if ~isempty(keep_idx)
        id_list = id_list(keep_idx);
        mrn_list = mrn_list(keep_idx);
        fx_dates = fx_dates(keep_idx, :);
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

% Load clinical outcome data (local failure status and immunotherapy use)
% from a curated spreadsheet. Local failure (LF) is the primary endpoint
% for correlating diffusion biomarkers with treatment response. Immuno-
% therapy status is tracked as a potential confounder because checkpoint
% inhibitors can alter tumor microenvironment independently of RT.
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

% Fraction labels matching the clinical RT schedule. Pancreatic SBRT
% typically delivers 5 fractions (Fx1-Fx5) over 1-2 weeks, with a
% post-treatment MRI acquired 4-6 weeks after completion. Each fraction
% has a DWI acquisition for longitudinal biomarker tracking.
fx_search = {'Fx1','Fx2','Fx3','Fx4','Fx5','post'};

% Check available system memory before processing
[mem_gb_total, mem_gb_avail] = get_system_memory();
fprintf('System memory: %.1f GB total, %.1f GB available\n', mem_gb_total, mem_gb_avail);

% Initialise output struct arrays for GTVp (primary tumor) and GTVn
% (nodal disease). Pancreatic tumors may have both a primary pancreatic
% mass (GTVp) and involved lymph nodes (GTVn), which can show different
% diffusion characteristics and treatment response patterns.
data_vectors_gtvp = struct;
data_vectors_gtvn = struct;

% Pre-allocate summary metric arrays (patient x fraction x repeat).
% NaN initialization ensures that missing data (e.g., patient without Fx3
% scan) is naturally handled by nanmean/nanstd in downstream analysis
% rather than corrupting calculations with zeros.
adc_mean = nan(size(dwi_locations));
adc_kurtosis = nan(size(dwi_locations));

% IVIM true diffusion coefficient D — mean values per pipeline variant.
% Tracking D separately from ADC is important because ADC conflates true
% diffusion with pseudo-diffusion contributions, especially at low b-values.
% D isolates the tissue diffusivity component, which more directly reflects
% cellularity changes during treatment.
d_mean = nan(size(dwi_locations));         % standard IVIM fit
d_mean_dncnn = nan(size(dwi_locations));   % DnCNN-denoised IVIM fit
d_mean_ivimnet = nan(size(dwi_locations)); % IVIMnet deep-learning fit

d_kurtosis = nan(size(dwi_locations));

% DVH (Dose-Volume Histogram) parameters within GTVp and GTVn.
% These dose metrics characterize the spatial distribution of radiation
% dose within the tumor contour and are essential for correlating
% diffusion changes with delivered dose — the fundamental dose-response
% analysis in this study.
% Sized to 6 columns to accommodate on-treatment fractions + Post-RT scan
dmean_gtvp = nan(length(mrn_list), 6);  % mean dose in GTVp (Gy) — overall dose intensity
dmean_gtvn = nan(length(mrn_list), 6);  % mean dose in GTVn (Gy)

d95_gtvp = nan(length(mrn_list), 6);    % D95% — dose received by 95% of GTV (Gy)
d95_gtvn = nan(length(mrn_list), 6);    % D95 is a coverage metric: low D95 indicates
                                         % cold spots within the tumor that may correlate
                                         % with local failure

v50gy_gtvp = nan(length(mrn_list), 6);  % V50Gy — fraction of GTV receiving >=50 Gy
v50gy_gtvn = nan(length(mrn_list), 6);  % V50Gy captures high-dose coverage relevant
                                         % to dose escalation studies

% Clinical outcome flags (per patient)
lf = zeros(size(mrn_list));      % local failure (1 = yes) — primary endpoint
immuno = zeros(size(mrn_list));  % received immunotherapy (1 = yes) — potential confounder

% Track problematic DWI acquisitions for manual review. DWI artifacts
% (motion, geometric distortion, incomplete acquisitions) are common in
% pancreatic imaging due to respiratory motion and bowel peristalsis.
% Flagging these allows the physicist to review and decide whether to
% exclude the data or apply motion correction.
bad_dwi_locations_per_patient = cell(length(mrn_list), 1);

% Checkpoint directory setup. Per-patient checkpointing is critical because
% processing a full cohort (DICOM conversion + model fitting for ~30
% patients x 6 fractions x multiple repeats) can take many hours. If the
% pipeline is interrupted (crash, timeout, resource limit), completed
% patients are preserved and only unfinished patients are re-processed.
checkpoint_dir = fullfile(dataloc, 'processed_patients');
if ~isfolder(checkpoint_dir)
    mkdir(checkpoint_dir);
    % Write provenance sentinel so cache-clearing knows this directory was
    % created by the pipeline and is safe to delete.
    sent_fid = fopen(fullfile(checkpoint_dir, '.pipeline_created'), 'w');
    if sent_fid > 0, fprintf(sent_fid, 'Created by load_dwi_data\n'); fclose(sent_fid); end
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
% Uses shared normalize_patient_ids utility for Octave-compatible ID matching.
% Always normalize folder IDs (independent of spreadsheet column names)
[~, id_list_normalized] = normalize_patient_ids({}, id_list);
% Find the spreadsheet column containing patient IDs.  The column may be
% named 'Pat', 'Patient', 'PatientID', etc. depending on the spreadsheet
% version.  Try common names in priority order.
T_Pat_normalized = {};
pat_col_candidates = {'Pat', 'Patient', 'PatientID', 'Patient_ID'};
pat_col_found = '';
for pci = 1:numel(pat_col_candidates)
    if isfield(T, pat_col_candidates{pci})
        pat_col_found = pat_col_candidates{pci};
        break;
    end
end
if ~isempty(pat_col_found)
    [T_Pat_normalized, ~] = normalize_patient_ids(T.(pat_col_found), id_list);
else
    % Last resort: try the first column of the table
    col_names = fieldnames(T);
    if ~isempty(col_names)
        first_col = T.(col_names{1});
        if iscell(first_col) || iscategorical(first_col) || ischar(first_col)
            [T_Pat_normalized, ~] = normalize_patient_ids(first_col, id_list);
            fprintf('  ⚠️ No ''Pat'' column found; using first column ''%s'' for patient matching.\n', col_names{1});
        end
    end
    if isempty(T_Pat_normalized)
        error('load_dwi_data:noPatColumn', ...
            'Clinical spreadsheet has no recognizable patient ID column. Cannot match patients to clinical data. Available columns: %s', ...
            strjoin(fieldnames(T), ', '));
    end
end

% --- DEBUG: print spreadsheet vs folder patient IDs for matching diagnosis ---
fprintf('\n--- DEBUG: Spreadsheet columns: %s ---\n', strjoin(fieldnames(T), ', '));
if ~isempty(pat_col_found)
    fprintf('--- DEBUG: Using column ''%s'' for patient matching ---\n', pat_col_found);
end
fprintf('--- DEBUG: Clinical spreadsheet Pat column (first 5) ---\n');
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
    n_fx = size(dwi_locations, 2);
    n_rp = size(dwi_locations, 3);
    pat_data_vectors_gtvp = struct;
    pat_data_vectors_gtvn = struct;
    pat_dmean_gtvp = nan(1, n_fx);
    pat_dmean_gtvn = nan(1, n_fx);
    pat_d95_gtvp = nan(1, n_fx);
    pat_d95_gtvn = nan(1, n_fx);
    pat_v50gy_gtvp = nan(1, n_fx);
    pat_v50gy_gtvn = nan(1, n_fx);
    pat_adc_mean = nan(1, n_fx, n_rp);
    pat_adc_kurtosis = nan(1, n_fx, n_rp);
    pat_d_mean = nan(1, n_fx, n_rp);
    pat_d_kurtosis = nan(1, n_fx, n_rp);
    pat_d_mean_dncnn = nan(1, n_fx, n_rp);
    pat_d_mean_ivimnet = nan(1, n_fx, n_rp);

    % Initialize potential fields to ensure struct consistency
    [pat_data_vectors_gtvp, pat_data_vectors_gtvn] = init_scan_structs(n_fx, n_rp);

    % Per-patient DIR (Deformable Image Registration) reference volumes.
    % Populated at Fx1 (baseline) and reused for all subsequent fractions.
    % The Fx1 b=0 volume serves as the fixed (target) image for DIR because:
    % 1. Fx1 is acquired before any RT-induced changes, providing a clean
    %    anatomical reference
    % 2. All GTV contours are drawn on the Fx1 anatomy by the physician
    % 3. Warping Fx2-Fx5 to Fx1 space enables direct voxel-to-voxel
    %    longitudinal comparison of diffusion parameters
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
    if isfield(T, 'Immuno')
        pat_immuno = T.Immuno(i_pat(1));
    elseif isfield(T, 'IO')
        pat_immuno = T.IO(i_pat(1));
    else
        pat_immuno = 0;
    end
    if isfield(T, 'LF')
        pat_lf = T.LF(i_pat(1));
    elseif isfield(T, 'LocalFailure')
        pat_lf = T.LocalFailure(i_pat(1));
    else
        pat_lf = 0;
    end

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
        % Use regexp with word-boundary to prevent 'Fx1' matching 'Fx10'
        % Case-insensitive to handle 'fx1', 'FX1', 'Fx1', etc.
        fxtmp_idx = ~cellfun(@isempty, regexpi({basefolder_contents.name}, [fx_search{fi} '(\b|$)'], 'once'));
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

            % Check memory requirements before processing this scan
            if ~isempty(dicomloc) && exist(dicomloc, 'dir')