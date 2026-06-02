function [data_vectors_gtvp, data_vectors_gtvn, adc_mean, adc_kurtosis, ...
    d_mean, d_mean_dncnn, d_mean_ivimnet, d_kurtosis, dmean_gtvp, ...
    dmean_gtvn, d95_gtvp, d95_gtvn, v50gy_gtvp, v50gy_gtvn, lf, immuno, ...
    bad_dwi_locations_per_patient, checkpoint_dir, patient_completed] = ...
    load_dwi_data_init_arrays(dwi_locations, mrn_list, id_list, dataloc)
    % LOAD_DWI_DATA_INIT_ARRAYS Pre-allocate cohort output arrays and set up checkpoints.
    %   Extracted verbatim from load_dwi_data.m Section 2. Allocates NaN/zero output
    %   arrays sized to the cohort, creates the per-patient checkpoint directory (with
    %   provenance sentinel), and scans existing checkpoints to mark completed patients
    %   while clearing any stale .lock / partial-checkpoint pairs from interrupted runs.

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
% Sentinel backfill for legacy directories happens inside
% clear_pipeline_cache (called earlier in prepare_pipeline_session) so
% the same clear_cache:true run can sweep them.

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
end
