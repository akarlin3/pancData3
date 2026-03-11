% execute_all_workflows.m
% Automates the /run_data workflow for sequential modeling pipelines
%
% [ANALYTICAL DESIGN RATIONALE — WHY THREE SEQUENTIAL DWI TYPES]:
% Pancreatic DWI analysis compares three signal processing approaches to
% determine which yields the most clinically predictive diffusion parameters:
%
%   1. Standard  — Conventional voxel-wise fitting of raw DWI signal.
%                  This is the baseline reference method used in most published
%                  pancreatic DWI studies. Signal is noisy (low SNR in the
%                  pancreas due to respiratory motion and small organ size),
%                  but results are directly interpretable.
%
%   2. dnCNN     — Deep learning denoising (DnCNN) applied to raw DWI images
%                  BEFORE model fitting. By reducing thermal and physiological
%                  noise, the subsequent IVIM/ADC fit converges more reliably,
%                  particularly for the perfusion fraction (f) and pseudo-
%                  diffusion (D*) parameters, which are notoriously unstable
%                  at low SNR. The hypothesis is that denoised fits produce
%                  more reproducible and prognostically discriminative metrics.
%
%   3. IVIMnet   — Neural network-based IVIM parameter estimation that learns
%                  the signal-to-parameter mapping from training data rather
%                  than iterative nonlinear least-squares fitting. This avoids
%                  local minima in the bi-exponential IVIM cost function and
%                  may better separate D, f, and D* in noisy pancreatic data.
%
% All three types share the same raw DICOM data and GTV masks. The only
% difference is HOW the voxel-wise diffusion parameters are estimated.
% Running them sequentially through the same analytical pipeline (sanity
% checks, metrics, survival analysis) enables a controlled comparison of
% processing methodology impact on clinical outcome prediction.
%
% [WHY SEQUENTIAL AND NOT PARALLEL]:
% The three pipelines share disk I/O to the same patient data directory and
% write to the same config.json. Running them in parallel would create race
% conditions on config.json and competing disk access. Sequential execution
% with skip_to_reload=true on the 2nd and 3rd runs avoids redundant DICOM
% loading (the most expensive I/O step), since the raw images are identical
% across all three methods — only the fitting algorithm differs.

% Resolve the repository root from this script's location on disk, not from
% pwd(), so the script works correctly regardless of the user's current
% working directory (important for cluster/batch job submission).
% This script lives in pipeline/, so go up one level to reach repo root.
pipeline_root = fileparts(mfilename('fullpath'));
repo_root = fileparts(pipeline_root);
config_file = fullfile(repo_root, 'config.json');

% --- 1. SET UP ENVIRONMENT AND POOL (max 2 workers) ---
% [ANALYTICAL RATIONALE — PARALLEL POOL CONFIGURATION]:
% The load_dwi_data step uses parfor to process patients in parallel.
% Each patient's DWI data involves DICOM-to-NIFTI conversion followed by
% voxel-wise IVIM fitting — a compute-intensive nonlinear optimization for
% every voxel within the GTV mask across multiple b-values and fractions.
% Parallelism at the patient level (not voxel level) is chosen because:
%   (a) Each patient's data is independent — no inter-patient dependencies.
%   (b) Per-patient granularity matches the checkpointing unit, so a crashed
%       worker loses at most one patient's results.
%   (c) Only 2 workers are used to avoid memory exhaustion — each worker
%       holds an entire patient's 4D DWI volume (potentially >1 GB) in RAM
%       during IVIM fitting. More workers would exceed typical workstation
%       memory on a 60+ patient pancreatic cohort.
% Parallel pool setup is skipped under GNU Octave, which does not support
% MATLAB's Parallel Computing Toolbox (parpool, parfor, pctRunOnAll).
% Octave compatibility mode runs the pipeline serially instead.
if ~exist('OCTAVE_VERSION', 'builtin')
    % Delete any stale parallel jobs before creating a new pool.
    % Stale jobs can occur when a previous pipeline run was killed (kill -9)
    % or MATLAB crashed during parfor execution, leaving zombie job entries
    % that prevent a new pool from starting cleanly.
    try
        allJobs = matlab.internal.parallel.getJobs();
        if ~isempty(allJobs)
            delete(allJobs);
            fprintf('⚙️ Deleted %d stale parallel job(s).\n', numel(allJobs));
        end
    catch
        % Fallback: use the job manager from the default cluster profile
        try
            c = parcluster;
            if ~isempty(c.Jobs)
                delete(c.Jobs);
                fprintf('⚙️ Deleted stale parallel jobs from cluster profile.\n');
            end
        catch
            % No jobs to clean up
        end
    end

    % Destroy any existing pool to ensure a clean state with exactly 2
    % workers and infinite idle timeout (the pipeline may take hours).
    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);
    end
    p = parpool('Processes', 2, 'IdleTimeout', Inf);
    % Attach the main data-loading function so workers can access it without
    % relying on path resolution, which can differ across workers.
    addAttachedFiles(p, {fullfile(pipeline_root, 'core', 'load_dwi_data.m')});
    % Replicate the path setup on all workers so that utility functions
    % (parse_config, safe_load_mask, escape_shell_arg, etc.) and third-party
    % dependencies (IVIM fitting, DVH tools) are available inside parfor.
    pctRunOnAll addpath(fullfile(pipeline_root, 'core'));
    pctRunOnAll addpath(fullfile(pipeline_root, 'utils'));
    pctRunOnAll addpath(fullfile(pipeline_root, 'dependencies'));
    % Suppress noisy warnings on workers but keep critical ones.
    % Specifically, NIFTI file-not-found warnings are expected when a patient
    % has not yet been converted from DICOM, and DELETE:FileNotFound is
    % benign during checkpoint cleanup. We intentionally do NOT use
    % warning('off','all') because that would silence IVIM convergence
    % warnings and numerical precision warnings that indicate fitting
    % problems — those are diagnostically valuable.
    w_state_before_pct = warning;  % save client warning state before pctRunOnAll clobbers it
    pctRunOnAll warning('off', 'MATLAB:imagesci:niftiinfo:fileDoesNotExist');
    pctRunOnAll warning('off', 'MATLAB:DELETE:FileNotFound');
    warning(w_state_before_pct);   % restore client warnings (pctRunOnAll also executes on client)
end

% Reset persistent output folder in run_dwi_pipeline so each full workflow
% sequence starts with a fresh timestamped output directory.
% run_dwi_pipeline uses a persistent variable (MASTER_OUTPUT_FOLDER) to
% track the output folder across multiple calls within the same session.
% Clearing the function resets that persistent state, ensuring we do not
% accidentally append results to a previous run's output directory.
clear run_dwi_pipeline;

% Create the timestamped output folder now so the diary can live there.
% The timestamp ensures each pipeline execution produces an isolated,
% reproducible output directory — critical for comparing results across
% different analysis runs (e.g., before vs. after reprocessing with
% updated GTV contours or different IVIM fitting parameters).
% run_dwi_pipeline will reuse this folder for all three DWI types.
timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
eaw_output_folder = fullfile(repo_root, sprintf('saved_files_%s', timestamp_str));
if ~exist(eaw_output_folder, 'dir'), mkdir(eaw_output_folder); end
% Write provenance sentinel so cleanup tools know this directory was
% created by the pipeline and is safe to delete.
sentinel = fullfile(eaw_output_folder, '.pipeline_created');
fid = fopen(sentinel, 'w');
if fid > 0, fprintf(fid, 'Created by execute_all_workflows at %s\n', timestamp_str); fclose(fid); end

% --- Workflow Progress GUI (GUI environments only) ---
% Show a top-level progress bar tracking the 3 DWI types (Standard, dnCNN,
% IVIMnet).  On headless/cluster environments, isDisplayAvailable() returns
% false and wfGUI stays empty — all progress is console-only.
% The onCleanup guard ensures the GUI window is closed even on error/Ctrl-C.
wfGUI = [];
if ProgressGUI.isDisplayAvailable()
    wfGUI = ProgressGUI('DWI Analysis Workflow', 3);
end
cleanup_wf_gui = onCleanup(@() closeWfGUI(wfGUI));

% Master diary for execute_all_workflows console output.
% MATLAB only supports one active diary at a time. This master diary
% captures the orchestrator-level output (which DWI type is running,
% test results, errors). Each core module will temporarily override this
% diary with its own file, and the orchestrator restarts it after each
% module returns. This layered diary architecture ensures that both
% high-level progress and per-module detail are captured to separate files.
eaw_diary_file = fullfile(eaw_output_folder, 'execute_all_workflows.log');
diary(eaw_diary_file);

% --- 1.5 RUN TEST SUITE BEFORE PIPELINE ---
% [ANALYTICAL RATIONALE — PRE-PIPELINE TESTING]:
% Running the full test suite before a multi-hour pipeline execution is a
% safeguard against producing scientifically invalid results. The tests
% verify critical invariants including:
%   - Patient-stratified CV folds (no data leakage across patients)
%   - KNN imputation temporal bounds (no future timepoint contamination)
%   - IVIM parameter range plausibility (D ~ 1e-3 mm^2/s, f in [0,1])
%   - Competing risks model correctness (IPCW weights, cause-specific hazards)
% A single broken invariant could silently inflate predictive accuracy or
% produce biologically implausible diffusion parameter maps. It is far
% cheaper to catch this before committing to the full cohort computation.
% Check config for skip_tests option before running the test suite
% Check if the user has opted out of pre-pipeline testing via config.json.
% This is useful during iterative development when tests have already been
% verified and the researcher wants to skip the ~2-minute test suite overhead.
eaw_skip_tests = false;
try
    eaw_raw = fileread(config_file);
    eaw_cfg = jsondecode(eaw_raw);
    if isfield(eaw_cfg, 'skip_tests') && eaw_cfg.skip_tests
        eaw_skip_tests = true;
    end
catch
    % If config can't be read, proceed with tests enabled (safe default)
end

if eaw_skip_tests
    disp('⏭️ Skipping test suite (skip_tests = true in config.json).');
else
    disp('====== RUNNING TEST SUITE BEFORE PIPELINE ======');
    test_diary_file = fullfile(eaw_output_folder, 'test_suite_output.log');
    diary(test_diary_file);
    try
        run(fullfile(pipeline_root, 'tests', 'run_all_tests.m'));
        disp('====== ALL TESTS PASSED — PROCEEDING WITH PIPELINE ======');
        diary(eaw_diary_file);  % switch back to master diary
    catch ME
        try
            diary(eaw_diary_file);  % switch back before error
        catch
            diary off;  % if master diary file is inaccessible, just turn diary off
        end
        fprintf('❌ Test failure: %s\n', ME.message);
        error('PipelineAborted:TestFailure', ...
            'Pipeline aborted: test suite did not pass.');
    end
end

% [CONFIG.JSON MUTATION AND RECOVERY STRATEGY]:
% This script modifies config.json in-place between DWI type runs to change
% dwi_type and skip_to_reload. This is necessary because run_dwi_pipeline
% reads config.json via parse_config, which expects a file path — not an
% in-memory struct. The original string is preserved in memory and restored
% via onCleanup on normal exit, errors, or Ctrl-C.

% Load the raw JSON string.  We modify individual fields via regex
% (json_set_field) instead of jsondecode/jsonencode so that the on-disk
% formatting (indentation, field order, numeric precision) is preserved
% during intermediate writes.  The original string is kept for rollback.
% Read config.json as a raw character string (not a decoded struct) so we
% can modify individual fields with json_set_field (regex-based) while
% preserving the original formatting, comments, and field ordering.
fid = fopen(config_file, 'r');
if fid < 0
    error('execute_all_workflows:fileOpenFailed', ...
        'Cannot open %s for reading. Check file exists and permissions.', config_file);
end
raw = fread(fid, inf);       % Read as uint8 byte vector
str = char(raw');             % Transpose to row vector and convert to char
fclose(fid);
original_config_str = str;    % Preserve original for rollback on exit/error
config_json = str;            % Mutable working copy for field-level edits

% Register an onCleanup handler that restores the original config.json
% content when this script exits — whether by normal completion, error,
% or Ctrl-C.  This is critical because the script mutates config.json
% between DWI type runs (changing dwi_type and skip_to_reload), and the
% user's original settings must be preserved.
restore_config = onCleanup(@() restore_config_file(config_file, original_config_str));

% Execute all 9 discrete target modules.
% [ANALYTICAL STEP ORDERING RATIONALE]:
% The step sequence reflects the scientific analysis dependency chain:
%   load                    -> Raw DICOM to fitted IVIM/ADC parameter maps
%   sanity                  -> Validate parameter plausibility (convergence QC)
%   visualize               -> Generate parameter maps and distribution plots
%   metrics_baseline        -> Extract pre-treatment (Fx1) summary metrics per patient
%   metrics_longitudinal    -> Compute intra-treatment changes (Fx2-Fx5 vs Fx1)
%   metrics_dosimetry       -> Correlate RT dose (D95, V50) with diffusion sub-volumes
%   metrics_stats_comparisons -> Wilcoxon/Mann-Whitney group comparisons (responders vs non)
%   metrics_stats_predictive  -> Predictive modeling (logistic regression, feature selection)
%   metrics_survival        -> Cox regression, competing risks, Kaplan-Meier
%
% Each downstream step depends on outputs from upstream steps. For example,
% metrics_longitudinal needs baseline values from metrics_baseline to compute
% percent change; metrics_stats_predictive needs dose-response features from
% metrics_dosimetry as candidate predictors; metrics_survival needs risk
% scores from metrics_stats_predictive for stratified Kaplan-Meier curves.
steps = {'load', 'sanity', 'visualize', 'metrics_baseline', ...
         'metrics_longitudinal', 'metrics_dosimetry', ...
         'metrics_stats_comparisons', 'metrics_stats_predictive', 'metrics_survival'};

% Conditionally inject compare_cores after metrics_baseline when enabled.
% compare_cores runs all 11 tumor core delineation methods and computes
% pairwise Dice/Hausdorff agreement — a computationally expensive step
% that is off by default.  It must come after metrics_baseline because
% it requires the validated voxel data and summary metrics.
% Note: eaw_cfg was decoded earlier (~line 163) for the skip_tests check.
if exist('eaw_cfg', 'var') && isfield(eaw_cfg, 'run_compare_cores') && eaw_cfg.run_compare_cores
    cc_idx = find(strcmp(steps, 'metrics_baseline'));
    if ~isempty(cc_idx)
        steps = [steps(1:cc_idx), {'compare_cores'}, steps(cc_idx+1:end)];
    else
        steps{end+1} = 'compare_cores';  % Append if metrics_baseline not in steps
    end
end

% --- 2. RUN STANDARD PIPELINE ---
% Standard processing: conventional voxel-wise nonlinear least-squares
% fitting of the raw (undenoised) DWI signal to IVIM and ADC models.
% skip_to_reload = false means this first run performs the full DICOM
% conversion + model fitting pipeline from scratch.
disp('====== STARTING STANDARD PIPELINE ======');
if ~isempty(wfGUI) && wfGUI.isValid()
    counts = struct('completed', 0, 'total', 3, 'stepName', 'DWI Type 1/3: Standard');
    wfGUI.update(0, counts, 'Running Standard pipeline...', 'running');
end
% Set dwi_type to Standard and skip_to_reload to false — this first run
% performs full DICOM conversion + model fitting from raw DWI images.
config_json = json_set_field(config_json, 'dwi_type', 'Standard');
config_json = json_set_field(config_json, 'skip_to_reload', false);
% Write the modified config to disk so run_dwi_pipeline can read it.
fid = fopen(config_file, 'w');
if fid < 0
    error('execute_all_workflows:configWriteFailed', ...
        'Cannot write config.json for Standard run. Check file permissions.');
end
fwrite(fid, config_json); fclose(fid);
% Pass the shared output folder so all 3 DWI types write to the same
% timestamped directory, enabling cross-type comparison by analysis scripts.
run_dwi_pipeline(config_file, steps, eaw_output_folder);
diary(eaw_diary_file);  % restart master diary (module diaries override during run)

% --- 3. RUN dnCNN PIPELINE ---
% DnCNN (Denoising Convolutional Neural Network) processing: the same raw
% DWI images are first passed through a pre-trained DnCNN that removes
% noise while preserving diffusion-weighted contrast. The denoised images
% are then fit with the same IVIM/ADC models as Standard. The hypothesis
% is that noise reduction in the pancreas (inherently low-SNR due to
% respiratory motion artifacts and small tumor volumes) will yield more
% stable and reproducible IVIM parameter estimates, particularly for the
% perfusion fraction f and pseudo-diffusion coefficient D*.
%
% skip_to_reload = true skips DICOM conversion (already done in Standard
% run) and jumps directly to loading cached NIFTI volumes, applying DnCNN
% denoising, and refitting the diffusion models. This saves significant
% time since dcm2niix conversion is identical across all three methods.
disp('====== STARTING dnCNN PIPELINE ======');
if ~isempty(wfGUI) && wfGUI.isValid()
    counts = struct('completed', 1, 'total', 3, 'stepName', 'DWI Type 2/3: dnCNN');
    wfGUI.update(1/3, counts, 'Running dnCNN pipeline...', 'running');
end
config_json = json_set_field(config_json, 'dwi_type', 'dnCNN');
config_json = json_set_field(config_json, 'skip_to_reload', true);
fid = fopen(config_file, 'w');
if fid < 0
    error('execute_all_workflows:configWriteFailed', ...
        'Cannot write config.json for dnCNN run. Check file permissions.');
end
fwrite(fid, config_json); fclose(fid);
run_dwi_pipeline(config_file, steps, eaw_output_folder);
diary(eaw_diary_file);  % restart after pipeline run (module diaries override this)

% --- 4. RUN IVIMnet PIPELINE ---
% IVIMnet processing: instead of iterative nonlinear least-squares fitting,
% a trained neural network directly maps the multi-b-value DWI signal to
% IVIM parameters (D, f, D*). This approach:
%   (a) Avoids local minima in the bi-exponential IVIM cost function — a
%       major issue with conventional fitting where the D* and f parameters
%       are poorly conditioned and highly correlated.
%   (b) Provides orders-of-magnitude faster inference (forward pass vs.
%       iterative optimization per voxel).
%   (c) May regularize parameter estimates through the learned prior from
%       training data, acting as an implicit Bayesian constraint.
%
% The provenance of the IVIMnet training set is tracked by
% load_dl_provenance.m to ensure no analysis patients were in the training
% set — this would constitute a severe form of data leakage that would
% inflate apparent predictive accuracy.
disp('====== STARTING IVIMnet PIPELINE ======');
if ~isempty(wfGUI) && wfGUI.isValid()
    counts = struct('completed', 2, 'total', 3, 'stepName', 'DWI Type 3/3: IVIMnet');
    wfGUI.update(2/3, counts, 'Running IVIMnet pipeline...', 'running');
end
config_json = json_set_field(config_json, 'dwi_type', 'IVIMnet');
config_json = json_set_field(config_json, 'skip_to_reload', true);
fid = fopen(config_file, 'w');
if fid < 0
    error('execute_all_workflows:configWriteFailed', ...
        'Cannot write config.json for IVIMnet run. Check file permissions.');
end
fwrite(fid, config_json); fclose(fid);
run_dwi_pipeline(config_file, steps, eaw_output_folder);
diary(eaw_diary_file);  % restart after pipeline run (module diaries override this)

if ~isempty(wfGUI) && wfGUI.isValid()
    counts = struct('completed', 3, 'total', 3, 'stepName', 'All workflows complete');
    wfGUI.update(1, counts, 'All 3 DWI types processed', 'success');
end
disp('====== ALL WORKFLOWS COMPLETED ======');
diary off;

% Explicitly trigger the onCleanup handlers NOW rather than waiting for
% garbage collection.  In a MATLAB script (not function), onCleanup objects
% persist in the base workspace after the script ends and the cleanup
% callback never fires automatically.  This explicit delete() invokes
% the cleanup callbacks immediately, restoring config.json to its original
% state and closing the GUI window.  Without this, config.json would be
% left with dwi_type="IVIMnet" and skip_to_reload=true until the user
% clears the workspace or exits MATLAB.
delete(restore_config);
delete(cleanup_wf_gui);

function closeWfGUI(gui)
%CLOSEWFGUI  Close the workflow progress GUI if valid.
    if ~isempty(gui)
        try gui.close(); catch, end
    end
end

function restore_config_file(config_path, original_str)
% Restore config.json to its original state when the script exits (whether
% by normal completion or error).  The user's config.json must be left
% exactly as they wrote it after the workflow completes.
    try
        fid = fopen(config_path, 'w');
        if fid == -1
            error('restore:openFailed', 'Cannot open config.json for writing.');
        end
        fwrite(fid, original_str);
        fclose(fid);
    catch
        fprintf('❌ CRITICAL: config.json could not be restored to its original state.\n');
    end
end
