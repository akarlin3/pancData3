function fold_id = make_grouped_folds(id_list_cell, y, n_folds)
% make_grouped_folds  Build a stratified k-fold ID array grouped by patient ID.
%
%   fold_id = make_grouped_folds(id_list_cell, y, n_folds) assigns an integer
%   fold label (1..k) to every row so that ALL rows sharing the same
%   patient ID receive the same label. It strictly supports Stratified
%   Grouped K-Fold Cross-Validation, ensuring the ratio of events to
%   non-events is balanced across folds based on the patient-level max(y).
%
%   Inputs
%   ------
%   id_list_cell - (n_rows x 1) cell array of char vectors (patient IDs).
%   y            - (n_rows x 1) numeric array of event indicators (1=event, 0=censor).
%   n_folds      - requested number of folds k; clamped only to
%                  n_unique_patients if k exceeds the patient count.
%
%   Output
%   ------
%   fold_id - (n_rows x 1) integer column vector, values in 1..k.
%
%   Analytical Rationale — Why Patient-Grouped Folds are Essential:
%   ---------------------------------------------------------------
%   In longitudinal DWI studies, each patient contributes multiple rows
%   to the time-dependent panel (one per scan interval from build_td_panel).
%   If rows from the same patient appear in both train and test folds, the
%   model can memorize patient-specific patterns (e.g., baseline ADC) and
%   appear to generalize when it is actually recognizing the same patient.
%   This is intra-patient data leakage — the most common and most
%   damaging form of leakage in longitudinal oncology studies.
%
%   Patient-grouped folds guarantee that ALL rows from a given patient are
%   in the same fold, so the model is always evaluated on truly unseen
%   patients.  This yields honest estimates of how the model would perform
%   on a new patient walking into the clinic.
%
%   Stratification preserves the event rate across folds.  In pancreatic
%   cancer with ~30-40% local failure rate, unstratified random splits
%   could produce folds with 0 events, making it impossible to estimate
%   hazard ratios or concordance indices for that fold.

% Map each row to its unique patient index.  ic(row) gives the index
% into unique_ids, allowing us to propagate the patient-level fold
% assignment back to every row belonging to that patient.
[unique_ids, ~, ic] = unique(id_list_cell);
n_unique   = numel(unique_ids);

% Derive patient-level event status for stratification.  A patient is
% classified as an "event" if ANY of their longitudinal rows contains a
% local failure (y == 1).  This reflects the clinical reality that a
% patient either experienced local failure or did not — the event is a
% patient-level outcome, not a row-level outcome.
%
% Competing risk patients (y == 2, e.g., death from non-cancer causes)
% are NOT counted as events for stratification purposes.  In the
% Cause-Specific Hazard framework used downstream, competing events are
% treated as censored for the cause of interest (local failure).
% Grouping them with truly censored patients (y == 0) for fold balancing
% ensures each fold has a representative mix of local failures vs.
% non-failures, which is what the Cox model needs to estimate stable
% hazard ratios.
%
% accumarray with @any efficiently computes "did this patient have any
% local failure?" across all their rows without an explicit loop.
pt_y = double(accumarray(ic, double(y == 1), [n_unique, 1], @any));

% Safety check: don't request more folds than we have unique patients.
% With k > n_unique, some folds would be empty, and leave-one-patient-out
% CV (k = n_unique) is the maximum meaningful granularity for grouped CV.
k = min(n_folds, n_unique);

if k <= 1
    % Degenerate case: fewer than 2 unique patients
    fold_id = ones(numel(id_list_cell), 1);
    return;
end

% Attempt stratified partition; fall back to unstratified if the
% minority class is too small for k folds.  cvpartition requires at
% least one member of each class per fold; with rare events (e.g., 3
% local failures out of 50 patients), 5-fold stratified CV is impossible
% and we gracefully degrade to fewer folds or unstratified partitioning.
%
% When the minority class count is between 1 and k-1, cvpartition
% succeeds but emits warning stats:cvpartition:KFoldMissingGrp.  We
% suppress the warning AND clear lastwarn so it does not propagate to
% the pipeline error log via run_dwi_pipeline's lastwarn check.
% Reducing k to match the minority class size would eliminate the
% warning but would also halve the training set size with small event
% counts (e.g., 2-fold instead of 5-fold), hurting model stability.
try
    warnState = warning('off', 'stats:cvpartition:KFoldMissingGrp');
    restoreWarn = onCleanup(@() warning(warnState));
    [prev_msg, prev_id] = lastwarn;
    cvp = cvpartition(pt_y, 'KFold', k);
    % Clear lastwarn if cvpartition updated it with the suppressed warning.
    % Suppressing display via warning('off',...) does not prevent lastwarn
    % from being updated in all MATLAB versions.
    [~, cur_id] = lastwarn;
    if strcmp(cur_id, 'stats:cvpartition:KFoldMissingGrp')
        lastwarn(prev_msg, prev_id);
    end
catch
    % Retry with fewer folds to preserve stratification before falling back
    % to unstratified partitioning.
    % Ensure at least 2-fold stratification when any events exist, to
    % prevent all events from clustering in a single fold.
    k_try = max(2, min([k, sum(pt_y == 0), sum(pt_y == 1)]));
    if k_try >= 2
        try
            cvp = cvpartition(pt_y, 'KFold', k_try);
            fprintf('  💡 make_grouped_folds: Stratified CV with %d folds failed; using %d folds.\n', k, k_try);
            k = k_try;  % Sync loop bound with actual partition size
        catch
            fprintf('  💡 make_grouped_folds: Falling back to unstratified folds (minority class too small for %d folds).\n', k);
            cvp = cvpartition(n_unique, 'KFold', k);
        end
    else
        fprintf('  💡 make_grouped_folds: Falling back to unstratified folds (minority class too small for %d folds).\n', k);
        cvp = cvpartition(n_unique, 'KFold', k);
    end
end

% Convert cvpartition's test-set indicator back to fold labels.
% cvpartition stores fold membership as binary test indicators; we need
% integer fold labels (1..k) for downstream indexing.
pt_fold = zeros(n_unique, 1);
for f = 1:k
    pt_idx = find(test(cvp, f));
    pt_fold(pt_idx) = f;
end
% Propagate patient-level fold labels to every row.  ic maps each row
% to its patient index, so pt_fold(ic) gives each row the fold label
% of its parent patient — guaranteeing grouped assignment.
fold_id = pt_fold(ic);
end
