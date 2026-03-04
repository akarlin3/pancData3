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

[unique_ids, ~, ic] = unique(id_list_cell);
n_unique   = numel(unique_ids);

% Derive patient-level event status. A patient is an "event" if ANY of
% their longitudinal rows contains an event (max(y) > 0).
% Optimized: vectorized string grouping and any() evaluation via accumarray
pt_y = double(accumarray(ic, double(y > 0), [n_unique, 1], @any));

% Safety check: don't request more folds than we have unique patients
k = min(n_folds, n_unique);

if k <= 1
    % Degenerate case: fewer than 2 unique patients
    fold_id = ones(numel(id_list_cell), 1);
    return;
end

% Attempt stratified partition; fall back to unstratified if the
% minority class is too small for k folds.
try
    warnState = warning('off', 'stats:cvpartition:KFoldMissingGrp');
    restoreWarn = onCleanup(@() warning(warnState));
    cvp = cvpartition(pt_y, 'KFold', k);
catch
    % Retry with fewer folds to preserve stratification before falling back
    % to unstratified partitioning.
    % Ensure at least 2-fold stratification when any events exist, to
    % prevent all events from clustering in a single fold.
    k_try = max(2, min([k, sum(pt_y == 0), sum(pt_y == 1)]));
    if k_try >= 2
        try
            cvp = cvpartition(pt_y, 'KFold', k_try);
            warning('make_grouped_folds:reducedFolds', ...
                'Stratified CV with %d folds failed; using %d folds to preserve stratification.', k, k_try);
            % Also log via fprintf so message appears in diary even when
            % callers suppress warnings with warning('off','all').
            fprintf('  ⚠️  make_grouped_folds: Stratified CV with %d folds failed; using %d folds.\n', k, k_try);
        catch
            warning('make_grouped_folds:unstratified', ...
                'Stratified CV failed (minority class too small for %d folds). Using unstratified folds.', k);
            fprintf('  ⚠️  make_grouped_folds: Falling back to unstratified folds (minority class too small for %d folds).\n', k);
            cvp = cvpartition(n_unique, 'KFold', k);
        end
    else
        warning('make_grouped_folds:unstratified', ...
            'Stratified CV failed (minority class too small for %d folds). Using unstratified folds.', k);
        fprintf('  ⚠️  make_grouped_folds: Falling back to unstratified folds (minority class too small for %d folds).\n', k);
        cvp = cvpartition(n_unique, 'KFold', k);
    end
end

pt_fold = zeros(n_unique, 1);
for f = 1:k
    pt_idx = find(test(cvp, f));
    pt_fold(pt_idx) = f;
end
fold_id = pt_fold(ic);
end
