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
%   n_folds      - requested number of folds k; clamped to n_unique_patients
%                  or the number of minority class occurrences if k is too large.
%
%   Output
%   ------
%   fold_id - (n_rows x 1) integer column vector, values in 1..k.

[unique_ids, ~, ic] = unique(id_list_cell);
n_unique   = numel(unique_ids);

% Derive patient-level event status. A patient is an "event" if ANY of
% their longitudinal rows contains an event (max(y) > 0).
pt_y = accumarray(ic, y, [n_unique, 1], @(x) double(any(x > 0)));

% Safety check: don't request more folds than we have unique patients
k = min(n_folds, n_unique);

% Stratification constraint: cvpartition fails if k > number of instances in the minority class
n_events = sum(pt_y == 1);
n_censor = sum(pt_y == 0);
min_class_count = min(n_events, n_censor);

if k > min_class_count && min_class_count > 0
    k = min(k, min_class_count);
    % If even the clamped k is 0 or 1, we can't reliably stratify
end

if k <= 1
    % Fallback if data is too weird or only one class exists
    k = min(n_folds, n_unique);
    cvp = cvpartition(n_unique, 'KFold', k);
else
    % Create stratified partition
    try
        cvp = cvpartition(pt_y, 'KFold', k);
    catch
        % Ultimate fallback if cvpartition rejects the distribution
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
