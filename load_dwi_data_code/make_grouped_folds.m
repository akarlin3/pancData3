function fold_id = make_grouped_folds(id_list_cell, n_folds)
% make_grouped_folds  Build a k-fold ID array grouped by patient ID.
%
%   fold_id = make_grouped_folds(id_list_cell, n_folds) assigns an integer
%   fold label (1..k) to every row so that ALL rows sharing the same
%   patient ID receive the same label.  This prevents intra-patient data
%   leakage when a feature matrix contains multiple rows per patient (e.g.,
%   counting-process / start-stop survival panels where one patient
%   contributes a Week-1 row and a Week-3 row).
%
%   Inputs
%   ------
%   id_list_cell - (n_rows x 1) cell array of char vectors (patient IDs).
%   n_folds      - requested number of folds k; clamped to n_unique_patients
%                  so the function is safe even for very small datasets.
%
%   Output
%   ------
%   fold_id - (n_rows x 1) integer column vector, values in 1..k.
%
%   Algorithm
%   ---------
%   1. Collect the sorted list of unique patient IDs.
%   2. Partition the unique IDs into k folds using cvpartition (random,
%      without stratification — the caller controls the rng seed).
%   3. Propagate each patient's fold assignment to every row that belongs
%      to that patient via strcmp, so all rows for a patient land in the
%      same fold regardless of how many intervals that patient contributes.

unique_ids = unique(id_list_cell);
n_unique   = numel(unique_ids);
k          = min(n_folds, n_unique);

cvp = cvpartition(n_unique, 'KFold', k);

fold_id = zeros(numel(id_list_cell), 1);
for f = 1:k
    pt_idx = find(test(cvp, f));
    fold_id(ismember(id_list_cell, unique_ids(pt_idx))) = f;
end
end
