%% test_knn_temporal_leakage.m
% Verification script for patient-aware KNN imputation.
%
% In longitudinal DWI studies, each patient has multiple timepoints. When a
% feature value is missing at one timepoint, KNN imputation finds the k
% nearest rows to fill it in. Without patient-awareness, KNN will select
% other timepoints from the SAME patient as nearest neighbours (since their
% features are nearly identical), effectively leaking the patient's own
% temporal trajectory into the imputed value.
%
% knn_impute_train_test prevents this by setting infinite distances between
% rows belonging to the same patient, forcing the algorithm to borrow
% information only from other patients.
%
% Tests:
%   1. Train-set imputation excludes same-patient rows
%   2. Test-set imputation excludes same-patient rows
%
% Usage:
%   test_knn_temporal_leakage   % from MATLAB command window (run as script)

% -----------------------------------------------------------------------
%  1. Setup Mock Data
% -----------------------------------------------------------------------
% 2 patients, 3 timepoints each.
% Feature 1 is constant within each patient (1 for Pat 1, 5 for Pat 2),
% making intra-patient rows appear extremely close in feature space.
% This is an intentional leakage trap: without the patient-exclusion fix,
% KNN will always pick same-patient rows as nearest neighbours.
pat_ids = [1; 1; 1; 2; 2; 2];
X = [
    1, 10; % Pat 1, T1 — both features present
    1, NaN; % Pat 1, T2 — feature 2 missing (the value to impute)
    1, 12; % Pat 1, T3 — both features present
    5, 50; % Pat 2, T1
    5, 52; % Pat 2, T2
    5, 54; % Pat 2, T3
];

fprintf('Running Patient-Aware KNN Imputation Verification...\n');

% -----------------------------------------------------------------------
%  2. Run Train-Set Imputation
% -----------------------------------------------------------------------
% k=2 nearest neighbours. Empty test set ([]). Patient IDs enforce the
% exclusion constraint: rows from the same patient get infinite distance.
k = 2;
[X_imp, ~] = knn_impute_train_test(X, [], k, pat_ids, []);

% -----------------------------------------------------------------------
%  3. Check Result for Row 2 (Pat 1, T2)
% -----------------------------------------------------------------------
% Distance calculation uses only mutually-present features (feature 1 here,
% since feature 2 is NaN in the target row).
%
% Without leakage fix:
%   Dist to Row 1 (same patient): sqrt((1-1)^2) = 0   --> selected as neighbour
%   Dist to Row 3 (same patient): sqrt((1-1)^2) = 0   --> selected as neighbour
%   Imputed value = mean(10, 12) = 11  (LEAKED from own temporal trajectory)
%
% With leakage fix:
%   Same-patient rows get Inf distance, so KNN is forced to use Pat 2 rows.
%   Dist to Rows 4,5,6: sqrt((1-5)^2) = 4  (all equidistant)
%   Imputed value ~ mean of 2 nearest from {50, 52, 54} ~ 51 (CORRECT)

val_row2 = X_imp(2, 2);
fprintf('Imputed value for Patient 1, T2: %.2f\n', val_row2);

% The threshold of 40 cleanly separates the leakage case (~11) from the
% correct case (~51). Any value above 40 confirms cross-patient imputation.
if val_row2 > 40
    fprintf('SUCCESS: Temporal leakage prevented. Neighbors selected from independent patients.\n');
else
    fprintf('FAILURE: Temporal leakage detected. Neighbors likely selected from same patient.\n');
    error('KNN Imputation failed to exclude same-patient rows.');
end

% -----------------------------------------------------------------------
%  4. Verify Test-Set Imputation
% -----------------------------------------------------------------------
% Simulate a held-out test row for Patient 1 with feature 2 missing.
% The train set still contains Patient 1's rows, so leakage would occur
% if the function does not exclude train rows belonging to the same patient.
X_te = [1, NaN]; % New test row for Pat 1 (feature 1 = 1, feature 2 missing)
pat_id_te = 1;   % Same patient ID as train rows 1-3
[~, X_te_imp] = knn_impute_train_test(X, X_te, k, pat_ids, pat_id_te);
val_te = X_te_imp(1, 2);
fprintf('Imputed value for Test Patient 1: %.2f\n', val_te);

% Same logic: if value > 40, neighbours came from Pat 2 (correct).
% If value ~ 11, neighbours came from Pat 1 train rows (leakage).
if val_te > 40
    fprintf('SUCCESS: Test leakage prevented.\n');
else
    fprintf('FAILURE: Test leakage detected.\n');
    error('KNN Imputation failed to exclude same-patient rows for test set.');
end
