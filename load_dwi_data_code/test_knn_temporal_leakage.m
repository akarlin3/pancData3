%% test_knn_temporal_leakage.m
% Verification script for patient-aware KNN imputation

% 1. Setup Mock Data
% 2 patients, 3 timepoints each.
% Baseline features are identical for the same patient (temporal leakage trap).
pat_ids = [1; 1; 1; 2; 2; 2];
X = [
    1, 10; % Pat 1, T1
    1, NaN; % Pat 1, T2 (Missing feat 2)
    1, 12; % Pat 1, T3
    5, 50; % Pat 2, T1
    5, 52; % Pat 2, T2
    5, 54; % Pat 2, T3
];

fprintf('Running Patient-Aware KNN Imputation Verification...\n');

% 2. Run Imputation
k = 2;
[X_imp, ~] = knn_impute_train_test(X, [], k, pat_ids, []);

% 3. Check Result for Row 2 (Pat 1, T2)
% Euclidean distance to Row 1 and Row 3 (same patient) is 0 + (feat 2 diff).
% Since feat 2 is missing in Row 2, the current logic uses mutual features (feat 1).
% Dist to Row 1: sqrt((1-1)^2 / 1) = 0
% Dist to Row 3: sqrt((1-1)^2 / 1) = 0
% Dist to Row 4,5,6: sqrt((1-5)^2 / 1) = 4
%
% If leakage is NOT FIXED, Row 2 should be imputed using Row 1 and 3 (mean = 11).
% If leakage IS FIXED, Row 2 should be imputed using neighbors from Pat 2 (Row 4, 5, 6).

val_row2 = X_imp(2, 2);
fprintf('Imputed value for Patient 1, T2: %.2f\n', val_row2);

% Neighbors from Pat 2 are 50, 52, 54. Mean of 2 closest (or any 2) would be 51, 53, or 52.
% Definitely NOT ~11.
if val_row2 > 40
    fprintf('SUCCESS: Temporal leakage prevented. Neighbors selected from independent patients.\n');
else
    fprintf('FAILURE: Temporal leakage detected. Neighbors likely selected from same patient.\n');
    error('KNN Imputation failed to exclude same-patient rows.');
end

% 4. Verify Test Imputation
X_te = [1, NaN]; % New row for Pat 1
pat_id_te = 1;
[~, X_te_imp] = knn_impute_train_test(X, X_te, k, pat_ids, pat_id_te);
val_te = X_te_imp(1, 2);
fprintf('Imputed value for Test Patient 1: %.2f\n', val_te);

if val_te > 40
    fprintf('SUCCESS: Test leakage prevented.\n');
else
    fprintf('FAILURE: Test leakage detected.\n');
    error('KNN Imputation failed to exclude same-patient rows for test set.');
end
