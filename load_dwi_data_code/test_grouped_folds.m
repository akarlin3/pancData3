%% test_grouped_folds.m
% Validates make_grouped_folds: grouped k-fold ID assignment for patient-aware
% cross-validation to prevent intra-patient data leakage in Elastic Net CV.
%
% Usage:
%   test_grouped_folds     % from MATLAB command window (run as script)

n_pass = 0;
n_fail = 0;

fprintf('Running test_grouped_folds...\n\n');

% -----------------------------------------------------------------------
%  1. All rows for the same patient land in the same fold
% -----------------------------------------------------------------------
try
    % 3 patients, each with 2 rows (simulating start-stop format)
    ids = {'P1'; 'P1'; 'P2'; 'P2'; 'P3'; 'P3'};
    rng(0);
    fold_id = make_grouped_folds(ids, 3);

    % Each patient contributes 2 rows; both rows must share the same fold.
    assert(fold_id(1) == fold_id(2), 'P1 rows must share a fold');
    assert(fold_id(3) == fold_id(4), 'P2 rows must share a fold');
    assert(fold_id(5) == fold_id(6), 'P3 rows must share a fold');
    fprintf('[PASS] test_same_patient_rows_same_fold\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] test_same_patient_rows_same_fold: %s\n', e.message);
    n_fail = n_fail + 1;
end

% -----------------------------------------------------------------------
%  2. Fold IDs cover exactly 1..k (no empty or out-of-range folds)
% -----------------------------------------------------------------------
try
    ids = {'A'; 'A'; 'B'; 'C'; 'C'; 'D'; 'D'; 'D'; 'E'};
    rng(1);
    fold_id = make_grouped_folds(ids, 5);

    assert(min(fold_id) == 1,           'Minimum fold ID must be 1');
    assert(max(fold_id) == 5,           'Maximum fold ID must equal k');
    assert(all(fold_id >= 1),           'All fold IDs must be >= 1');
    assert(all(floor(fold_id) == fold_id), 'Fold IDs must be integers');
    fprintf('[PASS] test_fold_range\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] test_fold_range: %s\n', e.message);
    n_fail = n_fail + 1;
end

% -----------------------------------------------------------------------
%  3. k is clamped when n_folds > n_unique_patients
% -----------------------------------------------------------------------
try
    % 2 unique patients, request 5 folds -> should give 2 folds at most
    ids = {'X'; 'X'; 'Y'; 'Y'};
    rng(2);
    fold_id = make_grouped_folds(ids, 5);

    actual_k = max(fold_id);
    assert(actual_k <= 2, sprintf('Actual folds (%d) must be <= n_unique_patients (2)', actual_k));
    fprintf('[PASS] test_fold_clamping\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] test_fold_clamping: %s\n', e.message);
    n_fail = n_fail + 1;
end

% -----------------------------------------------------------------------
%  4. No row is left unassigned (fold_id == 0)
% -----------------------------------------------------------------------
try
    rng(3);
    n_patients = 10;
    % Each patient contributes between 1 and 4 rows
    ids = {};
    for p = 1:n_patients
        n_rows = randi(4);
        for r = 1:n_rows
            ids{end+1} = sprintf('PAT%02d', p); %#ok<SAGROW>
        end
    end
    ids = ids(:);  % column cell array

    fold_id = make_grouped_folds(ids, 5);
    assert(all(fold_id ~= 0), 'Every row must receive a valid fold ID');
    assert(numel(fold_id) == numel(ids), 'Output length must match input length');
    fprintf('[PASS] test_no_unassigned_rows\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] test_no_unassigned_rows: %s\n', e.message);
    n_fail = n_fail + 1;
end

% -----------------------------------------------------------------------
%  5. Cross-sectional data (one row per patient) — verify grouping matches
%     a standard k-fold partition exactly in terms of patient isolation
% -----------------------------------------------------------------------
try
    rng(4);
    n = 20;
    ids = arrayfun(@(i) sprintf('S%02d', i), 1:n, 'UniformOutput', false)';
    fold_id = make_grouped_folds(ids, 5);

    % With one row per patient, each fold should have ~n/5 = 4 patients
    for f = 1:5
        n_in_fold = sum(fold_id == f);
        assert(n_in_fold >= 1, sprintf('Fold %d must be non-empty', f));
    end
    fprintf('[PASS] test_cross_sectional_one_row_per_patient\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] test_cross_sectional_one_row_per_patient: %s\n', e.message);
    n_fail = n_fail + 1;
end

% -----------------------------------------------------------------------
%  Summary
% -----------------------------------------------------------------------
fprintf('\n--- test_grouped_folds results: %d passed, %d failed ---\n', n_pass, n_fail);
if n_fail > 0
    error('test_grouped_folds: %d test(s) FAILED', n_fail);
end
