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
    y = [1; 0; 0; 0; 1; 1];
    rng(0);
    fold_id = make_grouped_folds(ids, y, 3);

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
    y = zeros(length(ids), 1);
    rng(1);
    fold_id = make_grouped_folds(ids, y, 5);

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
%  3. k is clamped when n_folds > min_class_count or n_unique
% -----------------------------------------------------------------------
try
    % 2 unique patients, request 5 folds -> should give 2 folds at most
    ids = {'X'; 'X'; 'Y'; 'Y'};
    y = [1; 0; 0; 0];
    rng(2);
    fold_id = make_grouped_folds(ids, y, 5);

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
    y = randi([0 1], length(ids), 1);

    fold_id = make_grouped_folds(ids, y, 5);
    assert(all(fold_id ~= 0), 'Every row must receive a valid fold ID');
    assert(numel(fold_id) == numel(ids), 'Output length must match input length');
    fprintf('[PASS] test_no_unassigned_rows\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] test_no_unassigned_rows: %s\n', e.message);
    n_fail = n_fail + 1;
end

% -----------------------------------------------------------------------
%  5. Stratification test (minority class spread)
% -----------------------------------------------------------------------
try
    rng(4);
    n = 10;
    ids = cell((n*2), 1);
    y = zeros((n*2), 1);
    
    for i=1:n
        ids{2*i - 1} = sprintf('S%02d', i);
        ids{2*i} = sprintf('S%02d', i);
    end
    
    % Only patient S01 and S02 have an event
    y(1:2) = 1; % S01
    y(3:4) = 1; % S02
    
    % We have 10 patients, 2 events. Limit to 2 folds so stratification places 1 event per fold.
    fold_id = make_grouped_folds(ids, y, 2);
    
    s01_fold = fold_id(1);
    s02_fold = fold_id(3);
    
    assert(s01_fold ~= s02_fold, 'Events were clustered in the same fold rather than stratified');
    assert(max(fold_id) == 2, 'K should have clamped/succeeded down to 2');
    
    fprintf('[PASS] test_stratified_spread\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] test_stratified_spread: %s\n', e.message);
    n_fail = n_fail + 1;
end

% -----------------------------------------------------------------------
%  Summary
% -----------------------------------------------------------------------
fprintf('\n--- test_grouped_folds results: %d passed, %d failed ---\n', n_pass, n_fail);
if n_fail > 0
    error('test_grouped_folds: %d test(s) FAILED', n_fail);
end
