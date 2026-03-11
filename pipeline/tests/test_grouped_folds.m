%% test_grouped_folds.m
% Validates make_grouped_folds: grouped k-fold ID assignment for patient-aware
% cross-validation to prevent intra-patient data leakage in Elastic Net CV.
%
% In survival analysis with start-stop (interval) data, each patient may
% contribute multiple rows. If rows from the same patient end up in different
% CV folds, the model can memorise within-patient patterns (temporal leakage).
% make_grouped_folds prevents this by assigning all rows for a given patient
% to the same fold, while attempting to stratify events across folds.
%
% Tests cover:
%   1. Same-patient rows land in the same fold
%   2. Fold IDs span exactly 1..k with no gaps
%   3. k is clamped when fewer unique patients than requested folds
%   3b. Class imbalance does NOT reduce k below n_unique_patients
%   4. Every row receives a valid (non-zero) fold assignment
%   5. Minority-class events are stratified across folds
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
    % 3 patients, each with 2 rows (simulating start-stop survival format).
    % y encodes event status: P1 has an event, P2 does not, P3 does.
    ids = {'P1'; 'P1'; 'P2'; 'P2'; 'P3'; 'P3'};
    y = [1; 0; 0; 0; 1; 1];
    rng(0); % Fix seed so fold assignment is deterministic
    fold_id = make_grouped_folds(ids, y, 3);

    % Each patient contributes 2 rows; both rows must share the same fold.
    % This is the fundamental patient-grouping guarantee.
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
    % 5 unique patients with varying row counts (A=2, B=1, C=2, D=3, E=1).
    % All events set to 0 (no stratification pressure).
    ids = {'A'; 'A'; 'B'; 'C'; 'C'; 'D'; 'D'; 'D'; 'E'};
    y = zeros(length(ids), 1);
    rng(1);
    fold_id = make_grouped_folds(ids, y, 5);

    % With exactly 5 unique patients and k=5, every fold ID from 1..5 must appear.
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
    % Only 2 unique patients but 5 folds requested. Cannot split 2 patients
    % into more than 2 groups, so k must be clamped to at most 2.
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
%  3b. Class imbalance must NOT reduce k below n_unique_patients
% -----------------------------------------------------------------------
try
    % 20 patients, only 2 with events — severe class imbalance.
    % Request 5 folds. Old code clamped k to min_class_count (2),
    % destroying training set size. k must stay at 5.
    rng(5);
    n_pts = 20;
    ids = {};
    for p = 1:n_pts
        ids{end+1} = sprintf('IMB%02d', p); %#ok<SAGROW>
        ids{end+1} = sprintf('IMB%02d', p); %#ok<SAGROW>
    end
    ids = ids(:);
    y = zeros(numel(ids), 1);
    y(1:2) = 1;  % IMB01 is an event
    y(3:4) = 1;  % IMB02 is an event

    fold_id = make_grouped_folds(ids, y, 5);
    actual_k = max(fold_id);
    assert(actual_k == 5, ...
        sprintf('With 20 patients and 2 events, k should be 5, got %d', actual_k));
    assert(all(fold_id >= 1), 'All fold IDs must be >= 1');
    fprintf('[PASS] test_imbalance_preserves_k\n');
    n_pass = n_pass + 1;
catch e
    fprintf('[FAIL] test_imbalance_preserves_k: %s\n', e.message);
    n_fail = n_fail + 1;
end

% -----------------------------------------------------------------------
%  4. No row is left unassigned (fold_id == 0)
% -----------------------------------------------------------------------
try
    rng(3);
    n_patients = 10;
    % Each patient contributes between 1 and 4 rows (random).
    % This tests the common case of variable-length patient records.
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
    % A fold_id of 0 would indicate an unassigned row (bug).
    assert(all(fold_id ~= 0), 'Every row must receive a valid fold ID');
    % Output vector must be the same length as the input row count.
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

    % Create 10 patients, each with 2 rows (start-stop format)
    for i=1:n
        ids{2*i - 1} = sprintf('S%02d', i);
        ids{2*i} = sprintf('S%02d', i);
    end

    % Only patient S01 and S02 have an event (severe class imbalance: 2/10)
    y(1:2) = 1; % S01
    y(3:4) = 1; % S02

    % With 2 folds and 2 event-patients, stratification should place
    % exactly 1 event-patient in each fold so that every fold has at
    % least one positive example for model training.
    fold_id = make_grouped_folds(ids, y, 2);

    s01_fold = fold_id(1);
    s02_fold = fold_id(3);

    % If both event-patients land in the same fold, one training fold has
    % zero events, which breaks survival model fitting.
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
