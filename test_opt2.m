n_rows = 100000;
id_list_cell = cell(n_rows, 1);
k = 5;
for i=1:n_rows
    id_list_cell{i} = sprintf('PT%04d', randi(1000));
end
[unique_ids, ~, ic] = unique(id_list_cell);
n_unique = numel(unique_ids);

% simulate cvpartition 'test' function
% randomly assign each unique id to a fold
pt_fold_assignments = randi(k, n_unique, 1);
test_mock = @(f) pt_fold_assignments == f;

tic;
fold_id1 = zeros(n_rows, 1);
for f = 1:k
    pt_idx = find(test_mock(f));
    fold_id1(ismember(id_list_cell, unique_ids(pt_idx))) = f;
end
t1 = toc;

tic;
pt_fold = zeros(n_unique, 1);
for f = 1:k
    pt_idx = find(test_mock(f));
    pt_fold(pt_idx) = f;
end
fold_id2 = pt_fold(ic);
t2 = toc;

fprintf('Baseline part 2: %f s\n', t1);
fprintf('Optimized part 2: %f s\n', t2);
fprintf('Same output? %d\n', isequal(fold_id1, fold_id2));
