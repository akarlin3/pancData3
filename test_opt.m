n_rows = 100000;
id_list_cell = cell(n_rows, 1);
y = rand(n_rows, 1) > 0.9;
for i=1:n_rows
    id_list_cell{i} = sprintf('PT%04d', randi(1000));
end

tic;
unique_ids = unique(id_list_cell);
n_unique   = numel(unique_ids);
pt_y = zeros(n_unique, 1);
for i = 1:n_unique
    pt_y(i) = double(any(y(strcmp(id_list_cell, unique_ids{i})) > 0));
end
t1 = toc;

tic;
[unique_ids2, ~, ic] = unique(id_list_cell);
n_unique2 = numel(unique_ids2);
pt_y2 = zeros(n_unique2, 1);
for i = 1:n_unique2
    pt_y2(i) = double(any(y(ic == i) > 0));
end
t2 = toc;

tic;
[unique_ids3, ~, ic3] = unique(id_list_cell);
n_unique3 = numel(unique_ids3);
pt_y3 = accumarray(ic3, y, [n_unique3, 1], @(x) double(any(x > 0)));
t3 = toc;

fprintf('Baseline: %f s\n', t1);
fprintf('ic loop: %f s\n', t2);
fprintf('accumarray: %f s\n', t3);
