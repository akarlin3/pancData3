% bench_metrics_string.m

% Parameters
n_patients = 1000; % Number of patients in id_list
n_clinical = 5000; % Number of rows in clinical table T

% Mock Data Generation
disp('Generating mock data...');
% Create random patient IDs
chars = ['A':'Z', '0':'9', '_'];
id_len = 10;

% Generate id_list (cell array of char vectors)
id_list = cell(n_patients, 1);
for i = 1:n_patients
    id_list{i} = chars(randi(length(chars), 1, id_len));
end

% Generate T.Pat (cell array of char vectors)
% Ensure some overlap
T_Pat = cell(n_clinical, 1);
for i = 1:n_clinical
    if i <= n_patients
        % Make some match (with potential underscores)
        T_Pat{i} = id_list{i};
    else
        T_Pat{i} = chars(randi(length(chars), 1, id_len));
    end
end

% Mock T table (struct for simplicity in benchmark)
T = struct();
T.Pat = T_Pat;
T.LocalOrRegionalFailure = randi([0 1], n_clinical, 1);

% Pre-allocation for results
lf_original = nan(n_patients, 1);
lf_optimized = nan(n_patients, 1);

% --- Benchmark Original Approach ---
disp('Benchmarking Original Approach...');
tic;
for j = 1:length(id_list)
    % Original code logic
    % Note: strrep on a cell array returns a cell array
    i_find = find(contains(strrep(T.Pat,'_','-'), strrep(id_list{j},'_','-')));
    if ~isempty(i_find)
        i_find = i_find(1);
        lf_original(j) = T.LocalOrRegionalFailure(i_find);
    end
end
time_original = toc;
fprintf('Original Approach Time: %.4f seconds\n', time_original);

% --- Benchmark Optimized Approach ---
disp('Benchmarking Optimized Approach...');
tic;
% Optimization: Pre-process strings outside the loop
pat_normalized = strrep(T.Pat, '_', '-');
id_list_normalized = strrep(id_list, '_', '-');

for j = 1:length(id_list)
    % Optimized logic
    i_find = find(contains(pat_normalized, id_list_normalized{j}));
    if ~isempty(i_find)
        i_find = i_find(1);
        lf_optimized(j) = T.LocalOrRegionalFailure(i_find);
    end
end
time_optimized = toc;
fprintf('Optimized Approach Time: %.4f seconds\n', time_optimized);

% Verification
if isequaln(lf_original, lf_optimized)
    disp('Verification Passed: Results match.');
    fprintf('Speedup: %.2fx\n', time_original / time_optimized);
else
    disp('Verification FAILED: Results do not match!');
end
