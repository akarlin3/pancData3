patterns = {'*GTV_MR', '*GTVp', '*GTV_panc*'};
index = 1;

% Mock gtv names to run the loop many times
gtv_names = cell(1, 10000);
for i = 1:10000
    gtv_names{i} = sprintf('prefix_GTVp1_token3_%d.mat', i);
end

% Original implementation
tic;
gtv_search_result_old = zeros(1, 10000);
for gi = 1:length(gtv_names)
    gtmp = strsplit(gtv_names{gi}, '_');
    if length(gtmp) >= 2
        gtmp_tok = gtmp{2};
        for p = 1:length(patterns)
            pat = patterns{p};
            regex_pattern = regexptranslate('wildcard', [pat int2str(index)]);
            isfound = regexp(gtmp_tok, regex_pattern);
            if ~isempty(isfound) && isfound(1) == 1
                gtv_search_result_old(gi) = 1;
                break;
            end
        end
    end
end
time_old = toc;

% Optimized implementation
tic;
precomputed_regex_patterns = cell(1, length(patterns));
for p = 1:length(patterns)
    pat = patterns{p};
    precomputed_regex_patterns{p} = regexptranslate('wildcard', [pat int2str(index)]);
end

gtv_search_result_new = zeros(1, 10000);
for gi = 1:length(gtv_names)
    gtmp = strsplit(gtv_names{gi}, '_');
    if length(gtmp) >= 2
        gtmp_tok = gtmp{2};
        for p = 1:length(patterns)
            regex_pattern = precomputed_regex_patterns{p};
            isfound = regexp(gtmp_tok, regex_pattern);
            if ~isempty(isfound) && isfound(1) == 1
                gtv_search_result_new(gi) = 1;
                break;
            end
        end
    end
end
time_new = toc;

fprintf('Original time: %.4f s\n', time_old);
fprintf('Optimized time: %.4f s\n', time_new);
fprintf('Speedup: %.2fx\n', time_old / time_new);

assert(isequal(gtv_search_result_old, gtv_search_result_new), 'Results do not match!');
disp('Results match!');
