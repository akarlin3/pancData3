classdef test_perf < matlab.unittest.TestCase
% TEST_PERF  Benchmark for GTV file search pattern matching used in
%   discover_gtv_file and find_gtv_files.
%
%   GTV mask files are located by matching their names against wildcard
%   patterns (e.g., '*GTV_MR', '*GTVp'). This benchmark compares two
%   approaches over 10,000 filenames:
%     1. Original: calls regexptranslate('wildcard', ...) inside the inner
%        loop, recompiling the regex on every file.
%     2. Optimized: pre-compiles regex patterns once before the loop.
%
%   Both must produce identical match results.

    methods (Test)
        function testOptimizedMatchesOriginal(testCase)
        %TESTOPTIMIZEDMATCHESORIGINAL Verify pre-compiled regex matches
        %   produce the same results as in-loop compilation, and measure
        %   the speedup.

            % GTV wildcard search patterns used by the pipeline
            patterns = {'*GTV_MR', '*GTVp', '*GTV_panc*'};
            index = 1;  % Patient index appended to pattern

            % Generate 10,000 synthetic GTV-like filenames
            gtv_names = cell(1, 10000);
            for i = 1:10000
                gtv_names{i} = sprintf('prefix_GTVp1_token3_%d.mat', i);
            end

            % Approach 1 (Original): compile regex inside the inner loop.
            % regexptranslate is called patterns x files times.
            tic;
            gtv_search_result_old = zeros(1, 10000);
            for gi = 1:length(gtv_names)
                % Extract the second underscore-delimited token for matching
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

            % Approach 2 (Optimized): pre-compile regex patterns once,
            % then reuse them across all files.
            tic;
            precomputed_regex_patterns = cell(1, length(patterns));
            for p = 1:length(patterns)
                precomputed_regex_patterns{p} = regexptranslate('wildcard', [patterns{p} int2str(index)]);
            end
            gtv_search_result_new = zeros(1, 10000);
            for gi = 1:length(gtv_names)
                gtmp = strsplit(gtv_names{gi}, '_');
                if length(gtmp) >= 2
                    gtmp_tok = gtmp{2};
                    for p = 1:length(patterns)
                        isfound = regexp(gtmp_tok, precomputed_regex_patterns{p});
                        if ~isempty(isfound) && isfound(1) == 1
                            gtv_search_result_new(gi) = 1;
                            break;
                        end
                    end
                end
            end
            time_new = toc;

            fprintf('Original: %.4f s, Optimized: %.4f s, Speedup: %.2fx\n', time_old, time_new, time_old / time_new);
            % Both approaches must produce identical match vectors
            testCase.verifyEqual(gtv_search_result_old, gtv_search_result_new);
        end
    end
end
