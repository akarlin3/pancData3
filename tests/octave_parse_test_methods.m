function methods_list = octave_parse_test_methods(txt, attr)
    % OCTAVE_PARSE_TEST_METHODS  Extract function names from a methods(attr) block.
    methods_list = {};
    pat = ['methods\s*\(\s*' attr '[\s,\)]+'];
    idx = regexp(txt, pat);
    if isempty(idx); return; end
    for s = 1:numel(idx)
        depth = 0;
        lines = strsplit(txt(idx(s):end), '\n');
        for i = 1:numel(lines)
            ln = strtrim(lines{i});
            if ~isempty(regexp(ln, '^\s*methods', 'once'))
                depth = depth + 1;
            end
            if strcmp(ln, 'end')
                depth = depth - 1;
                if depth <= 0; break; end
            end
            tokens = regexp(ln, '^\s*function\s+(?:\[?\w+[,\s\w]*\]?\s*=\s*)?(\w+)\s*\(', 'tokens');
            if ~isempty(tokens) && depth >= 1
                methods_list{end+1} = tokens{1}{1};
            end
        end
    end
end
