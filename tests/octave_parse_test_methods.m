function methods_list = octave_parse_test_methods(txt, attr)
    % OCTAVE_PARSE_TEST_METHODS  Extract function names from a methods(attr) block.
    %   Parses the raw text of a classdef test file and returns a cell array
    %   of function names found inside methods blocks that match the given
    %   attribute (e.g., 'Test', 'TestMethodSetup', 'TestMethodTeardown').
    %
    %   This is needed because GNU Octave does not fully support MATLAB's
    %   unittest metaclass introspection, so we parse the source text directly.
    %
    %   Inputs:
    %     txt  - Full text content of a .m classdef file (char vector).
    %     attr - The methods block attribute to search for (e.g., 'Test').
    %
    %   Outputs:
    %     methods_list - Cell array of function name strings found in matching
    %                    methods blocks.

    methods_list = {};

    % Build a regex to find 'methods(attr...)' block headers
    pat = ['methods\s*\(\s*' attr '[\s,\)]+'];
    idx = regexp(txt, pat);
    if isempty(idx); return; end

    % Process each matching methods block
    for s = 1:numel(idx)
        depth = 0;
        % Split the text from the block header onward into lines
        lines = strsplit(txt(idx(s):end), '\n');
        for i = 1:numel(lines)
            ln = strtrim(lines{i});
            % Track nesting depth via 'methods' and 'end' keywords
            if ~isempty(regexp(ln, '^\s*methods', 'once'))
                depth = depth + 1;
            end
            if strcmp(ln, 'end')
                depth = depth - 1;
                if depth <= 0; break; end  % Exited the target methods block
            end
            % Match function definitions, capturing just the function name.
            % Handles both 'function name(...)' and 'function ret = name(...)'.
            tokens = regexp(ln, '^\s*function\s+(?:\[?\w+[,\s\w]*\]?\s*=\s*)?(\w+)\s*\(', 'tokens');
            if ~isempty(tokens) && depth >= 1
                methods_list{end+1} = tokens{1}{1};
            end
        end
    end
end
