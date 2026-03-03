classdef TestSuite
    % TestSuite  Minimal shim for matlab.unittest.TestSuite under Octave.
    properties
        Name = '';
        TestClass = '';
        TestMethod = '';
    end

    methods (Static)
        function suite = fromFolder(folder, varargin)
            % Discover test classes in folder (and optionally subfolders).
            includeSubfolders = false;
            for k = 1:2:numel(varargin)
                if strcmpi(varargin{k}, 'IncludingSubfolders')
                    includeSubfolders = varargin{k+1};
                end
            end

            if includeSubfolders
                files = dir(fullfile(folder, '**', 'test_*.m'));
                % Also pick up benchmark_*.m from subfolders
                files2 = dir(fullfile(folder, '**', 'benchmark_*.m'));
                files = [files; files2];
            else
                files = dir(fullfile(folder, 'test_*.m'));
            end

            suite = [];
            for i = 1:numel(files)
                [~, className, ~] = fileparts(files(i).name);
                % Skip non-classdef files (script-based tests) by peeking at first line
                fpath = fullfile(files(i).folder, files(i).name);
                fid = fopen(fpath, 'r');
                if fid == -1; continue; end
                firstLine = fgetl(fid);
                fclose(fid);
                if isempty(firstLine) || ~contains(firstLine, 'classdef')
                    continue;
                end

                % Find all Test methods by parsing the file
                methods_list = TestSuite.parseTestMethods(fpath);
                for j = 1:numel(methods_list)
                    entry = struct();
                    entry.Name = [className '/' methods_list{j}];
                    entry.TestClass = className;
                    entry.TestMethod = methods_list{j};
                    entry.FilePath = fpath;
                    entry.Folder = files(i).folder;
                    if isempty(suite)
                        suite = entry;
                    else
                        suite(end+1) = entry;
                    end
                end
            end
        end
    end

    methods (Static, Access = private)
        function methods_list = parseTestMethods(filepath)
            % Parse a classdef file and extract method names from methods(Test) blocks
            methods_list = {};
            fid = fopen(filepath, 'r');
            if fid == -1; return; end
            txt = fread(fid, '*char')';
            fclose(fid);

            in_test_block = false;
            brace_depth = 0;
            lines = strsplit(txt, '\n');
            for i = 1:numel(lines)
                ln = strtrim(lines{i});
                % Detect start of methods (Test) block
                if ~isempty(regexp(ln, '^\s*methods\s*\(\s*Test\s*\)', 'once'))
                    in_test_block = true;
                    brace_depth = 0;
                    continue;
                end
                if in_test_block
                    % Track nesting by counting function/end or looking for next methods/end
                    if ~isempty(regexp(ln, '^\s*end\s*$', 'once'))
                        if brace_depth <= 0
                            in_test_block = false;
                            continue;
                        else
                            brace_depth = brace_depth - 1;
                            continue;
                        end
                    end
                    % Match function definitions
                    tokens = regexp(ln, '^\s*function\s+(?:\w+\s*=\s*)?(\w+)\s*\(', 'tokens');
                    if ~isempty(tokens)
                        methods_list{end+1} = tokens{1}{1};
                        brace_depth = brace_depth + 1;
                    end
                end
            end
        end
    end
end
