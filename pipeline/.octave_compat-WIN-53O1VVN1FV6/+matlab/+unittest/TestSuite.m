% TESTSUITE  Octave-compatible shim for MATLAB's matlab.unittest.TestSuite.
%
%   MATLAB's TestSuite class (part of the unittest framework, R2013a+)
%   provides test discovery and organization. Octave does not ship with
%   this framework, so this shim reimplements the fromFolder() discovery
%   method by parsing classdef files on disk.
%
%   Behavioral differences from MATLAB's TestSuite:
%   - Only fromFolder() is implemented (no fromClass, fromMethod, etc.).
%   - Test discovery is file-based: finds test_*.m and benchmark_*.m files,
%     verifies they are classdef files (not scripts), and parses method
%     names from methods(Test) blocks via regex.
%   - Returns a struct array instead of a TestSuite object array; the
%     TestRunner shim consumes this struct format.
%   - Does not support test parameterization or shared fixtures.
classdef TestSuite
    % TestSuite  Minimal shim for matlab.unittest.TestSuite under Octave.
    properties
        Name = '';
        TestClass = '';
        TestMethod = '';
    end

    methods (Static)
        function suite = fromFolder(folder, varargin)
            % FROMFOLDER  Discover test classes in a folder.
            %   suite = TestSuite.fromFolder(folder) scans for test_*.m files.
            %   suite = TestSuite.fromFolder(folder, 'IncludingSubfolders', true)
            %   also recurses into subdirectories and picks up benchmark_*.m.
            %
            %   Returns a struct array with fields: Name, TestClass,
            %   TestMethod, FilePath, Folder -- one entry per test method.
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

            % Build the suite by iterating over discovered files.
            suite = [];
            for i = 1:numel(files)
                [~, className, ~] = fileparts(files(i).name);
                % Skip non-classdef files (script-based tests) by peeking at first line.
                % MATLAB's fromFolder also skips scripts, so this matches that behavior.
                fpath = fullfile(files(i).folder, files(i).name);
                fid = fopen(fpath, 'r');
                if fid == -1; continue; end
                firstLine = fgetl(fid);
                fclose(fid);
                if isempty(firstLine) || ~contains(firstLine, 'classdef')
                    continue;
                end

                % Find all Test methods by parsing the file's source text.
                % MATLAB introspects class metadata; we parse regex instead.
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
            % PARSETESTMETHODS  Extract test method names from a classdef file.
            %   Reads the file as text and uses a simple state machine to find
            %   methods(Test) blocks, then extracts function names within them.
            %   Tracks nesting depth to correctly handle the closing 'end' of
            %   the methods block vs. individual function 'end' keywords.
            methods_list = {};
            fid = fopen(filepath, 'r');
            if fid == -1; return; end
            txt = fread(fid, '*char')';
            fclose(fid);

            % State machine: in_test_block tracks whether we are inside a
            % methods(Test) block; brace_depth counts nested function/end pairs
            % so we know when the methods block itself closes.
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
