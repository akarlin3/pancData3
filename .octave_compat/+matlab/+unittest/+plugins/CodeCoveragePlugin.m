% CODECOVERAGEPLUGIN  Octave-compatible no-op shim for MATLAB's CodeCoveragePlugin.
%
%   MATLAB's CodeCoveragePlugin (R2014b+) instruments source files during
%   test execution and produces a code coverage report. Octave has no
%   equivalent profiling/coverage mechanism that integrates with the unittest
%   framework. This shim accepts the same API calls (forFolder) but performs
%   no coverage instrumentation or reporting -- it exists solely to prevent
%   errors when test runner code tries to add the plugin.
classdef CodeCoveragePlugin
    % CodeCoveragePlugin  No-op shim for Octave compatibility.
    methods (Static)
        function obj = forFolder(folders)
            % FORFOLDER  Create a CodeCoveragePlugin for the given source folders.
            %   Returns an empty plugin object. No coverage data is collected.
            obj = matlab.unittest.plugins.CodeCoveragePlugin();
        end
    end
end
