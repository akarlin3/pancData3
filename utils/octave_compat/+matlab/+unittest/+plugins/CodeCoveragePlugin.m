classdef CodeCoveragePlugin
    % CodeCoveragePlugin  No-op shim for Octave compatibility.
    methods (Static)
        function obj = forFolder(folders)
            obj = matlab.unittest.plugins.CodeCoveragePlugin();
        end
    end
end
