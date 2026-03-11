% PATHFIXTURE  Octave-compatible shim for matlab.unittest.fixtures.PathFixture.
%
%   MATLAB's PathFixture temporarily adds a folder to the MATLAB path for
%   the duration of a test, then removes it on teardown. This shim only
%   performs the addpath step; it does not remove the path on teardown,
%   because Octave's old-style classdef does not support destructor-based
%   cleanup. In practice this is acceptable for test runs where the path
%   additions are harmless after the test completes.
classdef PathFixture
    % PathFixture  Shim that adds a folder to the path (Octave compat).
    properties
        Folder = '';  % The folder path to add to the MATLAB/Octave search path.
    end
    methods
        function obj = PathFixture(folder)
            % PATHFIXTURE  Constructor: stores the folder and adds it to the path.
            %   Unlike MATLAB's version, this does not track the path for later
            %   removal -- the folder remains on the path after the fixture is
            %   disposed.
            obj.Folder = folder;
            addpath(folder);
        end
    end
end
