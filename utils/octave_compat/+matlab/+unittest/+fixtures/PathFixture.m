classdef PathFixture
    % PathFixture  Shim that adds a folder to the path (Octave compat).
    properties
        Folder = '';
    end
    methods
        function obj = PathFixture(folder)
            obj.Folder = folder;
            addpath(folder);
        end
    end
end
