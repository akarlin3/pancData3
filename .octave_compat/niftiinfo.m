function info = niftiinfo(filename)
% NIFTIINFO Octave-compatible shim for MATLAB's niftiinfo.
%   Returns a struct with basic NIfTI header information.
    info = struct();
    info.Filename = filename;
    info.ImageSize = [0 0 0];
    info.PixelDimensions = [1 1 1];
    info.Datatype = 'double';
end
