% NIFTIINFO  Octave-compatible shim for MATLAB's niftiinfo (R2017b+).
%
%   info = niftiinfo(filename) returns a struct describing a NIfTI file's
%   header metadata (image size, voxel dimensions, data type, etc.).
%
%   MATLAB's niftiinfo reads the actual .nii header. This shim cannot parse
%   real NIfTI files, so it returns a struct with placeholder values. The
%   Filename field is populated so that niftiread(info) works correctly with
%   the companion niftiread shim.
%
%   Behavioral differences from MATLAB's niftiinfo:
%   - ImageSize is [0 0 0] (unknown) rather than the real volume dimensions.
%   - PixelDimensions defaults to [1 1 1] (isotropic 1mm).
%   - Datatype is always 'double' regardless of the actual file.
%   - Fields like TransformMatrix, SpaceUnits, TimeUnits are absent.
function info = niftiinfo(filename)
    % Return a struct with the minimum fields needed by pipeline code.
    info = struct();
    info.Filename = filename;
    info.ImageSize = [0 0 0];
    info.PixelDimensions = [1 1 1];
    info.Datatype = 'double';
end
