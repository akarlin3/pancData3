% NIFTIWRITE  Octave-compatible shim for MATLAB's niftiwrite (R2017b+).
%
%   niftiwrite(data, filename) saves image data to disk.
%   niftiwrite(data, filename, info) accepts an optional niftiinfo struct
%   (ignored in this shim).
%
%   MATLAB's niftiwrite produces a proper .nii binary file. This shim
%   instead saves the data as a MATLAB .mat file with the naming convention
%   '<filename>_nifti.mat', which the niftiread shim can load back.
%
%   Behavioral differences from MATLAB's niftiwrite:
%   - Output is a .mat file, not a .nii file; not interoperable with
%     external NIfTI readers (FSL, FreeSurfer, etc.).
%   - Header/spatial metadata (voxel size, orientation) is discarded.
%   - Additional name-value arguments (e.g., 'Compressed') are accepted
%     via varargin but silently ignored.
function niftiwrite(data, filename, varargin)
    % Save the data matrix into a .mat proxy file.
    save([filename '_nifti.mat'], 'data');
end
