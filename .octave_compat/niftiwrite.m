function niftiwrite(data, filename, varargin)
% NIFTIWRITE Octave-compatible shim for MATLAB's niftiwrite.
%   Saves data as a .mat file with _nifti suffix as a fallback.
    save([filename '_nifti.mat'], 'data');
end
