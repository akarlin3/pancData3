% NIFTIREAD  Octave-compatible shim for MATLAB's niftiread (R2017b+).
%
%   data = niftiread(filename) reads a NIfTI image volume.
%   data = niftiread(info) reads using a niftiinfo struct's .Filename field.
%
%   MATLAB's niftiread parses the .nii/.nii.gz binary format directly.
%   Octave does not include NIfTI I/O. This shim uses a .mat file fallback:
%   it looks for a file named '<filename>_nifti.mat' containing a variable
%   'data'. This convention pairs with the niftiwrite shim which saves data
%   in this format.
%
%   Behavioral differences from MATLAB's niftiread:
%   - Cannot read actual .nii or .nii.gz files; only the _nifti.mat proxy.
%   - The pipeline's DICOM-to-NIfTI conversion (via dcm2niix) produces real
%     .nii files; this shim is only useful for test mocks or pre-converted data.
function data = niftiread(filename_or_info)
    % Accept either a filename string or a niftiinfo struct.
    if isstruct(filename_or_info)
        filename = filename_or_info.Filename;
    else
        filename = filename_or_info;
    end
    % Look for the companion .mat file written by the niftiwrite shim.
    mat_file = [filename '_nifti.mat'];
    if exist(mat_file, 'file')
        tmp = load(mat_file, 'data');
        data = tmp.data;
    else
        error('niftiread:FileNotFound', 'Cannot read NIfTI file: %s', filename);
    end
end
