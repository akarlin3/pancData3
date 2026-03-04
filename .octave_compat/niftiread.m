function data = niftiread(filename_or_info)
% NIFTIREAD Octave-compatible shim for MATLAB's niftiread.
%   Accepts a filename string or a niftiinfo struct (with .Filename field).
    if isstruct(filename_or_info)
        filename = filename_or_info.Filename;
    else
        filename = filename_or_info;
    end
    mat_file = [filename '_nifti.mat'];
    if exist(mat_file, 'file')
        tmp = load(mat_file, 'data');
        data = tmp.data;
    else
        error('niftiread:FileNotFound', 'Cannot read NIfTI file: %s', filename);
    end
end
