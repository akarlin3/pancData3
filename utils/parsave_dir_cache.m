function parsave_dir_cache(filename, gtv_mask_warped, D_forward, ref3d)
% PARSAVE_DIR_CACHE Helper function to save DIR results inside a parfor loop.
%   MATLAB's 'save' function is not transparent inside parfor loops because
%   variable names in the string arguments are not visible to the static workspace
%   analyzer. This function wraps the save call, taking the variables as inputs.
%
%   Usage:
%       parsave_dir_cache(filename, gtv_mask_warped, D_forward, ref3d);

    save(filename, 'gtv_mask_warped', 'D_forward', 'ref3d');
end
