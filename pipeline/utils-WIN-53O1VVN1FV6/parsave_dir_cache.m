function parsave_dir_cache(filename, gtv_mask_warped, D_forward, ref3d)
% PARSAVE_DIR_CACHE Helper function to save DIR results inside a parfor loop.
%   MATLAB's 'save' function is not transparent inside parfor loops because
%   variable names in the string arguments are not visible to the static workspace
%   analyzer. This function wraps the save call, taking the variables as inputs.
%
%   Usage:
%       parsave_dir_cache(filename, gtv_mask_warped, D_forward, ref3d);
%
%   Inputs:
%       filename        - String/char describing absolute path of the output .mat file
%       gtv_mask_warped - Follow-up scan mask deformed to original coordinate space
%       D_forward       - Deformation map matrix outputted by demons algorithm
%       ref3d           - Reference coordinate map specifying the original frame
%
%   Outputs:
%       None
%
% --- Analytical Rationale ---
% In longitudinal pancreatic DWI studies, each patient has multiple MRI scans
% acquired at different treatment fractions. To compare tumor regions across
% timepoints, Deformable Image Registration (DIR) aligns follow-up GTV masks
% to the baseline coordinate frame. This is computationally expensive, so
% results are cached to .mat files to avoid redundant recalculation.
%
% MATLAB's parfor loop restricts direct calls to 'save' because the static
% workspace analyzer cannot resolve variable names passed as strings. By
% wrapping save() in an explicit function that receives variables as arguments,
% we satisfy the parfor transparency requirement while enabling parallel
% processing of the patient cohort.
%
% The three cached variables capture the full DIR result:
%   - gtv_mask_warped: The GTV contour from a follow-up scan, warped into the
%     baseline scan's coordinate space. This allows voxel-by-voxel comparison
%     of diffusion parameters (ADC, D, f, D*) within consistent anatomy.
%   - D_forward: The dense displacement field from the demons registration
%     algorithm, mapping every voxel from baseline to follow-up space. Stored
%     for potential re-application to other volumes (e.g., dose maps) without
%     re-running registration.
%   - ref3d: The spatial reference object (imref3d) defining the baseline
%     image's world-coordinate system (voxel size, origin, dimensions). Required
%     to correctly interpret the warped mask and displacement field in physical
%     coordinates (mm), which is essential for accurate dosimetric correlation.
%
    % --- Perform the Cached Save ---
    % By receiving each variable as a function argument, MATLAB's parfor
    % static analyzer can verify that no workspace transparency violations
    % occur. The save call inside this wrapper is functionally identical to
    % calling save() directly, but the function boundary satisfies the
    % requirement that all variables referenced by save() are explicitly
    % present in the local scope.
    %
    % The '-v7.3' flag is intentionally omitted to use the default MAT-file
    % version (v7), which is faster for the relatively small 3D mask arrays
    % typical of pancreatic GTV contours (~64x64x30 voxels). Version 7.3
    % (HDF5-based) would add unnecessary overhead for arrays well under 2 GB.
    save(filename, 'gtv_mask_warped', 'D_forward', 'ref3d');

end
