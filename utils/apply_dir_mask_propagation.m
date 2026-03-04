function [gtv_mask_warped, D_forward, ref3d] = apply_dir_mask_propagation(b0_fixed, b0_moving, gtv_mask_fixed)
% apply_dir_mask_propagation  Propagates a GTV mask from a reference scan to
%   a follow-up scan using Deformable Image Registration (Demons algorithm).
%
%   gtv_mask_warped = apply_dir_mask_propagation(b0_fixed, b0_moving, gtv_mask_fixed)
%
%   Inputs:
%     b0_fixed       - 3D double array: b=0 DWI image at baseline (Fx1). Used
%                      as the fixed (reference) image for registration.
%     b0_moving      - 3D double array: b=0 DWI image at the target fraction.
%                      Must be the SAME spatial dimensions as b0_fixed.
%     gtv_mask_fixed - 3D logical/binary array: GTV mask at baseline (Fx1),
%                      same dimensions as b0_fixed.
%
%   Output:
%     gtv_mask_warped - 3D logical array: GTV mask propagated to the moving
%                       image space via the estimated deformation field. Empty
%                       ([]) if registration fails or image sizes mismatch.
%     D_forward       - Demons displacement field (same size as b0_moving,
%                       4th dim = 3 for XYZ). Maps baseline to current-fraction
%                       coordinates. Returned so callers can reuse the field to
%                       warp other volumes (e.g., dose maps) without re-running
%                       Demons. To warp a current-fraction volume back into
%                       baseline geometry (e.g., quantitative parameter maps
%                       produced by DnCNN-IVIM fitting in native space), use
%                       the approximate inverse field: -D_forward.
%     ref3d           - imref3d spatial reference object matching b0_moving.
%
%   Algorithm:
%     1. Normalise both b=0 images to [0, 1] for numerically stable Demons.
%     2. Run imregdemons with a 3-level multi-resolution pyramid (coarse-to-
%        fine), which handles large inter-fraction deformations.
%     3. Apply the displacement field to the baseline mask using imwarp with
%        linear interpolation to generate a probability map.
%     4. Threshold the continuous map at 0.5 to binarize, ensuring smooth
%        scaling of the mask boundary.
%
%   Requires: Image Processing Toolbox (imregdemons, imwarp, imref3d).

    gtv_mask_warped = [];
    D_forward       = [];
    ref3d           = [];

    % --- Input validation ---
    if isempty(b0_fixed) || isempty(b0_moving) || isempty(gtv_mask_fixed)
        warning('apply_dir_mask_propagation:emptyInput', ...
            'One or more inputs are empty. Skipping.');
        return;
    end

    if ~isequal(size(b0_fixed), size(b0_moving))
        warning('apply_dir_mask_propagation:sizeMismatch', ...
            'b0_fixed (%s) and b0_moving (%s) have different sizes. Cannot register.', ...
            mat2str(size(b0_fixed)), mat2str(size(b0_moving)));
        return;
    end

    % Warn if voxel spacing information is available and differs.
    % imregdemons operates in voxel coordinates, so anisotropic or
    % mismatched spacing can distort the displacement field.
    if isfield(b0_fixed, 'PixelDimensions') && isfield(b0_moving, 'PixelDimensions')
        if ~isequal(b0_fixed.PixelDimensions, b0_moving.PixelDimensions)
            warning('apply_dir_mask_propagation:spacingMismatch', ...
                'Fixed and moving images have different voxel spacing. Registration accuracy may be degraded.');
        end
    end

    if ~isequal(size(b0_fixed), size(gtv_mask_fixed))
        warning('apply_dir_mask_propagation:sizeMismatch', ...
            'b0_fixed and gtv_mask_fixed have different sizes. Cannot apply mask.');
        return;
    end

    % --- Robust percentile-based normalisation to [0, 1] ---
    % mat2gray uses global min/max, which is highly vulnerable to contrast
    % compression from isolated bright artifacts in DWI data. Instead we
    % clip to the [1st, 99th] percentile of non-zero voxels and linearly
    % rescale, preserving tissue contrast for imregdemons.
    fixed_norm = robust_percentile_norm(double(b0_fixed));
    moving_norm = robust_percentile_norm(double(b0_moving));

    % --- Symmetric Diffeomorphic Registration (Halfway-Space / Midpoint) ---
    has_demons = exist('imregdemons', 'file') || exist('imregdemons', 'builtin');
    if has_demons
        try
            mid_img = (fixed_norm + moving_norm) / 2;
            D_fixed_to_mid = imregdemons(fixed_norm, mid_img, [100 50 25], ...
                'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);
            D_moving_to_mid = imregdemons(moving_norm, mid_img, [100 50 25], ...
                'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);
            fixed_warped = imwarp(fixed_norm, D_fixed_to_mid);
            moving_warped = imwarp(moving_norm, D_moving_to_mid);
            mid_img_refined = (fixed_warped + moving_warped) / 2;
            D_forward_mid = imregdemons(fixed_norm, mid_img_refined, [100 50 25], ...
                'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);
            D_backward_mid = imregdemons(moving_norm, mid_img_refined, [100 50 25], ...
                'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);
            % Compose displacement fields: D_forward(x) = D_fwd(x) + (-D_bwd)(x + D_fwd(x))
            % Simple subtraction is a linear approximation that degrades for
            % the several-mm inter-fraction deformations common in pancreas.
            D_forward = compose_displacement_fields(D_forward_mid, -D_backward_mid);
            ref3d = imref3d(size(b0_moving));
            mask_warped_float = imwarp(double(gtv_mask_fixed), D_forward, ...
                'Interp', 'linear', 'FillValues', 0);
        catch ME
            warning('apply_dir_mask_propagation:imregdemonsFailed', '%s', ME.message);
            return;
        end
    else
        warning('apply_dir_mask_propagation:noDemons', ...
            'imregdemons not available. Using identity transform (no deformation).');
        sz = size(b0_moving);
        D_forward = zeros([sz, 3]);
        ref3d = struct('ImageSize', sz);
        mask_warped_float = double(gtv_mask_fixed);
    end


    % Threshold and return as logical
    gtv_mask_warped = logical(mask_warped_float >= 0.5);
end

function img_norm = robust_percentile_norm(img)
% robust_percentile_norm  Normalise a 3-D image to [0, 1] using the 1st and
%   99th percentiles of non-zero voxels. This avoids contrast compression
%   caused by isolated bright artifacts that would dominate a simple
%   min/max (mat2gray) rescaling.

    nz = img(img ~= 0);
    if isempty(nz)
        img_norm = zeros(size(img));
        return;
    end

    p = prctile(nz, [1 99]);
    lo = p(1);
    hi = p(2);

    if hi <= lo
        % Degenerate case: all non-zero voxels have the same value.
        img_norm = zeros(size(img));
        return;
    end

    img_norm = (img - lo) / (hi - lo);
    img_norm = max(min(img_norm, 1), 0);  % clip to [0, 1]
end

function D_out = compose_displacement_fields(D1, D2)
% compose_displacement_fields  Compose two 3-D displacement fields.
%   D_out(x) = D1(x) + D2(x + D1(x))
%   where D2 is evaluated at the displaced coordinates via trilinear
%   interpolation.  This is more accurate than simple addition or
%   subtraction for deformations larger than ~1 voxel.
%
%   imregdemons returns displacement fields in meshgrid convention:
%     D(:,:,:,1) = column (X) displacement
%     D(:,:,:,2) = row    (Y) displacement
%     D(:,:,:,3) = slice  (Z) displacement
%   interpn indexes arrays in ndgrid order (row, col, slice), so we must
%   add D(:,:,:,2) to the row coordinate and D(:,:,:,1) to the column
%   coordinate.

    sz = size(D1);
    [R, C, S] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));

    % Displaced coordinates after applying D1 (meshgrid → ndgrid mapping)
    R1 = R + D1(:,:,:,2);   % rows displaced by Y component
    C1 = C + D1(:,:,:,1);   % cols displaced by X component
    S1 = S + D1(:,:,:,3);   % slices displaced by Z component

    D_out = D1;
    for dim = 1:3
        D_out(:,:,:,dim) = D1(:,:,:,dim) + ...
            interpn(D2(:,:,:,dim), R1, C1, S1, 'linear', 0);
    end
end
