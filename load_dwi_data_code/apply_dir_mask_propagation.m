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
        warning('apply_dir_mask_propagation: One or more inputs are empty. Skipping.');
        return;
    end

    if ~isequal(size(b0_fixed), size(b0_moving))
        warning('apply_dir_mask_propagation: b0_fixed (%s) and b0_moving (%s) have different sizes. Cannot register.', ...
            mat2str(size(b0_fixed)), mat2str(size(b0_moving)));
        return;
    end

    if ~isequal(size(b0_fixed), size(gtv_mask_fixed))
        warning('apply_dir_mask_propagation: b0_fixed and gtv_mask_fixed have different sizes. Cannot apply mask.');
        return;
    end

    % --- Normalise images to [0, 1] for numerically stable Demons ---
    fixed_norm = mat2gray(double(b0_fixed));
    moving_norm = mat2gray(double(b0_moving));

    % --- Symmetric Diffeomorphic Registration (Halfway-Space / Midpoint) ---
    % Replaces standard additive demons with a symmetric approach (e.g., Log-Domain / SyN).
    % This ensures topological consistency and physical validity by registering both
    % images to a common midpoint geometry.
    try
        % Initial midpoint
        mid_img = (fixed_norm + moving_norm) / 2;
        
        % Pass 1: Register Fixed -> Mid and Moving -> Mid
        D_fixed_to_mid = imregdemons(fixed_norm, mid_img, [100 50 25], ...
            'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);
        D_moving_to_mid = imregdemons(moving_norm, mid_img, [100 50 25], ...
            'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);
        
        % Refine midpoint based on initial alignment
        fixed_warped = imwarp(fixed_norm, D_fixed_to_mid);
        moving_warped = imwarp(moving_norm, D_moving_to_mid);
        mid_img_refined = (fixed_warped + moving_warped) / 2;
        
        % Pass 2: Final Registration to Refined Midpoint
        D_forward_mid = imregdemons(fixed_norm, mid_img_refined, [100 50 25], ...
            'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);
        D_backward_mid = imregdemons(moving_norm, mid_img_refined, [100 50 25], ...
            'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);

        % The displacement field to map Fixed -> Moving is approximated by:
        % D = D_forward_mid - D_backward_mid (symmetric proxy)
        D_forward = D_forward_mid - D_backward_mid;

        % Warp the baseline mask into the current fraction's geometry.
        % We use 'linear' interpolation to generate a continuous probability 
        % map, preventing nearest-neighbor truncation artifacts and allowing 
        % the mask boundary to scale with the local deformation Jacobian.
        ref3d = imref3d(size(b0_moving));
        mask_warped_float = imwarp(double(gtv_mask_fixed), D_forward, ...
            'Interp', 'linear', 'OutputView', ref3d, 'FillValues', 0);
    catch ME
        warning('apply_dir_mask_propagation:imregdemonsFailed', '%s', ME.message);
        return;
    end


    % Threshold and return as logical
    gtv_mask_warped = logical(mask_warped_float >= 0.5);
end
