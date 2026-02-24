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
%                       4th dim = 3 for XYZ). Returned so callers can reuse
%                       the field to warp other volumes (e.g., dose maps)
%                       without re-running Demons.
%     ref3d           - imref3d spatial reference object matching b0_moving.
%
%   Algorithm:
%     1. Normalise both b=0 images to [0, 1] for numerically stable Demons.
%     2. Run imregdemons with a 3-level multi-resolution pyramid (coarse-to-
%        fine), which handles large inter-fraction deformations.
%     3. Apply the displacement field to the baseline mask using imwarp with
%        nearest-neighbour interpolation to preserve binary values.
%     4. Threshold at 0.5 and return as logical array.
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

    % --- Run Demons non-rigid registration ---
    % 3-level pyramid: iterations = [200 100 50] (coarse → fine)
    % AccumulatedFieldSmoothing: regularises the field to prevent folding.
    try
        % Calculate deformation from fixed (Fx1) to moving (current Fx)
        [D_forward, ~] = imregdemons(fixed_norm, moving_norm, [200 100 50], ...
            'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);

        % Warp the baseline mask into the current fraction's geometry
        ref3d = imref3d(size(b0_moving));
        mask_warped_float = imwarp(double(gtv_mask_fixed), D_forward, ...
            'Interp', 'nearest', 'OutputView', ref3d, 'FillValues', 0);
    catch ME
        warning('apply_dir_mask_propagation:imregdemonsFailed', '%s', ME.message);
        return;
    end


    % Threshold and return as logical
    gtv_mask_warped = logical(mask_warped_float >= 0.5);
end
