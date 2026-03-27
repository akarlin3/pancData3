function [gtv_mask_warped, D_forward, ref3d] = apply_dir_mask_propagation(b0_fixed, b0_moving, gtv_mask_fixed, varargin)
% apply_dir_mask_propagation  Propagates a GTV mask from a reference scan to
%   a follow-up scan using Deformable Image Registration (Demons algorithm).
%
%   gtv_mask_warped = apply_dir_mask_propagation(b0_fixed, b0_moving, gtv_mask_fixed)
%   [gtv_mask_warped, D_forward, ref3d] = apply_dir_mask_propagation(b0_fixed, b0_moving, gtv_mask_fixed, 'Name', Value, ...)
%
%   Inputs:
%     b0_fixed       - 3D double array: b=0 DWI image at baseline (Fx1). Used
%                      as the fixed (reference) image for registration.
%     b0_moving      - 3D double array: b=0 DWI image at the target fraction.
%                      Must be the SAME spatial dimensions as b0_fixed.
%     gtv_mask_fixed - 3D logical/binary array: GTV mask at baseline (Fx1),
%                      same dimensions as b0_fixed.
%
%   Optional Name-Value Arguments:
%     'FixedPixelSpacing'    - [1x3] double: voxel spacing of fixed image [x y z] in mm
%     'MovingPixelSpacing'   - [1x3] double: voxel spacing of moving image [x y z] in mm
%     'FixedOrientation'     - [3x3] double: orientation matrix of fixed image
%     'MovingOrientation'    - [3x3] double: orientation matrix of moving image
%     'SpacingTolerance'     - double: tolerance for voxel spacing compatibility (default: 0.1 mm)
%     'OrientationTolerance' - double: tolerance for orientation matrix compatibility (default: 0.1)
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
%     1. Validate voxel spacing and orientation compatibility between images.
%     2. Normalise both b=0 images to [0, 1] for numerically stable Demons.
%     3. Run imregdemons with a 3-level multi-resolution pyramid (coarse-to-
%        fine), which handles large inter-fraction deformations.
%     4. Apply the displacement field to the baseline mask using imwarp with
%        linear interpolation to generate a probability map.
%     5. Threshold the continuous map at 0.5 to binarize, ensuring smooth
%        scaling of the mask boundary.
%
%   Requires: Image Processing Toolbox (imregdemons, imwarp, imref3d).

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'b0_fixed', @(x) isnumeric(x) && (isempty(x) || ndims(x) == 3));
    addRequired(p, 'b0_moving', @(x) isnumeric(x) && (isempty(x) || ndims(x) == 3));
    addRequired(p, 'gtv_mask_fixed', @(x) (isnumeric(x) || islogical(x)) && (isempty(x) || ndims(x) == 3));
    addParameter(p, 'FixedPixelSpacing', [], @(x) isnumeric(x) && length(x) == 3);
    addParameter(p, 'MovingPixelSpacing', [], @(x) isnumeric(x) && length(x) == 3);
    addParameter(p, 'FixedOrientation', [], @(x) isnumeric(x) && isequal(size(x), [3 3]));
    addParameter(p, 'MovingOrientation', [], @(x) isnumeric(x) && isequal(size(x), [3 3]));
    addParameter(p, 'SpacingTolerance', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'OrientationTolerance', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    
    parse(p, b0_fixed, b0_moving, gtv_mask_fixed, varargin{:});
    
    fixed_spacing = p.Results.FixedPixelSpacing;
    moving_spacing = p.Results.MovingPixelSpacing;
    fixed_orientation = p.Results.FixedOrientation;
    moving_orientation = p.Results.MovingOrientation;
    spacing_tol = p.Results.SpacingTolerance;
    orientation_tol = p.Results.OrientationTolerance;

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

    if ~isequal(size(b0_fixed), size(gtv_mask_fixed))
        warning('apply_dir_mask_propagation:sizeMismatch', ...
            'b0_fixed and gtv_mask_fixed have different sizes. Cannot apply mask.');
        return;
    end

    % --- Voxel spacing and orientation validation ---
    if ~isempty(fixed_spacing) && ~isempty(moving_spacing)
        spacing_diff = abs(fixed_spacing - moving_spacing);
        if any(spacing_diff > spacing_tol)
            warning('apply_dir_mask_propagation:spacingMismatch', ...
                'Voxel spacing mismatch exceeds tolerance. Fixed: [%.3f %.3f %.3f] mm, Moving: [%.3f %.3f %.3f] mm, Max diff: %.3f mm > %.3f mm tolerance.', ...
                fixed_spacing(1), fixed_spacing(2), fixed_spacing(3), ...
                moving_spacing(1), moving_spacing(2), moving_spacing(3), ...
                max(spacing_diff), spacing_tol);
            return;
        end
    end

    if ~isempty(fixed_orientation) && ~isempty(moving_orientation)
        orientation_diff = abs(fixed_orientation - moving_orientation);
        if any(orientation_diff(:) > orientation_tol)
            warning('apply_dir_mask_propagation:orientationMismatch', ...
                'Orientation matrix mismatch exceeds tolerance. Max difference: %.3f > %.3f tolerance.', ...
                max(orientation_diff(:)), orientation_tol);
            return;
        end
        
        % Additional check for orthogonality of orientation matrices
        if ~is_orthogonal_matrix(fixed_orientation, orientation_tol) || ...
           ~is_orthogonal_matrix(moving_orientation, orientation_tol)
            warning('apply_dir_mask_propagation:nonOrthogonalOrientation', ...
                'One or both orientation matrices are not orthogonal within tolerance.');
            return;
        end
    end

    % --- Robust percentile-based normalisation to [0, 1] ---
    % DWI b=0 images can contain isolated bright voxels from susceptibility
    % artifacts (common near the air-tissue interface in the pancreas/
    % duodenum region) or coil sensitivity hotspots.  mat2gray uses global
    % min/max, which would compress the tissue signal range to a narrow
    % band, degrading the intensity-based cost function in imregdemons.
    % Instead we clip to the [1st, 99th] percentile of non-zero voxels
    % and linearly rescale, preserving the tissue contrast that drives
    % accurate deformation estimation.
    fixed_norm = robust_percentile_norm(double(b0_fixed));
    moving_norm = robust_percentile_norm(double(b0_moving));

    % --- Symmetric Diffeomorphic Registration (Halfway-Space / Midpoint) ---
    % The pancreas is a highly mobile organ that undergoes significant
    % inter-fraction deformation due to respiratory motion, stomach/
    % duodenal filling, and treatment-induced tumor shrinkage.  Standard
    % (asymmetric) Demons registration is biased: the result depends on
    % which image is designated as fixed vs. moving.  Symmetric registration
    % via a halfway-space (midpoint image) eliminates this bias by treating
    % both images equally, producing a more physically plausible deformation
    % field for propagating the GTV contour.
    %
    % The algorithm:
    %   1. Estimate an initial midpoint image as the voxel-wise average.
    %   2. Register both images to the midpoint to refine it.
    %   3. Re-register both images to the refined midpoint.
    %   4. Compose the forward and inverse (negated backward) fields to get
    %      a single fixed-to-moving transformation.
    %
    % The multi-resolution pyramid [100 50 25] iterations at 3 scales
    % handles large deformations at coarse resolution first (capturing
    % bulk organ motion), then refines at finer scales (capturing local
    % tumor boundary changes).  AccumulatedFieldSmoothing = 1.5 regularizes
    % the deformation field to prevent folding (non-invertible transforms),
    % which would produce physically impossible tissue transformations.
    has_demons = exist('imregdemons', 'file') || exist('imregdemons', 'builtin');
    if has_demons
        try
            % Step 1: Initial midpoint estimate (simple average of the two images)
            mid_img = (fixed_norm + moving_norm) / 2;
            % Step 2: Register both images toward the midpoint
            D_fixed_to_mid = imregdemons(fixed_norm, mid_img, [100 50 25], ...
                'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);
            D_moving_to_mid = imregdemons(moving_norm, mid_img, [100 50 25], ...
                'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);
                
            % Validate transformation matrices
            if ~validate_displacement_field(D_fixed_to_mid) || ~validate_displacement_field(D_moving_to_mid)
                warning('apply_dir_mask_propagation:invalidTransformation', ...
                    'Registration to midpoint produced invalid displacement fields.');
                return;
            end
            
            % Step 3: Refine the midpoint using the warped images
            fixed_warped = imwarp(fixed_norm, D_fixed_to_mid);
            moving_warped = imwarp(moving_norm, D_moving_to_mid);
            mid_img_refined = (fixed_warped + moving_warped) / 2;
            % Step 4: Final registration to the refined midpoint
            D_forward_mid = imregdemons(fixed_norm, mid_img_refined, [100 50 25], ...
                'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);
            D_backward_mid = imregdemons(moving_norm, mid_img_refined, [100 50 25], ...
                'AccumulatedFieldSmoothing', 1.5, 'DisplayWaitBar', false);
                
            % Validate final transformation matrices
            if ~validate_displacement_field(D_forward_mid) || ~validate_displacement_field(D_backward_mid)
                warning('apply_dir_mask_propagation:invalidTransformation', ...
                    'Final registration produced invalid displacement fields.');
                return;
            end
                
            % Compose displacement fields: D_forward(x) = D_fwd(x) + (-D_bwd)(x + D_fwd(x))
            % Simple subtraction is a linear approximation that degrades for
            % the several-mm inter-fraction deformations common in pancreas.
            % Proper composition (via trilinear interpolation of D2 at
            % displaced coordinates) is more accurate for deformations > 1 voxel.
            D_forward = compose_displacement_fields(D_forward_mid, -D_backward_mid);
            
            % Validate composed displacement field
            if ~validate_displacement_field(D_forward)
                warning('apply_dir_mask_propagation:invalidTransformation', ...
                    'Composed displacement field is invalid.');
                return;
            end
            
            ref3d = imref3d(size(b0_moving));
            % Warp the GTV mask using linear interpolation (not nearest-
            % neighbor) to produce a continuous probability map.  This
            % avoids jagged mask boundaries that nearest-neighbor would
            % create, especially for the irregularly-shaped pancreatic GTV.
            mask_warped_float = imwarp(double(gtv_mask_fixed), D_forward, ...
                'Interp', 'linear', 'FillValues', 0);
        catch ME
            warning('apply_dir_mask_propagation:imregdemonsFailed', '%s', ME.message);
            return;
        end
    else
        warning('apply_dir_mask_propagation:noDemons', ...
            'imregdemons not available. Returning NaN mask to prevent downstream analysis on misaligned data.');
        gtv_mask_warped = [];
        D_forward = [];
        return;
    end


    % Threshold the continuous probability map at 0.5 to produce a binary
    % mask.  The 0.5 threshold is the natural decision boundary: voxels
    % with > 50% probability of being inside the GTV (based on the
    % deformed contour) are included.  This approach preserves smooth mask
    % boundaries and prevents systematic over- or under-estimation of GTV
    % volume that would occur with nearest-neighbor interpolation.
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

function is_valid = validate_displacement_field(D)
% validate_displacement_field  Validate that a displacement field represents
%   an invertible transformation by checking local Jacobian determinants and
%   condition numbers.
%
%   Input:
%     D - 4D displacement field array (sz(1) x sz(2) x sz(3) x 3)
%
%   Output:
%     is_valid - logical scalar, true if transformation is valid

    is_valid = false;
    
    if isempty(D) || size(D, 4) ~= 3
        return;
    end
    
    sz = size(D);
    
    % Sample points for validation (every 4th voxel to balance accuracy vs. speed)
    sample_step = 4;
    [r_idx, c_idx, s_idx] = ndgrid(2:sample_step:sz(1)-1, ...
                                   2:sample_step:sz(2)-1, ...
                                   2:sample_step:sz(3)-1);
    
    n_samples = numel(r_idx);
    determinants = zeros(n_samples, 1);
    condition_numbers = zeros(n_samples, 1);
    
    for i = 1:n_samples
        r = r_idx(i);
        c = c_idx(i);
        s = s_idx(i);
        
        % Compute local Jacobian matrix using central differences
        % J = I + grad(D), where I is identity and grad(D) is displacement gradient
        J = eye(3);
        
        % X component gradients
        J(1,1) = J(1,1) + (D(r, c+1, s, 1) - D(r, c-1, s, 1)) / 2;  % dDx/dx
        J(2,1) = J(2,1) + (D(r+1, c, s, 1) - D(r-1, c, s, 1)) / 2;  % dDx/dy
        J(3,1) = J(3,1) + (D(r, c, s+1, 1) - D(r, c, s-1, 1)) / 2;  % dDx/dz
        
        % Y component gradients
        J(1,2) = J(1,2) + (D(r, c+1, s, 2) - D(r, c-1, s, 2)) / 2;  % dDy/dx
        J(2,2) = J(2,2) + (D(r+1, c, s, 2) - D(r-1, c, s, 2)) / 2;  % dDy/dy
        J(3,2) = J(3,2) + (D(r, c, s+1, 2) - D(r, c, s-1, 2)) / 2;  % dDy/dz
        
        % Z component gradients
        J(1,3) = J(1,3) + (D(r, c+1, s, 3) - D(r, c-1, s, 3)) / 2;  % dDz/dx
        J(2,3) = J(2,3) + (D(r+1, c, s, 3) - D(r-1, c, s, 3)) / 2;  % dDz/dy
        J(3,3) = J(3,3) + (D(r, c, s+1, 3) - D(r, c, s-1, 3)) / 2;  % dDz/dz
        
        determinants(i) = det(J);
        condition_numbers(i) = cond(J);
    end
    
    % Check validation criteria
    min_det = min(determinants);
    max_cond = max(condition_numbers);
    
    % Transformation is invalid if:
    % 1. Any Jacobian determinant <= 0 (indicates folding/non-invertible regions)
    % 2. Any condition number > 1000 (indicates severe ill-conditioning)
    % 3. More than 5% of sampled points have determinant < 0.1 (excessive compression)
    
    fold_threshold = 0.0;
    compression_threshold = 0.1;
    condition_threshold = 1000;
    compression_fraction_limit = 0.05;
    
    has_folding = min_det <= fold_threshold;
    is_poorly_conditioned = max_cond > condition_threshold;
    compression_fraction = sum(determinants < compression_threshold) / n_samples;
    has_excessive_compression = compression_fraction > compression_fraction_limit;
    
    if has_folding
        warning('apply_dir_mask_propagation:transformation_folding', ...
            'Displacement field contains folding (min determinant: %.3f)', min_det);
        return;
    end
    
    if is_poorly_conditioned
        warning('apply_dir_mask_propagation:transformation_conditioning', ...
            'Displacement field is poorly conditioned (max condition number: %.1f)', max_cond);
        return;
    end
    
    if has_excessive_compression
        warning('apply_dir_mask_propagation:transformation_compression', ...
            'Displacement field has excessive compression (%.1f%% of points with det < %.2f)', ...
            compression_fraction * 100, compression_threshold);
        return;
    end
    
    is_valid = true;
end

function is_orthogonal = is_orthogonal_matrix(M, tolerance)
% is_orthogonal_matrix  Check if a matrix is orthogonal within tolerance.
%   An orthogonal matrix satisfies M * M' = I (identity matrix).
%
%   Input:
%     M         - [3x3] double: matrix to test
%     tolerance - double: tolerance for identity check
%
%   Output:
%     is_orthogonal - logical scalar, true if matrix is orthogonal

    if ~isequal(size(M), [3 3])
        is_orthogonal = false;
        return;
    end
    
    % Check if M * M' is close to identity matrix
    product = M * M';
    identity_diff = abs(product - eye(3));
    
    is_orthogonal = all(identity_diff(:) <= tolerance);
end