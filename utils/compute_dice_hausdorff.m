function [dice, hd_max, hd95] = compute_dice_hausdorff(mask_a, mask_b, vox_dims)
% COMPUTE_DICE_HAUSDORFF Compute Dice coefficient and Hausdorff distances
%   between two 3D binary masks with anisotropic voxel support.
%
%   [dice, hd_max, hd95] = compute_dice_hausdorff(mask_a, mask_b, vox_dims)
%
%   Inputs:
%       mask_a    - 3D logical mask (first sub-volume)
%       mask_b    - 3D logical mask (second sub-volume, same size as mask_a)
%       vox_dims  - [dx, dy, dz] voxel dimensions in mm (default: [1 1 1])
%
%   Outputs:
%       dice      - Dice similarity coefficient [0, 1].
%                   Both empty -> NaN; one empty, other non-empty -> 0.
%       hd_max    - Maximum (symmetric) Hausdorff distance in mm.
%                   Both empty -> NaN; one empty -> Inf.
%       hd95      - 95th percentile Hausdorff distance in mm.
%                   Both empty -> NaN; one empty -> Inf.
%
%   ANALYTICAL RATIONALE:
%   Dice measures volumetric overlap: Dice = 2|A & B| / (|A| + |B|).
%   Values near 1 indicate high spatial agreement between repeat scan
%   sub-volumes; values near 0 indicate the sub-volume definition is
%   spatially unstable across same-session acquisitions.
%
%   Hausdorff distance measures the worst-case surface discrepancy in mm.
%   The 95th percentile variant (HD95) is robust to single outlier voxels
%   and is the standard metric for contour comparison in radiation oncology.
%   Both are computed as symmetric (bidirectional) distances:
%     HD(A,B) = max(directed_HD(A->B), directed_HD(B->A))

    if nargin < 3 || isempty(vox_dims)
        vox_dims = [1 1 1];
    end

    mask_a = logical(mask_a);
    mask_b = logical(mask_b);

    sum_a = sum(mask_a(:));
    sum_b = sum(mask_b(:));

    % --- Edge cases ---
    if sum_a == 0 && sum_b == 0
        dice = NaN;
        hd_max = NaN;
        hd95 = NaN;
        return;
    end
    if sum_a == 0 || sum_b == 0
        dice = 0;
        hd_max = Inf;
        hd95 = Inf;
        return;
    end

    % --- Dice coefficient ---
    intersection = sum(mask_a(:) & mask_b(:));
    dice = 2 * intersection / (sum_a + sum_b);

    % --- Hausdorff distance ---
    % Extract surface voxels. bwperim returns the perimeter (boundary)
    % voxels of each mask.  For very small masks where bwperim returns
    % all-false (e.g., a single voxel), fall back to the full mask.
    if exist('bwperim', 'file') || exist('bwperim', 'builtin')
        surf_a = bwperim(mask_a);
        surf_b = bwperim(mask_b);
        if ~any(surf_a(:)), surf_a = mask_a; end
        if ~any(surf_b(:)), surf_b = mask_b; end
    else
        % Octave fallback: use all mask voxels as surface points
        surf_a = mask_a;
        surf_b = mask_b;
    end

    % Convert surface voxel indices to physical coordinates (mm)
    [ra, ca, sa] = ind2sub(size(mask_a), find(surf_a));
    pts_a = [ra(:), ca(:), sa(:)] .* vox_dims(:)';

    [rb, cb, sb] = ind2sub(size(mask_b), find(surf_b));
    pts_b = [rb(:), cb(:), sb(:)] .* vox_dims(:)';

    % Compute pairwise Euclidean distances between surface point sets
    if exist('pdist2', 'file') || exist('pdist2', 'builtin')
        D_mat = pdist2(pts_a, pts_b);
    else
        % Octave fallback: manual pairwise distance computation
        n_a = size(pts_a, 1);
        n_b = size(pts_b, 1);
        D_mat = zeros(n_a, n_b);
        for i = 1:n_a
            diff = pts_b - pts_a(i, :);
            D_mat(i, :) = sqrt(sum(diff .^ 2, 2))';
        end
    end

    % Directed Hausdorff distances:
    %   For each point in A, find minimum distance to any point in B
    %   (and vice versa). The symmetric Hausdorff is the max of both.
    min_a_to_b = min(D_mat, [], 2);   % n_a x 1
    min_b_to_a = min(D_mat, [], 1)';  % n_b x 1

    hd_max = max(max(min_a_to_b), max(min_b_to_a));
    hd95 = max(prctile(min_a_to_b, 95), prctile(min_b_to_a, 95));
end
