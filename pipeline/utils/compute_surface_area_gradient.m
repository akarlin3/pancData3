function surface_area = compute_surface_area_gradient(mask, voxel_spacing)
%COMPUTE_SURFACE_AREA_GRADIENT  Gradient-based surface area for 3D binary masks.
%   Uses the magnitude of the gradient of the distance transform at the
%   boundary to estimate surface area, weighted by voxel face areas.

    mask = logical(mask);
    dx = voxel_spacing(1);
    dy = voxel_spacing(2);
    dz = voxel_spacing(3);

    % Detect boundary voxels using 6-connectivity erosion
    se = false(3,3,3);
    se(2,2,1) = true; se(2,2,3) = true;
    se(2,1,2) = true; se(2,3,2) = true;
    se(1,2,2) = true; se(3,2,2) = true;
    se(2,2,2) = true;

    eroded = imerode(mask, se);
    boundary = mask & ~eroded;

    % Each boundary voxel contributes face area for each exposed face
    surface_area = 0;
    [rows, cols, slices] = ind2sub(size(mask), find(boundary));
    for idx = 1:length(rows)
        r = rows(idx); c = cols(idx); s = slices(idx);
        % Check each of 6 neighbors; exposed face contributes area
        if r == 1 || ~mask(r-1, c, s), surface_area = surface_area + dy * dz; end
        if r == size(mask,1) || ~mask(r+1, c, s), surface_area = surface_area + dy * dz; end
        if c == 1 || ~mask(r, c-1, s), surface_area = surface_area + dx * dz; end
        if c == size(mask,2) || ~mask(r, c+1, s), surface_area = surface_area + dx * dz; end
        if s == 1 || ~mask(r, c, s-1), surface_area = surface_area + dx * dy; end
        if s == size(mask,3) || ~mask(r, c, s+1), surface_area = surface_area + dx * dy; end
    end
end
