function elongation = compute_elongation_optimized(mask, voxel_spacing)
%COMPUTE_ELONGATION_OPTIMIZED  Elongation from covariance matrix eigenvalues.
%   Elongation = 1 - sqrt(lambda_min / lambda_max), where lambda are the
%   eigenvalues of the spatial covariance matrix of mask voxel coordinates.

    elongation = NaN;
    mask = logical(mask);
    if ~any(mask(:))
        return;
    end

    is_3d = ndims(mask) == 3 && size(mask, 3) > 1;

    if is_3d
        [r, c, s] = ind2sub(size(mask), find(mask));
        coords = [r * voxel_spacing(1), c * voxel_spacing(2), s * voxel_spacing(3)];
    else
        [r, c] = find(mask);
        coords = [r * voxel_spacing(1), c * voxel_spacing(2)];
    end

    if size(coords, 1) < 2
        elongation = 0;
        return;
    end

    C = cov(coords);
    eigenvalues = sort(eig(C), 'descend');
    eigenvalues = eigenvalues(eigenvalues > 0);

    if isempty(eigenvalues) || eigenvalues(1) == 0
        elongation = 0;
    else
        elongation = 1 - sqrt(eigenvalues(end) / eigenvalues(1));
    end
end
